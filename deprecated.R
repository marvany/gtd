



    # ========================================= #
    #       BUILD RANDOM FOREST CLASSIFIER      #
    # ========================================= #








    # ===================== #
    #       FINAL TABLE     #
    # ===================== #


text = paste('\nROWS WITH -9 VALUES (N =', dt[doubtterr == -9, .N],'), WERE FILTERED OUT.')
writeReport(text)
dt <- dt[!is.na(doubtterr) & (doubtterr != -9)]

dt <- 
text = paste('\nAFTER FILTERING OUT FOR ALL MISSING VALUES ')

## -------- 2.  UNIVARIATE SCREENING  --------------------------------------
binary_y <- as.factor(dt[[target]])
pred_cols <- setdiff(names(dt), target)

cat_vars <- pred_cols[sapply(dt[, ..pred_cols], function(x) is.character(x) || is.factor(x))]
num_vars <- setdiff(pred_cols, cat_vars)

# χ² / Fisher for categoricals
chi_p <- sapply(cat_vars, function(col) {
  tbl <- table(dt[[col]], binary_y)
  if (min(dim(tbl)) < 2) return(NA_real_)
  suppressWarnings(chisq.test(tbl)$p.value)
})
chi_rank <- sort(chi_p, na.last = NA)

# Wilcoxon for numerics
wil_p <- sapply(num_vars, function(col) {
  x <- dt[[col]]
  if (all(is.na(x))) return(NA_real_)
  suppressWarnings(wilcox.test(x ~ binary_y)$p.value)
})
wil_rank <- sort(wil_p, na.last = NA)

cat("\nTop-10 categorical (χ² p-value):\n"); print(head(chi_rank, 10))
cat("\nTop-10 numeric (Wilcoxon p-value):\n"); print(head(wil_rank, 10))

## -------- 3.  MODEL-BASED IMPORTANCE (Random Forest)  --------------------
prep <- copy(dt)

# replace GTD sentinel –99 with NA in numerics
sentinel <- -99
num_cols <- names(prep)[sapply(prep, is.numeric)]
for (c in num_cols) set(prep, which(prep[[c]] == sentinel), c, NA)

# collapse ultra-rare levels for big categoricals
hi_card <- cat_vars[sapply(prep[, ..cat_vars], function(x) uniqueN(x, na.rm = TRUE) > 50)]
for (c in hi_card) prep[[c]] <- fct_other(prep[[c]], keep = names(sort(table(prep[[c]]), TRUE)[1:30]))

# one-hot encode categoricals
X_cat <- one_hot(as.data.table(prep[, ..cat_vars]), sparsifyNAs = TRUE)
X_num <- as.matrix(prep[, ..num_vars])
X_num[is.na(X_num)] <- 0
X <- cbind(X_cat, X_num)

# drop rows with NA target
keep <- !is.na(prep[[target]])
X <- X[keep, ]
y <- prep[[target]][keep]

rf <- ranger(
  x = X,
  y = as.factor(y),
  num.trees = 400,
  importance = "permutation",
  classification = TRUE,
  respect.unordered.factors = "partition",
  num.threads = parallel::detectCores()
)

cat(sprintf("\nOOB error rate: %.3f\n", rf$prediction.error))

imp <- sort(rf$variable.importance, TRUE)
top20 <- head(imp, 20)
print(round(top20, 4))

vip(rf, num_features = 20) +
  ggtitle("Random-Forest Permutation Importance")

## -------- 4.  QUICK REFINEMENT LOOP  -------------------------------------
drop_low <- names(imp)[imp < quantile(imp, 0.10)]   # bottom 10 %
keep_mat <- X[, setdiff(colnames(X), drop_low)]

rf2 <- ranger(
  x = keep_mat,
  y = as.factor(y),
  num.trees = 400,
  importance = "permutation",
  classification = TRUE,
  respect.unordered.factors = "partition",
  num.threads = parallel::detectCores()
)

cat(sprintf("\nRefined OOB error: %.3f\n", rf2$prediction.error))

vip(rf2, num_features = 20) +
  ggtitle("Refined Feature Set Importance")

## -------- 5.  EXPORT RANKINGS  ------------------------------------------
fwrite(data.table(feature = names(chi_rank), chi_p = chi_rank), "chi2_rank.csv")
fwrite(data.table(feature = names(wil_rank), wilcoxon_p = wil_rank), "wilcox_rank.csv")
fwrite(data.table(feature = names(imp), rf_perm_imp = imp), "rf_importance.csv")










####################################################################################
####################################################################################

prepare_dt <- function(dt){
  dt <- remove_neg9_features(dt)
  dt <- drop_sparse_strings(dt)
  char_cols <- names(which(sapply(dt, is.character)))
  dt[, (char_cols) := lapply(.SD, function(x) as.integer(factor(x, levels = unique(x)))), 
    .SDcols = char_cols]

  dt <- na.omit(dt)
  dt <- drop_bad_rows(dt)
  dt <- dt[, (names(dt)) := lapply(.SD, as.numeric), .SDcols = names(dt)]
  dt
}


set.seed(123)





ord_vars <- c(
  "weapsubtype1",  "targtype1_txt", "targsubtype1", "natlty1",
  "property", "nkill", "nwound", "success", "ishostkid", 'suicide',
  "country",       "guncertain1", 'natlty1_txt', 'weaptype1' 
)
ord_vars <- c(
  "weapsubtype1",  "targtype1_txt", "targsubtype1", "natlty1",
  "property"
)

ord_vars <- c(
  "weapsubtype1",  "targtype1_txt", "targsubtype1", "natlty1",
  "property", "country", "guncertain1", 'natlty1_txt', 'weaptype1' 
)

dt <- copy(backup)
dt <- prepare_dt(dt)
names(dt)

y         <- "doubtterr"                 # 0 / 1
pred_cols <- setdiff(ord_vars, y)
dt[, (pred_cols) := lapply(.SD, as.numeric), .SDcols = pred_cols]
dt <- dt[, c(y, pred_cols), with = FALSE]
names(dt)


## ------------------------------------------------------------------------
## 1.  stratified k-fold indices  ------------------------------------------
k     <- 10
set.seed(1)

dt[, fold := {                # creates a “fold” column
     idx <- sample(.N)        # random shuffle *inside* class
     rep(seq_len(k), length.out = .N)[order(idx)]
   }, by = y]

## ------------------------------------------------------------------------
## 2.  cross-validated fitting & evaluation  -------------------------------
thresh <- 0.5        # decision threshold, *will revisit later*

dt[, doubtterr := as.integer(doubtterr)]

results <- rbindlist(lapply(seq_len(k), \(i) {

  train <- dt[ fold != i ]
  test  <- dt[ fold == i ]

  ## ---- scale predictors *inside* each fold ------------------------------
  mu  <- train[, lapply(.SD, mean), .SDcols = pred_cols]
  s   <- train[, lapply(.SD, sd  ), .SDcols = pred_cols]

  train[, (pred_cols) := Map(\(x,m,s) (x - m)/s, .SD, mu, s),
        .SDcols = pred_cols]
  test [, (pred_cols) := Map(\(x,m,s) (x - m)/s, .SD, mu, s),
        .SDcols = pred_cols]

  ## ---- optional class-weighting  ----------------------------------------
 # w <- ifelse(train[[y]] == 1,
 #             0.5 / mean(train[[y]] == 1),   # ↑ weight for minority class
 #             0.5 / mean(train[[y]] == 0))   # ↓ weight for majority class

  fit <- glm(
    as.formula(paste(y, "~ .")),
    data    = as.data.frame(train[, c(y, pred_cols), with = FALSE]),
  #  weights = w                              # ← ***delete this line***
    family  = binomial()
  )
  

  ## ---- predict & score  --------------------------------------------------
  pr   <- predict(fit, newdata = as.data.frame(test), type = "response")
  tru  <- test[[y]]
  pred <- as.integer(pr >= thresh)

  data.table(
    fold      = i,
    auc_roc   = as.numeric(auc(tru, pr)),
    auc_pr    = pr.curve(scores.class0 = pr[tru == 1],
                         scores.class1 = pr[tru == 0])$auc.integral,
    logloss   = -mean(tru * log(pr) + (1 - tru) * log(1 - pr)),
    prevalence= mean(tru),                    # sanity-check: ≈0.17 in every fold
    accuracy  = mean(pred == tru)             # included only for reference
  )
}))

## ------------------------------------------------------------------------
## 3.  aggregate CV performance  -------------------------------------------
summary <- results[, .(
  mean_auc_roc = mean(auc_roc),   sd_auc_roc = sd(auc_roc),
  mean_auc_pr  = mean(auc_pr ),   sd_auc_pr  = sd(auc_pr ),
  mean_logloss = mean(logloss),   sd_logloss = sd(logloss),
  mean_acc     = mean(accuracy),  sd_acc     = sd(accuracy)
)]
print(summary)



