###########################################################
#  GTD feature-triage template
#  ---------------------------------------------------------
#  • Assumes your data.frame/data.table is called  dt
#  • Pick your target in the  target  variable
#  • Requires: data.table, ranger, vip, mltools, forcats
############################################################
setwd("/sc/arion/projects/va-biobank/PROJECTS/ma_extras/personal/gtd/gtd")
source("classification_helper.R")

## -------- 0.  SET-UP  ----------------------------------------------------
setDTthreads(max(1L, parallel::detectCores() - 1))
dt = fread('Resources/globalterrorismdb.csv')
backup <- copy(dt)


## -------- 1.  TRIAGE / QUICK PROFILE  ------------------------------------
profile <- rbindlist(
  lapply(names(dt), function(col) {
    x <- dt[[col]]
    data.table(
      variable     = col,
      class        = class(x)[1],
      n_unique     = uniqueN(x, na.rm = TRUE),
      pct_missing  = round(mean(is.na(x))*100, 1)
    )
  })
)[order(-pct_missing, -n_unique)]
fwrite(profile, 'results/summary_sizes.tsv', sep = '\t')
print(head(profile, 20))  # glance at worst offenders



    # ================================== #
    #       FIND DEPENDENT VARIABLE      #
    # ================================== #

dt <- backup
# BASED ON COLUMN pct_missing WE FIND THAT WE WANT TO EXCLUDE VARIABLES WITH >30% DATA, SINCE THESE 
# MAY SIGNIFICANTLY IMPACT THE POWER FOR OUR ANALYSIS
# 1. Identify columns with >30 % missing
drop_cols <- profile[pct_missing > 50, variable]

#-------- 2. Remove them from dt  (data.table syntax)
dt[, (drop_cols) := NULL]

#report
writeLines(paste("FOR THE DEPENDENT VARIABLE SELECTION WE FIRST FILTERED OUT FEATURES WITH >50% MISSING VALUES (N=",length(drop_cols),")."), "report.txt")

#-------- 3. Remove variables with >3 unique values
candidate_dependents <- profile[n_unique <= 3, variable]
overlapping = intersect(candidate_dependents, names(dt))
dt <- dt[, ..overlapping]

#report
text = (paste("\nWE THEN PROCEEDED IN FILTERING OUT VARIABLES THAT HAD ENLISTED MORE THAN 3 OUTCOMES LEAVING US IN TOTAL 
WITH N = ",length(overlapping)," FEATURES. THE FULL LIST OF VARIABLES CONSIDERED MAY BE FOUND IN OUR SUPPLEMENTARY TABLE"))
writeReport(text)
text = paste("\nDEPENDENT VARIABLES CONSIDERED:\n",
                paste(overlapping, collapse = "\n"))
writeReport(text, "supplementary.txt") # for supplementary

#-------- 4. Manually inspect
sort(names(dt))
dt[, .N, by = doubtterr] # BINGO

#-------- 5. doubtterr will be the independent variable
text = paste('\nAFTER MANUALLY INSPECTING THE REMAINING VARIABLES (N =',ncol(dt),') WE FOUND THAT VARIABLE doubterr WAS VERY INTERESTING FOR OUR DOWNSTREAM ANALYSES. ')
writeReport(text)
text = paste('\nBRIEFLY doubtterr EXPRESSES THE CERTAINTY ON WHETHER THIS WAS A TERRORIST ATTACK, WITH 0 MEANING THERE IS NO DOUBT THAT THIS WAS INDEED A TERRORIST ATTACK, 1 EXPRESSING UNCERTAINTY AND -9 SUGGESTING UNAVAILABLE INFORMATION.')
writeReport(text)


    # ===================================== #
    #       FIND INDEPENDENT VARIABLES      #
    # ===================================== #

# BACK TO ORIGINAL TABLE WITH FILTERED OUT IN THE > 30 % OF MISSING VALUES
dt <- copy(backup)


#-------- 2. Remove columns with >30 % missing
drop_cols <- profile[pct_missing > 50, variable]
dt[, (drop_cols) := NULL]
#text = paste("SIMILARLY TO WHAT WE DID BEFORE WE FILTER OUT VARIABLES WHOSE MISSING VALUES ARE >50% OF OUR TOTAL DATASET.")
#writeReport(text)


#--------- 3. HANDLE -9 VALUES and ''; a. remove features with many -9
dt <- remove_neg9_features(dt)
dt <- drop_sparse_strings(dt)

#---------- 4. remove rows with NA, -9, '' values
dt <- na.omit(dt)
dt <- drop_bad_rows(dt)

#--------- 5. Turn text to numbers
char_cols <- names(which(sapply(dt, is.character)))
dt[, (char_cols) := lapply(.SD, function(x) as.integer(factor(x, levels = unique(x)))), 
   .SDcols = char_cols]




#--------- 6. Filter out columns with near 0 variance
dt <- remove_zero_var(dt)


#--------- 7. Pairwise-correlation filter (fast, numeric columns only)
dt <- correlation_filter(dt)

if(FALSE){
#--------- 8. Filter out high VIF 
# make sure all are numeric
dt[, (names(dt)) := lapply(.SD, as.numeric), .SDcols = names(dt)]

# outcome name
y <- "doubtterr"

# predictors = every other column
preds <- setdiff(names(dt), y)

# y ~ x1 + x2 + x3 + ...
form <- as.formula(paste(y, "~", paste(preds, collapse = " + ")))

m <- scale(dt)
mod <- lm(form, data = as.data.frame(na.omit(m)), family = 'poisson')   # target = your outcome

vif_vals <- vif(mod)                              # named numeric vector
too_high  <- names(vif_vals[vif_vals > 10])
too_high
}




#-------- 9. Regularize
# First we restore all 
cols_to_keep <- names(dt)
id_to_keep <- dt$eventid
fwrite(data.table(eventid = id_to_keep), 'data_inspection/ids_of_rows_considered.csv')


temp <- copy(backup)
temp <- temp[, ..cols_to_keep]
temp <- temp[eventid %in% id_to_keep]
dim(temp)
dim(dt)
dt <- copy(temp)

class_dt <- data.table(
  variable = names(dt),
  class    = vapply(dt, function(x) paste(class(x), collapse = ", "), character(1L))
)
class_dt


char_cols <- names(which(sapply(dt, is.character)))
dt[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]

# 
class_dt <- data.table(
  variable = names(dt),
  class    = vapply(dt, function(x) paste(class(x), collapse = ", "), character(1L))
)
class_dt

insp.dt = data.table(
  variable = names(dt),
  n_unique = vapply(dt, uniqueN, integer(1L)),
  class = vapply(dt, function(x) paste(class(x), sep = ", "), character(1L))
)
insp.dt

# WE FILTER OUT VARIABLES WITH MORE THAN 40 UNIQUE VALUES FOR COMPUATIONAL SIMPLICITY
to_keep = c(insp.dt[n_unique < 40, variable])
dt <- dt[, ..to_keep]

data.table(
  variable = names(dt),
  n_unique = vapply(dt, uniqueN, integer(1L)),
  class = vapply(dt, class, character(1L))
)

# Build xmat
xmat <- build_xmat(dt)

# Build yvec
yvec <- scale(dt[[ycol]])
yvec <- -yvec

# Build lambda
lambda_grid <- 10^seq(1, -6, length.out = 400)     # 10¹ … 10⁻⁶
fwrite(cbind(yvec, xmat), 'data_inspection/final_predictors_considered.csv') # V1 is the doubterr variable


fit <- glmnet(
  x               = xmat,
  y               = yvec,
  family          = "binomial",
  alpha           = 1,               # pure LASSO
  standardize     = TRUE,
  lambda          = lambda_grid
)

k <- 20

# count how many β ≠ 0 at each λ (skip the intercept row 1)
nz_per_lambda <- colSums(coef(fit)[-1, ] != 0)

# get all λ indices where 1…k variables survive
eligible_idx <- which(nz_per_lambda > 0 & nz_per_lambda <= k)

# if none exist, extend the λ grid automatically
if (length(eligible_idx) == 0) {
  message("No λ produced ≤ ", k, " variables; extending the grid …")
  lambda_grid <- 10^seq(1, -8, length.out = 600)   # go looser
  fit <- glmnet(xmat, yvec, family = "binomial", alpha = 1,
                standardize = TRUE, lambda = lambda_grid)
  nz_per_lambda <- colSums(coef(fit)[-1, ] != 0)
  eligible_idx  <- which(nz_per_lambda > 0 & nz_per_lambda <= k)
  if (length(eligible_idx) == 0)
    stop("Even the extended grid couldn’t satisfy the sparsity target.")
}

# pick the λ with *largest* index among eligible  (i.e. smallest λ, least penalty)
best_idx   <- max(eligible_idx)
best_lambda<- fit$lambda[best_idx]

#── 3.  EXTRACT FINAL COEFFICIENTS & SELECTED FEATURES ────────────────────────
beta_full  <- drop(coef(fit, s = best_lambda))      # includes intercept
beta       <- beta_full[-1]                         # remove intercept
sel_mask   <- beta != 0
sel_vars   <- names(beta)[sel_mask]
sel_beta   <- beta[sel_mask]

cat("\nSelected λ  :", best_lambda,
    "\n# Non-zero β :", length(sel_beta), "\n\nVariables kept:\n")
print(abs(sel_beta[order(-abs(sel_beta))]))

# crit3 is filtered out because it is used to create doubtterr
ord_vars <- abs(sel_beta[order(-abs(sel_beta))])
ord_vars <- setdiff(names(ord_vars), c('crit3', 'gname')) # we remove these because they will introduce bias
ord_vars

writeLines(ord_vars, 'data_inspection/top_20_predictors.csv')
#ord_vars <- ord_vars[1:10] # we keep only the first 10

"after inspecting the ordered variables we see that targtype1_txt , region_txt , dbsource , weapsubtype1_txt , INT_MISC ,
multiple"

# these are the broader variables categories
vars_selected <- c('targtype1_txt' , 'region_txt' , 'dbsource' , 'weapsubtype1_txt' , 'INT_MISC' ,'multiple')



####################################################################################
####################################################################################


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







