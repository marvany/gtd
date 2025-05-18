    # ==================================== #
    #       BUILD LOGISTIC CLASSIFIER      #
    # ==================================== #
source('/sc/arion/projects/va-biobank/PROJECTS/ma_extras/personal/gtd/gtd/classification_helper.R')
set.seed(123)

dt <- fread('Resources/globalterrorismdb.csv')
backup <- copy(dt)

y       = "doubtterr"     # <- name of outcome column (0/1 or FALSE/TRUE)

ord_vars <- c(
  "weapsubtype1", "targtype1_txt", "targsubtype1",
  "natlty1", "property"
)

ord_vars <- c("weapsubtype1", "targtype1_txt", "targsubtype1", "natlty1",
"property", "multiple", "dbsource", "success", "country", "guncertain1")


dt <- dt[, c(y, ord_vars), with = FALSE]

dt <- remove_neg9_features(dt)
dt <- drop_sparse_strings(dt)
char_cols <- names(which(sapply(dt, is.character)))
dt[, (char_cols) := lapply(.SD, function(x) as.integer(factor(x, levels = unique(x)))), 
   .SDcols = char_cols]

dt <- na.omit(dt)
dt <- drop_bad_rows(dt)
dt <- dt[, (names(dt)) := lapply(.SD, as.numeric), .SDcols = names(dt)]

#ord_vars <- c(
#  "weapsubtype1",  "targtype1_txt", "targsubtype1", "natlty1",
#  "property",      "multiple",      "dbsource",     "success",
#  "country",       "guncertain1"
#)

k       = 10           # <- # folds (use 5 or 10 by convention)
seed    = 1L           # <- reproducible shuffling
thresh  = 0.5        # <- decision threshold

form <- as.formula(paste(y, "~", paste(setdiff(names(dt), y), collapse = " + ")))

folds <- sample(rep_len(seq_len(k), nrow(dt)))

all_lg_results <- rbindlist(
  pbmclapply(
    seq_len(k), 
    mc.cores = parallel::detectCores() - 2,
    function(i) {
      
      train <- scale(dt[ folds != i])
      train <- restore_doubtterr(train)
      test  <- scale(dt[ folds == i ])
      test <- restore_doubtterr(test)
      
      fit <- glm(
        formula = form,
        data    = as.data.frame(train),
        family  = binomial()
      )
      
      ## predicted probabilities & class labels
      pr   <- predict(fit, newdata = as.data.frame(test), type = "response")
      pred <- as.integer(pr >= thresh)
      tru  <- as.integer(as.data.frame(test)[[y]])
      
      data.table(
                  fold_id = i,
                  #row_id  = row_ids,
                  tru     = tru,
                  pred    = pred
                )
      }
    )
  )


## ──────────────────────────────────────────────────────────────
## 5.  calculate performance ------------------------------------
performance_results <- all_lg_results[,.(
    auc_roc   = as.numeric(pROC::auc(tru, pred)),
    auc_pr    = PRROC::pr.curve(scores.class0 = pred[tru == 1],
                                scores.class1 = pred[tru == 0])$auc.integral,
    logloss   = -mean(tru * log(pred) + (1 - tru) * log(1 - pred)),
    accuracy  = mean(pred == tru)),
    by = fold_id]

performance_results
fwrite(performance_results, 'lg/performance_results.csv')

## ──────────────────────────────────────────────────────────────
## 5.  aggregate performance -----------------------------------
summary <- performance_results[, .(auc_roc, auc_pr,logloss,accuracy)][, .(
  mean_auc_roc = mean(auc_roc), sd_auc_roc = sd(auc_roc),
  mean_auc_pr  = mean(auc_pr),  sd_auc_pr  = sd(auc_pr),
  mean_logloss = mean(logloss), sd_logloss = sd(logloss),
  mean_acc     = mean(accuracy),sd_acc     = sd(accuracy)
)]
summary
fwrite(summary, 'lg/final_model_summary.csv')

## ──────────────────────────────────────────────────────────────
## 6.   calculate F1
f1 <- f1_from_dt(all_lg_results[, .(tru,pred), by= fold_id])
f1
fwrite(f1, 'lg/f1_score.csv')

