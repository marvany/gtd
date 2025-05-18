    # ========================================= #
    #       BUILD A RANDOM FOREST CLASSIFER     #
    # ========================================= #
setwd("/sc/arion/projects/va-biobank/PROJECTS/ma_extras/personal/gtd/gtd")
source("classification_helper.R")
library(data.table)
library(ranger)         # fast RF
library(pbmcapply)      # parallel lapply
library(pROC)           # AUC-ROC
library(PRROC)          # AUC-PR
dt <- fread('Resources/globalterrorismdb.csv')
backup <- copy(dt)

## ──────────────────────────────────────────────────────────────
## 0.  data prep  (your existing helpers) -----------------------
dt <- copy(backup)          # full GTD slice
y         <- "doubtterr"                 
pred_cols <- c(
  "weapsubtype1", "targtype1_txt", "targsubtype1",
  "natlty1", "property"
)

dt <- dt[, c(y, pred_cols), with = FALSE]
dt <- prepare_dt(dt)



## ──────────────────────────────────────────────────────────────
## 1.  outer 10-fold, stratified on y ---------------------------
k_outer <- 10
set.seed(1)

dt[, fold := { idx <- sample(.N)
               rep(seq_len(k_outer), length.out = .N)[order(idx)] },
   by = doubtterr]

## ──────────────────────────────────────────────────────────────
## 2.  hyper-parameter grid (pruned automatically if mtry > p) --
p <- length(pred_cols)
grid <- CJ(
  mtry            = pmax(1L, floor(c(0.3, 0.6, 0.9) * p)),  # 30/60/90 % of p
  min.node.size   = c(1L, 5L, 10L),
  sample.fraction = c(0.6, 0.8, 1.0),
  num.trees       = c(300L, 600L)
)
grid <- unique(grid[mtry <= p])

## ──────────────────────────────────────────────────────────────
## 3.  outer-fold worker ---------------------------------------


## ──────────────────────────────────────────────────────────────
## 4.  run outer CV in parallel --------------------------------
unique(dt$doubtterr)

all_results <- rbindlist(
  lapply(seq_len(k_outer), eval_fold)
)
fwrite(all_results, 'rf/final_results.csv')
results[fold_id == 1,]
nrow(all_results)


## ──────────────────────────────────────────────────────────────
## 5.  calculate performance ------------------------------------
performance_results <- all_results[,.(
    auc_roc   = as.numeric(pROC::auc(tru, pred)),
    auc_pr    = PRROC::pr.curve(scores.class0 = pred[tru == 1],
                                scores.class1 = pred[tru == 0])$auc.integral,
    logloss   = -mean(tru * log(pred) + (1 - tru) * log(1 - pred)),
    accuracy  = mean(pred == tru)),
    by = fold_id]
fwrite(performance_results, 'rf/performance_results.csv')


## ──────────────────────────────────────────────────────────────
## 6.  aggregate performance -----------------------------------
summary <- performance_results[, .(auc_roc, auc_pr,logloss,accuracy)][, .(
  mean_auc_roc = mean(auc_roc), sd_auc_roc = sd(auc_roc),
  mean_auc_pr  = mean(auc_pr),  sd_auc_pr  = sd(auc_pr),
  mean_logloss = mean(logloss), sd_logloss = sd(logloss),
  mean_acc     = mean(accuracy),sd_acc     = sd(accuracy)
)]
print(summary)
fwrite(summary, 'rf/final_model_summary.csv')


## ──────────────────────────────────────────────────────────────
## 7.   calculate F1
f1 <- f1_from_dt(all_results[, .(tru,pred), by= fold_id])

fwrite(f1, 'rf/f1_score.csv')


