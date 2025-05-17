    # ========================================= #
    #       BUILD A RANDOM FOREST CLASSIFER     #
    # ========================================= #
library(data.table)
library(ranger)         # fast RF
library(pbmcapply)      # parallel lapply
library(pROC)           # AUC-ROC
library(PRROC)          # AUC-PR

## ──────────────────────────────────────────────────────────────
## 0.  data prep  (your existing helpers) -----------------------
dt <- copy(backup)          # full GTD slice
dt <- prepare_dt(dt)
dt
y         <- "doubtterr"                 # 0 / 1
pred_cols <- c(
  "weapsubtype1", "targtype1_txt", "targsubtype1",
  "natlty1", "property"
)

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
results <- rbindlist(
  pbmclapply(seq_len(k_outer), eval_fold,
             mc.cores = parallel::detectCores() - 1)
)
fwrite(results, 'final_results.csv')

## ──────────────────────────────────────────────────────────────
## 5.  aggregate performance -----------------------------------
summary <- results[, .(
  mean_auc_roc = mean(auc_roc), sd_auc_roc = sd(auc_roc),
  mean_auc_pr  = mean(auc_pr),  sd_auc_pr  = sd(auc_pr),
  mean_logloss = mean(logloss), sd_logloss = sd(logloss),
  mean_acc     = mean(accuracy),sd_acc     = sd(accuracy)
)]
print(summary)
fwrite(summary, 'final_model_summary.csv')


   mean_auc_roc sd_auc_roc mean_auc_pr   sd_auc_pr mean_logloss sd_logloss
          <num>      <num>       <num>       <num>        <num>      <num>
1:    0.9533981 0.00338798    0.869247 0.006087762          NaN         NA
    mean_acc      sd_acc
       <num>       <num>
1: 0.9404722 0.002079409




