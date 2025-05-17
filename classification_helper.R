
library(car)   # install.packages("car")
library(data.table)
library(caret)
library(data.table)
library(forcats)
library(ranger)
library(vip)
library(mltools)
library(Matrix)
library(ggplot2)
library(text2vec)      # fast tokenization + vectorizers
library(reticulate) 
library(glmnet)
library(pbmcapply)
library(pROC)
library(PRROC)

writeReport <- function(txt, output = "report.txt"){
    ## character vector you want to add
    ## open a connection in append-text mode
    con <- file(output, open = "at")   # "a" = append, "t" = text
    writeLines(txt, con, sep = "\n")         # write + newline
    close(con)                               # always close the connection
}


# REMOVES FEATURES WHOSE NUMBER OF -9 VALUES EXCEEDS A CERTAIN THRESHOLD
remove_neg9_features <- function(dt){
    val        <- -9
    thr        <- 0.30                       # 30 %
    n_rows     <- nrow(dt)
    ## 1. find columns where at least 30 % of entries equal -9
    num_cols <- names(which(sapply(dt, is.numeric)))   # restrict to numeric vars

    # proportion of rows that equal -9 for each numeric column
    prop_neg9 <- dt[ ,
    sapply(.SD, function(x) mean(x == val, na.rm = TRUE)),
    .SDcols = num_cols]

    # offending columns
    cols_bad <- names(prop_neg9[prop_neg9 > thr])

    ## drop them in-place
    dt[ , (cols_bad) := NULL ]
}

# DROPS COLUMNS
drop_sparse_strings <- function(dt, thresh = 0.30) {
  stopifnot(is.data.table(dt))

  n <- nrow(dt)

  # find offending columns
  sparse_cols <- dt[, which(
    vapply(.SD, function(col)
      mean((col == "" | is.na(col))) > thresh,
      logical(1L)
    )
  )]

  # drop them (nothing happens if none are found)
  if (length(sparse_cols))
    dt[, (sparse_cols) := NULL]

  invisible(dt)            # return the modified data.table invisibly
}

# DROPS ROWS
drop_bad_rows <- function(DT, bad_str = "", bad_int = -9L) {
  stopifnot(is.data.table(DT))

  # TRUE for rows that contain a bad value in *any* column
  bad_row <- DT[, Reduce(`|`, lapply(.SD, function(col) {
    if      (is.character(col))        col == bad_str
    else if (is.integer(col))          col == bad_int
    else if (is.numeric(col))          col == bad_int
    else                               FALSE           # other types untouched
  }))]

  DT[!bad_row]        # return the filtered table (original DT unchanged)
}


# restores the outcome to be non-scaled
restore_doubtterr <- function(mat, tgt = "doubtterr") {
  stopifnot(is.matrix(mat), tgt %in% colnames(mat))
  j <- match(tgt, colnames(mat))
  mat[, j] <- as.integer(mat[, j] > 0)   # −1 ↦ 0, +1 ↦ 1   (robust to small FP noise)
  invisible(mat)                         # return the same matrix invisibly
}



# BUILDS RANDOM FOREST
eval_fold <- function(i) {
  train <- dt[fold != i]
  test  <- dt[fold == i]

  ## ---- inner 3-fold tuning via OOB/CV ------------------------
  k_inner <- 3
  set.seed(i)
  train[, inner := sample(rep(seq_len(k_inner), length.out = .N)),
        by = doubtterr]

  grid$score <- pbmclapply(
    seq_len(nrow(grid)),
    function(g) {
      hp <- grid[g]
      auc_inner <- numeric(k_inner)

      for (j in seq_len(k_inner)) {
        tr_in  <- train[inner != j]
        val_in <- train[inner == j]

        rf <- ranger(
          formula   = doubtterr ~ .,
          data      = tr_in[, c(y, pred_cols), with = FALSE],
          probability            = TRUE,
          num.trees              = hp$num.trees,
          mtry                   = hp$mtry,
          min.node.size          = hp$min.node.size,
          sample.fraction        = hp$sample.fraction,
          respect.unordered.factors = "order",
          seed = i * 100 + j
        )

        pr <- predict(rf,
          data = val_in[, c(y, pred_cols), with = FALSE])$predictions[, 2]
        auc_inner[j] <- as.numeric(pROC::auc(val_in[[y]], pr))
      }
      mean(auc_inner)
    },
    mc.cores = parallel::detectCores() - 1
  )

  best_hp <- grid[which.max(unlist(grid$score))]

  ## ---- final RF on full training fold ------------------------
  rf_final <- ranger(
    doubtterr ~ .,
    data      = train[, c(y, pred_cols), with = FALSE],
    probability            = TRUE,
    num.trees              = best_hp$num.trees,
    mtry                   = best_hp$mtry,
    min.node.size          = best_hp$min.node.size,
    sample.fraction        = best_hp$sample.fraction,
    respect.unordered.factors = "order",
    seed = 10 + i
  )

  ## ---- predict & metrics -------------------------------------
  pr  <- predict(rf_final,
          data = test[, c(y, pred_cols), with = FALSE])$predictions[, 2]
  tru <- test[[y]]
  pred <- as.integer(pr >= 0.5)

  data.table(
    fold      = i,
    auc_roc   = as.numeric(pROC::auc(tru, pr)),
    auc_pr    = PRROC::pr.curve(scores.class0 = pr[tru == 1],
                                scores.class1 = pr[tru == 0])$auc.integral,
    logloss   = -mean(tru * log(pr) + (1 - tru) * log(1 - pr)),
    accuracy  = mean(pred == tru)
  )
}

# PREPARE DATA FOR MODELING
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