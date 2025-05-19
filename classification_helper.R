
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
library(pbmcapply)
library(sparcl)
library(cluster)
library(fpc)


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
    fold_id = i,           # identifies which outer fold produced the row
    row_id  = test[, .I],  # row index **within the original dt** (drop this if not needed)
    tru     = tru,
    pred    = pred
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

# PREPARE DATA FOR CLUSTER
prepare_dt_clustering <- function(dt){
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

# CALCULATE F1
f1_from_dt <- function(dt) {
  # fast 1-pass counts
  cm <- dt[, .(
    TP = sum(tru == 1 & pred == 1),
    FP = sum(tru == 0 & pred == 1),
    FN = sum(tru == 1 & pred == 0)
  )]

  with(cm, {
    precision <- TP / (TP + FP)
    recall    <- TP / (TP + FN)
    f1        <- if (precision + recall == 0) NA_real_
                 else 2 * precision * recall / (precision + recall)

    data.table(precision, recall, f1)
  })
}


# REMOVE ZERO VARIANCE COLUMNS
remove_zero_var <- function(dt){
  nzv_idx   <- nearZeroVar(dt)          # dt is your data.table / data.frame
  nzv_names <- names(dt)[nzv_idx]       # column names with near-zero variance
  dt[,(nzv_names):=NULL]
}


# CORRELATION MATRIX FILTER
correlation_filter <- function(dt){

  #dt[, (eventid):=NULL]
  num_cols <- names(which(sapply(dt, is.numeric)))        # keep only numerics
  #num_cols <- setdiff(num_cols, 'eventid') # leave eventid out cause it's highly correlated with year
  corr_mat <- cor(dt[, ..num_cols], use = "pairwise.complete.obs")

  # Drop one variable from every pair with |r| ≥ 0.90
  high_corr_idx   <- findCorrelation(corr_mat, cutoff = 0.90, verbose = TRUE)
  high_corr_names <- num_cols[high_corr_idx]

  dt[, !high_corr_names, with = FALSE]              # in-place removal
}


prep_for_clustering <- function(dt,
                          nzv_cut   = .95,   # drop if >95 % identical values
                          corr_cut  = .90) { # drop one of any |r| > .9 pair
  dt <- remove_neg9_features(dt)
  dt <- drop_sparse_strings(dt)
  
  # keep only numeric for distance-based clustering
  num_cols <- names(which(sapply(dt, is.numeric)))
  dt_num   <- dt[, ..num_cols]
  
  ## 1A. near-zero variance
  idx_nzv <- nearZeroVar(dt_num, freqCut = 95/5, uniqueCut = nzv_cut * 100)
  if (length(idx_nzv)) dt_num <- dt_num[, -idx_nzv, with = FALSE]
  
  ## 1B. high pairwise correlation
  R       <- cor(dt_num, use = "pairwise.complete.obs")
  idx_cor <- findCorrelation(R, cutoff = corr_cut)
  if (length(idx_cor)) dt_num <- dt_num[, -idx_cor, with = FALSE]

  dt <- na.omit(dt)
  dt <- drop_bad_rows(dt)
  
  ## 1C. scale for clustering
  as.matrix(scale(dt_num))
}
library(data.table)


prep_numeric_matrix <- function(dt) {
  # 1. keep only numeric columns
  num_cols <- names(which(sapply(dt, is.numeric)))
  Xdt      <- copy(dt)[, ..num_cols]

  # 2. drop columns that are *all* NA  ------------------------------
  all_na   <- names(which(colSums(!is.finite(as.matrix(Xdt))) == nrow(Xdt)))
  if (length(all_na)) Xdt[, (all_na) := NULL]

  # 3. drop zero-variance columns *before* scaling -----------------
  zero_sd  <- names(which(sapply(Xdt, function(x)
                       var(x, na.rm = TRUE) == 0 | all(is.na(x)))))
  if (length(zero_sd)) Xdt[, (zero_sd) := NULL]

  # 4. optional: impute remaining NAs (rf-based missRanger is fast)
  #   library(missRanger)
  #   Xdt <- missRanger(Xdt, pmm.k = 3, seed = 1)

  # 5. scale & final sanity check ----------------------------------
  X <- scale(as.matrix(Xdt))


  X
}






permute_sparsek <- function(X, K,
                            nperms   = 30,
                            wbounds  = NULL,   # let sparcl choose a grid
                            nstart   = 20,
                            seed     = 1,
                            cores    = detectCores() - 2) {

  ## (a) parallel permutation fits
  perm_gap <- pbmclapply(
    X = seq_len(nperms),
    mc.cores = cores,
    FUN = function(i) {
      set.seed(seed + i)
      Xperm <- apply(X, 2, sample)                    # column-wise shuffle
      fit   <- KMeansSparseCluster(Xperm, K = K,
                                   wbounds = wbounds,
                                   nstart  = nstart,
                                   silent  = TRUE)
      max(fit[[1]]$gaps)                              # return gap for this λ
    }
  )

  ## (b) best λ = wbounds value that maximises mean gap over permutations
  gap_vec <- unlist(perm_gap, use.names = FALSE)
  best_w  <- wbounds[ which.max(gap_vec) ]

  list(bestw = best_w, gaps = gap_vec)
}
