
# PREPARE DATA FOR CLUSTER
prepare_dt_clustering <- function(dt){
  dt <- remove_neg9_features(dt)
  dt <- drop_sparse_strings(dt)
  char_cols <- names(which(sapply(dt, is.character)))
  #dt[, (char_cols) := lapply(.SD, function(x) as.integer(factor(x, levels = unique(x)))), 
  #  .SDcols = char_cols]

  dt <- na.omit(dt)
  dt <- drop_bad_rows(dt)
  #dt <- dt[, (names(dt)) := lapply(.SD, as.numeric), .SDcols = names(dt)]
  dt
}

# Hierarchical clustering & stability utilities -------------------------------------------------
# 1. prepare_matrix(dt)       → scaled numeric matrix X
# 2. evaluate_internal(X, ks) → internal indices (silhouette, Dunn, CH, entropy)
# 3. evaluate_stability(X, ks, B) → bootstrap Jaccard stability per k
# ------------------------------------------------------------------------------

# ------------------------- package loader --------------------------------------
required_pkgs <- c("data.table", "cluster", "fpc", "pbmcapply")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)
lapply(required_pkgs, library, character.only = TRUE)

# ------------------------- 1. matrix prep --------------------------------------
prepare_matrix <- function(dt) {
  num_cols <- names(which(sapply(dt, is.numeric)))
  if (length(num_cols) == 0) stop("No numeric columns to cluster.")
  scale(as.matrix(dt[, ..num_cols]))
}

# ------------------------- 2. internal indices ---------------------------------
evaluate_internal <- function(X, ks = 2:5, method = "ward.D2",
                              ncores = max(1L, parallel::detectCores() - 1L)) {
  d  <- dist(X, method = "euclidean")
  hc <- hclust(d, method = method)

  data.table::rbindlist(
    pbmclapply(
      ks,
      function(k) {
        cl  <- cutree(hc, k)
        sil <- cluster::silhouette(cl, d)
        cs  <- fpc::cluster.stats(d, cl)
        data.table::data.table(
          k                 = k,
          avg_silhouette    = mean(sil[, 3]),
          dunn              = cs$dunn,
          calinski_harabaz  = cs$ch
        )
      },
      mc.cores = ncores
    )
  )
}

# ------------------------- 3. progressive bootstrap stability -------------------
# For each bootstrap b, we draw a subset of size n * frac[b] *without replacement*.
# frac[1] = 1, frac[B] = min_frac — producing a progressively smaller sample.
# Each subset is evaluated via fpc::clusterboot (bootmethod = "subset").

evaluate_stability <- function(X, ks = 2:5, B = 10, min_frac = 0.01,
                               method = "ward.D2",
                               ncores = max(1L, parallel::detectCores() - 1L)) {
  if (min_frac <= 0 || min_frac > 1) stop("min_frac must be in (0,1].")
  d      <- dist(X, method = "euclidean")
  fracs  <- seq(1, min_frac, length.out = B)  # descending sample sizes

  data.table::rbindlist(
    lapply(
      ks,
      function(k) {
        jacc <- unlist(pbmclapply(fracs, function(f) {
          cb <- fpc::clusterboot(
            d,
            B             = 1,              # one resample per fraction
            bootmethod    = "subset",      # without‑replacement subset
            subtuning = round(nrow(X)*f),
            bootpar       = f,              # fraction of rows to keep
            distances     = TRUE,
            clustermethod = disthclustCBI,
            method        = method,
            k             = k,
            seed          = 123,
            count = FALSE
          )
          mean(cb$subsetmean)   # average Jaccard for this subset size
        }, mc.cores = ncores))
        
        data.table::data.table(
          k              = k,
          frac           = fracs,
          avg_jaccard    = jacc
        )
      } #,mc.cores = ncores
    )
  )
}

if(F){
  # FOR DEBUGGING
  evaluate_stability <- function(X, ks = 2:5, B = 10, min_frac = 0.5,
                                method = "ward.D2",
                                ncores = max(1L, parallel::detectCores() - 1L)) {
    if (min_frac <= 0 || min_frac > 1) stop("min_frac must be in (0,1].")
    d      <- dist(X, method = "euclidean")
    fracs  <- seq(1, min_frac, length.out = B)  # descending sample sizes

    data.table::rbindlist(
      lapply(
        ks,
        function(k) {
          jacc <- unlist(lapply(fracs, function(f) {
            cb <- fpc::clusterboot(
              d,
              B             = 2,              # one resample per fraction
              bootmethod    = "subset",      # without‑replacement subset
              bootpar       = f,              # fraction of rows to keep
              distances     = TRUE,
              clustermethod = disthclustCBI,
              method        = method,
              k             = k,
              seed          = 123
            )
            mean(cb$subsetmean)   # average Jaccard for this subset size
          }))#, mc.cores = ncores))
          data.table::data.table(
            k              = k,
            frac           = fracs,
            avg_jaccard    = jacc
          )
        } #,mc.cores = ncores
      )
    )
  }
}
if(F){
# ------------------------- 3. bootstrap stability ------------------------------
# Uses fpc::clusterboot + disthclustCBI (Ward linkage) to compute Jaccard scores.
evaluate_stability <- function(X, ks = 2:5, B = 50, method = "ward.D2",
                               ncores = max(1L, parallel::detectCores() - 1L)) {
  d <- dist(X, method = "euclidean")

  data.table::rbindlist(
    lapply(
      ks,
      function(k) {
        cb <- fpc::clusterboot(
          d,
          B             = B,
          bootmethod    = "boot",    # simple bootstrap
          distances     = TRUE,
          clustermethod = disthclustCBI,
          method        = method,
          k             = k,
          seed          = 123
        )
        data.table::data.table(
          k                    = k,
          avg_jaccard          = mean(cb$bootmean),
          clusterwise_jaccard  = list(cb$bootmean)  # store per‑cluster
        )
      }#,  mc.cores = ncores
    )
  )[order(-avg_jaccard)]
}
}
# ------------------------- example workflow ------------------------------------
# dt <- fread("data.csv")
# X  <- prepare_matrix(dt)
# int_metrics  <- evaluate_internal(X, ks = 2:5)
# stab_metrics <- evaluate_stability(X, ks = 2:5, B = 100)
# print(int_metrics)
# print(stab_metrics)



# ------------------------- 4. quick feature‑distribution plots ------------------
# Visualise selected numeric variables across clusters using ggplot2.

plot_feature_distributions <- function(dt, cl, features = NULL, max_feats = 20,
                                    output_path = NULL,
                                    width = 12, height = 8, dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  library(ggplot2)
  library(data.table)
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  library(ggplot2)
  library(data.table)

  dt <- copy(dt)
  dt[, cluster := factor(cl)]

  num_cols <- names(which(sapply(dt, is.numeric)))[!names(dt) %in% "cluster"]
  if (is.null(features)) features <- num_cols
  features <- head(features, max_feats)   # keep plot manageable

  long <- melt(dt[, c("cluster", features), with = FALSE],
               id.vars = "cluster",
               variable.name = "feature",
               value.name = "value")

  p <- ggplot(long, aes(x = cluster, y = value, fill = cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~ feature, scales = "free_y", ncol = 4) +
    guides(fill = FALSE) +
    labs(x = "Cluster", y = "Value (original scale)") +
    theme_bw(base_size = 12)

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = output_path, plot = p, width = width, height = height, dpi = dpi)
  }
  p
}





# FULLY DEPRECATED
if(FALSE){
    ## ---------- 1. specify outcome and column groups --------------------
  ycol       <- "doubtterr"

  char_cols  <- names(which(sapply(dt, is.factor)))
  num_cols   <- setdiff(names(which(sapply(dt, is.numeric))), ycol)   # drop outcome

  ## ---------- 2. scale numeric predictors ----------------------------
  X_num   <- scale(dt[ , ..num_cols])          # matrix of z-scored numerics

  ## ---------- 3. one-hot–encode character predictors -----------------
  if (length(char_cols)) {
    dt[ , (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]   # chars → factors
    X_char <- model.matrix(~ . - 1, data = dt[ , ..char_cols])          # -1 → no intercept
    # If thousands of dummy columns, prefer:  X_char <- sparse.model.matrix(~ . - 1, data = dt[ , ..char_cols])
    xmat   <- cbind(X_num, X_char)             # final design matrix
  } else {
    xmat   <- X_num                            # char_cols absent
  }
}

