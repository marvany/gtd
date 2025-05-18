# ================================ #
#                        #
# ================================ #

source('/sc/arion/projects/va-biobank/PROJECTS/ma_extras/personal/gtd/gtd/classification_helper.R')
set.seed(123)

dt <- fread('Resources/globalterrorismdb.csv')
backup <- copy(dt)

#remove_neg9_features(dt)
#drop_sparse_strings(dt)
#drop_bad_rows(dt)
dt <- copy(backup)
X <- prep_for_clustering(dt)
#X <- prep_numeric_matrix(dt)
## dt = your original data.table -------------------------------------------

good_idx <- rowSums(!is.finite(X)) == 0   #  find rows with infinite counts
X_clean  <- X[good_idx, , drop = FALSE] # remove inf counts
nf_counts <- colSums(!is.finite(X_clean)) # test that there are no inf counts
X <- X_clean



Ks <- 2:6
cores_outer <- 1                       # avoid nested forks
perm_list <- lapply(
  Ks,
  function(k) permute_sparsek(X, k, nperms = 30, cores = detectCores() - 2)
)

set.seed(1)
Ks        <- 2:6
perm_list <- lapply(
  Ks,
  #mc.cores = parallel::detectCores() - 2, 
  function(k) {
  set.seed(1)
  KMeansSparseCluster.permute(X, K = k, nperms = 30)
  browser()
})

## pick the K with the largest gap statistic
gap_vals <- vapply(perm_list, function(p) max(p$gaps), numeric(1))
best_idx <- which.max(gap_vals)

best_K   <- Ks[best_idx]
best_w   <- perm_list[[best_idx]]$bestw

skm <- KMeansSparseCluster(X, K = best_K, wbounds = best_w)

clusters      <- skm[[1]]$Cs
selected_vars <- colnames(X)[skm[[1]]$ws != 0]



#############################################################################

perm     <- KMeansSparseCluster.permute(X, K = 2:6, nperms = 30) # tune λ
best_K   <- perm$bestK        # #clusters chosen by the gap statistic
best_w   <- perm$bestw        # λ (wbounds) that maximised the gap

skm      <- KMeansSparseCluster(X, K = best_K, wbounds = best_w)

clusters      <- skm[[1]]$Cs               # cluster labels for each row
selected_vars <- colnames(X)[skm[[1]]$ws != 0]   # non-zero feature weights