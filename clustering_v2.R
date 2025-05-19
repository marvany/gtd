    # ========================== #
    #         CLUSTERING         #
    # ========================== #

source('/sc/arion/projects/va-biobank/PROJECTS/ma_extras/personal/gtd/gtd/classification_helper.R')
source('/sc/arion/projects/va-biobank/PROJECTS/ma_extras/personal/gtd/gtd/clustering_helper.R')
set.seed(123)

dt <- fread('data_inspection/final_predictors_considered.csv')
backup <- copy(dt)
dt <- copy(backup)

# Variables based on which we will be clustering
pred_cols <- c(
  "weapsubtype1", "targtype1_txt", "targsubtype1",
  "natlty1", "property"
)

"Because targtype1_txt and targsubtype1 refer to the same thing conceptualy, we 
dropped the latter for simplicity"

pred_cols <- c("weapsubtype1", "targtype1_txt", "natlty1", "property")

# ---------------------------- 1. prepare the matrix ----------------------------
# Keep only numeric columns; convert to matrix & z‑score each variable
dt <- dt[, ..pred_cols]
dt <- prepare_dt_clustering(dt)
fwrite(head(dt), 'Resources/preprocessing_clusters_example.csv')


dt <- dt[sample(.N, 10000)] # subset to reduce computational burden
num_cols <- names(which(sapply(dt, is.numeric)))
X        <- scale(as.matrix(dt[, ..num_cols]))
dim(X) # 132k , 5

# -------------- Determination of cluster number based on HC

## Separation metrics
int_metrics  <- evaluate_internal(X, ks = 2:8) # both show that n = 3 is the best way to go
"Silhouette climbs steadily, but the incremental gain flattens after k ≈ 6 (from 0.30 to 0.33).
Dunn—your strictest separation metric—peaks at k = 5 and 6 (0.033 4) then falls by a third at k = 7 and halves again at k = 8.
Calinski–Harabasz keeps increasing, as expected, but the rise between k = 6 → 8 is modest compared with earlier jumps."

fwrite(int_metrics, 'clustering/optN_sep.csv')

## Stability metrics
stab_metrics <- evaluate_stability(X, ks = 2:8, B = 100) 
fwrite(stab_metrics, 'clustering/optN_stab.csv')

## Summary of stability metrics
# avg_jaccard is the average for the clusters in the given iteration
summ_stab_metrics <- stab_metrics[
  , .(median_avg_jaccard = median(avg_jaccard)), 
  by = k
]
fwrite(summ_stab_metrics, 'clustering/optN_summary_stab.csv')
"We can see that for clustering solution 5 and 3 the greatest stability is demonstrated, when we view these results
in conjunction to the above we reach the conclusion the k = 5 is the optimal number for our cluster solution."


    # ======================================= #
    #     CONSTRUCT CLUSTERING SOLUTIONS      #
    # ======================================= #

# Set the number of our clustering solution
k = 5

# Euclidean distance + Ward’s linkage (works well for compact, spherical clusters)
dt <- dt[sample(.N, 10000)] # subset to reduce computational burden
num_cols <- names(which(sapply(dt, is.numeric)))
X        <- scale(as.matrix(dt[, ..num_cols]))

d   <- dist(X, method = "euclidean")
hc  <- hclust(d, method = "ward.D2")
cl <- cutree(hc, k)
plot_feature_distributions(
  dt, cl,
  output_path = "clustering/hc_results/feature_distribution.png"
)

"Παρατηρούμε πως στην 5η ομαδα ενω η εθνικότητα (natlty1) εμφανίζει μεγάλη ομοιογένεια μεταξύ των ομαδων 1-4, ειναι μοναδική 
για την 5η ομάδα. Στην 3η, 4η ομάδα, η επιλογή
του στόχου (targtype1_txt, targsubtype1) ήταν σχεδόν μοναδικές για τις ομάδες αυτές. Απ' την αλλη στην 2η ομαδα "

# 
sizes = table(cl)
sizes



# ---------------------------- 4. report results -------------------------------
# Show metrics sorted by average silhouette (feel free to sort by another index)
print(metrics[order(-avg_silhouette)])

# ---------------------------- 5. visualisation --------------------------------
best_k <- metrics[order(-avg_silhouette)]$k[1]

plot(hc, labels = FALSE,
     main = sprintf("Hierarchical clustering dendrogram (best k = %d)", best_k))
rect.hclust(hc, k = best_k, border = 2:6) # outline clusters on dendrogram

# ---------------------------- 6. return value ---------------------------------
# Optionally return the metrics table for downstream use
metrics
