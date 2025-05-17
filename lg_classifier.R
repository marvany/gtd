
    # ==================================== #
    #       BUILD LOGISTIC CLASSIFIER      #
    # ==================================== #
library(pbmcapply)
library(pROC)
set.seed(123)
dt <- copy(backup)
dt <- remove_neg9_features(dt)
dt <- drop_sparse_strings(dt)
char_cols <- names(which(sapply(dt, is.character)))
dt[, (char_cols) := lapply(.SD, function(x) as.integer(factor(x, levels = unique(x)))), 
   .SDcols = char_cols]

dt <- na.omit(dt)
dt <- drop_bad_rows(dt)
dt <- dt[, (names(dt)) := lapply(.SD, as.numeric), .SDcols = names(dt)]

ord_vars <- c(
  "weapsubtype1",  "targtype1_txt", "targsubtype1", "natlty1",
  "property",      "multiple",      "dbsource",     "success",
  "country",       "guncertain1"
)

y       = "doubtterr",     # <- name of outcome column (0/1 or FALSE/TRUE)
k       = 10           # <- # folds (use 5 or 10 by convention)
seed    = 1L           # <- reproducible shuffling
thresh  = 0.5        # <- decision threshold

dt <- dt[, c(y, ord_vars), with = FALSE]
form <- as.formula(paste(y, "~", paste(setdiff(names(dt), y), collapse = " + ")))

set.seed(seed)
folds <- sample(rep_len(seq_len(k), nrow(dt)))

results <- rbindlist(lapply(seq_len(k), function(i) {
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
      fold      = i,
      auc       = auc(tru, pr)[1],
      accuracy  = mean(pred == tru),
      n_test    = .N                              # size of hold-out fold
    )

}))

## aggregate performance
summary <- results[,
.(mean_auc = mean(auc),
    sd_auc   = sd(auc),
    mean_acc = mean(accuracy),
    sd_acc   = sd(accuracy))
]
summary
