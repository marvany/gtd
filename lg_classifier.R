    # ==================================== #
    #       BUILD LOGISTIC CLASSIFIER      #
    # ==================================== #
source('/sc/arion/projects/va-biobank/PROJECTS/ma_extras/personal/gtd/gtd/classification_helper.R')
set.seed(123)

dt <- fread('Resources/globalterrorismdb.csv')
backup <- copy(dt)

# these are the top 20 informative variables
top_20_predictors <- c(
  "targtype1_txtTerrorists/Non-State Militia",
  "targtype1_txtViolent Political Party",
  "region_txtNorth America",
  "dbsourceCETIS",
  "targtype1_txtMilitary",
  "region_txtEastern Europe",
  "dbsourceHewitt Project",
  "attacktype1_txtBombing/Explosion",
  "dbsourceUMD Schmid 2012",
  "targtype1_txtPrivate Citizens & Property",
  "dbsourceUMD Algeria 2010-2012",
  "dbsourceSTART Primary Collection",
  "targtype1_txtPolice",
  "weapsubtype1_txtUnknown Explosive Type",
  "region_txtMiddle East & North Africa",
  "region_txtSouth America",
  "weapsubtype1",
  "INT_MISC",
  "multiple"
)


# ---------------------------------------------
#   REBUILD THE TABLE USING SELECTED PREDICTORS
y       = "doubtterr"     # <- name of outcome column (0/1 or FALSE/TRUE)
vars_selected <- c('targtype1_txt' , 'region_txt' , 'dbsource' , 'weapsubtype1', 'weapsubtype1_txt' , 'INT_MISC' ,'multiple', 'attacktype1_txt')
ord_vars <- vars_selected # for convenience


# Inner join with rows_to_keep ; created in prepare_data.R
rows_to_keep = fread('data_inspection/ids_of_rows_considered.csv')
dt <- rows_to_keep[dt, on = 'eventid', nomatch = 0]
dt <- dt[, c(y, ord_vars), with = FALSE]

dim(dt) # 87k , 7


# --------------------- PREPARE FOR REGRESSION
# TURN CHARACTERS TO FACTORS
char_cols <- names(which(sapply(dt, is.character)))
dt[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]

insp.dt = data.table(
  variable = names(dt),
  n_unique = vapply(dt, uniqueN, integer(1L)),
  class = vapply(dt, function(x) paste(class(x), sep = ", "), character(1L))
)
insp.dt
binary_vars = unlist(insp.dt[n_unique == 2, 'variable'])
fwrite(insp.dt, 'vars_considered_summary.csv')

# Build xmat - does the encoding
xmat <- build_xmat(dt)

# Build yvec
yvec <- scale(dt[[y]])
yvec <- -yvec

dt <- cbind(
  data.table(
    doubtterr = yvec
  ),
  xmat
)
dt <- data.table(doubtterr = as.numeric(yvec))
dt[, (colnames(xmat)) := as.data.frame(xmat)]

tokeep <- c(y, safe_names(top_20_predictors))
setnames(dt, safe_names(names(dt)))
setdiff(tokeep, names(dt)) # SHOULD BE 0
dt <- dt[, ..tokeep]

dt[, (binary_vars) := lapply(.SD, function(x) ifelse(x > 0, 1, 0)), .SDcols = binary_vars]


head(dt)
# --------------------- BUILD MODEL
k       = 10           # <- # folds (use 5 or 10 by convention)
seed    = 1L           # <- reproducible shuffling
thresh  = 0.5        # <- decision threshold

form <- as.formula(paste(y, "~", paste(setdiff(names(dt), y), collapse = " + ")))

folds <- sample(rep_len(seq_len(k), nrow(dt)))


fwrite(dt, 'lg/model_table.csv')
all_lg_results <- rbindlist(
  pbmclapply(
    seq_len(k), 
    mc.cores = parallel::detectCores() - 2,
    function(i) {
      
      train <- dt[ folds != i]
      test <- dt[ folds == i ]
      
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
    accuracy  = mean(pred == tru)),
    by = fold_id]

performance_results
fwrite(performance_results, 'lg/performance_results.csv')

## ──────────────────────────────────────────────────────────────
## 5.  aggregate performance -----------------------------------
summary <- performance_results[, .(auc_roc, auc_pr,accuracy)][, .(
  mean_auc_roc = mean(auc_roc), sd_auc_roc = sd(auc_roc),
  mean_auc_pr  = mean(auc_pr),  sd_auc_pr  = sd(auc_pr),
  mean_acc     = mean(accuracy),sd_acc     = sd(accuracy)
)]
summary
fwrite(summary, 'lg/final_model_summary.csv')

## ──────────────────────────────────────────────────────────────
## 6.   calculate F1
f1 <- f1_from_dt(all_lg_results[, .(tru,pred), by= fold_id])
f1
fwrite(f1, 'lg/f1_score.csv')

