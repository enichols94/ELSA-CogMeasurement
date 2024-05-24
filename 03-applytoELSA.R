##########################################################################
### Author: Emma Nichols
### Project: Development and assessment of analytic methods to improve the 
###          measurement of cognition in longitudinal studies of aging through 
###          the use of sub-studies with comprehensive neuropsychological testing
### Purpose: Apply regression and CFA measures to the core ELSA data
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate,
               simputation, VIM, MplusAutomation, glmnet)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# SET OBJECTS -------------------------------------------------------------

dir <- paste0(dropbox_dir, "projects/ips_cog_prediction/")
rawdata_dir <- paste0(dir, "data/source/")
derived_dir <- paste0(dir, "data/derived/")
model_dir <- "C:/Users/emmanich/mplus_models/elsa_cogimp/"
plot_dir <- paste0(dir, "plots/")

model_date <- "2023_11_14"
cfa_dir <- paste0(model_dir, model_date, "/fulldata_BAYES/")

regression_date <- "2023_11_14"

# GET DATA ----------------------------------------------------------------

elsa_dt <- read_rds(paste0(derived_dir, "elsacore_long.rds"))

hcap_dt <- read_rds(paste0(derived_dir, "elsaprocessed.rds"))

varmap_dt <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))

# FORMAT AND APPEND HCAP DATA ---------------------------------------------

## format HCAP data and append this to the core ELSA data

cog_vars <- varmap_dt[!is.na(cfa_label), label]
cog_cfalabels <- varmap_dt[!is.na(cfa_label), cfa_label]

format_cfadata <- function(data){
  dt <- copy(data)
  cv_map <- cog_vars %in% names(dt)
  setnames(dt, cog_vars[cv_map], cog_cfalabels[cv_map])
  dt <- dt[, c("id", cog_cfalabels[cv_map]), with = F]
  dt[, wave := 0]
  return(dt)
}

hcap_dt <- format_cfadata(hcap_dt)

# GET CFA MODEL PARAMETERS ------------------------------------------------

## read in CFA model parameters from saved model

cfa <- readModels(cfa_dir)

params <- as.data.table(cfa$parameters$unstandardized)

get_statement <- function(row_num){
  row <- params[row_num]
  if (grepl("BY", row[, paramHeader]) | grepl("WITH", row[, paramHeader])){
    statement <- paste0(tolower(str_extract(row[, paramHeader], "^[A-Z|0-9]*")), " ", 
                        str_extract(row[, paramHeader], "[A-Z]*$"), " ",
                        tolower(row[, param]), "@", row[, est], ";\n")
  } else if (row[, paramHeader] %in% c("Means", "Intercepts", "Thresholds")){
    statement <- paste0("[", tolower(row[, param]), "@", row[, est], "];\n")
  } else if (row[, paramHeader] %in% c("Variances", "Residual.Variances")){
    statement <- paste0(tolower(row[, param]), "@", row[, est], ";\n")
  }
  return(statement)
}

model_statement <- paste(unlist(lapply(1:nrow(params), get_statement)), collapse = " ")

# RUN MODEL ON STACKED DATA ------------------------------------------

## CFA model with all parameters fixed to estimate scores in the full ELSA sample

runmodel_fixed <- function(data, estimator = "MLR", statement = model_statement, sens = ""){
  dt <- copy(data)
  
  cfa_vars <- varmap_dt[!is.na(cfa_label) & cfa_label %in% names(dt), cfa_label]
  cat_vars <- varmap_dt[cfa_label %in% cfa_vars & is.na(continuous), cfa_label]
  dt <- dt[, c("id", "wave", cfa_vars), with = F]
  dt <- rbind(dt, hcap_dt) ## append on hcap data to have at least one observation of all variables
  dt[, new_id := as.numeric(paste0(id, wave))]
  
  ## get categorical items and use in variable command
  categorical_items <- c()
  for(i in 1:length(cat_vars)){
    if (i%%4==0) {
      unit <- paste0("\n", cat_vars[i])
    } else {
      unit <- cat_vars[i]
    }
    categorical_items <- c(categorical_items, unit)
  }
  variable_command <- paste0("\nCATEGORICAL = ", paste0(categorical_items, collapse = " "), ";\nIDVARIABLE = new_id;")
  
  ## create directories
  if (sens == ""){
    dir.create(paste0(model_dir, model_date, "/FIXED_", estimator), recursive = T)
    mdir <- paste0(model_dir, model_date, "/FIXED_", estimator, "/")
  } else {
    dir.create(paste0(model_dir, model_date, "/FIXED_", estimator, "_", sens), recursive = T)
    mdir <- paste0(model_dir, model_date, "/FIXED_", estimator, "_", sens, "/")
  }
  
  ## toggle commands by estimator
  if (estimator == "BAYES"){
    analysis_command <- paste0("ESTIMATOR = ", estimator, "; COVERAGE = 0; PROCESSORS = 4;")
    save_command <- "FILE = output.sav; SAVE = FSCORES(200);"
  } else if (estimator == "MLR"){
    analysis_command <- paste0("ESTIMATOR = ", estimator, "; LINK = PROBIT; COVERAGE = 0; PROCESSORS = 4;")
    save_command <- "FILE = output.sav; SAVE = FSCORES;"
  }
  
  ## run model
  object <- mplusObject(
    TITLE = paste0("ELSA fixed CFA - longitudinal"),
    VARIABLE = variable_command,
    ANALYSIS = analysis_command, 
    MODEL = model_statement,
    OUTPUT = "STDYX; MODINDICES; TECH1; TECH8; TECH10;",
    SAVEDATA = save_command,
    usevariables = c(cfa_vars, "new_id"),
    rdata = as.data.frame(dt)
  )
  
  model <- mplusModeler(
    object = object,
    modelout = paste0(mdir, "model.inp"),
    hashfilename = F,
    run = 1L
  )
  
}

if(file.exists(paste0(model_dir, model_date, "/FIXED_BAYES/model.out"))){ ## read in old model if it exists, otherwise re-run 
  message("Reading model")
  cfa_model <- readModels(paste0(model_dir, model_date, "/FIXED_BAYES"))
  fscore_dt <- as.data.table(cfa_model$savedata)
} else {
  message("Running model")
  cfa_model <- runmodel_fixed(data = elsa_dt, estimator = "BAYES", statement = model_statement)
  fscore_dt <- as.data.table(cfa_model$results$savedata)
}

fscore_dt[, `:=` (id = as.numeric(gsub("[0-9]$", "", NEW_ID)), wave = as.numeric(str_extract(NEW_ID, "[0-9]$")))]
fscore_dt <- fscore_dt[!wave == 0, .(id, wave, cfa_score = G.Mean)]

# RUN SENSITIVITY WITHOUT WAVE 7-9 MEASURES -------------------------------

sens1_dt <- copy(elsa_dt)

## variable names
cog_vars <- varmap_dt[!is.na(cfa_label), label]
cog_cfalabels <- varmap_dt[!is.na(cfa_label), cfa_label]

## get rid of extra variables not in pred19 - set to missing
preds79 <- c("scissors", "cactus", "serial7", "primeminister", "monarch", "uspresident", "backcount20")
sens1_dt[, c(preds79, paste0(preds79, "I")) := NULL]
sens1_dt[, (preds79) := NA]

## get rid of imputed versions of 1-9 variables
preds19 <- c("recall10imm", "recall10del", "month", "year", "dayweek", "date", "animals")
sens1_dt[, paste0(preds19, "I") := NULL]

## get rid of old cfa vars
sens1_dt[, (cog_cfalabels) := NULL]

## change all variable names
for (var in cog_vars){
  sens1_dt[, c(cog_cfalabels[cog_vars == var]) := get(var)]
}

if(file.exists(paste0(model_dir, model_date, "/FIXED_BAYES_s1/model.out"))){ ## read in old model if it exists, otherwise re-run
  message("Reading model")
  cfa_model_s1 <- readModels(paste0(model_dir, model_date, "/FIXED_BAYES_s1"))
  fscore_dt_s1 <- as.data.table(cfa_model_s1$savedata)
} else {
  message("Running model")
  cfa_model_s1 <- runmodel_fixed(data = sens1_dt, estimator = "BAYES", statement = model_statement, sens = "s1")
  fscore_dt_s1 <- as.data.table(cfa_model_s1$results$savedata)
}

fscore_dt_s1[, `:=` (id = as.numeric(gsub("[0-9]$", "", NEW_ID)), wave = as.numeric(str_extract(NEW_ID, "[0-9]$")))]
fscore_dt_s1 <- fscore_dt_s1[!wave == 0, .(id, wave, cfa_score_s1 = G.Mean)]

# GET REGRESSION MODELS ---------------------------------------------------

regression_models <- read_rds(paste0(derived_dir, "regressionmodels_", regression_date, ".rds"))

# GET REGRESSION PREDICTIONS ----------------------------------------------

## WITH IMPUTATIONS

regression_dt <- copy(elsa_dt)
impute_vars <- names(regression_dt)[grepl("I$", names(regression_dt))]
regression_dt <- regression_dt[, c("id", "wave", impute_vars), with = F]
setnames(regression_dt, impute_vars, gsub("I$", "", impute_vars))

## get data subsets
regression_dt19 <- copy(regression_dt)
regression_dt79 <- copy(regression_dt[wave > 6])
regression_dt79[, serial7 := as.factor(serial7)]

## get predictions
regression_dt19[, pred19 := predict(regression_models$model18, newdata = regression_dt19)]
regression_dt79[, pred79 := predict(regression_models$model78, newdata = regression_dt79)]

## WITHOUT IMPUTATIONS

regression_noimp_dt <- copy(elsa_dt)
coefvars <- attr(regression_models$model78$terms, "term.labels")
regression_noimp_dt <- regression_noimp_dt[, c("id", "wave", coefvars), with = F]

## get data subsets
regression_noimp_dt19 <- copy(regression_noimp_dt)
regression_noimp_dt79 <- copy(regression_noimp_dt[wave > 6])
regression_noimp_dt79[, serial7 := as.factor(serial7)]

## get predictions
regression_noimp_dt19[, pred19_noimp := predict(regression_models$model18, newdata = regression_noimp_dt19)]
regression_noimp_dt79[, pred79_noimp := predict(regression_models$model78, newdata = regression_noimp_dt79)]

## combine predictions
regressionpred_dt <- merge(regression_dt19[, .(id, wave, pred19)], regression_dt79[, .(id, wave, pred79)], by = c("id", "wave"), all = T)
regressionpred_noimp_dt <- merge(regression_noimp_dt19[, .(id, wave, pred19_noimp)], regression_noimp_dt79[, .(id, wave, pred79_noimp)], by = c("id", "wave"), all = T)

# COMBINE ALL DATA --------------------------------------------------------

elsa_dt <- Reduce(function(x,y) merge(x,y, by=c("id","wave"), all.x=T),
                  list(elsa_dt, fscore_dt, fscore_dt_s1,
                       regressionpred_dt, regressionpred_noimp_dt))

# SAVE DATA ---------------------------------------------------------------

write_rds(elsa_dt, paste0(derived_dir, "coredata_withpreds_", date, ".rds"))
