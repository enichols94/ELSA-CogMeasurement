##########################################################################
### Author: Emma Nichols
### Project: Development and assessment of analytic methods to improve the 
###          measurement of cognition in longitudinal studies of aging through 
###          the use of sub-studies with comprehensive neuropsychological testing
### Purpose: Estimate CFA models to create gold standard cognitive factor scores 
###          and the full CFA for prediction to the ELSA sample
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate,
               simputation, VIM, MplusAutomation)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# SET OBJECTS -------------------------------------------------------------

dir <- paste0(dropbox_dir, "projects/ips_cog_prediction/")
rawdata_dir <- paste0(dir, "data/source/")
derived_dir <- paste0(dir, "data/derived/")
model_dir <- "C:/Users/emmanich/mplus_models/elsa_cogimp/"
plot_dir <- paste0(dir, "plots/")

## update model_dir with date
model_dir <- paste0(model_dir, date, "/")

# GET DATA ----------------------------------------------------------------

elsa_dt <- read_rds(paste0(derived_dir, "elsaprocessed.rds"))

varmap_dt <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))

# CREATE MISSINGNESS ------------------------------------------------------

nm1 <- c("recall10imm", "recall10del", "month", "year", "dayweek", "date", "animals")
nm2 <- c(nm1, "scissors", "cactus", "serial7", "primeminister", "monarch", "uspresident", "backcount20"); tm2 <- c("monarch", "uspresident", "backcount20")

## nm1 = items in core survey all waves
## nm2 = items in core survey waves 7-8
## tm2 = items in waves 7-8 but not HCAP

cog_vars <- varmap_dt[!is.na(cfa_label), label]

elsacfa_dt <- copy(elsa_dt) ## version with no missingness to identify bifactors and misfit

elsacfa_gs <- copy(elsa_dt)
elsacfa_gs[, c(tm2) := NULL] ## only HCAP items for gold standard comparison

elsacfa_dt1 <- copy(elsa_dt)
elsacfa_dt1[train_sample == 0, c(setdiff(cog_vars, nm1)) := NA]
elsacfa_dt1[, c(tm2) := NULL] ## exclude because introduced in waves 7-8

elsacfa_dt2 <- copy(elsa_dt)
elsacfa_dt2[train_sample == 0, c(setdiff(cog_vars, nm2)) := NA][train_sample == 1, c(tm2) := NA]

# FORMAT DATA FOR CFA -----------------------------------------------------

## rename all variables to mplus friendly names, subset datasets to only include these 
## variables

cog_vars <- varmap_dt[!is.na(cfa_label), label]
cog_cfalabels <- varmap_dt[!is.na(cfa_label), cfa_label]

format_cfadata <- function(data){
  dt <- copy(data)
  cv_map <- cog_vars %in% names(dt)
  setnames(dt, cog_vars[cv_map], cog_cfalabels[cv_map])
  dt <- dt[, c("id", cog_cfalabels[cv_map]), with = F]
  return(dt)
}

elsacfa_dt <- format_cfadata(elsacfa_dt)
elsacfa_gs <- format_cfadata(elsacfa_gs)
elsacfa_dt1 <- format_cfadata(elsacfa_dt1)
elsacfa_dt2 <- format_cfadata(elsacfa_dt2)

# FULL CFA - WLSMV ESTIMATOR ---------------------------------------------

## use model with wlsmv estimator to get model fit

run_fullcfa <- function(data, name, estimator = "WLSMV",
                        bifactor_statements = 
                          c("sp2 BY m7* m8 (sp2); sp2@1;\n sp3 BY m9* m10 m11; sp3@1;\n sp4 BY m1* m2 (sp4); sp4@1; sp5 BY m4* m5 m6; sp5@1;\n", ## updated bifactors
                            "g WITH sp2@0 sp3@0 sp4@0 sp5@0;\n", "sp2 WITH sp3@0 sp4@0 sp5@0;\n", "sp3 WITH sp4@0 sp5@0;\n", "sp4 WITH sp5@0;\n")){
  dt <- copy(data)
  
  cfa_vars <- varmap_dt[!is.na(cfa_label) & cfa_label %in% names(dt), cfa_label]
  cat_vars <- varmap_dt[cfa_label %in% cfa_vars & is.na(continuous), cfa_label]
  dt <- dt[, c("id", cfa_vars), with = F]
  
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
  variable_command <- paste0("\nCATEGORICAL = ", paste0(categorical_items, collapse = " "), ";\nIDVARIABLE = id;")
  
  ## get all factor items and make model command
  factor_items <- c()
  for(i in 1:length(cfa_vars)){
    if (i == 1) {
      unit <- paste0(cfa_vars[i], "*")
    } else if (i%%4==0) {
      unit <- paste0("\n", cfa_vars[i])
    } else {
      unit <- cfa_vars[i]
    }
    factor_items <- c(factor_items, unit)
  }
  model_command <- paste0("g BY ", paste(factor_items, collapse = " "), ";\n g@1; [g@0];")
  
  ## add bifactors
  model_command <- paste0(model_command,
                          paste0(bifactor_statements, collapse = ""),
                          collapse = "")
  
  ## create directories
  dir.create(paste0(model_dir, name, "_", estimator), recursive = T)
  mdir <- paste0(model_dir, name, "_", estimator, "/")
  
  ## run model
  object <- mplusObject(
    TITLE = paste0("ELSA full CFA - not hierarchical"),
    VARIABLE = variable_command,
    ANALYSIS = paste0("ESTIMATOR = ", estimator, "; PARAMETERIZATION=THETA; COVERAGE = 0;"), 
    MODEL = model_command,
    OUTPUT = "STDYX; MODINDICES; TECH1; TECH10;",
    SAVEDATA = "FILE = output.sav; SAVE = FSCORES;",
    usevariables = c(cfa_vars, "id"),
    rdata = as.data.frame(dt)
  )
  
  model <- mplusModeler(
    object = object,
    modelout = paste0(mdir, "model.inp"),
    hashfilename = F,
    run = 1L
  )
  
}

full <- run_fullcfa(elsacfa_dt, "fulldata",
                    bifactor_statements = c("sp2 BY m7* m8 (sp2); sp2@1;\n sp3 BY m9* m10 m11; sp3@1;\n sp4 BY m1* m2 (sp4); sp4@1; sp5 BY m4* m5 m6; sp5@1;\n", ## updated bifactors
                                            "g WITH sp2@0 sp3@0 sp4@0 sp5@0;\n", "sp2 WITH sp3@0 sp4@0 sp5@0;\n", "sp3 WITH sp4@0 sp5@0;\n", "sp4 WITH sp5@0;\n"))

# BAYESIAN ESTIMATOR ------------------------------------------------------

## bayesian estimator has beetter handling of missing data

run_fullcfa_bayes <- function(data, name, estimator = "BAYES", bifactor_statements = 
                                c("sp2 BY m7* m8 (sp2); sp2@1;\n sp3 BY m9* m10 m11; sp3@1;\n sp4 BY m1* m2 (sp4); sp4@1; sp5 BY m4* m5 m6; sp5@1;\n", ## updated bifactors
                                  "g WITH sp2@0 sp3@0 sp4@0 sp5@0;\n", "sp2 WITH sp3@0 sp4@0 sp5@0;\n", "sp3 WITH sp4@0 sp5@0;\n", "sp4 WITH sp5@0;\n")){
  dt <- copy(data)
  
  cfa_vars <- varmap_dt[!is.na(cfa_label) & cfa_label %in% names(dt), cfa_label]
  cat_vars <- varmap_dt[cfa_label %in% cfa_vars & is.na(continuous), cfa_label]
  dt <- dt[, c("id", cfa_vars), with = F]
  
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
  variable_command <- paste0("\nCATEGORICAL = ", paste0(categorical_items, collapse = " "), ";\nIDVARIABLE = id;")
  
  ## get all factor items and make model command
  factor_items <- c()
  for(i in 1:length(cfa_vars)){
    if (i == 1) {
      unit <- paste0(cfa_vars[i], "*")
    } else if (i%%4==0) {
      unit <- paste0("\n", cfa_vars[i])
    } else {
      unit <- cfa_vars[i]
    }
    factor_items <- c(factor_items, unit)
  }
  model_command <- paste0("g BY ", paste(factor_items, collapse = " "), ";\n g@1; [g@0];")
  
  ## add bifactors
  model_command <- paste0(model_command,
                          paste0(bifactor_statements, collapse = ""),
                          collapse = "")
  
  ## create directories
  dir.create(paste0(model_dir, name, "_", estimator), recursive = T)
  mdir <- paste0(model_dir, name, "_", estimator, "/")
  
  ## run model
  object <- mplusObject(
    TITLE = paste0("ELSA full CFA - not hierarchical"),
    VARIABLE = variable_command,
    ANALYSIS = paste0("ESTIMATOR = ", estimator, "; COVERAGE = 0; PROCESSORS = 4;"), 
    MODEL = model_command,
    OUTPUT = "STDYX; MODINDICES; TECH1; TECH8; TECH10;",
    SAVEDATA = "FILE = output.sav; SAVE = FSCORES(200);",
    usevariables = c(cfa_vars, "id"),
    rdata = as.data.frame(dt)
  )
  
  model <- mplusModeler(
    object = object,
    modelout = paste0(mdir, "model.inp"),
    hashfilename = F,
    run = 1L
  )
  
}

fullbayes <- run_fullcfa_bayes(elsacfa_dt, "fulldata")

# GET GOLD STANDARD MODEL AND SCORES ----------------------------------------

goldstandard <- run_fullcfa_bayes(elsacfa_gs, "gold", 
                                  bifactor_statements = 
                                    c("sp2 BY m7* m8 (sp2); sp2@1;\n sp3 BY m9* m10 m11; sp3@1;\n sp4 BY m1* m2 (sp4); sp4@1; sp5 BY m4* m5 m6; sp5@1;\n", ## updated bifactors
                                      "g WITH sp2@0 sp3@0 sp4@0 sp5@0;\n", "sp2 WITH sp3@0 sp4@0 sp5@0;\n", "sp3 WITH sp4@0 sp5@0;\n", "sp4 WITH sp5@0;\n"))

scoredt_gs <- data.table(goldstandard$results$savedata)
savescore_dt <- copy(scoredt_gs[, .(id = ID, gs = G.Mean, gs_se = G.Standard.Deviation)])
savescore_dt <- merge(savescore_dt, elsa_dt[, .(id, train_sample)], by = "id")

write_rds(savescore_dt, paste0(derived_dir, "fscores_", date, ".rds"))

