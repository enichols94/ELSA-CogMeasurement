##########################################################################
### Author: Emma Nichols
### Project: Development and assessment of analytic methods to improve the 
###          measurement of cognition in longitudinal studies of aging through 
###          the use of sub-studies with comprehensive neuropsychological testing
### Purpose: Test measurement approaches in bootstrapped versions of the ELSA-HCAP
###          sample with different train/test splits
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate,
               simputation, VIM, MplusAutomation, parallel, arm)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# SET OBJECTS -------------------------------------------------------------

dir <- paste0(dropbox_dir, "projects/ips_cog_prediction/")
rawdata_dir <- paste0(dir, "data/source/")
derived_dir <- paste0(dir, "data/derived/")
model_dir <- "C:/Users/emmanich/mplus_models/elsa_cogimp/"
plot_dir <- paste0(dir, "plots/")

## create save dir
dir.create(paste0("C:/Users/emmanich/mplus_models/elsa_cogimp/bootstrap_save/", date))

# GET DATA ----------------------------------------------------------------

elsa_dt <- read_rds(paste0(derived_dir, "elsaprocessed.rds"))

fscore_dt <- read_rds(paste0(derived_dir, "fscores_", date, ".rds"))
elsa_dt <- merge(elsa_dt, fscore_dt[, .(id, gs)], by = "id")

varmap_dt <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))

# STEP 1 SAVE DATASETS ----------------------------------------------------

## to do once - save datasets to file path

create_resampled_data <- function(iteration_num, data = elsa_dt){
  
  print(iteration_num)
  
  analysis_data <- copy(data)
  
  ## create new random test subset
  sample_ids <- sample(x = analysis_data[, id], size = round(nrow(analysis_data)*2/3), replace = F)
  analysis_data[, train_sample := ifelse(id %in% sample_ids, 1, 0)]
  
  ## save data
  save_dir <- paste0("C:/Users/emmanich/mplus_models/elsa_cogimp/bootstrap_save/", date, "/")
  write_rds(analysis_data, paste0(save_dir, "data_", iteration_num, ".rds"))
  
}

datasets <- lapply(1:500, create_resampled_data) ## only do this once before you run everything

# GET ARGUMENTS -----------------------------------------------------------

## to be able to run iterations simultaneously - to launch via command line

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
num_start <- as.numeric(args[1]) ## replication number to start
num_end <- as.numeric(args[2]) ## replication number to end
print(paste0(num_start, " to ", num_end))

# HELPER FUNCTIONS --------------------------------------------------------

run_fullcfa_bayes <- function(data, map, inum = iteration_num, estimator = "BAYES", bifactor_statements = 
                                c("sp2 BY m7* m8 (sp2); sp2@1;\n sp3 BY m9* m10 m11; sp3@1;\n sp4 BY m1* m2 (sp4); sp4@1; sp5 BY m4* m5 m6; sp5@1;\n", ## updated bifactors
                                  "g WITH sp2@0 sp3@0 sp4@0 sp5@0;\n", "sp2 WITH sp3@0 sp4@0 sp5@0;\n", "sp3 WITH sp4@0 sp5@0;\n", "sp4 WITH sp5@0;\n")){
  dt <- copy(data)
  
  cfa_vars <- map[!is.na(cfa_label) & cfa_label %in% names(dt), cfa_label]
  cat_vars <- map[cfa_label %in% cfa_vars & is.na(continuous), cfa_label]
  dt <- dt[, c("id", cfa_vars), with = F]
  
  ## drop from bifactor statements if there are items that are dropped for low contingency cells
  ## this won't work if the item is the first item in a bifactor or one of only two items, but I don't think this should be too much of an issue
  bifactor_items <- gsub(" ", "", (unlist(str_extract_all(bifactor_statements, " [a-z][0-9]")))) 
  delete_items <- setdiff(bifactor_items, cfa_vars)
  if (length(delete_items) > 0) bifactor_statements <- gsub(paste0(" ", delete_items), "", bifactor_statements)
  
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
  variable_command <- paste0("\nCATEGORICAL = ", paste0(categorical_items, collapse = " "), ";\nAUXILIARY = id;")
  
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
  
  ## set directories
  mdir <- "C:/Users/emmanich/mplus_models/elsa_cogimp/bootstrap_models/"
  
  ## run model
  object <- mplusObject(
    TITLE = paste0("ELSA full CFA - not hierarchical"),
    VARIABLE = variable_command,
    ANALYSIS = paste0("ESTIMATOR = ", estimator, "; COVERAGE = 0;"), 
    MODEL = model_command,
    OUTPUT = "STDYX; MODINDICES; TECH1; TECH10;",
    SAVEDATA = paste0("FILE = output", inum, ".sav; SAVE = FSCORES(200);"),
    usevariables = c(cfa_vars, "id"),
    rdata = as.data.frame(dt)
  )
  
  model <- mplusModeler(
    object = object,
    modelout = paste0(mdir, "model", inum, ".inp"),
    hashfilename = F,
    run = 1L
  )
  
  return(model)
}

prep_cfa <- function(data, map = variable_map, wave = "18"){
  
  dt <- copy(data)
  
  ## remove cognitive variables that shouldn't be in various dataset
  nm1 <- c("recall10imm", "recall10del", "month", "year", "dayweek", "date", "animals")
  nm2 <- c(nm1, "scissors", "cactus", "serial7", "primeminister", "monarch", "uspresident", "backcount20"); tm2 <- c("monarch", "uspresident", "backcount20")
  cog_vars <- map[!is.na(cfa_label), label]
  
  if (wave == "18"){
    dt[train_sample == 0, c(setdiff(cog_vars, nm1)) := NA]
    dt[, c(tm2) := NULL] ## exclude because introduced in waves 7-8
  } else if (wave == "78"){
    dt[train_sample == 0, c(setdiff(cog_vars, nm2)) := NA][train_sample == 1, c(tm2) := NA]
  }
  
  ## format data with correct names
  cog_vars <- map[!is.na(cfa_label), label]
  cog_cfalabels <- map[!is.na(cfa_label), cfa_label]
  
  format_cfadata <- function(data){
    dt <- copy(data)
    cv_map <- cog_vars %in% names(dt)
    setnames(dt, cog_vars[cv_map], cog_cfalabels[cv_map])
    dt <- dt[, c("id", cog_cfalabels[cv_map]), with = F]
    return(dt)
  }
  
  dt <- format_cfadata(dt)
  
  ## remove variables with low contingency because cause estimation issues
  deletevars <- c()
  for (v in names(dt)){
    if (dt[!is.na(get(v)), length(unique(get(v)))] == 2){
      if (min(as.numeric(dt[, table(get(v))])) < 10){
        deletevars <- c(deletevars, v)
      } 
    } else if (dt[!is.na(get(v)), length(unique(get(v)))] == 1){
      deletevars <- c(deletevars, v)
    }
  }
  dt[, c(deletevars) := NULL]
  
  return(dt)
}

run_prediction <- function(data, map){
  
  dt <- copy(data)
  
  preds18 <- c("recall10imm", "recall10del", "month", "year", "dayweek", "date", "animals")
  preds78 <- c(preds18, "scissors", "cactus", "serial7", "primeminister", "monarch", "uspresident", "backcount20")
    
  coefsets <- c("preds18", "preds78")
  
  ## make predictions
  make_predictions <- function(predictor_set){
    model_preds <- get(predictor_set)
    model <- lm(as.formula(paste0("gs ~ ", paste0(model_preds, collapse = " + "))), data = dt[train_sample == 1])
    pred_dt <- data.table(id = dt[, id], newvar = predict(model, newdata = dt))
    pred_sims <- arm::sim(model, 1)
    pred_mult <- data.matrix(cbind(rep(1, nrow(dt)), dt[, c(model_preds), with = FALSE])) %*% pred_sims@coef[1,]
    pred_dt[, newvar_sim := rnorm(nrow(dt), pred_mult, 0)] ## don't add data uncertainty into simulations
    setnames(pred_dt, c("newvar", "newvar_sim"), c(predictor_set, paste0(predictor_set, "_sim")))
    return(list(pred_dt, model))
  }
  
  all_predictions <- lapply(coefsets, make_predictions)
  preds <- Reduce(function(x,y) merge(x,y,by="id",all=T,sort=F), lapply(1:length(all_predictions), function(x) as.data.table(all_predictions[[x]][1])))
  return(preds)
}

# BOOTSTRAP FUNCTION ------------------------------------------------------

run_full_analysis <- function(iteration_num, variable_map = varmap_dt){
  
  print(iteration_num)
  
  ## read in data
  save_dir <- paste0("C:/Users/emmanich/mplus_models/elsa_cogimp/bootstrap_save/", date, "/")
  analysis_data <- read_rds(paste0(save_dir, "data_", iteration_num, ".rds"))
  
  ## get word recall for comparisons
  analysis_data[, wordrecall := recall10imm*10+recall10del*10]
  analysis_data[, wordrecall := (wordrecall-mean(wordrecall))/sd(wordrecall)]
  
  ## get simple sum
  analysis_data[, simplesum := (simplesum-mean(simplesum))/sd(simplesum)]
  
  ## prep CFA data
  cfa18_data <- prep_cfa(data = analysis_data, map = variable_map)
  cfa78_data <- prep_cfa(data = analysis_data, wave = "78", map = variable_map)
  
  ## run CFAs
  cfa_w18 <- run_fullcfa_bayes(data = cfa18_data, inum = iteration_num, map = variable_map,
                           bifactor_statements = 
                             c("sp2 BY m7* m8 (sp2); sp2@1;\n sp3 BY m9* m10 m11; sp3@1;\n sp4 BY m1* m2 (sp4); sp4@1; sp5 BY m4* m5 m6; sp5@1;\n", ## updated bifactors
                               "g WITH sp2@0 sp3@0 sp4@0 sp5@0;\n", "sp2 WITH sp3@0 sp4@0 sp5@0;\n", "sp3 WITH sp4@0 sp5@0;\n", "sp4 WITH sp5@0;\n"))
  cfa_w78 <- run_fullcfa_bayes(data = cfa78_data, inum = iteration_num, map = variable_map, bifactor_statements = 
                             c("sp2 BY m7* m8 (sp2); sp2@1;\n sp3 BY m9* m10 m11; sp3@1;\n sp4 BY m1* m2 (sp4); sp4@1; sp5 BY m4* m5 m6; sp5@1;\n", ## updated bifactors
                               "g WITH sp2@0 sp3@0 sp4@0 sp5@0;\n", "sp2 WITH sp3@0 sp4@0 sp5@0;\n", "sp3 WITH sp4@0 sp5@0;\n", "sp4 WITH sp5@0;\n"))
  
  ## get predictions
  predictions <- run_prediction(data = analysis_data, map = variable_map)
  
  ## add on predictions from cfas
  scoredt_w18 <- data.table(cfa_w18$results$savedata)
  scoredt_w78 <- data.table(cfa_w78$results$savedata)
  predictions <- merge(predictions,
                        scoredt_w18[, .(id = ID, w18 = G.Mean, w18_se = G.Standard.Deviation)], by = "id")
  predictions <- merge(predictions, 
                       scoredt_w78[, .(id = ID, w78 = G.Mean, w78_se = G.Standard.Deviation)], by = "id")
  predictions <- merge(predictions, analysis_data[, .(id, gs, wordrecall, simplesum, train_sample)], by = "id")
  
  ## save this iteration
  save_dir <- paste0("C:/Users/emmanich/mplus_models/elsa_cogimp/bootstrap_save/", date, "/")
  write_rds(predictions, paste0(save_dir, iteration_num, ".rds"))
  
}

lapply(num_start:num_end, function(x) run_full_analysis(iteration_num = x))

