##########################################################################
### Author: Emma Nichols
### Project: Development and assessment of analytic methods to improve the 
###          measurement of cognition in longitudinal studies of aging through 
###          the use of sub-studies with comprehensive neuropsychological testing
### Purpose: Estimate regression models to use for prediction to the ELSA sample
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

# GET DATA ----------------------------------------------------------------

elsa_dt <- read_rds(paste0(derived_dir, "elsaprocessed.rds"))

varmap_dt <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))

## merge factor scores into ELSA data
fscore_dt <- read_rds(paste0(derived_dir, "fscores_2023_11_14.rds")) ## need to change the date here when have new GS factor score
elsa_dt <- merge(elsa_dt, fscore_dt[, .(id, gs)], by = "id")

# GET PREDICTORS ----------------------------------------------------------

preds18 <- c("recall10imm", "recall10del", "month", "year", "dayweek", "date", "animals")
preds78 <- c(preds18, "scissors", "cactus", "serial7", "primeminister", "monarch", "uspresident", "backcount20")

coefsets <- c("preds18", "preds78")

# SAVE FULLDATA MODELS ----------------------------------------------------

model18 <- lm(as.formula(paste0("gs ~ ", paste0(preds18, collapse = " + "))), data = elsa_dt)
model78 <- lm(as.formula(paste0("gs ~ ", paste0(preds78, collapse = " + "))), data = elsa_dt)

write_rds(list(model18 = model18, model78 = model78), paste0(derived_dir, "regressionmodels_", date, ".rds"))
