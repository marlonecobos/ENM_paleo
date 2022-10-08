# ------------------------------------------------------------------------------
# Course: ECOLOGICAL MODELS APPLIED TO FOSSIL DATA: MODEL EVALUATION
# Author: Marlon E. Cobos, Hannah L. Owens
# Date modified: 07/10/2022
# ------------------------------------------------------------------------------


# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses for the practical section of 
# the course dedicated to model evaluation. Basic plots will be produced 
# to check results. 
# 
# All data required can be obtained using the code in this script. Initial data 
# was obtained and prepared using the scrip "data_preparation.R". Previous 
# results were obtained using the script ""
#
# Note: some results will be written in your working directory.
# ------------------------------------------------------------------------------


# R packages needed ------------------------------------------------------------
# the next lines of code load packages 

## if needed, install packages (remotes, raster) using 
## install.packages("package_name")
## kuenm is on GitHub so install it as 
## remotes::install_github("marlonecobos/kuenm") 
## https://github.com/marlonecobos/kuenm#installing-the-package

library(raster)
library(kuenm)
# ------------------------------------------------------------------------------


# Working directory and sub-directories -----------------------------------------
# working directory
setwd("YOUR/DIRECTORY")

## new sub-directories
thres_dir <- "ENM/Thresholded"

if (!dir.exists(thres_dir)) {
  dir.create(thres_dir)
}

ieval_dir <- "ENM/Independent_evaluation"

if (!dir.exists(ieval_dir)) {
  dir.create(ieval_dir)
}
# ------------------------------------------------------------------------------


# Variable importance (contribution to models) ---------------------------------
# after model calibration we created a final model with selected parameter 
# settings and variables, the information generated there is used to learn about
# how variables contributed to models

# variable contribution analysis
var_imp <- model_var_contrib(fmod.dir = "ENM/Final_models", project = TRUE, 
                             ext.type = "E")

var_imp

# plot of variable contribution
plot_contribution(contribution_list = var_imp[[1]][[1]])
# ------------------------------------------------------------------------------



# Thresholding results ---------------------------------------------------------
# we calibrated models using data and a background from HS1
# models were projected to the whole accessible area for periods HS1 and MSI19

# data
## records used in final model
occ_hs1 <- read.csv("ENM/hs1_joint.csv")

## raster layers (median projections of replicated final model)
mod_hs1 <- raster(list.files(path = "ENM/Final_models", 
                             pattern = "hs1_median.asc", 
                             recursive = TRUE, full.names = TRUE))
mod_msi19 <- raster(list.files(path = "ENM/Final_models", 
                               pattern = "msi19_median.asc", 
                               recursive = TRUE, full.names = TRUE))

## suitability values in records
occ_suit <- data.frame(occ_hs1, 
                       Cloglog.prediction = extract(mod_hs1, occ_hs1[, 2:3]))


# determine threshold values
## minimum training presence
mtp_thres <- min(occ_suit$Cloglog.prediction)

## modified threshold (E = 5% omission)
e5_n <- nrow(occ_suit) * 0.05
e5_n <- ceiling(e5_n) + 1# many people do not sum 1 
e5_thres <- sort(occ_suit$Cloglog.prediction)[e5_n]

# threshold suitability layers
## minimum training presence threshold
mod_hs1_mtp <- mod_hs1 >= mtp_thres
mod_msi9_mtp <- mod_msi19 >= mtp_thres

## modified threshold (E = 5% omission)
mod_hs1_e5 <- mod_hs1 >= e5_thres
mod_msi9_e5 <- mod_msi19 >= e5_thres


# write results as raster layers
writeRaster(mod_hs1_mtp, filename = paste0(thres_dir, "/med_hs1_mtp.tif"), 
            format = "GTiff", overwrite = TRUE)
writeRaster(mod_hs1_e5, filename = paste0(thres_dir, "/med_hs1_e5%.tif"), 
            format = "GTiff", overwrite = TRUE)
writeRaster(mod_msi9_mtp, filename = paste0(thres_dir, "/med_msi9_mtp.tif"), 
            format = "GTiff", overwrite = TRUE)
writeRaster(mod_msi9_e5, filename = paste0(thres_dir, "/med_msi9_e5%.tif"), 
            format = "GTiff", overwrite = TRUE)


# simple plot of suitable areas
par(mfrow = c(2, 2))

plot(mod_msi9_mtp, main = "Suitable areas (MTP)\nMIS (ca. 787 ka)")
plot(mod_msi9_e5, main = "Suitable areas (5%)\nMIS (ca. 787 ka")

plot(mod_hs1_mtp, main = "HS 1 (17.0-14.7 ka)")
plot(mod_hs1_e5, main = "HS 1 (17.0-14.7 ka)")
# ------------------------------------------------------------------------------



# Model evaluation -------------------------------------------------------------
# for this section we are going to assume that the records from ~0.13 mya can be
# used as independent testing records for our final model calibrated and 
# projected to the period HS1 ~0.015 mya (not necessarily true)

# independent testing records
test_oc <- read.csv("Occurrence_data/mammut_americanum_0.13mya.csv") 


# statistical significance (partial ROC) 
p_roc <- kuenm_proc(occ.test = test_oc[, 2:3], model = mod_hs1, threshold = 5)

p_roc$pROC_summary


# predicting ability (omission rates)
om_rates <- kuenm_omrat(model = mod_hs1, threshold = c(0, 5, 10), 
                        occ.tra = occ_hs1[, 2:3], occ.test = test_oc[, 2:3])

om_rates

# fitting and complexity (using candidate models)
## selected model (find it ENM/Calibration_results/selected_models.csv)
sel_mod <- "M_0.1_F_lq_Set_4"

## model constructed with raw outputs
raw_mx <- paste0("ENM/Candidate_models/", sel_mod, "_all/Mammut_americanum.csv")
raw_model_output <- read.csv(raw_mx)

## number of parameters in the model
ld <- paste0("ENM/Candidate_models/", sel_mod, "_all/Mammut_americanum.lambdas")
lambdas <- readLines(ld)
npar <- n_par(x = lambdas)

## aicc calculation (we have this, just for practice, no test records needed)
model_aicc <- aicc(occ = occ_hs1[, 2:3], prediction = raw_model_output,
                   npar = npar)

model_aicc


# write all results together
## all results
all_eval <- data.frame(Model = sel_mod, 
                       t(data.frame(p_roc$pROC_summary)), 
                       t(data.frame(om_rates)), model_aicc)

all_eval

## write all evaluation results
evalname <- paste0(ieval_dir, "/evaluation_results.csv")

if (!file.exists(evalname)) {
  write.csv(all_eval, evalname, row.names = FALSE)
}
# ------------------------------------------------------------------------------