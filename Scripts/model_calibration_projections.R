# ------------------------------------------------------------------------------
# Course: ECOLOGICAL MODELS APPLIED TO FOSIL DATA: MODEL CALIBRATION AND
#         PROJECTIONS
# Author: Marlon E. Cobos, Hannah L. Owens
# Date: 04/10/2022
# ------------------------------------------------------------------------------


# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses for the practical section of 
# the course dedicated to model calibration and projections. Basic plots will be 
# produced to check data. 
# 
# All data required can be obtained using the code in this script. Initial data 
# was obtained and prepared using the scrip "data_preparation.R"
#
# Note: several results will be written in your working directory.
# ------------------------------------------------------------------------------


# R packages needed ------------------------------------------------------------
# these lines load packages (if needed packages are installed first)
if (!require("remotes")) {
  install.packages("remotes")
  Sys.sleep(2)
  library(remotes)
}
if (!require("kuenm")) { # it may be tricky (requires compilation tools)
  remotes::install_github("marlonecobos/kuenm")
  Sys.sleep(2)
  library(kuenm)
}

library(raster)
# ------------------------------------------------------------------------------


# Preparing data for modeling --------------------------------------------------
# reading occurrence data
oc_msi19 <- read.csv("Occurrence_data/oc_msi19_thin1.csv")
oc_hs1 <- read.csv("Occurrence_data/oc_hs1_thin1.csv")

# reading raster layers
m_msi19 <- stack(list.files("Environmental_data/accessible_area_mis19", 
                            pattern = ".tif$", full.names = TRUE))
m_hs1 <- stack(list.files("Environmental_data/accessible_area_hs1", 
                          pattern = ".tif$", full.names = TRUE))

# reducing environmental dimensions using 
# the overall variance of all variables is explained with fewer variables
# PCs for msi19 are created using the loadings from the pca created with hs1
pcs <- kuenm_rpca(variables = m_hs1, var.scale = TRUE, project = TRUE, 
                  proj.vars = m_msi19, write.result = TRUE, n.pcs = 3,
                  out.format = "GTiff", out.dir = "Environmental_data/PCA")

# preparing data for model calibration
## new directory
dir.create("ENM")

data_prep <- prepare_swd(occ = oc_hs1, species = "taxon_name", longitude = "lng",
                         latitude = "lat", data.split.method = "random", 
                         train.proportion = 0.7, 
                         raster.layers = pcs$PCRasters_initial, 
                         sample.size = 10000, var.sets = "all_comb", 
                         min.number = 2, save = TRUE, 
                         name.occ = "ENM/hs1", 
                         back.folder = "ENM/background")
# ------------------------------------------------------------------------------



# Model calibration ------------------------------------------------------------
# define arguments for simplicity
oj <- "ENM/hs1_joint.csv"
otr <- "ENM/hs1_train.csv"
ote <- "ENM/hs1_test.csv"
back_dir <- "ENM/background"
bat <- "ENM/bash"
dir_candidates <- "ENM/Candidate_models"
rgmul <- c(0.1, 0.5, 1, 2.5, 5)
fclas <- c("lq", "lp", "lqp")
mx <- "/mnt/backup/Maxent"
sel_type <- "OR_AICc"
error <- 5
dir_eval <- "ENM/Calibration_results"

# run calibration (create candidate models and selection of best models)
cal <- kuenm_cal_swd(occ.joint = oj, occ.tra = otr, occ.test = ote, 
                     back.dir = back_dir, batch = bat, 
                     out.dir.models = dir_candidates, reg.mult = rgmul, 
                     f.clas = fclas, maxent.path = mx, selection = sel_type, 
                     threshold = error, out.dir.eval = dir_eval)

# explore results
## summary
cal$calibration_statistics

## all calibration results
View(cal$calibration_results)

## selected model(s)
cal$selected_models

## what variables are in set 4
data_prep$sets$Set_4
# ------------------------------------------------------------------------------


# Model projections--------------------------------------------------------------
# prepare data for projections
## directories
proj_var <- "ENM/Proj_variables"
dir.create(proj_var) # general directory for projection layers 

dir.create("ENM/Proj_variables/Set_4") # directory for Set of variable selected
dir.create("ENM/Proj_variables/Set_4/hs1") # calibration scenario
dir.create("ENM/Proj_variables/Set_4/msi19") # projection scenario

## write PC raster layers created for msi19 scenario (all 3 PCs as in Set_4)
pcnames <- paste0(names(pcs$PCRasters_Projected_PCs), ".asc")

for (i in 1:length(pcnames)) {
  writeRaster(pcs$PCRasters_initial[[i]], 
              filename = paste0("ENM/Proj_variables/Set_4/hs1/", pcnames[i]), 
              format = "ascii")
  writeRaster(pcs$PCRasters_Projected_PCs[[i]], 
              filename = paste0("ENM/Proj_variables/Set_4/msi19/", pcnames[i]), 
              format = "ascii")
}


# project selected model
## arguments
bat_mod <- "ENM/bat_model"
reps <- 5
rep_type <- "Bootstrap"
jack <- TRUE
out_form <- "cloglog"
proj <- TRUE
extrap <- "ext"
mess <- TRUE
mod_dir <- "ENM/Final_models"
wait <- TRUE


## final model with projections
kuenm_mod_swd(occ.joint = oj, back.dir = back_dir, out.eval = dir_eval, 
              batch = bat_mod, rep.n = reps, rep.type = rep_type, 
              jackknife = jack, out.format = out_form, project = proj, 
              G.var.dir = proj_var, ext.type = extrap, write.mess = mess, 
              maxent.path = mx, out.dir = mod_dir, wait = wait)
# ------------------------------------------------------------------------------
