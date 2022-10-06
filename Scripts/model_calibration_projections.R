# ------------------------------------------------------------------------------
# Course: ECOLOGICAL MODELS APPLIED TO FOSIL DATA: MODEL CALIBRATION,
#         PROJECTIONS, AND UNCERTAINTY
# Author: Marlon E. Cobos, Hannah L. Owens
# Date modified: 05/10/2022
# ------------------------------------------------------------------------------


# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses for the practical section of 
# the course dedicated to model calibration, projections, and uncertainty. Basic 
# plots will be produced to check results. 
# 
# All data required can be obtained using the code in this script. Initial data 
# was obtained and prepared using the scrip "data_preparation.R"
#
# Note: several results will be written in your working directory.
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
if (!dir.exists("ENM")) {
  dir.create("ENM")
}

## preparing data to run Maxent in SWD format
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
oj <- "ENM/hs1_joint.csv" # name of file with all records
otr <- "ENM/hs1_train.csv" # name of file with training records
ote <- "ENM/hs1_test.csv" # name of file with testing records
back_dir <- "ENM/background" # name of directory with backgrounds per set of variables 
bat <- "ENM/bash"  # name of file to write Maxent java code
dir_candidates <- "ENM/Candidate_models"  # name of directory to write candidate models 
rgmul <- c(0.1, 0.5, 1, 2.5, 5) # regularization multipliers to test
fclas <- c("lq", "lp", "lqp") # response types to be used
mx <- "/mnt/backup/Maxent" # complete path to maxent.jar
sel_type <- "OR_AICc"  # type of selection to be done
error <- 5  # omission error allowed (corrected threshold)
dir_eval <- "ENM/Calibration_results"  # name of directory to write calibration results

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
proj_var <- "ENM/Proj_variables" # name of directory with projection scenarios
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
bat_mod <- "ENM/bat_model" # name of file to write Maxent java code 
reps <- 5 # number of replicates
rep_type <- "Bootstrap" # type of replicate to be done 
jack <- TRUE # to get jackknife results 
out_form <- "cloglog" # Maxent output format
proj <- TRUE # to project model to other scenarios 
extrap <- "ext" # type of Maxent extrapolation (free extrapolation)
mess <- TRUE # to get MESS results
mod_dir <- "ENM/Final_models" # name of directory to write final model and projections 
wait <- TRUE # to wait until Maxent finishes  


## final model with projections
kuenm_mod_swd(occ.joint = oj, back.dir = back_dir, out.eval = dir_eval, 
              batch = bat_mod, rep.n = reps, rep.type = rep_type, 
              jackknife = jack, out.format = out_form, project = proj, 
              G.var.dir = proj_var, ext.type = extrap, write.mess = mess, 
              maxent.path = mx, out.dir = mod_dir, wait = wait)

# simple plots of median results
## read raster layers
mod_hs1 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_hs1_median.asc")
mod_msi19 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_msi19_median.asc")

## plot
par(mfrow = c(1, 2))
plot(mod_msi19, main = "Suitability for Mammoth\nPleistocene MIS (ca. 787 ka)")
points(oc_msi19[, 2:3], pch = 16, cex = 0.4)
plot(mod_hs1, main = "Suitability for Mammoth\nPleistocene HS 1 (17.0-14.7 ka)")
points(oc_hs1[, 2:3], pch = 16, cex = 0.4)
# ------------------------------------------------------------------------------



# Extrapolation risk MOP -------------------------------------------------------
# MOP as in Owens et al. (2013)
## arguments
swd <- TRUE
var_set <- "Set_4"
mop_dir <- "ENM/MOP_results"
ref_percent <- 5

## run MOP
kuenm_mmop(G.var.dir = proj_var, M.var.dir = back_dir, is.swd = swd, 
           sets.var = var_set, out.mop = mop_dir, percent = ref_percent)


# A more detailed version of MOP
## bring function from GitHub
source("https://raw.githubusercontent.com/marlonecobos/new_MOP/main/Scripts/MOP_function.R")

## calibration background
back <- read.csv("ENM/background/Set_4.csv")
back <- as.matrix(back[, 4:6])

## variable stacks for projections were prepared in the first section

## other arguments
type_mop <- "detailed"
dist_type <- "euclidean"
scal <- TRUE
cent <- TRUE
rescal <- FALSE

## run MOP for projection to hs1
mop_hs1 <- mop(m = back, g = pcs$PCRasters_initial, mop_type = type_mop, 
               distance = dist_type, scale = scal, center = cent, 
               percent = ref_percent, rescale_mop = rescal)

## run MOP for projection to msi19
mop_msi19 <- mop(m = back, g = pcs$PCRasters_Projected_PCs, mop_type = type_mop, 
                 distance = dist_type, scale = scal, center = cent, 
                 percent = ref_percent, rescale_mop = rescal)
# ------------------------------------------------------------------------------



# Plots of extrapolation risk analyses -----------------------------------------
# just as a remainder, we calibrated models using data and a background from HS1
# models were projected to the whole accessible area for periods HS1 and MSI19

# MESS results
## dissimilarity
novel_hs1 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_0_hs1_novel.asc")
novel_msi19 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_0_msi19_novel.asc") 

## novel conditions
nov_lim_hs1 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_0_hs1_novel_limiting.asc")
nov_lim_msi19 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_0_msi19_novel_limiting.asc") 

# results from kuenm MOP (similarity and novel conditions)
kuenm_mop_hs1 <- raster("ENM/MOP_results/Set_4/MOP_5%_hs1.tif")
kuenm_mop_msi19 <- raster("ENM/MOP_results/Set_4/MOP_5%_msi19.tif") 

# plots of results uncertainty analyses
## dissimilarity results
par(mfrow = c(2, 3))

plot(novel_msi19, main = "MESS\nMIS (ca. 787 ka)")
plot(1 - kuenm_mop_msi19, main = "kuenm MOP\nMIS (ca. 787 ka)")
plot(mop_msi19$mop_basic, main = "New MOP\nMIS (ca. 787 ka)")

plot(novel_hs1, main = "HS 1 (17.0-14.7 ka)")
plot(1 - kuenm_mop_hs1, main = "HS 1 (17.0-14.7 ka)")
plot(mop_hs1$mop_basic, main = "HS 1 (17.0-14.7 ka)")

## novel conditions
par(mfrow = c(2, 3))

plot(nov_lim_msi19, col = c("orange", "red", "white"), 
     main = "MESS\nMIS (ca. 787 ka)")
plot(kuenm_mop_msi19 == 0, col = c("white", "red"), 
     main = "kuenm MOP\nMIS (ca. 787 ka)")
plot(mop_msi19$mop_simple, col = c("orange", "red"), 
     main = "New MOP\nMIS (ca. 787 ka)")

plot(nov_lim_hs1, col = c("orange", "red", "white"), 
     main = "HS 1 (17.0-14.7 ka)")
plot(kuenm_mop_hs1 == 0, col = c("white", "red"), main = "HS 1 (17.0-14.7 ka)")
plot(mop_hs1$mop_simple, col = c("orange", "red"), main = "HS 1 (17.0-14.7 ka)")

## novel conditions towards specific end of variables (derived from new MOP)
### MIS 
plot(mop_msi19$mop_detailed$towards_low_end, col = "blue")
plot(mop_msi19$mop_detailed$towards_high_end, col = "red")

### HS1
plot(mop_hs1$mop_detailed$towards_low_end, col = "blue")
plot(mop_hs1$mop_detailed$towards_high_end, col = "red")
# ------------------------------------------------------------------------------

