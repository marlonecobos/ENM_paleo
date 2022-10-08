# ------------------------------------------------------------------------------
# Course: ECOLOGICAL MODELS APPLIED TO FOSSIL DATA: NICHE COMPARISONS
# Author: Marlon E. Cobos, Hannah L. Owens
# Date modified: 07/10/2022
# ------------------------------------------------------------------------------


# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses for the practical section of 
# the course dedicated to niche comparisons. Basic plots will be 
# produced to check results. 
# 
# All data required can be obtained using the code in this script. Initial data 
# was obtained and prepared using the scrip "data_preparation.R"
#
# Note: some results will be written in your working directory.
# ------------------------------------------------------------------------------


# R packages needed ------------------------------------------------------------
# the next lines of code load packages 

## if needed, install packages (remotes, raster) using 
## install.packages("package_name")
## ellipsenm is on GitHub so install it as 
## remotes::install_github("marlonecobos/ellipsenm") 
## https://github.com/marlonecobos/ellipsenm#installing-the-package

library(raster)
library(ellipsenm)
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

over_dir <- "Niche_comparison"

if (!dir.exists(over_dir)) {
  dir.create(over_dir)
}
# ------------------------------------------------------------------------------



# Extent of suitable habitat in projections ------------------------------------
# we thresholded suitability layers for the periods HS1 and MSI19

# read thresholded suitability layers
## minimum training presence threshold
mod_hs1_mtp <- raster("ENM/Thresholded/med_hs1_mtp.tif")
mod_msi9_mtp <- raster("ENM/Thresholded/med_msi9_mtp.tif")

## modified threshold (E = 5% omission)
mod_hs1_e5 <- raster("ENM/Thresholded/med_hs1_e5%.tif")
mod_msi9_e5 <- raster("ENM/Thresholded/med_msi9_e5%.tif")

# area calculation
## minimum training presence threshold
### keep only suitable pixels
mod_hs1_mtp_suit <- mod_hs1_mtp
mod_msi9_mtp_suit <- mod_msi9_mtp

mod_hs1_mtp_suit[mod_hs1_mtp_suit[] == 0] <- NA
mod_msi9_mtp_suit[mod_msi9_mtp_suit[] == 0] <- NA

### area
area_hs1_mtp <- area(mod_hs1_mtp_suit)
area_msi9_mtp <- area(mod_msi9_mtp_suit)

area_hs1_mtp <- sum(area_hs1_mtp[][mod_hs1_mtp_suit[] == 1], na.rm = TRUE)
area_msi9_mtp <- sum(area_msi9_mtp[][mod_msi9_mtp_suit[] == 1], na.rm = TRUE)

## modified threshold (E = 5% omission)
### keep only suitable pixels
mod_hs1_e5_suit <- mod_hs1_e5
mod_msi9_e5_suit <- mod_msi9_e5

mod_hs1_e5_suit[mod_hs1_e5_suit[] == 0] <- NA
mod_msi9_e5_suit[mod_msi9_e5_suit[] == 0] <- NA

### area
area_hs1_e5 <- area(mod_hs1_e5_suit)
area_msi9_e5 <- area(mod_msi9_e5_suit)

area_hs1_e5 <- sum(area_hs1_e5[][mod_hs1_e5_suit[] == 1], na.rm = TRUE)
area_msi9_e5 <- sum(area_msi9_e5[][mod_msi9_e5_suit[] == 1], na.rm = TRUE)


# write results
# write raster layers
writeRaster(mod_hs1_mtp_suit, format = "GTiff", overwrite = TRUE,
            filename = paste0(thres_dir, "/med_hs1_mtp_suitable.tif"))
writeRaster(mod_hs1_e5_suit, format = "GTiff", overwrite = TRUE,
            filename = paste0(thres_dir, "/med_hs1_e5%_suitable.tif"))
writeRaster(mod_msi9_mtp_suit, format = "GTiff", overwrite = TRUE,
            filename = paste0(thres_dir, "/med_msi9_mtp_suitable.tif"))
writeRaster(mod_msi9_e5_suit, format = "GTiff", overwrite = TRUE,
            filename = paste0(thres_dir, "/med_msi9_e5%_suitable.tif"))

## table with areas
areas <- data.frame(HS1_mtp_area_km2 = area_hs1_mtp, 
                    MIS_mtp_area_km2 = area_msi9_mtp,
                    HS1_e5_area_km2 = area_hs1_e5, 
                    MIS_e5_area_km2 = area_msi9_e5)

areas

a_name <- "ENM/Thresholded/Suitable_area_km2.csv"

if (!file.exists(a_name)) {
  write.csv(areas, a_name, row.names = FALSE)
}


# simple plot of suitable areas
par(mfrow = c(2, 2))

plot(mod_msi9_mtp_suit, main = "Suitable areas (MTP)\nMIS (ca. 787 ka)")
plot(mod_msi9_e5_suit, main = "Suitable areas (5%)\nMIS (ca. 787 ka")

plot(mod_hs1_mtp_suit, main = "HS 1 (17.0-14.7 ka)")
plot(mod_hs1_e5_suit, main = "HS 1 (17.0-14.7 ka)")
# ------------------------------------------------------------------------------



# Predicting records with projections ------------------------------------------
# we are going to find out how well records from MIS (ca. 787 ka) are predicted
# by a projection of the model constructed with data from HS 1 (17.0-14.7 ka)

# records from 
oc_msi19 <- read.csv("Occurrence_data/oc_msi19_thin1.csv")


# MTP threshold
## extract suitable/unsuitable value from layers
suitable_mtp <- extract(mod_msi9_mtp, oc_msi19[, 2:3])

## find omitted records
om_oc_mtp <- oc_msi19[suitable_mtp == 0, ]

## omission rate
or_mtp <- sum(suitable_mtp == 0) / length(suitable_mtp)


# 5% threshold
## extract suitable/unsuitable value from layers
suitable_e5 <- extract(mod_msi9_e5, oc_msi19[, 2:3])

## find omitted records
om_oc_e5 <- oc_msi19[suitable_e5 == 0, ]

## omission rate
or_e5 <- sum(suitable_e5 == 0) / length(suitable_e5)


# write results
## write omission rates
ors <- data.frame(Total_records_MIS = nrow(oc_msi19), 
                  Omission_rate_mtp = or_mtp, Omission_rate_e5 = or_e5)

ors

or_name <- paste0(over_dir, "/MIS_omitted_in_predictions_from_HS1.csv")

if (!file.exists(or_name)) {
  write.csv(ors, or_name, row.names = FALSE)
}


# simple plot of suitable areas and records 
par(mfrow = c(1, 2))

plot(mod_msi9_mtp_suit, main = "Records and suitable areas (MTP)\nMIS (ca. 787 ka)")
points(oc_msi19[, 2:3], pch = 1, cex = 0.6)
plot(mod_msi9_e5_suit, main = "Records and suitable areas (5%)\nMIS (ca. 787 ka")
points(oc_msi19[, 2:3], pch = 1, cex = 0.6)
# ------------------------------------------------------------------------------



# Ecological niche overlap -----------------------------------------------------
# data
## records 
oc_msi19 <- read.csv("Occurrence_data/oc_msi19_thin1.csv")
oc_hs1 <- read.csv("Occurrence_data/oc_hs1_thin1.csv")

## raster layers masked to calibration areas representative of each period
mc_msi19 <- stack(list.files("Environmental_data/relevant_area_mis19", 
                             pattern = ".tif$", full.names = TRUE))
mc_hs1 <- stack(list.files("Environmental_data/relevant_area_hs1", 
                           pattern = ".tif$", full.names = TRUE))


# reducing environmental dimensions using PCA
# the overall variance of all variables is explained with fewer variables
# PCs for msi19 are created using the loadings from the PCA created with hs1
pcs_rel <- kuenm_rpca(variables = mc_hs1, var.scale = TRUE, project = TRUE, 
                      proj.vars = mc_msi19, write.result = TRUE, n.pcs = 3,
                      out.format = "GTiff", 
                      out.dir = "Environmental_data/PCA_relevant_area")

pc_msi19 <- pcs_rel$PCRasters_Projected_PCs
pc_hs1 <- pcs_rel$PCRasters_initial


# preparing data
## records with environmental data
oc_msi19_env <- data.frame(oc_msi19, extract(pc_msi19, oc_msi19[, 2:3]))
oc_hs1_env <- data.frame(oc_hs1, extract(pc_hs1, oc_hs1[, 2:3]))

## background for comparisons
back_msi19 <- na.omit(pc_msi19[])
back_hs1 <- na.omit(pc_hs1[])

## objects to perform overlap analysis
niche1_msi19 <- overlap_object(data = oc_msi19_env, species = "taxon_name",
                               longitude = "lng", latitude = "lat", 
                               method = "covmat", level = 95, 
                               variables = back_msi19)

niche2_hs1 <- overlap_object(data = oc_hs1_env, species = "taxon_name",
                             longitude = "lng", latitude = "lat", 
                             method = "covmat", level = 95, 
                             variables = back_hs1)


# overlap analysis
## simple analysis
overlap <- ellipsoid_overlap(niche1_msi19, niche2_hs1, # the two niches
                             overlap_type = "back_union") # type of overlap

summary(overlap)


## analysis with significance test
overlap_sig <- ellipsoid_overlap(niche1_msi19, niche2_hs1,
                                 overlap_type = "back_union", 
                                 significance_test = TRUE, # this runs the significance test
                                 replicates = 100) # number of randomization processes

summary(overlap_sig)


# write results
## save overlap results (easiest way to save multiple results)
o_name <- paste0(over_dir, "/overlap.RData")

if (!file.exists(o_name)) {
  save(niche2_hs1, niche1_msi19, overlap_sig, file = o_name)
}
  

# simple plots
# 3D plot of niches and their overlap
plot_overlap(overlap_sig)

# significance test results
par(mfrow = c(1, 1))
plot_overlap_sig(overlap_sig, 
                 main = "Overlap analysis\nSignificance test results")
legend(x = 0.30, y = 8, legend = c("Observed", "5% CL"), lty = c(1, 2), lwd = 2, 
       col = c("blue", "darkgreen"), bty = "n")
# ------------------------------------------------------------------------------
