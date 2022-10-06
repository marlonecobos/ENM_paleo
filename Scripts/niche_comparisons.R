# ------------------------------------------------------------------------------
# Course: ECOLOGICAL MODELS APPLIED TO FOSIL DATA: NICHE COMPARISONS
# Author: Marlon E. Cobos, Hannah L. Owens
# Date modified: 05/10/2022
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
# ------------------------------------------------------------------------------


# Extent of suitable habitat in projections ------------------------------------
# we calibrated models using data and a background from HS1
# models were projected to the whole accessible area for periods HS1 and MSI19

# data
## records with suitability values
occ_suit <- read.csv("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_0_samplePredictions.csv")

## raster layers (median projections of replicated final model)
mod_hs1 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_hs1_median.asc")
mod_msi19 <- raster("ENM/Final_models/M_5_F_lqp_Set_4_E/Mammut_americanum_msi19_median.asc")

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


# simple plot of suitable areas
par(mfrow = c(2, 2))

plot(mod_msi9_mtp, main = "Suitable areas (MTP)\nMIS (ca. 787 ka)")
plot(mod_msi9_e5, main = "Suitable areas (5%)\nMIS (ca. 787 ka")

plot(mod_hs1_mtp, main = "HS 1 (17.0-14.7 ka)")
plot(mod_hs1_e5, main = "HS 1 (17.0-14.7 ka)")
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
# ------------------------------------------------------------------------------



# Ecological niche overlap -----------------------------------------------------
# data
## records 
oc_msi19 <- read.csv("Occurrence_data/oc_msi19_thin1.csv")
oc_hs1 <- read.csv("Occurrence_data/oc_hs1_thin1.csv")

## raster layers masked to calibration areas representative of each period
pc_msi19 <- stack(list.files("Environmental_data/PCA/Projected_PCs", 
                             pattern = ".tif$", full.names = TRUE))
pc_hs1 <- stack(list.files("Environmental_data/PCA/Initial", 
                           pattern = ".tif$", full.names = TRUE))

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


# simple plots
# 3D plot of niches and their overlap
plot_overlap(overlap_sig)

# significance test results
par(mfrow = c(1, 1))
plot_overlap_sig(overlap_sig, 
                 main = "Overlap analysis\nSignificance test results")
legend(x = 0.22, y = 14, legend = c("Observed", "5% CL"), lty = c(1, 2), lwd = 2, 
       col = c("blue", "darkgreen"), bty = "n")
# ------------------------------------------------------------------------------