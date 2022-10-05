# ------------------------------------------------------------------------------
# Course: ECOLOGICAL MODELS APPLIED TO FOSIL DATA: DATA PREPARATION
# Author: Marlon E. Cobos, Hannah L. Owens
# Date: 04/10/2022
# ------------------------------------------------------------------------------


# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses for the practical section of 
# the course dedicated to data preparation. Basic plots will be produced to 
# check data. 
# 
# All data required can be obtained using the code in this script. Make sure to 
# have internet connection to download such data.
#
# Note: results will be written in your working directory.
# ------------------------------------------------------------------------------


# R packages needed ------------------------------------------------------------
# these lines load packages (if needed packages are installed first)
if (!require("remotes")) {
  install.packages("remotes")
  Sys.sleep(2)
  library(remotes)
}
if (!require("paleobioDB")) {
  install.packages("paleobioDB")
  Sys.sleep(2)
  library(paleobioDB)
}
if (!require("ellipsenm")) { # it may be tricky (requires compilation tools)
  remotes::install_github("marlonecobos/ellipsenm")
  Sys.sleep(2)
  library(ellipsenm)
}
if (!require("spThin")) {
  install.packages("spThin")
  Sys.sleep(2)
  library(spThin)
}

library(raster)
# ------------------------------------------------------------------------------


# Working directory and sub-directories -----------------------------------------
# working directory
setwd("YOUR/DIRECTORY")

## sub-directory for occurrence data
if(!file.exists("Occurrence_data")){
  dir.create("Occurrence_data")
}
## sub-directory for environmental data
if(!file.exists("Environmental_data")){
  dir.create("Environmental_data")
}
# ------------------------------------------------------------------------------


# Getting occurrence data ------------------------------------------------------
# search for occurrences
mammoth_oc <-  pbdb_occurrences (limit = "all", base_name = "Mammut americanum", 
                                vocab = "pbdb", interval = "Quaternary",
                                show = c("coords", "phylo", "ident"))

# info in downloaded data
colnames(mammoth_oc)

# filter records for unique occurrences and relevant columns
## info to keep 
selected <- c("taxon_name", "lng", "lat", "early_age", "late_age") 

## filter data to keep selected columns
mammoth_oc <- mammoth_oc[, selected]
# ------------------------------------------------------------------------------


# Occurrence data cleaning and filtering ---------------------------------------
# filter to keep only unique events in the data
mammoth_oc <- unique(mammoth_oc)

# check what kind of data we have
str(mammoth_oc)

# some columns that should contain numeric info are in character format 
## this fixes the issue of all columns with character info
numeric_info <- as.numeric(unlist(mammoth_oc[, c("lng", "lat", "early_age", 
                                                 "late_age")]))
## put corrected info back with in the main table
mammoth_oc[, c("lng", "lat", "early_age", "late_age")] <- numeric_info


# check taxon name 
## all names in the column of species name
table(mammoth_oc$taxon_name)

## “cf” basically means “I’m not sure, the preservation isn’t that great
## exclude that one
mammoth_oc <- mammoth_oc[mammoth_oc$taxon_name != "Mammut cf. americanum", ]

## all other names are synonyms 
## change all to "Mammut_americanum"
mammoth_oc$taxon_name <- "Mammut_americanum"
table(mammoth_oc$taxon_name)


# simple plot to check the data using environmental layers and records
## limits for plot
lims <- apply(mammoth_oc[, c("lng","lat")], 2, range)
lims

## making limits a little wider
lims <- lims * matrix(c(1.1, 0.9, 0.9, 1.1), nrow = 2)
lims

## raster layer 
x11()
plot(variables$pleistocene_mis19$bio_1, xlim = lims[, 1], ylim = lims[, 2], 
     main = "All Mammoth Points") 

## the records
points(mammoth_oc[, c("lng","lat")], pch = 20, col = "black") 


# split points into relevant time slices
## estimated age of records
mammoth_oc$estimated_age <- (mammoth_oc$early_age + mammoth_oc$late_age) / 2

## date ranges from paleoclim.org
time_slices <- data.frame (maxTime = c(3.3, 3.205, .787, 0.130, 0.017, 0.0129, 
                                       0.0117, 0.00826, 0.0042),
                           minTime = c(3.3, 3.205, .787, 0.130, 0.0147, 0.0117, 
                                       0.00826, 0.0042, 0.0003))

time_slices$meanAge <- (time_slices$maxTime + time_slices$minTime) / 2

## assign time slices to occurrences based on mean age
for (i in 1:nrow(mammoth_oc)) {
  time_diff <- abs(time_slices$meanAge - mammoth_oc[i, "estimated_age"])
  mammoth_oc$time_index[i] <- which.min(time_diff)
  mammoth_oc$time_slice[i] <- time_slices$meanAge[mammoth_oc$time_index[i]]
}


# plot to check how the data look like now
## define colors
indexes <- unique(mammoth_oc$time_index)
ncolors <- length(indexes)
col_pal <- colorRampPalette(c('red','yellow', 'blue'))
mammoth_oc$col <- col_pal(ncolors)[as.factor(mammoth_oc$time_index)]

## plot with legend
x11()
plot(variables$pleistocene_mis19$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "All Mammoth Points by Date", col = "#D6D5D5") 
points(mammoth_oc[, c("lng","lat")], col = mammoth_oc$col, pch = 20)
legend("bottomleft", title = "Date (My)", col = col_pal(ncolors), pch = 20,
       legend = time_slices$meanAge[sort(indexes)], bt = "n", cex = 0.8)


# writing data to directory
write.csv(mammoth_oc, file = "Occurrence_data/mammut_americanum_all.csv", 
          row.names = FALSE)
# ------------------------------------------------------------------------------


# Downloading environmental data -----------------------------------------------
# get data from paleoclim.org
## lets check how many records per time slice 
table(mammoth_oc$time_slice)

## relevant links (pleistocene_mis19 = 0.787 mya, pleistocene_hs1 = 0.01585 mya)
## check all data here http://www.paleoclim.org/
periods <- c(
  pleistocene_mis19 = "http://sdmtoolbox.org/paleoclim.org/data/MIS19/MIS19_v1_r10m.zip",
  pleistocene_hs1 = "http://sdmtoolbox.org/paleoclim.org/data/HS1/HS1_v1_10m.zip"
)

## download data
variables <- lapply(1:length(periods), function(x) {
  pname <- names(periods)[x] # name of periods
  dfile <- paste0("Environmental_data/", pname, "_bio.zip") # name downloaded file
  
  download.file(periods[x], destfile = dfile, method = "auto", # download
                quiet = FALSE, mode = "wb", cacheOK = TRUE)
  
  dfol <- paste0("Environmental_data/", pname, "_bio") # folder to unzip data
  dir.create(dfol)
  unzip(zipfile = dfile, exdir = dfol) # unzip
  
  files <- list.files(dfol, pattern = ".tif$", full.names = TRUE, 
                      recursive = TRUE) # raster layers downloaded
  stack(files) # stack of layers
})

names(variables) <- names(periods) # name for variable sets

# check data
## variables in each set of layers
sapply(variables, names)

## simple plot
x11()
par(mfrow = c(2, 1))
plot(variables$pleistocene_mis19$bio_1) 
plot(variables$pleistocene_hs1$bio_1) 
# ------------------------------------------------------------------------------


# Processing environmental data ------------------------------------------------
# keep variables that match between periods
## match variable names
match_var <- match(names(variables$pleistocene_mis19), 
                   names(variables$pleistocene_hs1))

## variables to keep and sort names
vars_keep <- names(variables$pleistocene_hs1)[match_var]
vars_keep

vars_keep <- vars_keep[c(1, 12:14, 2:11)]
vars_keep

## filtering variables in stacks
variables$pleistocene_mis19 <- variables$pleistocene_mis19[[vars_keep]] 
variables$pleistocene_hs1 <- variables$pleistocene_hs1[[vars_keep]] 


# region of interest for model calibration and niche comparisons
## area for model calibration (accessible area)
cal_area <- concave_area(data = mammoth_oc, longitude = "lng",
                         latitude = "lat", buffer_distance = 500)

## simple plot
x11()
plot(variables$pleistocene_hs1$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Calibration area (acccessible area)")
points(mammoth_oc[, 2:3], pch = 16, cex = 0.3)
plot(cal_area, border = "blue", lwd = 2, add = TRUE)


## areas relevant for each time period (niche comparison areas)
cal_area_787 <- concave_area(data = mammoth_oc_split$`0.787`, longitude = "lng",
                             latitude = "lat", buffer_distance = 500)

cal_area_15 <- concave_area(data = mammoth_oc_split$`0.01585`, longitude = "lng",
                            latitude = "lat", buffer_distance = 500)

## simple plots
x11()
plot(variables$pleistocene_mis19$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Calibration area ~787 mya")
points(mammoth_oc_split$`0.787`[, 2:3], pch = 16, cex = 0.3)
plot(cal_area_787, border = "blue", lwd = 2, add = TRUE)

x11()
plot(variables$pleistocene_hs1$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Calibration area ~0.015 mya")
points(mammoth_oc_split$`0.01585`[, 2:3], pch = 16, cex = 0.3)
plot(cal_area_15, border = "blue", lwd = 2, add = TRUE)


# masking raster layers to regions of interest
## to model calibration area
acces_mis19 <- mask(crop(variables$pleistocene_mis19, cal_area), cal_area)
acces_hs1 <- mask(crop(variables$pleistocene_hs1, cal_area), cal_area)

## to areas relevant for niche comparison
var_mis19 <- mask(crop(variables$pleistocene_mis19, cal_area_787), cal_area_787)
var_hs1 <- mask(crop(variables$pleistocene_hs1, cal_area_15), cal_area_15)


# writing layers in directory
## prepare directories and variable names 
dirs <- c("accessible_area_mis19", "accessible_area_hs1", "relevant_area_mis19",
          "relevant_area_hs1")

rnames <- lapply(dirs, function (x) {
  cal_dir <- paste0("Environmental_data/", x)
  dir.create(cal_dir)
  
  paste0(cal_dir, "/", vars_keep, ".tif")
}) 

## writing layers in loop
for (i in 1:length(vars_keep)) {
  writeRaster(acces_mis19[[i]], filename = rnames[[1]][i], format = "GTiff")
  writeRaster(acces_hs1[[i]], filename = rnames[[2]][i], format = "GTiff")
  writeRaster(var_mis19[[i]], filename = rnames[[3]][i], format = "GTiff")
  writeRaster(var_hs1[[i]], filename = rnames[[4]][i], format = "GTiff")
}
# ------------------------------------------------------------------------------


# More data cleaning and spatial thinning --------------------------------------
## data split by time slices 
mammoth_oc_split <- split(mammoth_oc, f = as.factor(mammoth_oc$time_slice))

# check records outside raster layers
## records 0.787 mya
outside_mis19 <- which(is.na(extract(var_mis19[[1]], 
                                     mammoth_oc_split$`0.787`[, 2:3])))

outside_mis19 # all OK

## records 0.015 mya
outside_hs1 <- which(is.na(extract(var_hs1[[1]], 
                                     mammoth_oc_split$`0.01585`[, 2:3])))

outside_hs1 # one is outside

## checking weird record
x11()
plot(variables$pleistocene_hs1$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Records with no data")
points(mammoth_oc_split$`0.01585`[outside_hs1, 2:3], pch = 16, cex = 1)

plot(variables$pleistocene_hs1$bio_1, xlim = c(-73, -72), ylim = c(39, 40),
     main = "Records with no data")
points(mammoth_oc_split$`0.01585`[outside_hs1, 2:3], pch = 16, cex = 1)

## correcting manually that record (if many records are outside, the process
## will be longer, but see ellipsenm::to_closest())
mammoth_oc_split$`0.01585`[outside_hs1, 2:3] <- c(-72.55, 39.55)

## plotting point again to check
points(mammoth_oc_split$`0.01585`[outside_hs1, 2:3], pch = 1, cex = 1)


# Spatial thinning (geographic rarefaction of records to reduce autocorrelation)
## records 0.787 mya
thin(mammoth_oc_split$`0.787`, lat.col = "lat", 
     long.col = "lng", spec.col = "taxon_name", 
     thin.par = 30, reps = 10, write.files = TRUE, 
     max.files = 1, out.dir = "Occurrence_data", 
     locs.thinned.list.return = FALSE, out.base = "oc_msi19")

## records 0.015 mya
thin(mammoth_oc_split$`0.01585`, lat.col = "lat", 
     long.col = "lng", spec.col = "taxon_name", 
     thin.par = 30, reps = 10, write.files = TRUE, 
     max.files = 1, out.dir = "Occurrence_data", 
     locs.thinned.list.return = FALSE, out.base = "oc_hs1")


# writing clean data
## all time slices before thinning
for (i in 1:length(mammoth_oc_split)) {
  slice <- names(mammoth_oc_split)[i]
  fname <- paste0("Occurrence_data/mammut_americanum_", slice, "mya.csv")
  write.csv(mammoth_oc_split[[i]], file = fname, row.names = FALSE)
}
# ------------------------------------------------------------------------------
