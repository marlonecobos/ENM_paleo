# ------------------------------------------------------------------------------
# Course: ECOLOGICAL MODELS APPLIED TO FOSIL DATA: DATA PREPARATION
# Author: Marlon E. Cobos, Hannah L. Owens
# Date modified: 05/10/2022
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
# the next lines of code load packages 

## if needed, install packages (remotes, paleobioDB, maps, spThin, raster) using 
## install.packages("package_name")
## ellipsenm is on GitHub so install it as 
## remotes::install_github("marlonecobos/ellipsenm") 
## https://github.com/marlonecobos/ellipsenm/#installing-the-package

library(paleobioDB)
library(maps)
library(ellipsenm)
library(spThin)
library(raster)
# ------------------------------------------------------------------------------


# Working directory and sub-directories -----------------------------------------
# working directory
setwd("YOUR/DIRECTORY")

## sub-directory for occurrence data
if (!dir.exists("Occurrence_data")) {
  dir.create("Occurrence_data")
}

## sub-directory for environmental data
if (!dir.exists("Environmental_data")) {
  dir.create("Environmental_data")
}
# ------------------------------------------------------------------------------


# Getting occurrence data ------------------------------------------------------
# search for occurrences
mastodon_oc <- pbdb_occurrences (limit = "all", base_name = "Mammut americanum", 
                                 vocab = "pbdb", interval = "Quaternary",
                                 show = c("coords", "phylo", "ident"))

# write raw downloaded data to directory
fname_raw <- "Occurrence_data/mammut_americanum_raw.csv"
if (!file.exists(fname_raw)) {
  write.csv(mastodon_oc, file = fname_raw, row.names = FALSE)  
}

# info in downloaded data
colnames(mastodon_oc)

# filter records for unique occurrences and relevant columns
## info to keep 
selected <- c("taxon_name", "lng", "lat", "early_age", "late_age") 

## filter data to keep selected columns
mastodon_oc <- mastodon_oc[, selected]
# ------------------------------------------------------------------------------


# Occurrence data cleaning and filtering ---------------------------------------
# filter to keep only unique events in the data
mastodon_oc <- unique(mastodon_oc)

# check what kind of data we have
str(mastodon_oc)

# some columns that should contain numeric info are in character format 
## this fixes the issue of all columns with character info
numeric_info <- as.numeric(unlist(mastodon_oc[, c("lng", "lat", "early_age", 
                                                  "late_age")]))
## put corrected info back within the main table
mastodon_oc[, c("lng", "lat", "early_age", "late_age")] <- numeric_info


# check taxon name 
## all names in the column of species name
table(mastodon_oc$taxon_name)

## “cf” basically means “I’m not sure, the preservation isn’t that great"
## exclude that one
mastodon_oc <- mastodon_oc[mastodon_oc$taxon_name != "Mammut cf. americanum", ]

## all other names are synonyms 
## change all to "Mammut_americanum"
mastodon_oc$taxon_name <- "Mammut_americanum"
table(mastodon_oc$taxon_name)


# simple plot to check the data using environmental layers and records
## limits for plot
lims <- apply(mastodon_oc[, c("lng","lat")], 2, range)
lims

## making limits a little wider
lims <- lims * matrix(c(1.1, 0.9, 0.9, 1.1), nrow = 2)
lims

## base map 
map(xlim = lims[, 1], ylim = lims[, 2], fill = TRUE, col = "gray75")

## the records
points(mastodon_oc[, c("lng","lat")], pch = 20, col = "black") 

title(main = "All Mastodon Points")
box()

# split points into relevant time slices
## estimated age of records
mastodon_oc$estimated_age <- (mastodon_oc$early_age + mastodon_oc$late_age) / 2

## date ranges from paleoclim.org
time_slices <- data.frame (maxTime = c(3.3, 3.205, .787, 0.130, 0.017, 0.0129, 
                                       0.0117, 0.00826, 0.0042),
                           minTime = c(3.3, 3.205, .787, 0.130, 0.0147, 0.0117, 
                                       0.00826, 0.0042, 0.0003))

time_slices$meanAge <- (time_slices$maxTime + time_slices$minTime) / 2

## assign time slices to occurrences based on mean age
for (i in 1:nrow(mastodon_oc)) {
  time_diff <- abs(time_slices$meanAge - mastodon_oc[i, "estimated_age"])
  mastodon_oc$time_index[i] <- which.min(time_diff)
  mastodon_oc$time_slice[i] <- time_slices$meanAge[mastodon_oc$time_index[i]]
}


# plot to visualize data
## define colors
indexes <- unique(mastodon_oc$time_index)
ncolors <- length(indexes)
col_pal <- colorRampPalette(c('red','yellow', 'blue'))
mastodon_oc$col <- col_pal(ncolors)[as.factor(mastodon_oc$time_index)]

## plot with legend
map(xlim = lims[, 1], ylim = lims[, 2], fill = TRUE, col = "gray75")
points(mastodon_oc[, c("lng","lat")], col = mastodon_oc$col, pch = 20)

title(main = "All Mastodon Points by Date")
legend("bottomleft", title = "Date (My)", col = col_pal(ncolors), pch = 20,
       legend = time_slices$meanAge[sort(indexes)], bt = "n", cex = 0.8)
box()

# write data to directory
fname_all <- "Occurrence_data/mammut_americanum_all.csv"
if (!file.exists(fname_all)) {
  write.csv(mastodon_oc, file = fname_all, row.names = FALSE)  
}


# data split by time slices 
mastodon_oc_split <- split(mastodon_oc, f = as.factor(mastodon_oc$time_slice))
# ------------------------------------------------------------------------------


# Downloading environmental data -----------------------------------------------
# get data from paleoclim.org
## lets check how many records per time slice 
table(mastodon_oc$time_slice)

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
  
  if(!file.exists(dfile)){ # Checks to see if data is already present before download
    download.file(periods[x], destfile = dfile, method = "auto", # download
                  quiet = FALSE, mode = "wb", cacheOK = TRUE)
  }
  
  dfol <- paste0("Environmental_data/", pname, "_bio") # folder to unzip data
  if (!dir.exists(dfol)) {
    dir.create(dfol)
  }
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
dev.off(); par(mfrow = c(2, 1))

plot(variables$pleistocene_mis19$bio_1, 
     main = "Annual Mean Temperature (C * 10),\nPleistocene MIS, ca. 787 ka") 
plot(variables$pleistocene_hs1$bio_1,
     main = "Annual Mean Temperature (C * 10),\nPleistocene Heinrich Stadial 1 (17.0-14.7 ka)")
# ------------------------------------------------------------------------------


# Processing environmental data ------------------------------------------------
# keep variables that match between periods
## match variable names
match_var <- match(names(variables$pleistocene_mis19), 
                   names(variables$pleistocene_hs1))

## variables to keep and sort names
vars_keep <- names(variables$pleistocene_hs1)[match_var]
vars_keep

vars_keep <- vars_keep[c(1, 12:14, 2:11)] # Sorts numerically instead of alphabetically
vars_keep

## filtering variables in stacks
variables$pleistocene_mis19 <- variables$pleistocene_mis19[[vars_keep]] 
variables$pleistocene_hs1 <- variables$pleistocene_hs1[[vars_keep]] 


# region of interest for model calibration and niche comparisons
## area for model calibration (accessible area)
cal_area <- concave_area(data = mastodon_oc, longitude = "lng",
                         latitude = "lat", buffer_distance = 500)

## simple plot
dev.off()
plot(variables$pleistocene_hs1$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Calibration area (acccessible area)")
points(mastodon_oc[, 2:3], pch = 16, cex = 0.3)
plot(cal_area, border = "blue", lwd = 2, add = TRUE)


## areas relevant for each time period (niche comparison areas)
cal_area_787 <- concave_area(data = mastodon_oc_split$`0.787`, longitude = "lng",
                             latitude = "lat", buffer_distance = 500)

cal_area_15 <- concave_area(data = mastodon_oc_split$`0.01585`, longitude = "lng",
                            latitude = "lat", buffer_distance = 500)

## simple plots
plot(variables$pleistocene_mis19$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Calibration area ~787 mya")
points(mastodon_oc_split$`0.787`[, 2:3], pch = 16, cex = 0.3)
plot(cal_area_787, border = "blue", lwd = 2, add = TRUE)


plot(variables$pleistocene_hs1$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Calibration area ~0.015 mya")
points(mastodon_oc_split$`0.01585`[, 2:3], pch = 16, cex = 0.3)
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
  if (!dir.exists(cal_dir)) {
    dir.create(cal_dir)
  }
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
# check records outside raster layers
## records 0.787 mya
outside_mis19 <- which(is.na(extract(var_mis19[[1]], 
                                     mastodon_oc_split$`0.787`[, 2:3])))

outside_mis19 # all OK

## records 0.015 mya
outside_hs1 <- which(is.na(extract(var_hs1[[1]], 
                                   mastodon_oc_split$`0.01585`[, 2:3])))

outside_hs1 # one is outside

## checking weird record
plot(variables$pleistocene_hs1$bio_1, xlim = lims[, 1], ylim = lims[, 2],
     main = "Records with no data")
points(mastodon_oc_split$`0.01585`[outside_hs1, 2:3], pch = 16, cex = 1)

plot(variables$pleistocene_hs1$bio_1, xlim = c(-73, -72), ylim = c(39, 40),
     main = "Records with no data")
points(mastodon_oc_split$`0.01585`[outside_hs1, 2:3], pch = 16, cex = 1)

## correcting manually that record (if many records are outside, the process
## will be longer, but see ellipsenm::to_closest())
mastodon_oc_split$`0.01585`[outside_hs1, 2:3] <- c(-72.55, 39.55)

## plotting point again to check
points(mastodon_oc_split$`0.01585`[outside_hs1, 2:3], pch = 1, cex = 1)


# Spatial thinning (geographic rarefaction of records to reduce autocorrelation)
## records 0.787 mya
thin(mastodon_oc_split$`0.787`, lat.col = "lat", 
     long.col = "lng", spec.col = "taxon_name", 
     thin.par = 30, reps = 10, write.files = TRUE, 
     max.files = 1, out.dir = "Occurrence_data", 
     locs.thinned.list.return = FALSE, out.base = "oc_msi19")

## records 0.015 mya
thin(mastodon_oc_split$`0.01585`, lat.col = "lat", 
     long.col = "lng", spec.col = "taxon_name", 
     thin.par = 30, reps = 10, write.files = TRUE, 
     max.files = 1, out.dir = "Occurrence_data", 
     locs.thinned.list.return = FALSE, out.base = "oc_hs1")


# writing clean data
## all time slices before thinning
for (i in 1:length(mastodon_oc_split)) {
  slice <- names(mastodon_oc_split)[i]
  fname <- paste0("Occurrence_data/mammut_americanum_", slice, "mya.csv")
  
  if (!file.exists(fname)) {
    write.csv(mastodon_oc_split[[i]], file = fname, row.names = FALSE)
  }
}
# ------------------------------------------------------------------------------