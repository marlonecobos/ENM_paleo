library(paleobioDB) # Access paleobiodb.org from R
library(dplyr) # For data cleaning
library(rangeBuilder) # Has a convenient continent shapefile for quick, easy mapping

# Search for occurrences
mammothOccurrences <-  pbdb_occurrences (limit="all",
                                         base_name="Mammut americanum", vocab="pbdb",
                                         interval="Quaternary",
                                         show=c("coords", "phylo", "ident"))

# Filter records for Unique occurrences and relevant columns
mammothOccurrencesClean <- mammothOccurrences %>% select("taxon_name", "lng", "lat", "early_age", "late_age") %>% unique()
mammothOccurrencesClean[,c("lng", "lat", "early_age", "late_age")] <- as.numeric(unlist(mammothOccurrencesClean[,c("lng", "lat", 
                                                                                                                   "early_age", 
                                                                                                                   "late_age")]))

# Quick plot to see what's there
plot(mammothOccurrencesClean[,c("lng","lat")], pch = 20, col = "red", main = "All Mammoth Points", 
     xlab = "Longitude", ylab = "Latitude")
plot(gshhs, add = T)

# Split points into relevant time slices ----
mammothOccurrencesClean$estimated_age <- (mammothOccurrencesClean$early_age + mammothOccurrencesClean$late_age)/2

# Assign time slices to occurrences based on date ranges from paleoclim.org
relevantTimeSlices <- data.frame (maxTime = c(3.3, 3.205, .787, 0.130, 0.017, 0.0129, 0.0117, 0.00826, 0.0042), 
                                  minTime = c(3.3, 3.205, .787, 0.130, 0.0147, 0.0117, 0.00826, 0.0042, 0.0003))
relevantTimeSlices$meanAge <- (relevantTimeSlices$maxTime + relevantTimeSlices$minTime)/2

for(i in 1:nrow(mammothOccurrencesClean)){
  mammothOccurrencesClean$timeIndex[i] <- which.min(abs(relevantTimeSlices$meanAge - mammothOccurrencesClean[i,"estimated_age"]))
  mammothOccurrencesClean$timeSlice[i] <- relevantTimeSlices$meanAge[mammothOccurrencesClean$timeIndex[i]]
}

# How do the data look now?
ncolors <- length(unique(mammothOccurrencesClean$timeIndex))
rbPal <- colorRampPalette(c("#FFFFD4", "#FED98E", "#FE9929","#CC4C02"))
mammothOccurrencesClean$Col <- rbPal(ncolors)[mammothOccurrencesClean$timeIndex]

plot(mammothOccurrencesClean[,c("lng","lat")], main = "All Mammoth Points by Date", 
     xlab = "Longitude", ylab = "Latitude", col = mammothOccurrencesClean$Col, pch = ".", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(gshhs, col = "gray", add = T)
points(mammothOccurrencesClean[,c("lng","lat")], bg = mammothOccurrencesClean$Col, 
     pch = 21, cex = 1.5)
legend("topleft",title="Date (My)", legend=relevantTimeSlices$meanAge[sort(unique(mammothOccurrencesClean$timeIndex))],
       fill = rbPal(ncolors), cex = 1.5, bg = "white")

for (i in 1:length(relevantTimeSlices$meanAge)){
  tmp <- mammothOccurrencesClean[mammothOccurrencesClean$timeSlice == relevantTimeSlices$meanAge[i],-8:-9]
  if(nrow(tmp) > 0){
    write.csv(tmp, file = paste0("Exercises/SplitOccurrences/mammutAmericanum_", tmp$timeSlice[1], "mya.csv"))
  }
}
