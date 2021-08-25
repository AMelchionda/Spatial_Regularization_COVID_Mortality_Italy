###########################################################
##################### COMPLETE DATASET ####################
###########################################################

#### Load the data -----------------------------------------------------------

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/01 Buildings/files', sep = ''))
info <- read.csv('info.txt', encoding = "latin1")
Complete <- read.csv('Complete.txt')

setwd('../../')

#### Add missing comune names ---------------

names <- info[c('pro_com', 'comune', 'cod_reg')]
colnames(names) <- c("COMUNE", "NAME", "COD_REG")
definitive <- merge(names, Complete, by = "COMUNE")
colnames(definitive) <- c("id", "name", "cod_reg", "y", "x", "pop", "year", 1:366)

#### Convert municipalities coordinates to human numbers ---------------------
# Bless this incredible french site which given two geographical coordinates
# shows where they point in the map according to any system possible:
# https://app.dogeo.fr/Projection/#/coords-to-points
#
# I tried it with Aglie and Milano, and there are three possible choices
# that are reasonable, all in the EPSG format.
# I choose (randomly) to consider our coordinates to be in EPSG:23032 format
# and proceed to convert them in human.
#
# NOTE: I'm using the old dataset for this test, but as long as the names
#       on the new datasets are unchainged this code should still work.

library(rgdal)

# Extract the coordinates in a new dataset
orig_coords <- data.frame(lat = definitive$x, lon = definitive$y)
coordinates(orig_coords) <- c('lat', 'lon')

# Determine the projection of the lat-long coordinates, by default it is EPSG:4326
proj4string(orig_coords) <- CRS("+init=epsg:23032")
print(summary(orig_coords))

# Convert the coordinates to the used metric system (EPSG:4326)
new_coords<-spTransform(orig_coords,CRS("+init=epsg:4326"))
print(summary(new_coords))    

# Convert back the geo-object to a data frame
internal_nodes <- as.data.frame(new_coords)
colnames(internal_nodes) <- c('x', 'y')

rm(orig_coords, new_coords)

#### Replace into the final dataframe --------

definitive$x <- internal_nodes$x
definitive$y <- internal_nodes$y

#### Extract the final dataset -------------

setwd('./Resources')
write.csv(definitive, file = 'raw_complete.csv', row.names = FALSE)
