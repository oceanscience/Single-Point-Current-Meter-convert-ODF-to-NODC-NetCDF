# write a csv file with list of ships, instrument type, model and serial numbesr from the ODF Cruise Header
# Review and edit this file with correct info before writing ODF to NetCDF
# create a plot of the mooring locations on a map with the cruise number as labels

# author Diana Cardoso
# Oct 2023-2026
# Fisheries and Oceans Canada,Bedford Institute of Oceanography, Dartmouth, N.S. Canada B2Y 4A2

library(oce)
library(ocedata)
library(stringr)
library(glue)
data("coastlineWorldFine")

setwd("")

file <- list.files('.', pattern = "MCM*...*.ODF", ignore.case = TRUE)

b=1
for ( f in file ){
obj <- read.odf(f) #read ODF file with oce package
#print(f)
names(obj@metadata$header$CRUISE_HEADER) <- gsub(" ", "", names(obj@metadata$header$CRUISE_HEADER))

cruisedetials1 <- obj@metadata$header$CRUISE_HEADER$PLATFORM
INST_TYPE <- obj@metadata$header$INSTRUMENT_HEADER$INST_TYPE_1
DESCRIPTION <- obj@metadata$header$INSTRUMENT_HEADER$DESCRIPTION_1
if (is.null(DESCRIPTION)){
  DESCRIPTION <- NA
}
SERIAL_NUMBER <- obj@metadata$header$INSTRUMENT_HEADER$SERIAL_NUMBER_1
MODEL <- obj@metadata$header$INSTRUMENT_HEADER$MODEL_1
LAT <- obj@metadata$header$EVENT_HEADER$INITIAL_LATITUDE_1
LON <- obj@metadata$header$EVENT_HEADER$INITIAL_LONGITUDE_1
cruisenumber <- trim(obj@metadata$header$CRUISE_HEADER$CRUISE_NUMBER_1)
#remove quotes from cruisenumber and trim whitespace
cruisenumber <- trimws(gsub("[\'\"]", "", cruisenumber))

mp <- data.frame(cruisedetials1,f,INST_TYPE,DESCRIPTION,SERIAL_NUMBER,MODEL,LAT,LON,cruisenumber)

if (b==1){
#add column names to mp data frame
names(mp) <- c("PLATFORM","File_name","INST_TYPE","DESCRIPTION","SERIAL
_NUMBER","MODEL","LAT","LON","CRUISE_NUMBER")

#write mp data frame to csv file with column names, if b=1, otherwise append to csv file without column names
write.table(mp, file = "mooringships4.csv", append = TRUE, sep = ",",
            row.names = FALSE,
            col.names = TRUE)
}else {
  names(mp) <- c("PLATFORM","File_name","INST_TYPE","DESCRIPTION","SERIAL
_NUMBER","MODEL","LAT","LON","CRUISE_NUMBER")
  write.table(mp, file = "mooringships.csv", append = TRUE, sep = ",",
              row.names = FALSE,
              col.names = FALSE)
}

b=b+1
}
#read the csv file and plot the locations of the moorings on a map, with the cruise number as labels
mpall <- read.csv('mooringships.csv', header = TRUE)

min_lat <- (min(as.numeric(mpall$LAT), na.rm = TRUE)-1) 
min_lon <- (max(as.numeric(mpall$LON), na.rm = TRUE)+1)
max_lat <- (max(as.numeric(mpall$LAT), na.rm = TRUE)+1)
max_lon <- (min(as.numeric(mpall$LON), na.rm = TRUE)-1)

lonlim <- c(max_lon, min_lon)
latlim <- c(min_lat, max_lat)

#write a loop to replace "LON" and "LAT" with na values in the mpall$LON and mpall$LAT columns, if the values are not numeric
for (i in 1:nrow(mpall)){
  if (!is.numeric(mpall$LON[i])){
    mpall$LON[i] <- -30
    print(paste("Replacing LON value with -30 for row", i))
  }
  if (!is.numeric(mpall$LAT[i])){
    mpall$LAT[i] <- 60
    print(paste("Replacing LAT value with 60 for row", i))
  }
}

jpeg("mooring_map.jpeg", width = 1150, height = 750)
mapPlot(coastlineWorldFine, 
        projection="+proj=merc",
        col="lightgray", 
        longitudelim=lonlim, 
        latitudelim=latlim, lonlabels = TRUE,latlabels = TRUE,  )


mapPoints(longitude = mpall$LON,
          latitude = mpall$LAT,
          pch = 1,cex = 1.5, col="blue")

mapText(longitude = mpall$LON,
        latitude = mpall$LAT,mpall$CRUISE_NUMBER, cex = 1, col = "brown", pos=1)

mtext("MCM Mooring locations")

dev.off()
graphics.off()
