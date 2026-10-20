# write a csv file with CRUISE_DESCRIPTION and CRUISE_NAME  from the ODF Cruise Header
# Review and edit this file with correct info before writing ODF to NetCDF

# author Diana Cardoso
# Oct 2023-2026
# Fisheries and Oceans Canada,Bedford Institute of Oceanography, Dartmouth, N.S. Canada B2Y 4A2

rm(list=ls()) #remove everything in the working environment.
library(oce)

setwd("")

file <- list.files('.', pattern = "MCM*...*.ODF", ignore.case = TRUE)

for ( f in file ){
obj <- read.odf(f)
print(f)
cruisedetials1 <- obj@metadata$header$CRUISE_HEADER$CRUISE_DESCRIPTION
cruisedetials2 <- obj@metadata$header$CRUISE_HEADER$CRUISE_NAME
#cruisedetials2 <- obj@metadata$header$EVENT_HEADER$EVENT_COMMENTS_1
mp <- data.frame(cruisedetials1, cruisedetials2,f)
write.table(mp, file = "mooringprojects.csv", append = TRUE, sep = ",",
            row.names = FALSE,
            col.names = FALSE)
}

mpall <- read.csv('mooringprojects.csv', header = FALSE)
names(mpall) <- c("CRUISE_DESCRIPTION","CRUISE_NAME","File_name")
write.table(mpall, file = "mooringprojectsh.csv", append = FALSE, sep = ",",
            row.names = FALSE,col.names = TRUE)

