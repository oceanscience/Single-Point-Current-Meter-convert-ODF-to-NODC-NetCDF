# writing ODF to NetCDF for CIOOS - Microcat

# author Diana Cardoso
# Oct 2023-2026
# Fisheries and Oceans Canada,Bedford Institute of Oceanography, Dartmouth, N.S. Canada B2Y 4A2

rm(list=ls()) #remove everything in the working environment.
library(oce)
library(ncdf4)
library(stringr)
library(glue)

setwd("")
year <- 1983 #used as a folder name assuming data is organised in folders by year

# function file names
source("mcm_nc_edit.R") 
source("asP01.R")

file <- list.files(paste0('./',year), pattern = "MCM*...*.ODF", ignore.case = TRUE) # list of ODF files
gf32p01 <- read.csv('GF3 Code Map_final_v2.csv', header = TRUE) # GF3 to BODC mapping
#list of projects for given year, created using ODFmcm_title_program_project.R and manually editing
projectlist <- read.csv(paste0('./',year,'/mooringprojectsh.csv'), header = TRUE)
vessel_codes <- read.csv('Shared_models_vessel_BIO.csv', header = TRUE) # list of ship names and info

#  loop through each file, convert CRAT to conductivity, add program and project names, correct ship name
for ( f in file ){

obj <- read.odf(paste0('./',year,"/",f))

datanames <- names(obj@data)
unitnames <- names(obj@metadata$units)

#change the name of obj@metadata$header$ODF_HEADER$`FILE_SPECIFICATION _1` to obj@metadata$header$ODF_HEADER$FILE_SPECIFICATION_1
names(obj@metadata$header$CRUISE_HEADER) <- gsub(" ", "", names(obj@metadata$header$CRUISE_HEADER))
names(obj@metadata$header$ODF_HEADER)[names(obj@metadata$header$ODF_HEADER) == "FILE_SPECIFICATION _1"] <- "FILE_SPECIFICATION_1"
obj@metadata$header$CRUISE_HEADER$CRUISE_NAME <- obj@metadata$header$CRUISE_HEADER$CRUISE_NAME_1

# check if conductivity exists and is conductivity Ratio, if so, change name to conductivity ratio and add conductivity
if (exists("conductivity", where <- obj@metadata$dataNamesOriginal) && grepl("CRAT_01|CRAT|CRAT_1", obj@metadata$dataNamesOriginal$conductivity)){
#change data name of conductivity to conductivityratio and add conductivity
indexcond <- which(grepl("conductivity",datanames))
datanames <-sub("conductivity", "conductivityratio", datanames)
names(obj@data) <- c(datanames)
unitnames <-sub("conductivity", "conductivityratio", unitnames)
names(obj@metadata$units) <- c(unitnames)

#Calculate Electrical Conductivity(CNDC_01): c0S/m: conductivity [S/m] from
#conductivity ratio (S.Data.CNDC_01*10.0)./sw_c3515)
#add conductivity to obj
#calculate conductivity, 42.914 mS/cm is conductivity at a standard of 35 salinity, 15 °C temperature, 0 dbar pressure
#10 converts mS/cm to S/m

obj@data[[length(obj@data)+1]] <- (obj@data$conductivityratio*42.914)/10
names(obj@data)[[length(obj@data)]] <- c("conductivity")
obj@metadata$units[[length(obj@metadata$units)+1]] <- "xxx"
names(obj@metadata$units)[[length(obj@metadata$units)]] <- c("conductivity")

datanames <- names(obj@data)
#add units of conductivity

obj@metadata$units$conductivity <- list("unit", "scale")
names(obj@metadata$units$conductivity) <- c("unit", "scale")
obj@metadata$units$conductivity$unit <- "S/m"
obj@metadata$units$conductivityratio$unit <- "dimensionless"
#add gf3 name of conductivity
obj@metadata$dataNamesOriginal[[length(obj@metadata$dataNamesOriginal)+1]]<- "CNDC_01"
names(obj@metadata$dataNamesOriginal)[[length(obj@metadata$dataNamesOriginal)]] <- c("conductivity")
unitnames <- names(obj@metadata$units)

}

#find obj@metadata$dataNamesOriginal that contain noWMOcode in the name
indexnocode <- grep(pattern = "noWMOcode", names(obj@metadata$dataNamesOriginal), value = FALSE)
indexparcode <- grep(pattern = "PARAMETER_HEADER", names(obj@metadata$header), value = FALSE)
indexdatacode <- grep(pattern = "noWMOcode", names(obj@data), value = FALSE)
indexunitcode <- grep(pattern = "noWMOcode", names(obj@metadata$units), value = FALSE)

for (i in indexparcode){
  obj@metadata$header[[i]][[2]] <- trimws(gsub("[\'\"]", "",obj@metadata$header[[i]][[2]]))
  if (grepl("SREF", obj@metadata$header[[i]][[2]], ignore.case = TRUE)){
    SREF_index <- i
  }else if (grepl("DELTA-S", obj@metadata$header[[i]][[2]], ignore.case = TRUE)){
    DELTAS_index <- i     
  }else if (grepl("CCR", obj@metadata$header[[i]][[2]],ignore.case = TRUE)){
    CCR_index <- i
  }else {
    SREF_index <- NA
    DELTAS_index <- NA
    CCR_index <-  NA
  }
}

if (length(indexnocode) > 0){ 
  if (length(obj@metadata$units) > length(obj@data)){
  for (i in indexparcode){
    obj@metadata$header[[i]][[2]] <- trimws(gsub("[\'\"]", "",obj@metadata$header[[i]][[2]]))
    if (grepl("SREF", obj@metadata$header[[i]][[2]], ignore.case = TRUE)){
        obj@metadata$dataNamesOriginal[[which(SREF_index==indexparcode)-1]] <- "SREF"
        obj@metadata$units[[which(SREF_index==indexparcode)]]$unit <- "psu"
        names(obj@data)[[which(SREF_index==indexparcode)-1]]<- c("salinityreference")
        names(obj@metadata$dataNamesOriginal)[[which(SREF_index==indexparcode)-1]]<- c("salinityreference")
        names(obj@metadata$units)[[which(SREF_index==indexparcode)]]<- c("salinityreference")
    }else if (grepl("DELTA-S", obj@metadata$header[[i]][[2]], ignore.case = TRUE)){
      obj@metadata$dataNamesOriginal[[which(DELTAS_index==indexparcode)-1]] <- "DELTA-S"
      obj@metadata$units[[which(DELTAS_index==indexparcode)]]$unit <- "g kg-1"
      names(obj@data)[[which(DELTAS_index==indexparcode)-1]]<- c("salinitycorrection")
      names(obj@metadata$dataNamesOriginal)[[which(DELTAS_index==indexparcode)-1]]<- c("salinitycorrection")
      names(obj@metadata$units)[[which(DELTAS_index==indexparcode)]]<- c("salinitycorrection")      
    }else if (grepl("CCR", obj@metadata$header[[i]][[2]],ignore.case = TRUE)){
      obj@metadata$dataNamesOriginal[[which(CCR_index==indexparcode)-1]] <- "CCR"
      obj@metadata$units[[which(CCR_index==indexparcode)]]$unit <- "index dimensionless"
      names(obj@data)[[which(CCR_index==indexparcode)-1]]<- c("relcorrrelationConc")
      names(obj@metadata$dataNamesOriginal)[[which(CCR_index==indexparcode)-1]]<- c("relcorrrelationConc")
      names(obj@metadata$units)[[which(CCR_index==indexparcode)]]<- c("relcorrrelationConc")       
    }
  }
    
  } else {
    for (i in indexparcode){
      obj@metadata$header[[i]][[2]] <- trimws(gsub("[\'\"]", "",obj@metadata$header[[i]][[2]]))
      if (grepl("SREF", obj@metadata$header[[i]][[2]], ignore.case = TRUE)){
        obj@metadata$dataNamesOriginal[[which(SREF_index==indexparcode)]] <- "SREF"
        obj@metadata$units[[which(SREF_index==indexparcode)]]$unit <- "psu"
        names(obj@data)[[which(SREF_index==indexparcode)]]<- c("salinityreference")
        names(obj@metadata$dataNamesOriginal)[[which(SREF_index==indexparcode)]]<- c("salinityreference")
        names(obj@metadata$units)[[which(SREF_index==indexparcode)]]<- c("salinityreference")
      }else if (grepl("DELTA-S", obj@metadata$header[[i]][[2]], ignore.case = TRUE)){
        obj@metadata$dataNamesOriginal[[which(DELTAS_index==indexparcode)]] <- "DELTA-S"
        obj@metadata$units[[which(DELTAS_index==indexparcode)]]$unit <- "g kg-1"
        names(obj@data)[[which(DELTAS_index==indexparcode)]]<- c("salinitycorrection")
        names(obj@metadata$dataNamesOriginal)[[which(DELTAS_index==indexparcode)]]<- c("salinitycorrection")
        names(obj@metadata$units)[[which(DELTAS_index==indexparcode)]]<- c("salinitycorrection")      
      }else if (grepl("CCR", obj@metadata$header[[i]][[2]],ignore.case = TRUE)){
        obj@metadata$dataNamesOriginal[[which(CCR_index==indexparcode)]] <- "CCR"
        obj@metadata$units[[which(CCR_index==indexparcode)]]$unit <- "index dimensionless"
        names(obj@data)[[which(CCR_index==indexparcode)]]<- c("relcorrrelationConc")
        names(obj@metadata$dataNamesOriginal)[[which(CCR_index==indexparcode)]]<- c("relcorrrelationConc")
        names(obj@metadata$units)[[which(CCR_index==indexparcode)]]<- c("relcorrrelationConc")
}
}
  }
}
 
#get project, program and location description
index_projectlist <- grep(f, projectlist$File_name, value = FALSE)
obj@metadata$header$CRUISE_HEADER$CRUISE_DESCRIPTION <- projectlist$CRUISE_DESCRIPTION[index_projectlist]
obj@metadata$header$CRUISE_HEADER$CRUISE_NAME <- projectlist$CRUISE_NAME[index_projectlist]
project <- projectlist$project[index_projectlist]
program <- projectlist$program[index_projectlist]
location_description <- projectlist$location_description[index_projectlist]
cruise_name <- paste0("mooring deployment",", ",projectlist$CRUISE_NAME[index_projectlist]) #projectlist$CRUISE_NAME2[index_projectlist]
info_url <- projectlist$url[index_projectlist]

#get vessel ices code and NERC url
if (grepl("HUD", f)){
  obj@metadata$header$CRUISE_HEADER$PLATFORM <- 'hudson'
}

obj@metadata$header$CRUISE_HEADER$PLATFORM <- trimws(gsub("[\'\"]", "", obj@metadata$header$CRUISE_HEADER$PLATFORM))
obj@metadata$header$CRUISE_HEADER$PLATFORM<-gsub("'","",obj@metadata$header$CRUISE_HEADER$PLATFORM)
obj@metadata$header$CRUISE_HEADER$PLATFORM<-gsub("[()]","",obj@metadata$header$CRUISE_HEADER$PLATFORM)

if (nzchar(obj@metadata$header$CRUISE_HEADER$PLATFORM) && any(grepl(obj@metadata$header$CRUISE_HEADER$PLATFORM, vessel_codes$name, ignore.case = TRUE))){   
index_vessel_codes <- grep(obj@metadata$header$CRUISE_HEADER$PLATFORM, vessel_codes$name, value = FALSE,
                           ignore.case = TRUE, fixed = FALSE)
vessel_icescode <- vessel_codes$ICES_SHIPC_ship_codes[index_vessel_codes]
vessel_NERC <- vessel_codes$BODC.Vocabulary[index_vessel_codes]
print(obj@metadata$header$CRUISE_HEADER$PLATFORM)
}else {
  vessel_icescode <- NA
  vessel_NERC <- NA
  print("No vessel name found in file name or cruise header")
}
print(f)
obj[['institute']] <- 'DFO BIO'

mcm_nc(obj,location_description,info_url,project, program, cruise_name,vessel_icescode,vessel_NERC)

}
