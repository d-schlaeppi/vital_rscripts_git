rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Zone Cloner ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####
# compatible with Myrmidon 0.8.3


#### Libraries ####

library(FortMyrmidon)


#### Directories ####
USER <- "2A13_Office_Daniel"  # Replace with the desired USER option: Nath_office, 2A13_Office_Adriano, 2A13_Office_Daniel, AEL-laptop
HD <- "Nathalie" # alternative values "Daniel"
setUserAndHD <- function(USER, HD) {
  usr <- NULL  # Default value in case of an unrecognized USER option
  if (USER == "Nath_office") {
    usr <- "bzniks"
  } else if (USER == "2A13_Office_Adriano") {
    usr <- "cf19810"
  } else if (USER == "2A13_Office_Daniel") {
    usr <- "gw20248"
  } else if (USER == "AEL-laptop") {
    usr <- "ael"
  }
  if (!is.null(usr)) {print(usr)} else {print("define new user if necessary")}
  assign("usr", usr, envir = .GlobalEnv)
  hd <- NULL
  if (HD == "Nathalie") {
    hd <- "/DISK_B"
  } else if (HD == "Daniel") {
    hd <- "/gismo_hd5"
  }
  if (!is.null(hd)) {print(hd)} else {print("define new hd if necessary")}
  assign("hd", hd, envir = .GlobalEnv)  # Assign hd to the global environment
}
setUserAndHD(USER, HD)
DATADIR <- paste("/media/",usr, hd, "/vital/fc2",sep="")

#### Resources ####

#select the source file that contains the zones
source_file_zones  <- paste(DATADIR,"/f_zoned_c01.myrmidon",sep="")

# create a destination file list with myrmidon files that shall get the zones defined in the source file
dest_files_zones <- list.files(path=DATADIR, pattern="f_AntsCreated",full.names=T, recursive = T)

#### Extract zones ####
# extract zone information from source file and safe it in a list and get the 
zoned <- fmExperimentOpen(source_file_zones) 
zone_list <- list()
for (space in zoned$spaces){
  zone_list <- c(zone_list,space$zones)
}

# get information for each zone 
for (zone in zone_list){
  ID   <- zone$ID 
  name <- zone$name
  start <- zone$definitions[[1]]$start
  end   <- zone$definitions[[1]]$end
  shapes <- zone$definitions[[1]]$shapes
}


#### Clone zones ####
# clone the zones into their destination files.
for (unzoned_file in dest_files_zones){
    print(unzoned_file)
    unzoned <- fmExperimentOpen(unzoned_file) # open destination myrmidon file  
    for (space in unzoned$spaces){# check that there are no zones, or delete existing zones 
      if (length(space$zones)>0){
        for (zone in space$zones){
          space$deleteZone(zone$ID)
        }
      }
    }
    for (space in unzoned$spaces){ #duplicate the desired zones 
      for (zone in zone_list){
        space$createZone(zone$name) # create zone
        new_zone <- space$zones[[length(space$zones)]] # load new zone into object
        new_zone$addDefinition(zone$definitions[[1]]$shapes ,   zone$definitions[[1]]$start , zone$definitions[[1]]$end )# add zone definitions
      }
    }
  unzoned$save(paste0((DATADIR), "/f_zoned", substr(basename(unzoned_file), 14, nchar(basename(unzoned_file))-0))) 
  rm(list=c("unzoned"))
}

#### To do next ####
# Open each of the destination files and manually adjust position and size if necessary. 


