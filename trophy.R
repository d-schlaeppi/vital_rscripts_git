rm(list=ls())
### ### ### ###
### Trophy  ###
### ### ### ###

#### Background information | Read me ####
# this script will contain all steps required to automate the detection of trophallacting interactions in tracking data


# behaviour annotation will be done on the experimental data FC1 aka the first tracking experiment as well as a the curated trophallxis video sets
# first we set up the trophy myrmidon files for annotations. 

#### 1. prerequisites ####

# load necessary libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()

# set working directory (where the tracking data is)
directory <- '/media/gw20248/gismo_hd5/trophy_data/'
setwd(directory)
list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE)

# an empty myrmidon file "base_source.myrmidon" was copied manually into the directory
# get all folders in the directory and compile them as a list (only folders containing the tracking data)
data_list <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
data_list <- grep(data_list, pattern = "trophy", invert = FALSE, value = TRUE)
print(data_list)

data_collection <- NULL 
number_of_colonies <- 2
colony_identifiers <- c("trophy_01", "trophy_02")
tracking_system <- c("troj", "west")

for(i in 1:number_of_colonies) {
  # collect variables
  nr                    <- i
  colony_id             <- colony_identifiers[i]
  tracking_system_main  <- tracking_system[i]
  # combine variables to a data frame  
  data_collection <-  rbind(data_collection, data.frame(nr, 
                                                        colony_id,
                                                        tracking_system_main,
                                                        stringsAsFactors = F))
}
data_collection

#### create base myrmidon files and create ants ####
for (i in 1:nrow(data_collection)) {
      file_name <- paste0(data_collection[i,"colony_id"], '_base.myrmidon')
      tracking_data <- fmExperimentOpen("base_source.myrmidon")
      s <- tracking_data$createSpace(data_collection[i,"tracking_system_main"])
      printf("Space '%s' has ID: %d\n",s$name,s$ID)
      tracking_data$name <- paste0(data_collection[i,"colony_id"])
      # assign the tracking data
      for (folder_name in data_list) {
        if (grepl(data_collection[i,"tracking_system_main"], folder_name, fixed = TRUE) ) {
          tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name), fixCorruptedData = TRUE)
        } else {next}
      }
      tracking_data$save(paste0(directory, data_collection[i,"colony_id"], '.myrmidon'))
      print(paste0(data_collection[i,"colony_id"], ' base file complete'))
      data <- fmExperimentOpen(paste0(directory, data_collection[i,"colony_id"], '.myrmidon'))
      tag_statistics <- fmQueryComputeTagStatistics(data)
        for (j in 1:nrow(tag_statistics)) {  # loop over each tag
            if (tag_statistics[j,"count"] >= 0.01*max(tag_statistics[,"count"],na.rm=T) ) { # optional: create an antID only if the tag detection rate was more than 1/100 * the best tag detection rate reduce chance of false positives. Check in the files if this cut off works for you
              a <- data$createAnt(); # creates an antID, i.e. associates a decimal antID number to an tagID
              identification <- data$addIdentification(a$ID,tag_statistics[j,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
              print(identification)
            }
          }
      tracking_data$save(paste0(directory, data_collection[i,"colony_id"], 'ants_created.myrmidon'))
      paste0(data_collection[i,"colony_id"], ' ants created')
}


#### select ants randomly from the trophy set ####
#### select random focal ants from the experiment #### 

####
#Before doing any actual tests on the trophy data I will need to finish preprocessing it e.g. with queen and run the Extrapolation!
  
