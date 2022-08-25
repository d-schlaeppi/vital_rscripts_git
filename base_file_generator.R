rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      #### Base File Generator  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####

# automatically generate the myrmidon base files for all tracking files

#### prerequisites ####

# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()

setwd("/home/gw20248/Documents/vital_rscripts_git/")
dat <- read.csv("fc2_overview_data.csv", header = TRUE, stringsAsFactors = F)

# set working directory
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(directory)
list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE)


# idea: create a data frame containing information for each colony? and merge it to the data from the masterfile? 

data_collection <- NULL 

for(i in 1:nrow(dat)) {
  # collect variables
  nr                    <- i
  colony_id             <- dat[i, "colony_id"]
  block                 <- dat[i, "block"]
  colony_nr             <- paste0("c", sprintf("%02d", i))
  treatment             <- dat[i, "treatment"]
  food_position_1       <- dat[i, "food_position_1"]
  food_position_2       <- dat[i, "food_position_2"]
  tracking_system_main  <- dat[i, "tracking_system_main"]
  tracking_system_feeding <- dat[i, "tracking_system_feeding"] 
  # combine variables to a data frame  
  data_collection <-  rbind(data_collection, data.frame(nr, 
                                                        colony_id,
                                                        block, 
                                                        colony_nr,
                                                        treatment,
                                                        food_position_1,
                                                        food_position_2,
                                                        tracking_system_main,
                                                        tracking_system_feeding,
                                                        stringsAsFactors = F))
  }


# get all folders in the directory and compile them as a list (only folders containing the tracking data)
data_list <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
data_list <- grep(data_list, pattern = '.myrmidon', invert = TRUE, value = TRUE)
data_list <- grep(data_list, pattern = '.txt', invert = TRUE, value = TRUE)


# first an empty .myrmidon file is created manually as a source file using fort studio: base_source.myrmidon this is then used to create the rest of the files... 

#### Main Loop for file creation  ####
# first an empty .myrmidon file is created manually as a source file using fort studio: base_source.myrmidon this is then used to create the rest of the files... 

tracking_type <- c("main", "feeding")
for (i in 1:nrow(data_collection)) {
  for (type in tracking_type) {
    if (type == "main") {
      file_name <- paste0(data_collection[i,"colony_nr"], '_m_base.myrmidon') # file name will be of the following structure: c01_m_base.myrmidon (m = main tracking vs f = feeding)
      tracking_data <- fmExperimentOpen("base_source.myrmidon")
      s <- tracking_data$createSpace(data_collection[i,"tracking_system_main"])
      printf("Space '%s' has ID: %d\n",s$name,s$ID)
      tracking_data$name <- paste0("vital fc2 ",data_collection[i,"colony_nr"], " main")
      # assign the tracking data
      for (folder_name in data_list) {
        if (grepl(data_collection[i,"colony_nr"], folder_name, fixed = TRUE) & !grepl("feeding", folder_name, fixed = TRUE)) {
          tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name))
        } else {next}
      }
      # save the file base file with created ants 
      tracking_data$save(paste0(directory, data_collection[i,"colony_nr"], '_main.myrmidon'))
    } else {
      file_name <- paste0(data_collection[i,"colony_nr"], '_f_base.myrmidon') # file name will be of the following structure: c01_m_base.myrmidon (m = main tracking vs f = feeding)
      tracking_data <- fmExperimentOpen("base_source.myrmidon")
      s <- tracking_data$createSpace(data_collection[i,"tracking_system_feeding"])
      printf("Space '%s' has ID: %d\n",s$name,s$ID)
      tracking_data$name <- paste0("vital fc2 ",data_collection[i,"colony_nr"], " feeding")
      # assign the tracking data
      for (folder_name in data_list) {
        if (grepl(data_collection[i,"colony_nr"], folder_name, fixed = TRUE) & grepl("feeding", folder_name, fixed = TRUE)) {
          tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name))
        }
      }
      # save the file base file with created ants 
      tracking_data$save(paste0(directory, data_collection[i,"colony_nr"], '_feeding.myrmidon'))
    }
  }
}

### Error with addTrackingDataDirectory() --> reported on github 
# thus, the code does not work yet

#### once the issue with the tracking data valid method error is resolved The loop above should work to create all the base files..
# if it does not work the files can be created automatically, but the directories need to be added manually which is annoying.
# but if it works the ant creation code could be run within the same loop to save time. 







