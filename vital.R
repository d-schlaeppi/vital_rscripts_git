#### vital  R script - all the steps ####

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###                                                                                             ###
###   ###               ###  ###  ### ### ### ### ###           ###           ###               ###
###    ###             ###   ###  ### ### ### ### ###          ## ##          ###               ###
###     ###           ###    ###          ###                 ### ###         ###               ###
###      ###         ###     ###          ###                ###   ###        ###               ###
###       ###       ###      ###          ###               ###     ###       ###               ###
###        ###     ###       ###          ###              ###  ###  ###      ###               ###
###         ###   ###        ###          ###             ### ### ### ###     ###               ###
###          ### ###         ###          ###            ###           ###    ###               ###
###           ## ##          ###          ###           ###             ###   ### ### ### ###   ###
###            ###           ###          ###          ###               ###  ### ### ### ###   ###
###                                                                                             ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Information before tracking
# Tracking systems - don t change settings for different colonies (only small focus adjustments)
# If you use separate systems for main tracking and treatment tracking: Don't mix them. Always the same tracking systems for main tracking and different systems for treatment tracking (Unless you can use exactly the same setup) 
# Use dedicated queen tags (0x000) if possible!

# this script contains:

# auto creation of 



# Step by step processing your tracking data:
# Step 1: For each tracking system used in tracking select one exemplary colony get mean worker size for each of them (you might need antID matcher and pose cloner for the treatment systems)
#   Step 1.1 Create base myrmidon files
#   Step 1.2 Automatically generate the ants for the selected tracking files using the "ant_generator"
#   Step 1.3 Manually orient files in fort myrmidon
#   Step 1.4 Get the mean worker size per tracking system using the "body_builder“

# Step 2: Run data extrapolation for all the data based on the mean body size in each tracking system (Nathalies script)

# Step 3: Use the extrapolated data for all the colonies do all the necessary post processing

# Step 3.1: Create all base myrmidon files
# Step 3.2: Generate all ants 
# Step 3.3: Automatically orient all ants in all extrapolated data using the ant_orient express
# Step 3.4: Generate capsules for all ant
# Step 3.4: For each tracking file do the necessary post processing
# see XY for post processing tips
# https://uob.sharepoint.com/:w:/r/teams/grp-AntsEpidemiologyLab/_layouts/15/Doc.aspx?sourcedoc=%7B2562631B-A6E5-4289-907F-89502F6C27E6%7D&file=pre-processing_Adriano_June2022.docx&action=default&mobileredirect=true

# Step 4: Data analysis 

#### prerequisites ####


# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()

# read in table that contains a collection of information on the experiments and each colony
setwd("/home/gw20248/Documents/vital_rscripts_git/")
dat <- read.csv("fc2_overview_data.csv", header = TRUE, stringsAsFactors = F)

# create a dataframe containing usefull stuff for the myrmidon file
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



# set working directory for the myrmidon files 
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
# directory <- '/media/gw20248/gismo_hd3/vital/fc2/' 
setwd(directory)





#### 3.1 Create all base myrmidon files ####

# !!! Manually copy the file base_source.myrmidon into the directory containing the tracking files !!!

# get a list of all folders in the directory and compile them as a list containing only folders with the tracking data
list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE)
data_list <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
# exclusion of any myrmidon or text files 
data_list <- grep(data_list, pattern = '.myrmidon', invert = TRUE, value = TRUE) 
data_list <- grep(data_list, pattern = '.txt', invert = TRUE, value = TRUE)
data_list

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
      # assign the tracking data - does not currently work and needs to be done manually.
      #for (folder_name in data_list) {
      #  if (grepl(data_collection[i,"colony_nr"], folder_name, fixed = TRUE) & !grepl("feeding", folder_name, fixed = TRUE)) {
      #    tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name))
      #  } else {next}
      #}
      # save the file base file with created ants 
      tracking_data$save(paste0(directory, file_name))
    } else {
      file_name <- paste0(data_collection[i,"colony_nr"], '_f_base.myrmidon') # file name will be of the following structure: c01_m_base.myrmidon (m = main tracking vs f = feeding)
      tracking_data <- fmExperimentOpen("base_source.myrmidon")
      s <- tracking_data$createSpace(data_collection[i,"tracking_system_feeding"])
      printf("Space '%s' has ID: %d\n",s$name,s$ID)
      tracking_data$name <- paste0("vital fc2 ",data_collection[i,"colony_nr"], " feeding")
      # assign the tracking data
      # for (folder_name in data_list) {
      # if (grepl(data_collection[i,"colony_nr"], folder_name, fixed = TRUE) & grepl("feeding", folder_name, fixed = TRUE)) {
      #    tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name))
      #  }
      #}
      # save the file base file with created ants 
      tracking_data$save(paste0(directory, file_name))
    }
  }
}


#### 3.1.1 Extrastep manually assign to each of the base files the corresponding data in fort ####
 # needs to be done because addTrackingDataDirectory() does not work? otherwise this could all be done automatically and ant creation could go in the same loop


#### Create ant for all myrmidon files (main = m and feeding = f) ####

types <- c("m", "f")
for (i in 1:nrow(data_collection)) {
  for(type in types) {
    data <- fmExperimentOpen(paste0(data_collection[i,"colony_nr"], "_", type, "_base.myrmidon"))
    tag_statistics <- fmQueryComputeTagStatistics(data)
    for (j in 1:nrow(tag_statistics)) {  # loop over each tag
    if (tag_statistics[j,"count"] >= 0.01*max(tag_statistics[,"count"],na.rm=T) ) { # optional: create an antID only if the tag detection rate was more than 1/100 * the best tag detection rate reduce chance of false positives. Check in the files if this cut off works for you
      a <- data$createAnt(); # creates an antID, i.e. associates a decimal antID number to an tagID
      identification <- data$addIdentification(a$ID,tag_statistics[j,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
      print(identification)
    }
  }
  tracking_data$save(paste0(directory, substr(paste0(data_collection[i,"colony_nr"], "_", type, "_base.myrmidon"), 1, nchar(paste0(data_collection[i,"colony_nr"], "_",type, "_base.myrmidon"))-13), 'AntsCreated.myrmidon'))
  } #types
}#colonies


# Check on a few colonies if the ant creation has worked? 


### Check: print identifications
data <- fmExperimentOpen("c28_f_AntsCreated.myrmidon") # enter the name of a file you would like to check
ants <- data$ants
for (a in ants) {
  printf("Ant %s is identified by:\n", fmFormatAntID(a$ID))
  for (i in a$identifications){
    printf(" * %s\n", capture.output(i))
  }
}






