rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
     #### 1. ANT-ID-MATCHER  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####

# If ants are used in two separate tracking systems it is useful (or even essential) to make sure that based on the tagID the same AntID's are created in both files
# (e.g. tracking of a colony in one tracking system and then a few workers get sampled for treatment that is also tracked and then the workers get returned to their colony)
# when creating antID automatically based on tagID, antID are automatically created starting with 001. That means you will have e.g 200 in your main colony and you take a random sample of X workers from that
# in the new tracking file you also need to assign antID ant they will be automatically created from 001-0XX --> antID will not automatically be matched with the original antID 
# Because antID is read only and can not be overwritten we will create a new variable - meta_ID both for the main tracking and for the treatment tracking files
# For anything following this step we will need to call meta_ID instead of ant_ID if we refer to ants to avoid any confusion

# for more information on fort-myrmidon and fort-studio see: 
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"



# define/get the starting time of the experiment
# t <- fmTimeCreate(offset = 0) #SET TIME TO 1970 which is per definition way before the actual start of the experiment.
# t <- fmTimeCreate(offset = fmQueryGetDataInformations(main_tracking_data)$start)
# other useful syntax for other things when dealing with tracking system time: 
# from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start) # experiment start time
# to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start + 12*3600  ) # start time plus the duration in seconds
# loop to create create/assign each ant of the main set ID and use the tag 



#### prerequisites ####

# load libraries
library(FortMyrmidon) #R bindings

# set working directory
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(directory) # adjust the script the two files are not in the same folder

# select the files you want to work with (both need to have ants already created)
main_file_name <- "vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon" 
secondary_file_name <- "vital_fc2_esterhase_c02_feeding_DS_AntsCreated.myrmidon"




#### meta data creation ####

# define new key (metadata), including meta_ID that will be the one to be matched between the files

for (dataset_name in c(main_file_name, secondary_file_name)){
  # get data
  fort_data <- fmExperimentOpen(dataset_name)
  start_experiment <- fmTimeCreate(offset = fmQueryGetDataInformations(fmExperimentOpen(main_file_name))$start)
  # create key variable you want for your data sets (for now only the meta_ID is relevant)
  fort_data$setMetaDataKey(key = "meta_ID", default_Value = 001)
  fort_data$setMetaDataKey(key = "queen", default_Value = FALSE)
  fort_data$setMetaDataKey(key = "name", default_Value = "NA")
  fort_data$setMetaDataKey(key = "treated", default_Value = FALSE)
  # define meta_ID so it corresponds to antID (no overriding yet)
  for (i in 1:length(fort_data$ants)) {
    fort_data$ants[[i]]$setValue(key="meta_ID", value = c(fort_data$ants[[i]]$identifications[[1]]$targetAntID), time = start_experiment)
  }
  # save new version of myrmidon files
  fort_data$save(paste0(directory, substr(dataset_name, 1, nchar(dataset_name)-9),'_meta.myrmidon'))
  print("done")
}


#### Matching ant meta_ID in the secondary data (feeding tracking ####
# update filenames
main_file_name <- paste0(substr(main_file_name, 1, nchar(main_file_name)-9),'_meta.myrmidon')
secondary_file_name <- paste0(substr(secondary_file_name, 1, nchar(secondary_file_name)-9),'_meta.myrmidon')

# load the newly created myrmidon files
#main_data <- fmExperimentOpen(paste0(substr(main_file_name, 1, nchar(main_file_name)-9),'_meta.myrmidon'))
#secondary_data <- fmExperimentOpen(paste0(substr(secondary_file_name, 1, nchar(secondary_file_name)-9),'_meta.myrmidon')) 

# change metadata so that (i) the meta_ID of the ants in the treatment corresponds to the original antID based on tagValue and (ii) ants in the feeding tracking have metadata treated == TRUE

# recreate the loop so it matches the structure from the previous loop for variable creation 
# before running any further code check the last created file if there are two different antID variables assigned at different time points (start of main tracking and start of secondary tracking!

for(dataset_name in c(main_file_name, secondary_file_name)) {
  if (dataset_name == main_file_name) {print("no change required")} else {
    fort_data <- fmExperimentOpen(dataset_name)
    source_data <- fmExperimentOpen(main_file_name)
    start_experiment <- fmTimeCreate(offset=fmQueryGetDataInformations(source_data)$start)
    for (a in source_data$ants) {
      for (b in fort_data$ants) {
        if (a$identifications[[1]]$tagValue == b$identifications[[1]]$tagValue) {
          b$setValue("meta_ID",
                     value = as.numeric(source_data$ants[[a$getValue("meta_ID", time = start_experiment)]]$getValue(key="meta_ID", time = start_experiment)), 
                     time = start_experiment)
          } #if2
        }#for3
      }#for2
    }#else1
  fort_data$save(paste0(directory, substr(dataset_name, 1, nchar(dataset_name)-9),'_metaIDmatched.myrmidon'))
  print("done")
}#for1

    

                              
# last step, combine the two loops to one. so that we do not first create a pseudo meta_ID for the secondary file!
# next steps: include meta data treatment and think about making the previous steps as a loop for all tracking file simultaneously so I do not need to do it by hand
# either by (i) adjusting the file names or by (ii) running a special paste thing that creates a list with all the pairs of myrmidon files (main and treatment tracking)
# idea - create a list with all the trackings recucing the names to ColonyXX_main and ColonyXX_treatment


### combine the two loops to a single loop ###
 # include a line to define the queen. 


main_file_name <- "vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon" 
secondary_file_name <- "vital_fc2_esterhase_c02_feeding_DS_AntsCreated.myrmidon"





for (dataset_name in c(main_file_name, secondary_file_name)){
  # get data
  fort_data <- fmExperimentOpen(dataset_name)
  start_experiment <- fmTimeCreate(offset = fmQueryGetDataInformations(fmExperimentOpen(main_file_name))$start)
  # create key variable you want for your data sets (for now only the meta_ID is relevant)
  fort_data$setMetaDataKey(key = "meta_ID", default_Value = 001)
  fort_data$setMetaDataKey(key = "queen", default_Value = FALSE)
  fort_data$setMetaDataKey(key = "name", default_Value = "NA")
  fort_data$setMetaDataKey(key = "treated", default_Value = FALSE)
  # for the main file: define meta_ID so it corresponds to antID (no overriding yet)
  if (dataset_name == main_file_name) {
    for (i in 1:length(fort_data$ants)) {
      fort_data$ants[[i]]$setValue(key="meta_ID", value = c(fort_data$ants[[i]]$identifications[[1]]$targetAntID), time = start_experiment)
    }
    ants <- fort_data$ants
    for (i in ants) {
      if (i$identifications[[1]]$tagValue==0) {
        i$setValue("queen", TRUE, time = start_experiment)}
    } 
  } else { 
    source_data <- fmExperimentOpen(paste0(substr(main_file_name,1, nchar(main_file_name)-9),'_metaID.myrmidon'))
    for (a in source_data$ants) {
      for (b in fort_data$ants) {
        if (a$identifications[[1]]$tagValue == b$identifications[[1]]$tagValue) {
          b$setValue("meta_ID",
                     value = as.numeric(source_data$ants[[a$getValue(key="meta_ID", time=start_experiment)]]$getValue(key="meta_ID", time = start_experiment)), 
                     time = start_experiment)
        }
      }
    }
  }
  # save new version of myrmidon files
  fort_data$save(paste0(directory, substr(dataset_name, 1, nchar(dataset_name)-9),'_metaID.myrmidon'))
}




### next step copy and info on size and capsules





