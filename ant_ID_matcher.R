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
setwd(directory)

# select the two related files you want to work with (both need to have ants already created) 
main_file_name <- "vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon" # main tracking file
secondary_file_name <- "vital_fc2_esterhase_c02_feeding_DS_AntsCreated.myrmidon" # treatment tracking file 


#### loop ####

for (dataset_name in c(main_file_name, secondary_file_name)){
  # get data
  fort_data <- fmExperimentOpen(dataset_name)
  # create the key variables you want as metadata in your data sets 
  fort_data$setMetaDataKey(key = "meta_ID", default_Value = 001)
  fort_data$setMetaDataKey(key = "queen", default_Value = FALSE)
  fort_data$setMetaDataKey(key = "treated", default_Value = FALSE)
  # for the main file: define meta_ID so it corresponds to antID (no overriding yet) and set tagID OOO as the queen
  if (dataset_name == main_file_name) {
    for (i in 1:length(fort_data$ants)) {
      fort_data$ants[[i]]$setValue(key="meta_ID", value = c(fort_data$ants[[i]]$identifications[[1]]$targetAntID), time = fmTimeSinceEver())
    }
    # create vector of the treated ants
    treatment_data <- fmExperimentOpen(secondary_file_name)
    treated_ants <- treatment_data$ants
    tag_value_vector <- NULL
    tag_values <- NULL
    for (i in treated_ants) {
      tag_values <- i$identifications[[1]]$tagValue
      tag_value_vector <- rbind(tag_value_vector, data.frame(tag_values))
    }
    # for each ant adjust the meta data if it is the queen or a treated worker
    ants <- fort_data$ants
    for (i in ants) {
      if (i$identifications[[1]]$tagValue==0) {
        i$setValue("queen", TRUE, time = fmTimeSinceEver())}
      if(is.element(i$identifications[[1]]$tagValue, as.matrix(tag_value_vector))) {
        i$setValue(key="treated", value = TRUE, time = fmTimeSinceEver())}
    }
    # for the treatment data: define meta_ID so the id matches with the id generated for the main tracking file  
  } else {
    source_data <- fmExperimentOpen(paste0(substr(main_file_name,1, nchar(main_file_name)-9),'_metaID.myrmidon'))
    for (a in source_data$ants) {
      for (b in fort_data$ants) {
        if (a$identifications[[1]]$tagValue == b$identifications[[1]]$tagValue) {
          b$setValue("meta_ID",
                     value = as.numeric(source_data$ants[[a$getValue(key="meta_ID", time=fmTimeSinceEver())]]$getValue(key="meta_ID", time = fmTimeSinceEver())), 
                     time = fmTimeSinceEver())
        }
      }
    }
  }
  # save new version of myrmidon files
  fort_data$save(paste0(directory, substr(dataset_name, 1, nchar(dataset_name)-9),'_metaID.myrmidon'))
}


# think about making the previous steps as a loop for all tracking file simultaneously so I do not need to do it by hand
# either by (i) adjusting the file names or by (ii) running a special paste thing that creates a list with all the pairs of myrmidon files (main and treatment tracking)
# idea - create a list with all the trackings recucing the names to ColonyXX_main and ColonyXX_treatment











