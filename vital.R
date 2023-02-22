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

#### Background information | Read me ####

# Information before tracking
# Tracking systems - don t change settings for different colonies (only small focus adjustments)
# If you use separate systems for main tracking and treatment tracking: Don't mix them. Always the same tracking systems for main tracking and different systems for treatment tracking (Unless you can use exactly the same setup) 
# If you do short trackings (e.g. for treatments in addition to a main tracking increase the rate at which pictures are taken of each ant in the leto file)
# Use dedicated queen tags (0x000) if possible!

# used tracking systems Daniel: Maintracking - trojan, prideaux, guillam ; feeding sessions - guillam, esterhase 

# Step by step processing your tracking data:
# Step 1: For each tracking system setting (typically 1 per tracking system) used select one exemplary colony get mean worker size
#   Step 1.1 Create base myrmidon files
#   Step 1.2 Automatically generate the ants for the selected tracking files using the "ant_generator"
#   Step 1.3 Manually orient files in fort myrmidon
#   Step 1.4 Get the mean worker size per tracking system using the "ant_ruler“

# Step 2: Run data extrapolation for all the data based on the mean body size in each tracking system (Nathalies script)

# Step 3: Use the extrapolated data for all the colonies do all the necessary post processing

# Step 3.1: Create all base myrmidon files
# Step 3.2: Generate all ants
# Step 3.3: Automatically create the metadata variables needed
# Step 3.4: Adjust replaced or re-glued tags and do required manual post processing
# Step 3.5: Automatically orient all ants in all extrapolated data using the ant_orient express (Includes capsule generation)
# Step 3.6: Post processing of queen meta data (Manual)

# Step 4: Data analysis 

### ### ### ### ###
### Useful links ##
### ### ### ### ###

# for more information on fort-myrmidon and fort-studio see: 
# https://formicidae-tracker.github.io/myrmidon/latest/index.html

# Postprocessing tips
# https://uob.sharepoint.com/:w:/r/teams/grp-AntsEpidemiologyLab/_layouts/15/Doc.aspx?sourcedoc=%7B2562631B-A6E5-4289-907F-89502F6C27E6%7D&file=pre-processing_Adriano_June2022.docx&action=default&mobileredirect=true

# AEL Github repositories
# https://github.com/d-schlaeppi/vital_rscripts_git
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR
# https://github.com/Leckie-Bris/SICC
# https://github.com/EnricoGavagnin?tab=repositories



#### prerequisites ####
rm(list=ls())

# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()
library(Rcpp)         # contains sourceCpp (used for ant orientation)
library(circular)     # used for ant orientation
library(data.table)   # used to save files fwrite(list(myVector), file = "myFile.csv") 



### read in table that contains a collection of information on the experiments and each colony, and create a data frame containing all the useful things for the myrmidon files
setwd("/home/gw20248/Documents/vital_rscripts_git/")
dat <- read.csv("fc2_overview_data.csv", header = TRUE, stringsAsFactors = F)
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



# set working directory for the myrmidon files (before the extrapolation)
#directory <- "/home/gw20248/Documents/data_copy_for_trials/"  # preliminary processing on local computer
# or
#directory <- '/media/gw20248/gismo_hd3/vital/fc2/' # HD with raw data for preparations
# or 
directory <- '/media/gw20248/gismo_hd2/vital/fc2/'  # HD with extrapolated data
setwd(directory)

# use the below directory if you are already at the later steps
# set working directory to the new directory containing the extrapolated tracking data
# directory <- '/media/gw20248/gismo_hd2/vital/fc2/' 
# setwd(directory)


#### Step 1 ####
# Step 1: For each tracking system setting (typically 1 per tracking system) used select one exemplary colony get mean worker size

# select the tracking files you selected for size extraction 

#### 1.1 Create base myrmidon files ####  

# Was done manually so far, because there is an error with the base_file_generator
# Error with addTrackingDataDirectory() --> issue has been reported on github but no reply yet. 

#### Step 1.2 Automatically generate the ants for the selected tracking files using the "ant_generator" ####

files <- list(
  paste(directory,"vital_fc2_guillam_c03_DS_base.myrmidon",sep=''), 
  paste(directory,"vital_fc2_trojan_c27_DS_base.myrmidon",sep=''),
  paste(directory,"vital_fc2_prideaux_c02_DS_base.myrmidon",sep=''), 
  paste(directory,"vital_fc2_esterhase_c02_feeding_DS_base.myrmidon",sep=''),
  paste(directory,"vital_fc2_guillam_c12_feeding_DS_base.myrmidon",sep=''),
  paste(directory,"vital_fc2_guillam_c27_feeding_DS_base.myrmidon",sep='')
)

for (file in files) {
  tracking_data <- fmExperimentOpen(file) # files that need the ants created
  tag_statistics <- fmQueryComputeTagStatistics(tracking_data)   # extract the tag statistics to know how many times each tag was detected etz
  # create ants - decide on the cutoff: stronger cut off to reduce the chances of false positives if too many false ants are created | tag_statistics[,"count"] #check for the count numbers to see what range should be included (in this case the cutoff is reduced from 0.001 to 0,01)
  for ( i in 1:nrow(tag_statistics)) {  # loop over each tag
    if ( tag_statistics[i,"count"] >= 0.01*max(tag_statistics[,"count"],na.rm=T) ) { # antID only if the tag detection rate was more than 1/100 (adriano used 1/1000) of the best tag detection rate
      a <- tracking_data$createAnt(); # creates an antID, i.e. associates a decimal antID number to that particular tagID
      identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
      print(identification)
    }
  }
  tracking_data$save(paste0(substr(file, 1, nchar(file)-13), 'AntsCreated.myrmidon')) # Save to newly named myrmidon file: "previousname_AntsCreated"
}

#### 1.3 Manually orient files in fort myrmidon #### 

# still needs to be done manually... at least it is only one colony per tracking system / tracking system setting 

#### 1.4 Get the mean worker size per tracking system using the "ant_ruler“ ####

files <- list(
  paste(dir_data,"vital_fc2_guillam_c03_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''), 
  paste(dir_data,"vital_fc2_trojan_c27_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''), 
  paste(dir_data,"vital_fc2_esterhase_c02_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_guillam_c12_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_guillam_c27_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep='')
)

#### ANT-RULER ####

output_name <- file.path(paste0(dir_data,"Mean_ant_length_colonies_new.txt")) # define the name of the textfile containing your measurements

for (element in files) {
  # get tracking data
  ant_measurements <- NULL
  tracking_data <- fmExperimentOpen(element) 
  ants <- tracking_data$ants
  for (ant in ants){ # for each ant get the mean size and put it in a dataframe
    ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_px)
    ant_length_mm <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_mm)
    ant_measurements <- rbind(ant_measurements, data.frame(length_px = ant_length_px,
                                                           length_mm = ant_length_mm,
                                                           stringsAsFactors = F))
  }
  # queen exclusion from the dataframe
  interquartile_range <- quantile(ant_measurements$length_px,probs=c(0.25,0.75), na.rm =TRUE)
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  ant_measurements <- ant_measurements[which(ant_measurements$length_px>=outlier_bounds[1]&ant_measurements$length_px<=outlier_bounds[2]),]
  # printing and saving of the means from the dataframes without the queen
  print(element)  
  print(mean(ant_measurements$length_px, na.rm=TRUE))
  print(mean(ant_measurements$length_mm, na.rm=TRUE))
  table <- NULL
  table <- rbind(table, data.frame(mean(ant_measurements$length_px, na.rm=TRUE),
                                   mean(ant_measurements$length_mm, na.rm=TRUE),
                                   element, 
                                   stringsAsFactors = F))
  if (file.exists(output_name)){
    write.table(table,file=output_name,append=T,col.names=F,row.names=F,quote=T,sep=",")
  }else{
    write.table(table,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
  }
}

#### Data extrapolation ####
# Run data extrapolation for all the data based on the mean body size of each tracking system (Nathalie's script on the computer in her office)








#### 3.1 Create all base myrmidon files for extrapolated data  ####

# !!! Manually copy the file base_source.myrmidon into the directory containing the tracking files !!!

# set working directory to the new directory containing the extrapolated tracking data
directory <- '/media/gw20248/gismo_hd2/vital/fc2/' 
setwd(directory)


# get a list of all folders in the directory and compile them as a list containing only folders with the tracking data
list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE)
data_list <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
# exclusion of any myrmidon or text files 
data_list <- grep(data_list, pattern = '.myrmidon', invert = TRUE, value = TRUE) 
data_list <- grep(data_list, pattern = '.txt', invert = TRUE, value = TRUE)
data_list <- grep(data_list, pattern = 'climate', invert = TRUE, value = TRUE)
data_list <- grep(data_list, pattern = 'Data_Read_Me', invert = TRUE, value = TRUE)
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
      # how to fix the issue with adding the traking data directories
      # In the latest version of the R bindings, fmExperiment$addTrackingdataDirectory() takes an additional fixCorruptedData argument.
      # The way the rcpp package works (the package that allows a R program to interface with C++), if you fail to provide all needed arguments, it returns this very cryptic method, which does not tell you to add the missing value.
      # To get the old behavior please use FALSE, and the call will fail if there is a data corruption. Use TRUE to ask to not fail but try to fix any encountered error (will cause permanent data loss, but let you recover as much data as possible).
      tracking_data$save(paste0(directory, file_name)) # save the file base file with created ants
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


#### 3.1.1 Extra step manually assign to each of the base files the corresponding data in fort ####

 # needs to be done because addTrackingDataDirectory() does not work yet ? otherwise this could all be done automatically and ant creation could go in the same loop
 # there was a fix for this (see trophy data) but with the fix there seems to be a new error that causes fort to crash. No further investigations done yet.

#### 3.2 Create ants for all myrmidon files (main = m and feeding = f) ####

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
  data$save(paste0(directory, substr(paste0(data_collection[i,"colony_nr"], "_", type, "_base.myrmidon"), 1, nchar(paste0(data_collection[i,"colony_nr"], "_",type, "_base.myrmidon"))-13), 'AntsCreated.myrmidon'))
  } #types
}#colonies


# Check on a few colonies if the ant creation has worked? 




# now that we have all the ant identifications right and marked the times at which the tags were replaced and reoriented we can run the automated orientation script. 

#### 3.3 Create all the metadata keys ####

for (i in 1:nrow(data_collection)) {
    fort_data <- fmExperimentOpen(paste0(data_collection[i,"colony_nr"], "_m_AntsCreated.myrmidon"))
    fort_data$setMetaDataKey(key = "meta_ID",     default_Value = 001)  # create the key variables you want as metadata in your data sets 
    fort_data$setMetaDataKey(key = "IsQueen",     default_Value = FALSE)
    fort_data$setMetaDataKey(key = "IsTreated",   default_Value = FALSE)
    fort_data$setMetaDataKey(key = "IsAlive",     default_Value = TRUE)
    fort_data$setMetaDataKey(key = "treatment",   default_Value = "NA") # treated ants will get control or virus
    fort_data$setMetaDataKey(key = "glass_beads", default_Value = "NA") # treated ants will get yellow or blue
    fort_data$setMetaDataKey(key = "comment",     default_Value = "NA")
    fort_data$setMetaDataKey(key = "tag_reoriented", default_Value = FALSE)
    fort_data$spaces[[1]]$createZone(name = "nest") # create zones to be defined manually in the fort files 
    fort_data$spaces[[1]]$createZone(name = "arena")
    fort_data$spaces[[1]]$createZone(name = "water_left")
    fort_data$spaces[[1]]$createZone(name = "sugar_right")
    for (y in 1:length(fort_data$ants)) {
      fort_data$ants[[y]]$setValue(key="meta_ID", value = c(fort_data$ants[[y]]$identifications[[1]]$targetAntID), time = fmTimeSinceEver())
    }
    treatment_data <- fmExperimentOpen(paste0(data_collection[i,"colony_nr"], "_f_AntsCreated.myrmidon"))     # create vector of the treated ants
    treated_ants <- treatment_data$ants
    tag_value_vector <- NULL
    tag_values <- NULL
    for (z in treated_ants) {
      tag_values <- z$identifications[[1]]$tagValue
      tag_value_vector <- rbind(tag_value_vector, data.frame(tag_values))
    }
    ants <- fort_data$ants  # for each ant adjust the meta data if it is the queen or a treated worker
    for (x in ants) {
      if (x$identifications[[1]]$tagValue==0) {
        x$setValue("IsQueen", TRUE, time = fmTimeSinceEver())}
      if(is.element(x$identifications[[1]]$tagValue, as.matrix(tag_value_vector))) {
        x$setValue(key="IsTreated", value = TRUE, time = fmTimeSinceEver())}
    }
    fort_data$save(paste0(directory, substr(paste0(data_collection[i,"colony_nr"], "_m_base.myrmidon"), 1, nchar(paste0(data_collection[i,"colony_nr"], "_m_base.myrmidon"))-13), 'meta_keyed.myrmidon'))
}




#### 3.4 Adjust replaced or re-glued tags and other post processing steps ####
# For the tags that got lost and replaced or reglued with a new orientation (see experiment notes which ants had this)
# Manually adjust the myrmidon file according to adrianos post processing guide
# https://uob.sharepoint.com/:w:/r/teams/grp-AntsEpidemiologyLab/_layouts/15/Doc.aspx?sourcedoc=%7B2562631B-A6E5-4289-907F-89502F6C27E6%7D&file=pre-processing_Adriano_June2022.docx&action=default&mobileredirect=true&cid=5b8c1184-40b2-4f60-8bd7-e5801d42d6f5
# And do the rest of the manual post processing and meta data creation... e.g. Dead workers, Zones 
# for the ants that got a new tag also update the meta variable tag_reoriented = true so they will be skipped in the auto orientation just like the queen. 

#### Ant Orient Express ####
#### Ant Orient Express Part 1 ####
# get capsule information from a source file for which capsules were manually defined

# source C++ movement direction program
# https://uob.sharepoint.com/teams/grp-AntsEpidemiologyLab/Shared%20Documents/Forms/AllItems.aspx?id=%2Fteams%2Fgrp%2DAntsEpidemiologyLab%2FShared%20Documents%2FLinux%5Fand%5FR%5Fguides%2FR%5FScripts%2FAutomated%5FAnt%5FOrientation%5FNathalie&viewid=cbc49cd2%2D692a%2D42a7%2D9f4e%2D29be55b2f252
# copy the file from sharepoint into your directory and access it using the following code
sourceCpp(paste0(directory,"Get_Movement_Angle.cpp"))

# personal notes: 
# so far only main colonies included

# list a source myrmidon file containing the manually oriented data with capsules
# (choose 1 of the manually oriented colonies to define a capsule definition in fort studio for a medium sized ant and replicate the shape for all of the ants of the colony)
#source_data_list <- list("/home/gw20248/Documents/data_copy_for_trials/vital_fc2_trojan_c27_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon") #,
#                         "/home/gw20248/Documents/data_copy_for_trials/vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon")

source_data_list <- list(paste0(directory,"vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon"))

# get the information on the caps from the source file
oriented_metadata <- NULL
capsule_list <- list()
for (myrmidon_file in source_data_list){
  experiment_name <- unlist(strsplit(myrmidon_file,split="/"))[length(unlist(strsplit(myrmidon_file,split="/")))]
  oriented_data <- fmExperimentOpen(myrmidon_file)
  oriented_ants <- oriented_data$ants
  capsule_names <- oriented_data$antShapeTypeNames
  for (ant in oriented_ants){
    # extract ant length and capsules
    ant_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px)
    capsules      <- ant$capsules
    for (caps in 1:length(capsules)){
      capsule_name  <- capsule_names[[capsules[[caps]]$type]]
      capsule_coord <- capsules[[caps]]$capsule
      capsule_info <- data.frame(experiment = experiment_name,
                                 antID      = ant$ID,
                                 c1_ratio_x = capsule_coord$c1[1]/ant_length_px,
                                 c1_ratio_y = capsule_coord$c1[2]/ant_length_px,
                                 c2_ratio_x = capsule_coord$c2[1]/ant_length_px,
                                 c2_ratio_y = capsule_coord$c2[2]/ant_length_px,
                                 r1_ratio   = capsule_coord$r1[1]/ant_length_px,
                                 r2_ratio   = capsule_coord$r2[1]/ant_length_px
      )
      if (!capsule_name %in%names(capsule_list)){ # if this is the first time we encounter this capsule, add it to capsule list...
        capsule_list <- c(capsule_list,list(capsule_info)) 
        if(length(names(capsule_list))==0){
          names(capsule_list) <- capsule_name
        }else{
          names(capsule_list)[length(capsule_list)] <- capsule_name
        }
      }else{ #otherwise, add a line to the existing dataframe within capsule_list
        capsule_list[[capsule_name]] <- rbind(capsule_list[[capsule_name]] , capsule_info)
      }
    }
    # extract offset btewen tag centre and ant centre
    for (id in ant$identifications){
      oriented_metadata <- rbind(oriented_metadata,data.frame(experiment       = experiment_name,
                                                              antID            = ant$ID,
                                                              tagIDdecimal     = id$tagValue,
                                                              angle            = id$antAngle,
                                                              x_tag_coord      = id$antPosition[1], 
                                                              y_tag_coord      = id$antPosition[2],
                                                              x_ant_coord      = id$antPosition[1]*cos(-id$antAngle) - id$antPosition[2]*sin(-id$antAngle),
                                                              y_ant_coord      = id$antPosition[1]*sin(-id$antAngle) + id$antPosition[2]*cos(-id$antAngle),
                                                              length_px        = ant_length_px,
                                                              stringsAsFactors = F))
    }
  }
}

# Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen
# So we need to remove the queen from the computation by removing outliers in ant lenght measurements
interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75))
outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  # apply outlier exclusion to oriented_metadata and to capsule list
for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
}
# Once queen(s) has(have) been removed, get the mean coordinates of the offset between tag centre and ant centre
mean_x_ant_coord <- mean(oriented_metadata$x_ant_coord)
mean_y_ant_coord <- 0 ##set it to zero manually because we don't expect there to be a consistent bias in deviation / #mean_y_ant_coord <- mean(oriented_metadata$y_ant_coord) ###this is expected to be 0 or near 0 (check!) because the tag should be as likely to be to the right or to the left of the ant's bilateral symmetry line
mean_worker_length_px <- mean(oriented_metadata$length_px) # Get the average worker length from the data
for (caps in 1:length(capsule_list)){ # Finally, get information on each capsule
  capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
}

#### Ant Orient Express Part 2 ####

# open the tracking file to be auto oriented
# tracking_data <- fmExperimentOpen(paste0(dir_data,'c11_m_AntsCreated_cor.myrmidon'))

files <- list.files(directory)
files <- files[grep("tags_corrected.myrmidon",files)]
#files <- files[grep("pre_orientation",files)]

not_oriented <- NULL
to_orient_manually <- NULL

for (file in files) {
  tracking_data <- fmExperimentOpen(paste0(directory, file))
  for (caps in 1:length(capsule_list)){ #add the caps from the source file above
    tracking_data$createAntShapeType(names(capsule_list)[caps])
  }
  ants <- tracking_data$ants
  # get trajectory data to extract ant orientation:
  # short tracking using all data -> using start and end time from the experiment metadata
  #from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start) # experiment start time # from <- fmTimeSinceEver()
  #to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$end  ) # experiment end time   # to   <- fmTimeForever()
  from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$end -12*3600) # longer tracking data - computation takes very long -> only use a subset of the data e.g. last 12 hours
  to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$end)
  max_gap <- fmHour(24*365)  # use a  large value to make sure you get only one trajectory per ant (larger than the time difference between from to defined above)
  positions <- fmQueryComputeAntTrajectories(tracking_data, start = from, end = to, maximumGap = max_gap, computeZones = TRUE)
  #Then compute trajectories & hard-wire ant correspondence between trajectories_summary and trajectorues
  positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") # creates a ID string for each ant: ant1, ant2,...
  names(positions$trajectories)       <- positions$trajectories_summary$antID_str # and use the content of that column to rename the objects within trajectory list
  max_time_gap <- 0.5 # define a max temporal gap for which you are happy to calculate a movement angle; e.g. 0.5 s
  min_dist_moved <- 30 # define a minimum distance moved, as you don't want to use noise or small shifts in position in this calculation; e.g. 30 pix (to think about)
  for (i in 1:length(ants)){
    if(length(tracking_data$ants[[i]]$identifications) == 0) {
      print(paste(file, i, "no ant", sep = " -> "))
      not_oriented <- append(not_oriented, paste(file, i, "no ant", sep = " -> "))
      next
    }
    #skip ants that got retagged and create a table (to be saved in a file) with all the skipped ants that need manual orientation
    if (ants[[i]]$getValue("tag_reoriented", fmTimeForever())) {
      to_orient_manually <- append(to_orient_manually, paste(file, i, ants[[i]]$ID, "retagged", sep = " -> "))
      print(paste(file, i, "skipped because of retagging", sep = " -> "))
      next
    }
    if (tracking_data$ants[[i]]$identifications[[1]]$tagValue==0) {next} # skip the queen
    if (is.element(i,positions$trajectories_summary$antID)) {
      # to be fool proof, and be sure you extract the trajectory corresponding the correct ant, make sure you make use of the antID_str column!
      traj <- positions$trajectories [[   positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"antID_str"]    ]]
      traj <- cbind(traj,add_angles(traj,max_time_gap,min_dist_moved))   # feed traj to c++ program
      AntAngle <- as.numeric(- mean(circular(na.omit(traj$Tag_minus_Movement_Angle),units="radians",zero=0)))   #  get mean deviation angle between body and tag - the ant angle is equal to minus the Tag minus Movement angle output by C++ program
      x_tag_coord <- mean_x_ant_coord*cos(AntAngle) - mean_y_ant_coord*sin(AntAngle)  # now use trigonometry to calculate the pose, using AntAngle
      y_tag_coord <- mean_x_ant_coord*sin(AntAngle) + mean_y_ant_coord*cos(AntAngle)
      for (id in ants[[i]]$identifications){  # write this into ant metadata
        id$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)
      }
      # also add this to trajectories_summary
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"ant_angle"] <- AntAngle
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"x_tag_coord"] <- x_tag_coord
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"y_tag_coord"] <- y_tag_coord
      # finally, for each ant, add capsules using mean_ant_length and capsule_list
      for (caps in 1:length(capsule_list)){
        capsule_ratios <- capsule_list[[caps]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
        capsule_coords <- mean_worker_length_px*capsule_ratios
        ants[[i]]$addCapsule(caps, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
      }
    } else {
      print(paste(file, i, "FALSE", "exit before selected period", sep = " -> "))
      not_oriented <- append(not_oriented, paste(file, i, "FALSE", "exit before selected period", sep = " -> "))
    }
  }
  tracking_data$save(paste0(directory, substr(file, 1, nchar(file)-24),'oriented.myrmidon'))
  not_oriented <- append(not_oriented, paste("time", Sys.time() ,sep = " : "))
}
fwrite(list(not_oriented), file = paste(Sys.Date(), "",format(Sys.time(), "%H-%M-%S"), "not_oriented.txt", sep = "_"))
fwrite(list(to_orient_manually), file = paste(Sys.Date(), "",format(Sys.time(), "%H-%M-%S"), "to_orient_manually.txt", sep = "_"))


#### Queen modification and other checks ####
# Next check your files manually to see if orientations look all right (especially treatment ants and ants that got retagged or reoriented)
   # open each of the files and do the following 
   # 1 - For each queen manually assign at least 3 poses - done
   # 2 - manually adjust orientation of any re tagged ants (check file to_orient_manually.txt) - done
   # 3 - manually copy capsules for re tagged ants: open the myrmidon files for all colonies with re tagged ants. Select an ant with capsules, clone shape, select no scaling and no overwriting - done
   # 4 - run the code to copy queen capsules form a source file and adjust tag size of queens automatically

# Manually create a queen capsule source file from an average sized queen - > Select the most average sized queen from my data set. 

# step 1 use the ant ruler to select the most medium size queen. 
# get mean size of all queens
# select the queen with the smallest difference to the mean. 

# set working directory (where the tracking data is)
directory <- '/media/gw20248/gismo_hd2/vital/fc2/'  # HD with extrapolated data
setwd(directory)
# get a list with the latest myrmidon file with manually oriented queens
files <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
files <- grep(files, pattern = "queen_modified", invert = FALSE, value = TRUE)

# calculate the mean queen size and identify the colony with the most average sized queen for manual capsule definition. 
queen_size <- NULL
for (file in files) {
  tracking_data <- fmExperimentOpen(paste0(directory, file))
  queen_length_mm <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=1)$length_mm)
  queen_size <- rbind(queen_size, data.frame(colony = file,
                                             queen_length_mm = queen_length_mm,
                                             stringsAsFactors = F))
}
queen_size$diff <- abs(queen_size$queen_length_mm - mean(queen_size$queen_length_mm))
source_colony <-queen_size$colony[which.min(queen_size$diff)]
print(source_colony)

# open a source colony file to manually create a well fitting queen shape and safe the file as queen_shape_source.myrmidon
# Next: automatically assign the queen capsules and adjust queen tag size

### select queen shape source file
defined_capsule_file <- paste(directory,"queen_shape_source.myrmidon",sep='')
### select files that need the queen shape
no_capsule_list <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
no_capsule_list <- grep(files, pattern = "queen_modified", invert = FALSE, value = TRUE)


### create empty variables
data_list         <- list (defined_capsule_file)
oriented_metadata <- NULL
capsule_list <- list()
### run loop to get the capsule information from the source file
for (myrmidon_file in data_list){
  oriented_data <- fmExperimentOpen(myrmidon_file)
  experiment_name <- "FCII" 
  oriented_ants <- oriented_data$ants
  capsule_names <- oriented_data$antShapeTypeNames
  ### extract queen length and capsules (this script only works if all queens were done with a queen tag 0x000 and thus always with antID = 1)
  queen_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=1)$length_px)
  capsules      <- oriented_ants[[1]]$capsules
  for (caps in 1:length(capsules)){
    capsule_name  <- capsule_names[[capsules[[caps]]$type]]
    capsule_coord <- capsules[[caps]]$capsule
    capsule_info <- data.frame(experiment = experiment_name,
                               antID      = 1,
                               c1_ratio_x = capsule_coord$c1[1]/queen_length_px,
                               c1_ratio_y = capsule_coord$c1[2]/queen_length_px,
                               c2_ratio_x = capsule_coord$c2[1]/queen_length_px,
                               c2_ratio_y = capsule_coord$c2[2]/queen_length_px,
                               r1_ratio   = capsule_coord$r1[1]/queen_length_px,
                               r2_ratio   = capsule_coord$r2[1]/queen_length_px)
    if (!capsule_name %in%names(capsule_list)){ ###if this is the first time we encounter this capsule, add it to capsule list...
      capsule_list <- c(capsule_list,list(capsule_info)) 
      if(length(names(capsule_list))==0){
        names(capsule_list) <- capsule_name
      }else{
        names(capsule_list)[length(capsule_list)] <- capsule_name
      }
    }else{###otherwise, add a line to the existing dataframe within capsule_list
      capsule_list[[capsule_name]] <- rbind(capsule_list[[capsule_name]] , capsule_info)
    }
  }
  ###extract offset between tag center and ant center
  for (id in oriented_ants[[1]]$identifications){
    oriented_metadata <- rbind(oriented_metadata,data.frame(experiment       = experiment_name,
                                                            antID            = 1,
                                                            tagIDdecimal     = id$tagValue,
                                                            angle            = id$antAngle,
                                                            x_tag_coord      = id$antPosition[1], 
                                                            y_tag_coord      = id$antPosition[2],
                                                            x_ant_coord      = id$antPosition[1]*cos(-id$antAngle) - id$antPosition[2]*sin(-id$antAngle),
                                                            y_ant_coord      = id$antPosition[1]*sin(-id$antAngle) + id$antPosition[2]*cos(-id$antAngle),
                                                            length_px        = queen_length_px,
                                                            stringsAsFactors = F))
  }
}
for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
}

### Write capsule data to each file that requires the queen capsules
for (no_capsule_file in no_capsule_list) {
  # open tracking data which need new capsule
  tracking_data <- fmExperimentOpen(no_capsule_file) 
  ants <- tracking_data$ants
  ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=1)$length_px)
  ants[[1]]$clearCapsules()
  capsule_number <- 0   #assign capsule numbers that match the order of the looped capsule names positions
  for (capsule_name in unlist(tracking_data$antShapeTypeNames)) {
    capsule_number <- capsule_number +1
    capsule_ratios <- capsule_list[[capsule_name]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
    capsule_coords <- ant_length_px*capsule_ratios
    ants[[1]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
  }
  #tracking_data$save(no_capsule_file)
  try(ants[[1]]$identifications[[1]]$tagSize <- 1.52, silent = TRUE) # adjust tag size for the queen
  tracking_data$save(paste0(directory, substr(no_capsule_file, 1, nchar(no_capsule_file)-23),'queen_modified_test.myrmidon'))
  print(paste0("Added queen capsule to",no_capsule_file))
  #close experiment
  rm(list=(c("tracking_data")))
}

#### Next steps ####
  ### Manually apply the zones of the feeding area.

# Unfortunately the zone workspace does not allow us to have a screenshot of different times. Thus we have to modify the images used to display in the zone workspace. 
# First we copy the relevant images "frame_XXXX.png" into a new folder 
# Then we take a screenshot of the feeding that contains the food and use inkscape to overlay the new image on the original and safe it with the same original name in a new folder. 
# Last we copy all the files back into their respective folder
# the next time we open fort we should be able to see the newly created image in the zoning workspace

# get a list with the latest myrmidon file with manually oriented queens
files <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
files <- grep(files, pattern = "_feeding_DS.00", invert = FALSE, value = TRUE)
# Extract the replication numbers from the file names
rep_nums <- as.numeric(gsub(".*\\.(\\d+)$", "\\1", files))
# Identify and remove duplicates based on the file names without the replication numbers
unique_files <- files[!duplicated(gsub("\\.\\d+$", "", files), fromLast = TRUE)]

for (file in unique_files) {
  folderpath <- paste0(directory, file, "/ants/")
  files <- list.files(folderpath)
  files <- grep(files, pattern = "frame", invert = FALSE, value = TRUE)
  numbers <- as.integer(sub("^.*_(\\d+)\\.png$", "\\1", files))# Extract the numbers from the file names using regular expressions
  smallest_file <- files[which.min(numbers)]  # Find the file with the smallest number
  dest_path <- paste0(directory, "/screenshots/original_frame_pictures/", sub(".*(c[0-9]+).*", "\\1", file), "/")
  file.copy(file.path(folderpath, smallest_file), dest_path, overwrite = TRUE)
}




   ### Implement the ant pose cloner (check if scaling is needed) - Aim: get orientation and capsules for the feeding tracking from the already processed main tracking. 
   ### Extract the coordinates of the feeding zones
   ### Come up with a solution to that gives tells us which ants were feeding in the treatments and create a rule to define which colonies to include in the analyses. 



