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
# Step 3.3: Automatically orient all ants in all extrapolated data using the ant_orient express
# Step 3.4: Generate capsules for all ant
# Step 3.4: For each tracking file do the necessary post processing

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
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
# directory <- '/media/gw20248/gismo_hd3/vital/fc2/' 
setwd(directory)


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








#### 3.1 Create all base myrmidon files for extrapolated data ####

# !!! Manually copy the file base_source.myrmidon into the directory containing the tracking files !!!


# set working directory to the new directory containing the tracking data
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










#### 3.1.1 Extra step manually assign to each of the base files the corresponding data in fort ####
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
  data$save(paste0(directory, substr(paste0(data_collection[i,"colony_nr"], "_", type, "_base.myrmidon"), 1, nchar(paste0(data_collection[i,"colony_nr"], "_",type, "_base.myrmidon"))-13), 'AntsCreated.myrmidon'))
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



# if the code above runs, include the ant orientation code and the capsule code for a quick post processing, 
# then apply the ant pose cloner to synchronize feeding and main tracking files.
# (check if there is indeed no scaling required to make the capsule cloning work)
# then, run the ant orient_express!

# Information that might be included earlier: 
  
  # Steps to take first:
  # Create manually oriented base files
  # Orient 1 lanrge colony per tracking system (food flow experiment Daniel - 3 tracking systems used for main tracking and 2 tracking systems used to record feeding sessions)
  # used tracking systems Daniel: Maintracking - trojan, prideaux, guillam ; feeding sessions - guillam, esterhase 
  # choose 1 of the manually oriented colonies to define a capsule definition in fort studio for a medium sized ant and replicate the shape for all of the ants of the colony
  # apply this capsule definition for the remaining manually oriented colonies using the script below. The originals of these files have been stored as *.myrmidon.old
  

