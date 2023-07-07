# Metadata extraction file
rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### METADATA EXTRACTION FROM MYRMIDON FILES ADJUSTED FOR DANIEL ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
 
# Run this to extract ant metadata from all your myrmidon files of interest to create the metadata table required for the rest of the analysis
# It will create a file as an output that is saved into your data directory so it can be accessed in EXP1_base_analysis_DS.R

#### TO DO ####
# Currently c29 which is the trophallaxis colony is not included. 
# Maybe adjust its post processing according to the other files (i.e. the same meta data)
# After finishing this current run check the generated meta_data file and ants_to_check before continuing! !!!!!!!
# Safe ants to check somewhere. 
# continue exp1 base analysis script

#### LIBRARIES ####
{
library(FortMyrmidon) #R bindings
library(mallinfo)
library(MALDIquant)
library(dplyr)


#### Directories ####

# define user and hard drive
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


#### create the directories #### 

WORKDIR <- paste("/media/",usr, hd, "/vital/fc2",sep="")
DATADIR <- paste(WORKDIR, sep = "/")
SAVEDIR <- paste("/media/",usr, hd,"/vital/fc2/EXP1_base_analysis/EXP_summary_data",sep="") # where to save the interactions in the same structure as for Science 2018; check adrionos guide on how to copy the structure
INTDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/intermediary_analysis_steps",sep="")
BEHDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/processed_data/individual_behaviour",sep="")
SCRIPTDIR <- paste("/home/",usr,"/Documents/vital_rscripts_git",sep="") #SCRIPTDIR <- paste("/media/",usr, hd, "/vital/fc2/Documents/EXP1_base_analysis/EXP1_analysis_scripts", sep="")
BEH_FUNCTIONS <-  paste(SCRIPTDIR, "/Behavioural_Inference_DS",sep="")


#### Required ressources ####
FRAME_RATE <- 6
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # will add colony_metadata data frame to the environment so it can be accessed within this and other script
source(paste0(SCRIPTDIR,"/vital_metadata_feeders.R")) # will add metadata for the ants that were treated to the environment
source(paste(BEH_FUNCTIONS,"trajectory_extraction.R",sep="/")) # Ant Tasks to define functions
source(paste(SCRIPTDIR,"AntTasks_v082_DS.R",sep="/"))
Metadata_Exp1 <- file.path(DATADIR, paste0("Metadata_vital_", Sys.Date(), ".txt")) # define an output file path
ants_to_check <- data.frame(antID = character(), unique_ant_ID = character(), filename = character(), stringsAsFactors = FALSE)







#### List of myrmidon files ####

# create a list off all myrmidon files for which you want the worker metadata
  setwd(DATADIR)
# 2 get  list of all filenames (with path form directory) containting the capsule definitions
  add_directory <- function(filename) { # Function to add directory path to each filename
  paste0(DATADIR,"/", filename)
  }
  meta_files <- list.files()
  meta_files <- grep(meta_files, pattern = 'final', invert = FALSE, value = TRUE) 
  meta_files <- grep(meta_files, pattern = 'CapsuleDef', invert = TRUE, value = TRUE) 
  meta_files <- grep(meta_files, pattern = 'c29', invert = TRUE, value = TRUE)
  meta_files <- lapply(meta_files, add_directory)
  #meta_files <- meta_files[c(12, 17)]
}

#### Define GET METADATA FUNCTION ####
get_metadata <- function() {
  print(file)
  e <- fmExperimentOpen(file)  #open experiment
  print(paste0("Computing Ant Tasks and Zone Uses for", file))
  AntTasks <- AntTasks_DS(e=e, file=file) # COMPUTE THE ANT TASKS
  exp.Ants  <- e$ants
  # CREATE BASE FILE
  metadata <- NULL
  for (ant in exp.Ants){
    for (id in ant$identifications){
      metadata <- rbind(metadata,data.frame(colony_id        = colony_id,
                                            treatment        = treatment,
                                            treatment_simple = treatment_simple,
                                            antID            = ant$ID,
                                            tagIDdecimal     = id$tagValue,
                                            tagID            = sprintf("0x%03x", id$tagValue),
                                            identifStart     = capture.output(print(id$start)), 
                                            identifEnd       = capture.output(print(id$end)), 
                                            AntTask1perc     = AntTasks[which(AntTasks$antID==ant$ID),"AntTask"],
                                            prop_time_outside= round(AntTasks[which(AntTasks$antID==ant$ID),"prop_time_outside"],3),
                                            box              = e$spaces[[1]]$name, # tracking system
                                            stringsAsFactors = F))
      }
    }
  metadata[grepl("∞", metadata$identifEnd),"identifEnd"] <- NA
  metadata[grepl("∞", metadata$identifStart),"identifStart"] <- NA
  metadata$identifEnd <- as.POSIXct( metadata$identifEnd , format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )   #transform time in GMT
  #colony-wide metadata
  metadata$colony_size <- length(unique(metadata$antID))
  
  #empty metadata columns
  metadata$unique_ant_ID <- paste(metadata$colony_id, metadata$tagID, sep="_")
  metadata$comment     <- NA
  metadata$IsTreated   <- NA # if an ant got selected as a feeder
  metadata$IsAlive     <- NA 
  metadata$IsQueen     <- NA 
  metadata$surviv_time <- metadata$identifEnd # time of death or escape if applicanle, else end of exp
  metadata$position    <- NA # if an ant was a feeder on which position was it feeding  
  metadata$food_source <- NA # if an ant was a feeder - what food did it get "virus" or "control" (ants from control treatments will always have control and for the virus treatment half of the treated ants will have virus and the other half control)
  metadata$bead_colour <- NA # if an ant was a feeder - what color were the beads that were in its food source
  metadata$feeding_duration  <- NA # if an ant was a feeder - how long was it feeding for
  metadata$feeding_exclusion <- NA # some ants drowned themselves or similar and were no longer responsive afterwards and should thus be excluded for the Analyses "TRUE" if the ant is to be excluded, else set as "FALSE" which is the default for ants that are feeding
  metadata$status_ant <- NA
  time_treatment_start <- as.POSIXct(colony_metadata$time_treatment_start[colony_metadata$colony_id == colony_id],  format = "%Y-%m-%dT%H:%M:%OSZ", origin = "1970-01-01", tz = "GMT")       
  metadata$treatment_time <- time_treatment_start
  metadata$exp_stop_time  <- time_treatment_start+3*3600
  metadata$exp_start_time <- as.POSIXct(paste(format(time_treatment_start - 86400, "%Y-%m-%d"), "09:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  metadata$sampling_time  <- as.POSIXct(paste(format(time_treatment_start, "%Y-%m-%d"), "09:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

  # Fill empty columns: METADATA KEYS AND VALUES to be assigned to the empty columns
  list_keys <- list() 
  for (KEY in   1:length(e$metaDataKeys)) { # create a list that contains all the meta data keys from the fort-myrmidon file
    key <- names(e$metaDataKeys[KEY])
    list_keys <- c(list_keys,key)
  }
  # assign metadata values 
  for (ant in exp.Ants){
    individual  <- ant$ID
    #extract metadata info per key
    for (METADATA_KEY in list_keys) {
      #for more than one row, always last row will be picked (the relevant one if there is a timed change or the default one if there is not)
      for(ROW in 1:nrow(ant$getValues(METADATA_KEY))) {
        metadata[which(metadata$antID==ant$ID),METADATA_KEY] <- ant$getValues(METADATA_KEY)[ROW,"values"]        #assign metadata_key value when ID corresponds
      } # ROW
    } #METADATA_KEY
  }#ant

  ### check if ants in feeder data frame and ants in metadata that have IsTreated = TRUE correspond with each other!
  colony_feeders <- metadata_feeders[metadata_feeders$colony_id == colony_id, ] #subset feeder dataframe (sourced at the beginning) to only contain the colony in the current loop
  treated_ants <- metadata[metadata$IsTreated, "tagID"] 
  treated_ants_exist_in_feeders <- treated_ants %in% colony_feeders$tagID # check if treated ants exist in metadata 
  missing_ants <- treated_ants[!treated_ants_exist_in_feeders]
  print(paste0("For ", colony_id, ":"))
  
  if (length(missing_ants) > 0) {
    print("The following treated ants in metadata do not occur in metadata_feeders:")
    print(missing_ants)
    # ants_to_check <- c(ants_to_check, paste("Missing ant:", missing_ants, "in file:", file)) # break got replaced 
    for (i in 1:length(missing_ants)) {
      ants_to_check <<- rbind(ants_to_check, data.frame(antID = missing_ants[i], unique_ant_ID = paste(colony_id, missing_ants[i], sep="_") , filename = file)) 
    }
  } else {
    print(paste0("All treated ants in metadata occur in metadata_feeders ", "\U0001F44D"))
  }
  
  
  is_treated <- logical(nrow(colony_feeders)) # Create a logical vector to store the results
  missing_ants <- NULL
  
  for (i in 1:nrow(colony_feeders)) { # Iterate over each tagID in colony_feeders and check if IsTreated = TRUE in metadata 
    tagID <- colony_feeders$tagID[i]
    is_treated[i] <- any(metadata$IsTreated[metadata$tagID == tagID])
    if (!is_treated[i]) {
      missing_ants <- c(missing_ants, tagID)
    }
  }
  if (length(missing_ants) > 0) {    # Print the missing ants, if any
    cat("The following ants are missing in metadata or have IsTreated = FALSE:\n")
    print(missing_ants)
    # ants_to_check <- c(ants_to_check, paste("Missing ant:", missing_ants, "in file:", file)) # break got replaced
    for (i in 1:length(missing_ants)) {
      ants_to_check <<- rbind(ants_to_check, data.frame(antID = missing_ants[i], unique_ant_ID = paste(colony_id, missing_ants[i], sep="_") , filename = file)) 
    }
  } else {
    print(paste0("All ants in colony_feeders have IsTreated = TRUE in metadata. ", "\U0001F44D"))
  }
  
  # fill in remaining metavalues
  for (i in 1:nrow(metadata)) {  
    if (is.na(metadata$surviv_time[i])) { # set survive_time to the end of the experiment which is three hours after the treatment time. 
      metadata$surviv_time[i] <- as.POSIXct( time_treatment_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT")+3*3600
    }
    if (metadata$IsTreated[i] == TRUE) {
      if (metadata$unique_ant_ID[i] %in% ants_to_check$unique_ant_ID) {
        next  # Skip the current iteration of the loop
      }
      metadata$status_ant[i] <- "treated"
      metadata$position[i]    <- colony_feeders$position[which(colony_feeders$tagID == metadata$tagID[i])]
      metadata$food_source[i] <- colony_feeders$food[which(colony_feeders$tagID == metadata$tagID[i])]
      metadata$bead_colour[i] <- colony_feeders$beads[which(colony_feeders$tagID == metadata$tagID[i])]
      metadata$feeding_duration[i]  <- colony_feeders$total_duration[which(colony_feeders$tagID == metadata$tagID[i])]
      metadata$feeding_exclusion[i] <- colony_feeders$excluder[which(colony_feeders$tagID == metadata$tagID[i])]
    }
    if (metadata$IsTreated[i] == FALSE) {
      metadata$status_ant[i] <-"untreated"
    }
  }  
  # N_exposed in the colony
  metadata$N_treated <- sum(metadata$IsTreated) # eventually add number of feeders for each food source with and without excluder (but this can also be done at a later stage)
  # save
  if (file.exists(Metadata_Exp1)){
    write.table(metadata,file=Metadata_Exp1,append=T,col.names=F,row.names=F,quote=T,sep=",")
  }else{
    write.table(metadata,file=Metadata_Exp1,append=F,col.names=T,row.names=F,quote=T,sep=",")
  }
  gc()
  #clean up
  rm(list=ls()[which(!ls()%in%c("meta_files","colony_metadata","metadata_feeders","Metadata_Exp1","SCRIPTDIR","WORKDIR","DATADIR","AntTasks", "ants_to_check"))]) 
}




#### Run function ####

# loop it for each file... 

# for each file get some information replicate file
for (file in meta_files) {
  cat(crayon::yellow$bold("Start new colony\n"))
  print(file)
  colony_id <- sub(".*/final_(c\\d{2})\\.myrmidon", "\\1", file) # get unique colony identifier (Adriano had Rep_name)
  treatment <- colony_metadata[which (grepl(paste0(colony_id) ,colony_metadata$colony_id)), "treatment" ]# get treatment from meta data table (Adriano had Rep_treatment)
  treatment_simple <- colony_metadata[which (grepl(paste0(colony_id) ,colony_metadata$colony_id)), "treatment_simple" ] #Adriano had Rep_ts
  tracking_system_main <- colony_metadata[which (grepl(paste0(colony_id) ,colony_metadata$colony_id)), "tracking_system_main" ]# get tracking system from meta data 
  tracking_system_feeding <- colony_metadata[which (grepl(paste0(colony_id) ,colony_metadata$colony_id)), "tracking_system_feeding" ]
  time_treatment_start <- colony_metadata[which (grepl(paste0(colony_id) ,colony_metadata$colony_id)), "time_treatment_start" ]
  # check if the metadata file exists and if the colony has already been recorded: 
  if(file.exists(file.path(Metadata_Exp1))) { 
    metadata_present <- read.table(file.path(Metadata_Exp1),header=T,stringsAsFactors = F, sep=",")
    if (colony_id %in% unique(metadata_present$colony_id)) {
      print(paste0(colony_id," already present in metadata, skip"))
    } else {
      get_metadata() # if not present run the get_metadata() function 
      ants_to_check <- ants_to_check
    }
  } else {
    get_metadata() # if the file does not exist yet run the get_metadata() function
    ants_to_check <- ants_to_check
  }
} 

# ants_to_check <- c("c12_0x644",  "c17_0x25e", "c17_0x768", "c18_0x62c", "c25_0x103", "c25_0x013")


#### TEST #### 
# ants_to_check <- data.frame(antID = character(), unique_ant_ID = character(), filename = character(), stringsAsFactors = FALSE) 
# colony_ids <- c("c01", "c02")
# 
# test_function <- function() {
#   for (i in 1:3) {
#     ants_to_check <<- rbind(ants_to_check, data.frame(antID = paste("missing_ant", i,sep="_"), unique_ant_ID = paste(colony_id, "missing_ant", i, sep="_") , filename = colony_id)) 
#   }
# }
#       
# for (colony_id in colony_ids ) {
#   cat(crayon::yellow$bold("Start new colony\n"))
#   print(colony)
#   test_function()
# }
# 
# i <- 1
# colony_id <- colony_ids[1]





# > ants_to_check 
# antID unique_ant_ID                                           filename
# 1 0x644     c12_0x644 /media/gw20248/DISK_B/vital/fc2/final_c12.myrmidon
# 2 0x25e     c16_0x25e /media/gw20248/DISK_B/vital/fc2/final_c16.myrmidon
# 3 0x768     c16_0x768 /media/gw20248/DISK_B/vital/fc2/final_c16.myrmidon
# 4 0x62c     c18_0x62c /media/gw20248/DISK_B/vital/fc2/final_c18.myrmidon
# 5 0x103     c25_0x103 /media/gw20248/DISK_B/vital/fc2/final_c25.myrmidon
# 6 0x013     c25_0x013 /media/gw20248/DISK_B/vital/fc2/final_c25.myrmidon
