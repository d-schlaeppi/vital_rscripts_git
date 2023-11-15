rm(list = ls())
gc()
Sys.sleep(3)
mallinfo::malloc.trim(0L) # forced memory cleaning after every gc() | # install in terminal: > sudo apt-get install libtcmalloc-minimal4
                          # and in R: > install.packages("mallinfo", repos = "http://www.rforge.net/") | and run the following after gc(): > mallinfo::malloc.trim(0L)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### VITAL BASE ANALYSES Daniel  ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### READ ME ####
# This is Script is based on Adrianos, Tom & Nathalies Exp1_base_analysis script and got adjusted to the data structure and needs of Daniel
# Probably the whole lab contributed and if you use this you are probably aware of adrianos guides: https://github.com/AdrianoWanderlingh/Ant_Tracking/tree/main/Scripts

#### Index | Overview ####

#### To Do's ####

# When Calculating the ant tasks it might be better to do it over the full 24h of the pre treatment network?! Check with the script 19_facet_net 
# work through the whole script line by line
# get the exact end time for each of the experiments
# UPDATE all parts of the script that try to access BODYLENGTH_FILE by changing the code in a way that it refers to size in pixels based on colonies instead of tracking systems
# introduce all the fixes suggested via email:
# make sure everthing can run on laptop and uni computer (starting if statement)
# get a version of the myrmidon files with Nathalies 2018 capsule definitions for normal interactions. 
# What ever I do with the interaction network, make sure to get the same things with the trophallaxis networks
# Make sure to change frame rate to 6 (last in the last updated exp1 base analysis script)

# downstream: use 13_network_analysis.R and 14_summarise_interactions.R in the version linda sent around and include her fixes 
# In another script, 12_simulate_transmission.R (which is part of the 2018 Social Plasticity code and can be found here on Adriano's github), we discovered that the output of the Surv() function has changed the name for one of its columns, which is being referred to in lines 261 & 265: "rmean" (new) instead of "*rmean" (old); due to a library update perhaps? Please be aware of this since it does not give an error but just NAs for those specific cases.





# To correct when going through the code:

### From Adriano
#   in 13_network_analysis, the following block of code
# ###build NETWORK
# if (!grepl("survival",data_path)){
#   net <- graph.data.frame(interactions[c("Tag1","Tag2")],directed=F,vertices=actors)
#   ### add edge weights
#   if (edge_weights=="number"){
#     E(net)$weight <- interactions[,"number"]
#   }else if (edge_weights=="duration"){
#     E(net)$weight <- interactions[,"duration_min"]        
#   }
#   should be changed from "E(net)$weight <- interactions[,"number"]" to "E(net)$weight <- interactions[,"N"]" (the variable is called N).
#   
#   Also, when the code is run for the "grooming interactions", it should only be run for the "observed" folders (in theory, you should not have generated any randomisation for the grooming anyway), so to make sure your analysis of the network properties for grooming takes only "observed as input, you can place these lines before "for (input_folder in input_folders){":
# 
# #if grooming, perform analysis only on the observed cases
# if (grepl("grooming",input_path)) {
#   input_folders <- grep("observed", input_folders, value=TRUE)
# }
#
# for (input_folder in input_folders){



# P.S. sorry, already found another bug (I didn't test very thoroughly on Friday afternoon...) -  in 13_network_analysis add "edge_weights" to "to_keep". 
#to_keep <- c(ls(),"to_keep","input_folder","network_files","options","option","summary_collective","summary_individual","outputfolder","network_file","queenid", "edge_weights")






#### LIBRARIES ####
{
library(data.table)
library(lubridate)
library(pals)
library(dplyr)
library(adehabitatHR)
library(adehabitatLT)
library(changepoint)
library(e1071)
library(igraph)
library(gtools)
library(Rcpp)
library(survival)
}
#### starting parameters & directories ####
{
USER <- "2A13_Office_Daniel"  # Replace with the desired USER option: Nath_office, 2A13_Office_Adriano, 2A13_Office_Daniel, AEL-laptop
HD <- "/DISK_B" # One of the harddrives with the vital extrapolated data: possible values > Nathalies hd "/DISK_B" ; Daniels hds >  "/gismo_hd5" or  "/gismo_hd2" | just make sure to put the name of your HD in here
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
  hd <- NULL
  hd <- HD
  if (!is.null(hd)) {print(hd)} else {print("define new hd if necessary")}
  assign("hd", hd, envir = .GlobalEnv)  # Assign to the global environment
  assign("usr", usr, envir = .GlobalEnv)  # Assign to the global environment
}
setUserAndHD(USER, HD)

# set directories

WORKDIR <- paste("/media/",usr, hd, "/vital/fc2",sep="") 
DATADIR <- paste(WORKDIR, sep = "/") # working directory and directory where the tracking data is saved (was not the same for Adriano, but for me the same
SAVEDIR <- paste("/media/",usr, hd,"/vital/fc2/vital_experiment/summary_data",sep="") # where to save the interactions - I added a folder summary data to  vital_experiment within the science folder structure
INTDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/intermediary_analysis_steps",sep="") # remember to use the same folder structure as as for Science 2018
BEHDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/processed_data/individual_behaviour",sep="")
SCRIPTDIR <- paste("/home/",usr,"/Documents/vital_rscripts_git",sep="") # place where the needed r scripts are stored
BEH_FUNCTIONS <-  paste(SCRIPTDIR, "/Trophallaxis_Classifier",sep="") # check if and where that is needed...?

### 
# Note to myself keep updating the list of scripts used as part of this analysis pipeline.
# if not specified differently, the script version used in the pipeline is in the top git hub folder or in the folder with trophallacis_classifier. 
# previous version is either in the _old folder or in a folder referring to another persion e.g Adriano or Linda.
# At some point make a cleaner github version of all this... 

### source and set additional data (meta), functions, parameters and additional scripts ####
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # will add colony_metadata data frame to the environment so it can be accessed within this script (in my case containing bodylenght information)
metadata <- read.table(paste(DATADIR, "/Metadata_vital_2023-07-04.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",") # contains return time for each colony and end of experiment time!
metadata$status_ant <- NA
metadata$status_ant <-ifelse(metadata$IsTreated==TRUE,"treated","untreated")
metadata_colonies <- colony_metadata
  
str(metadata_colonies)

print("Loading functions and libraries...")
source(paste(SCRIPTDIR, "SpaceUse_v082.R", sep = "/")) # SPACE_USE
source(paste(SCRIPTDIR, "NetworkProperties_v082.R", sep = "/")) # NET_properties collective + individual
# suppressMessages(source(paste(SCRIPTDIR, "BEH_libraries_DS.R", sep = "/")))
source(paste(BEH_FUNCTIONS, "interaction_detection_ori.R", sep = "/"))
source(paste(BEH_FUNCTIONS, "trajectory_extraction.R", sep = "/"))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # required by merge_interactions.cpp
sourceCpp(paste(BEH_FUNCTIONS, "merge_interactions.cpp", sep = "/"))


# PARAMETERS
TimeWind <- 3600 ## 1h in seconds
N_DECIMALS <- 3
MAX_INTERACTION_GAP <- 10 # fmSecond(10) # the maximum gap in tracking before cutting the interactions in two different object.  Used in ComputeGraph function
MAX_DIST_GAP_MM <- 0.6
FRAME_RATE <- 6 # !!!! Usually we record at 8 frames
interactions_of_interest <- list(c("head", "body")) ### if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))


# Time dictionary DS (when time windows of interest:
# treatment was at 11:00 so we are interested what happens after and what happens always 24 h before. 
Time_dictionary <- data.frame(time_hours = c(-24,0), time_of_day = c(11,11)) # time_hours = time since exposure (negatives pre-exposure, positive post exposure), time_of_day = time hour slot of the day
Time_dictionary$period <- ifelse(Time_dictionary$time_hours<0, "pre", "post")

# Paths to save the files (Location according to Science2018 pipeline)
SPACE_USE     <-  file.path(BEHDIR,"pre_vs_post_treatment","individual_behavioural_data.txt")
SPACE_USE_PRE <-  file.path(BEHDIR,"pre_treatment","network_position_vs_time_outside.dat")


### FLAGS                                                                                                      # check if we still use the flags because the whole part on if = ture was #'d out. 
RUN_INTERACT     <- TRUE
RUN_SPACEUSE     <- TRUE
RUN_NETWORKS     <- TRUE 
warning(paste("RUN_INTERACT is set to:",RUN_INTERACT,
              "\nRUN_SPACEUSE is set to:",RUN_SPACEUSE,
              "\nRUN_NETWORKS is set to:",RUN_NETWORKS,sep="\t"))

# list all files for which the analysis should be run.
files_list <- list.files(DATADIR, pattern="CapsuleDef2018.myrmidon")
files_list <- files_list[which(!grepl("c29",files_list))]


### RUNNING TIME
loop_start_time <- Sys.time()

### create to_keep which will contain variables that should not be deleted when clearing memory between runs
to_keep <- c(ls(), c("to_keep"))
}




#### Loop over each colony to calculated ... ... ...####

for (REP.file in files_list) {                                                                        #for each colony...
  REP.file <- files_list[1]                                                                           # TO DELETE ONCE WE GET TO THE END OF THE LOOP!!!!
    colony_id <- unlist(strsplit(REP.file,split="_"))[grepl("c",unlist(strsplit(REP.file,split="_")))]  # extract colony id
    REP_colony <- colony_id                                                                             # just a safety measure because adriano used REP_colony as colony identifier and I might have missed an instance or it might be used in another script |  at some point try deleting to see if it can be deleted
    cat(paste("########################################\n",basename(REP.file)),sep="")                  # just a message indicating new colony. 
    Period_dataframe <- NULL                                                                          # checking time correspondances ???
    
    # start fresh data frames
    NetworkProp_collective <- data.frame()
    NetworkProp_collective_hour <- data.frame()
    NetworkProp_individual <- data.frame()
    NetworkProp_individual_hour <- data.frame()
    Interactions_REP_TREAT <- data.frame()
    
    # open experiment
    e <- fmExperimentOpen(paste(DATADIR,REP.file,sep="/"))
    print(paste0("Processing ", basename(REP.file)))
    
    # create variable containing mean worker size in mm and pixcel (used downstream)
    body_lengths <- NULL
    body_lengths <- metadata_colonies[metadata_colonies$colony_id == colony_id, c("colony_id", "mean_ant_lenght_px", "mean_ant_lenght_mm", "tracking_system_main")]
    colnames(body_lengths) <- c("colony_id", "body_length_px", "body_length_mm", "tracking_system")
    
    # get base info for this colony
    # REP_name <- unlist(strsplit(REP.file,split="/"))[length(unlist(strsplit(REP.file,split="/")))]                  to delete if everyting else works 
    treatment <- colony_metadata[which(colony_metadata$colony_id==colony_id),"treatment_simple"]
    tracking_system <-colony_metadata[which(colony_metadata$colony_id==colony_id),"tracking_system_main"]
    REP_TREAT <- paste(treatment, tracking_system, colony_id, sep="_")                                                # see whether we can run it without this.

    # here a big chunk was deleted 

    #### PERFORM THE FULL INTERACTION AND NETWORK ANALYSIS OUTSIDE OF THE HOURLY LOOP ####

    for (PERIOD in c("pre","post")) {                                                                                 # for the pre and the post treatment time block 
      PERIOD <-"pre"                                                                                                  # temporary - TO DELETE ONCE WE GET TO THE END OF THE LOOP!!!!
      Time_dictionary_PERIOD <- Time_dictionary[which(Time_dictionary$period==PERIOD),]
      colony            <- colony_id
      period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
      treatment_code    <- treatment
      # Select the Period of interest. (Usually it is the full 24 h, in my case it was 3h only because we froze the colonies to have snapshot of the food spread)
      cat(paste("#######","period:", PERIOD, sep= " "))
      treatment_start <-  as.POSIXct(metadata_colonies[    which   (   metadata_colonies$colony_id == colony_id   )     ,   "time_treatment_start"     ], format="%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
      ###
      ### Here maybe insert something to get the exact end of the experiment instead of the 3h block (but should make minimal difference) (And check if 24h is required for task definition)
      ###
      if (PERIOD=="pre"){
         time_start_GMT <- treatment_start - 24*3600
      }else{
        time_start_GMT <- treatment_start
        }
      time_end_GMT <-  time_start_GMT + 1*3600      # 3*3600    --!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # warning("using tiny time window for testing")

      time_start <- fmTimeCreate(time_start_GMT)       #get time into the right format
      time_end  <- fmTimeCreate(time_end_GMT)
      
      if (RUN_INTERACT) {
        INTERACT_loop_start_time <- Sys.time()

        INTERACTIONS_FULL   <-  file.path(INTDIR,"full_interaction_lists",period_code,"observed",            # DESTINATION FOLDER for the interactions 
                                 paste(colony,treatment_code,period_code,"interactions.txt",sep="_"))
   
        
     
        
        
        
      #### GET INTERACTIONS ####
      # Adriano had an additional thingy here to induce a skip should a file already be present (but I assume the default will be an overwrite)
      Interactions <- compute_Interactions(e = e, start = time_start, end = time_end, max_time_gap = MAX_INTERACTION_GAP)    # function in "NetworkProperties_v082.R" -  # INTERACTIONS IN THIS FUNCTION ARE CALCULATED ACCORDING TO STROEYMEYT ET AL, SCIENCE 2018

        
        
        
        
        
      # Remove interactions involving dead ants
      dead_by_REP <- metadata[which(metadata$IsAlive==FALSE & metadata$REP_treat==REP_TREAT),"antID"]
      Interactions <- Interactions[!(Interactions$Tag1 %in% dead_by_REP | Interactions$Tag2 %in% dead_by_REP), ]
      
      # create time vars
      Interactions$time_hours   <- NA
      Interactions$time_of_day  <- NA
      # TIME_HOURS zero is the moment of exposed ants return
      # warning("interaction binning loop should be fixed as done later for the SpaceUse, see TIME_HOURS ")
      for (TIME_HOURS in Time_dictionary_PERIOD$time_hours) { ## increments by 3 hours for duration of experiment
        # TIME_OF_DAY
        TIME_OF_DAY <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "time_of_day"]
        # PERIOD_DETAIL <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "period_detail"]
        # PERIOD_CIRCADIAN <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "period_circadian"]
        
        # if (TIME_HOURS==-24){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "pre_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   
        # }else if(TIME_HOURS==-21){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "pre_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   chunk_time_start <- chunk_time_start + 3*3600 
        #   
        # }else if(TIME_HOURS==-18){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "pre_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   chunk_time_start <- chunk_time_start + 6*3600 
        #   
        # }else if(TIME_HOURS==-12){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "middle_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   
        # }else if(TIME_HOURS==-9){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "middle_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   chunk_time_start <- chunk_time_start + 3*3600 
        #   
        # }else if (TIME_HOURS==-6){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "middle_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   chunk_time_start <- chunk_time_start + 6*3600 
        #   
        # }else if(TIME_HOURS==0){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "post_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   
        # }else if(TIME_HOURS==3){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "post_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   chunk_time_start <- chunk_time_start + 3*3600    
        #   
        # }else if(TIME_HOURS==6){
        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "post_start"     ]
        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
        #   chunk_time_start <- chunk_time_start + 6*3600 
        # }
        ###works for Daniel only!
        chunk_time_start <- time_start_GMT
        chunk_time_stop <- time_end_GMT
        
        #time windows
        From_TIME_HOURS <- as.numeric(chunk_time_start,3)
        To_TIME_HOURS <- as.numeric(chunk_time_stop,3)
        
        #adding time labels
        Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_hours"]   <- TIME_HOURS
        Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_of_day"]  <- TIME_OF_DAY
        # Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"period_detail"]   <- PERIOD_DETAIL
        # Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"period_circadian"]  <- PERIOD_CIRCADIAN
      }
      # LS: NAs in interactions time_hours and time_of_day because we also calculate interactions for the time between the chunks (which are not in the dictionary)
      unique(Interactions$time_hours)
      # Add metadata info
      Interactions <- cbind(data.frame(
        period = PERIOD,Interactions, REP_treat = REP_TREAT, colony=colony, treatment=treatment_code,
        treatment_detailed=metadata_colonies[which(metadata_colonies$colony_id==colony_id),"treatment"],
        
        stringsAsFactors = F
      ))
      
      # extras not present in Science files: (added at the end of the file output)
      #"REP_treat","period","ant1.zones","ant2.zones","duration"
      Interactions <- Interactions[,c("Tag1","Tag2","Startframe","Stopframe","Starttime","Stoptime","Box","Xcoor1","Ycoor1","Angle1","Xcoor2","Ycoor2","Angle2","Direction","Detections","time_hours","time_of_day","colony","treatment","treatment_detailed","REP_treat","period","ant1.zones","ant2.zones","duration")]
      # remove extra -3h gap leftovers (few mins)
      Interactions <- Interactions[which(Interactions$time_hours%in%Time_dictionary_PERIOD$time_hours),] 
      ## Interactions save (saved INSIDE the Network_analysis folder)
      #if (file.exists(INTERACTIONS_FULL)) {
      #  write.table(Interactions, file = INTERACTIONS_FULL, append = T, col.names = F, row.names = F, quote = F, sep = ",")
      #} else {
      # LS: add this to save one file each for PERIOD_DETAIL pre_1 and pre_2 bc I have two pre periods:
      # for (PERIOD_DETAIL in unique(Interactions$period_detail)){
      write.table(Interactions, file = INTERACTIONS_FULL, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
      # }

      #}
      
      ### split output into bins
        ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
        print(paste0("split files into 3-hours bins"))
        
        for (TIME_HOURS in unique(Interactions$time_hours)) {
          #labels for subsets
          TH <- paste0("TH",unique(Interactions[which(Interactions$time_hours==TIME_HOURS),"time_hours"]))
          TD <- paste0("TD",unique(Interactions[which(Interactions$time_hours==TIME_HOURS),"time_of_day"]))
          cat("\rTIME_HOURS", TH,"TIME_OF_DAY", TD,"PERIOD",PERIOD)
          
          INTERACTIONS_BINNED <-  file.path(INTDIR,"binned_interaction_lists",period_code,"observed", 
                                            paste(colony,treatment_code,period_code,TH,TD,"interactions.txt",sep="_"))
                    #save object by TH and TD
          write.table(Interactions[which(Interactions$time_hours==TIME_HOURS),], file = INTERACTIONS_BINNED, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
          }
     # } #TEMP
      INTERACT_loop_end_time <- Sys.time()
      print(paste("Interactions 3h chunk took ", round(as.numeric(difftime(INTERACT_loop_end_time, INTERACT_loop_start_time, units = "secs")),1), " sec to complete"))
      }# RUN_INTERACT
      
      if (file.exists(SPACE_USE)) {
        print("file exists")
        spaceUse_done <- read.table(SPACE_USE, header = T, stringsAsFactors = F, sep = "")
        spaceUse_done <- spaceUse_done[which(spaceUse_done$colony==colony),"time_hours"]
      }else{
        spaceUse_done <- c() # placeholder
      }
        
      
      if (RUN_SPACEUSE) {
        print("Computing SpaceUse based on full time-window pre AND post exposure")
        # warning("the Space Use Script can be substantially sped up by computing trajectories beforehand, then cutting the result by 3h chunks and performing operations on them.
        #         \n -This is because trajectories stitching is slow so may be better to perform it once. see inspiration from 8_process_trajectory_files.R
        #         \n -This will require saving the summary stats per ant in a safe place (possibly the trajectories summary as now but with clearing of the vars at each time cycle)
        #         \n -parallelisation of ant computation greatly improve speed. Before parallelisation of SpaceUse calculations: 2.1 sec per ant per 3h chunk (70h), after 0.21 sec with 10 cores. this includes only the summary stats computations (traj stitching excluded)")
        SpaceUse_loop_start_time <- Sys.time()
        # TIME_HOURS zero is the moment of exposed ants return
        for (TIME_HOURS in Time_dictionary_PERIOD$time_hours) { ## increments by 3 hours for duration of experiment
          # TIME_HOURS = -9 # temp
          # TIME_OF_DAY
          TIME_OF_DAY <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "time_of_day"]
            # has the time_hour been computed already? 
            if (TIME_HOURS %in% spaceUse_done) {
              cat("\rDONE UP TO: TIME_HOURS", TIME_HOURS,"TIME_OF_DAY", TIME_OF_DAY,"PERIOD",PERIOD,"| SKIP >>")
                      }else{
                        print(paste("TIME_HOURS", TIME_HOURS,"TIME_OF_DAY", TIME_OF_DAY,"PERIOD",PERIOD,sep=" "))
                        # PERIOD_DETAIL <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "period_detail"]
                        # PERIOD_CIRCADIAN <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "period_circadian"]
          #time windows LINDA REDEFINE AS ABOVE WITH THE NASTY SERIES OF IF STATEMENTS
          # chunk_time_start <- fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind)) # time_stop minus 48 hours plus incremental time
          # chunk_time_stop  <- fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 21) * TimeWind)) # time_stop minus 45 hours plus incremental time
                        chunk_time_start <- time_start_GMT
                        chunk_time_stop <- time_end_GMT
                        
                        # if (TIME_HOURS==-24){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "pre_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        #   
                        # }else if(TIME_HOURS==-21){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "pre_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        #   chunk_time_start <- chunk_time_start + 3*3600 
                        # 
                        # }else if(TIME_HOURS==-18){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "pre_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        #   chunk_time_start <- chunk_time_start + 6*3600 
                        # 
                        # }else if(TIME_HOURS==-12){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "middle_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        # 
                        # }else if(TIME_HOURS==-9){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "middle_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        #   chunk_time_start <- chunk_time_start + 3*3600 
                        # 
                        # }else if (TIME_HOURS==-6){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "middle_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        #   chunk_time_start <- chunk_time_start + 6*3600 
                        # 
                        # }else if(TIME_HOURS==0){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "post_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        # 
                        # }else if(TIME_HOURS==3){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "post_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        #   chunk_time_start <- chunk_time_start + 3*3600
                        # 
                        # }else if(TIME_HOURS==6){
                        #   chunk_time_start <- metadata_colonies[    which   (     grepl(   REP_TREAT ,  metadata_colonies$treat_TS_colony   )      )     ,   "post_start"     ]
                        #   chunk_time_start <- as.POSIXct(chunk_time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
                        #   chunk_time_start <- chunk_time_start + 6*3600
                        # }
                        # chunk_time_stop <- chunk_time_start + 3*3600     
                        chunk_time_start <- fmTimeCreate(chunk_time_start)
                        chunk_time_stop  <- fmTimeCreate(chunk_time_stop)
                        
                        
          ##TEMP TIME STOP
          # warning("using tiny time window for testing")
          #time_stop_h <-  fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind + 40*60) ) 
          
         
          ############ SPACE USE ##########################
          # SPACE USAGE
          SpaceUsage <- SpaceUse(e = e, start = chunk_time_start, end = chunk_time_stop)
          
          SpaceUsage <- dplyr::left_join(SpaceUsage, metadata[which(metadata$colony_id==colony_id),c("status_ant", "antID","IsAlive")], by = "antID") # Apply left_join dplyr function
          colnames(SpaceUsage)[which(colnames(SpaceUsage)=="antID")] <- "tag"
          colnames(SpaceUsage)[which(colnames(SpaceUsage)=="status_ant")] <- "status"
          
          #remove dead ants
          SpaceUsage <- SpaceUsage[which(SpaceUsage$IsAlive==T),]
          SpaceUsage$IsAlive <- NULL

          # Add metadata info
          SpaceUsage <- cbind(data.frame(
            colony = colony,
            colony_size = unique(metadata[which(metadata$colony_id==colony_id),"colony_size"]),
            treatment = treatment_code,
            treatment_detailed=metadata_colonies[which(metadata_colonies$colony_id==colony_id),"treatment"],
            age = NA,
            period = PERIOD,
            # period_detail = PERIOD_DETAIL,
            # period_circadian = PERIOD_CIRCADIAN,
            time_hours = TIME_HOURS,
            time_of_day = TIME_OF_DAY,
            SpaceUsage,
            REP_treat = REP_TREAT,
            stringsAsFactors = F
          ))

            ## Space Use save (saved INSIDE the Network_analysis folder)
            if (file.exists(SPACE_USE)) {
              # pre_vs_post_treatment/individual_behavioural_data.txt
              write.table(SpaceUsage, file = SPACE_USE, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
              # pre_treatment/network_position_vs_time_outside.dat
              write.table(SpaceUsage[which(SpaceUsage$period=="pre"),], file = SPACE_USE_PRE, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
            } else {
              # pre_vs_post_treatment/individual_behavioural_data.txt
              write.table(SpaceUsage, file = SPACE_USE, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
              # pre_treatment/network_position_vs_time_outside.dat
              write.table(SpaceUsage[which(SpaceUsage$period=="pre"),], file = SPACE_USE_PRE, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
              
            }
          } # end if TIME_HOURS not in spaceUse_done -> SpaceUse()
        } # end for all the TIME_HOURS loop
        SpaceUse_loop_end_time <- Sys.time()
        print(paste("spaceUse 3h chunk took ", as.numeric(difftime(SpaceUse_loop_end_time, SpaceUse_loop_start_time, units = "secs")), " sec to complete"))
        } # end if RUN_SPACEUSE
      } # end for PERIOD
   # }
    
  
    #### ADD EXTRA REPLICATE INFORMATIONS
  # LS: I changed the name size_treat to time_treat ####
  if (RUN_NETWORKS) {
    # add status (large, small) info to NetworkProp_collective
    NetworkProp_collective <- dplyr::left_join(NetworkProp_collective, unique(metadata[c("time_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function

    # add status (large, small) info to NetworkProp_individual
    NetworkProp_individual <- dplyr::left_join(NetworkProp_individual, unique(metadata[c("time_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function

    ########################################
    ##### SAVE FILES IN FOLDER #############
  
    ## Network properties Collective save (saved INSIDE the Network_analysis folder)
    if (file.exists(NET_properties_collective)) {
      write.table(NetworkProp_collective, file = NET_properties_collective, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    } else {
      write.table(NetworkProp_collective, file = NET_properties_collective, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
    }

    ## Network properties Individual save (saved INSIDE the Network_analysis folder)
    if (file.exists(NET_properties_individual)) {
      write.table(NetworkProp_individual, file = NET_properties_individual, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    } else {
      write.table(NetworkProp_individual, file = NET_properties_individual, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
    }
  } # end of RUN NETWORKS 

    # start fresh
    NetworkProp_collective <- data.frame()
    NetworkProp_individual <- data.frame()
    Interactions           <- data.frame()
    SpaceUsage             <- data.frame()

    # cleaning
    rm(list = ls()[which(!ls() %in% to_keep)])
    gc()
    mallinfo::malloc.trim(0L)
} #end of the very start to loop through the replicates

loop_end_time <- Sys.time()
print(paste("loop took ", as.numeric(difftime(loop_end_time, loop_start_time, units = "mins")), " minutes to complete"))
