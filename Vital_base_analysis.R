rm(list = ls())
gc()
Sys.sleep(3)
mallinfo::malloc.trim(0L) # forced memory cleaning after every gc() | # install in terminal: > sudo apt-get install libtcmalloc-minimal4
                          # and in R: > install.packages("mallinfo", repos = "http://www.rforge.net/") | and run the following after gc(): > mallinfo::malloc.trim(0L)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### VITAL BASE ANALYSIS Daniel  ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### READ ME ####
# This is Script is based on Adriano's, Tom's & Nathalie's Exp1_base_analysis script and got adjusted to the data structure and needs of Daniel
# Space use and interactions (as defined in Stroeymeyt et al 2018) will be calculated for the 3 hour post treatment window and the corresponding 3h pre exposure time window 24 h earlier
# results will be saved in folders based on the Stroeymeyt 2018 pipeline folder structure.


#### Index | Overview ####
# 1. Libraries
# 2. Starting parameters & directories
# 3. Loop over each colony to calculate basic interactions and space use
# 4. Define pre and post treatment time periods
# 5. Get basic interactions
# 6. Space use



#### 1. LIBRARIES ####
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
library(MALDIquant)
}
#### 2. Starting parameters & directories ####
{
USER <- "AEL-laptop"  # Replace with the desired USER option: Nath_office, 2A13_Office_Adriano, 2A13_Office_Daniel, AEL-laptop
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

# directories
WORKDIR <- paste("/media/",usr, hd, "/vital/fc2",sep="") 
DATADIR <- paste(WORKDIR, sep = "/") # working directory and directory where the tracking data is saved (was not the same for Adriano, but for me the same
SAVEDIR <- paste("/media/",usr, hd,"/vital/fc2/vital_experiment/summary_data",sep="") # where to save the interactions - I added a folder summary data to  vital_experiment within the science folder structure
INTDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/intermediary_analysis_steps",sep="") # remember to use the same folder structure as as for Science 2018
BEHDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/processed_data/individual_behaviour",sep="")
SCRIPTDIR <- paste("/home/",usr,"/Documents/vital_rscripts_git",sep="") # place where the needed r scripts are stored
BEH_FUNCTIONS <-  paste(SCRIPTDIR, "/Trophallaxis_Classifier",sep="") 

# source and set additional data (meta), functions, parameters and additional scripts
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # will add colony_metadata data frame to the environment so it can be accessed within this script (in my case containing bodylenght information)
metadata <- read.table(paste(DATADIR, "/Metadata_vital_2023-07-04.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",") # contains return time for each colony and end of experiment time!
metadata$status_ant <- NA
metadata$status_ant <-ifelse(metadata$IsTreated==TRUE,"treated","untreated")
metadata_colonies <- colony_metadata

print("Loading functions and libraries...")
source(paste(SCRIPTDIR, "SpaceUse_v082.R", sep = "/")) # SPACE_USE
source(paste(SCRIPTDIR, "NetworkProperties_v082.R", sep = "/")) # NET_properties collective + individual
# suppressMessages(source(paste(SCRIPTDIR, "BEH_libraries_DS.R", sep = "/")))                                    # a bunch of libraries from adriano, unclear if all of them are needed. 
source(paste(BEH_FUNCTIONS, "interaction_detection.R", sep = "/"))
source(paste(BEH_FUNCTIONS, "trajectory_extraction.R", sep = "/"))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # required by merge_interactions.cpp
sourceCpp(paste(BEH_FUNCTIONS, "merge_interactions.cpp", sep = "/"))


# PARAMETERS
TimeWind <- 3600                                         # 1h in seconds
N_DECIMALS <- 3
MAX_INTERACTION_GAP <- 10                                # fmSecond(10) # the maximum gap in tracking before cutting the interactions in two different object.  Used in ComputeGraph function
MAX_DIST_GAP_MM <- 0.6
FRAME_RATE <- 6                                          # Usually we record at 8 frames
interactions_of_interest <- list(c("head", "body"))      # if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))


# Time dictionary Daniel | treatment always around 11:00 -> we are interested what happens in the 3h window after at what happens at the same time 24h earlier 
Time_dictionary <- data.frame(time_hours = c(-24,0), time_of_day = c(11,11)) # time_hours = time since exposure (negatives pre-exposure, positive post exposure), time_of_day = time hour slot of the day
Time_dictionary$period <- ifelse(Time_dictionary$time_hours<0, "pre", "post")

# Paths to save the files (Location according to Science2018 pipeline)
SPACE_USE     <-  file.path(BEHDIR,"pre_vs_post_treatment","individual_behavioural_data.txt")
SPACE_USE_PRE <-  file.path(BEHDIR,"pre_treatment","network_position_vs_time_outside.dat")


# FLAGS                                                                                                     
RUN_INTERACT     <- FALSE 
RUN_SPACEUSE     <- TRUE
RUN_NETWORKS     <- FALSE  # this part of the code seems is not up to date and behind #s
warning(paste("RUN_INTERACT is set to:",RUN_INTERACT,
              "\nRUN_SPACEUSE is set to:",RUN_SPACEUSE,
              "\nRUN_NETWORKS is set to:",RUN_NETWORKS,sep="\t"))

# list all files for which the analysis should be run
files_list <- list.files(DATADIR, pattern="CapsuleDef2018.myrmidon")     
files_list <- files_list[which(!grepl("c29",files_list))]


# RUN TIME
loop_start_time <- Sys.time()
# create to_keep which will contain variables that should not be deleted when clearing memory between runs
to_keep <- c(ls(), c("to_keep"))
}


#### 3. Loop over each colony to calculate basic interactions and space use ####

  for (REP.file in files_list) {                                                                        
      # REP.file <- files_list[1]                                                                         # TO DELETE ONCE WE GET TO THE END OF THE LOOP!!!!
      colony_id <- unlist(strsplit(REP.file,split="_"))[grepl("c",unlist(strsplit(REP.file,split="_")))]  # extract colony id
      REP_colony <- colony_id                                                                             # just a safety measure because Adriano used REP_colony as colony identifier and I might have missed an instance or it might be used in another script |  at some point try deleting to see if it can be deleted
      cat(paste("########################################\n",basename(REP.file)),sep="")                  # just a message indicating new colony. 
      Period_dataframe <- NULL                                                                            
      
      # start fresh data frames
      NetworkProp_collective <- data.frame()
      NetworkProp_individual <- data.frame()
      Interactions           <- data.frame()
      SpaceUsage             <- data.frame()
      
      # open experiment
      e <- fmExperimentOpen(paste(DATADIR,REP.file,sep="/"))
      print(paste0("Processing ", basename(REP.file)))
      
      # create variable containing mean worker size in mm and pixcel (used downstream)
      body_lengths <- NULL
      body_lengths <- metadata_colonies[metadata_colonies$colony_id == colony_id, c("colony_id", "mean_ant_lenght_px", "mean_ant_lenght_mm", "tracking_system_main")]
      colnames(body_lengths) <- c("colony_id", "body_length_px", "body_length_mm", "tracking_system")
      
      # get base info for this colony
      treatment <- colony_metadata[which(colony_metadata$colony_id==colony_id),"treatment_simple"]
      tracking_system <-colony_metadata[which(colony_metadata$colony_id==colony_id),"tracking_system_main"]
      REP_TREAT <- paste(treatment, tracking_system, colony_id, sep="_")                                                
      
      #### 4. Define pre and post treatment time periods ####
      for (PERIOD in c("pre","post")) {                                                                                
        
        # PERIOD <-"pre"                                                                                                  # temporary - TO DELETE ONCE WE GET TO THE END OF THE LOOP!!!!
        Time_dictionary_PERIOD <- Time_dictionary[which(Time_dictionary$period==PERIOD),]
        colony            <- colony_id
        period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
        treatment_code    <- treatment
        
        # Select the Period of interest. (Usually it is the full 24 h, in my case it was 3h only because we froze the colonies to have snapshot of the food spread)
        cat(paste("#######","period:", PERIOD, sep= " "))
        treatment_start <-  as.POSIXct(metadata_colonies[    which   (   metadata_colonies$colony_id == colony_id   )     ,   "time_treatment_start"     ], format="%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  
        if (PERIOD=="pre"){
           time_start_GMT <- treatment_start - 24*3600
        }else{
          time_start_GMT <- treatment_start
        }
        
        # time_end_GMT <-  time_start_GMT + 1*3600      # 3*3600    (3h window)
        # warning("using tiny time window for testing")
        # new version for time end with the exact stop time of the experiment (still around 3 hours but sometimes a couple of minutes longer or shorter)
        
        treatment_end <-  as.POSIXct(metadata_colonies[    which   (   metadata_colonies$colony_id == colony_id   )     ,   "time_treatment_end"     ], format="%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
        if (PERIOD=="pre"){
          time_end_GMT <- treatment_end - 24*3600
        }else{
          time_end_GMT <- treatment_end
        }
        warning("using 3h time window, but not exactly 3h")
        
        # formatting times
        time_start <- fmTimeCreate(time_start_GMT)       
        time_end  <- fmTimeCreate(time_end_GMT)
        
        if (RUN_INTERACT) {
          INTERACT_loop_start_time <- Sys.time()
          INTERACTIONS_FULL   <-  file.path(INTDIR,"full_interaction_lists",period_code,"observed",            # DESTINATION FOLDER for the interactions 
                                   paste(colony,treatment_code,period_code,"interactions.txt",sep="_"))
     
        #### 5. GET INTERACTIONS ####
        Interactions <- compute_Interactions(e = e, start = time_start, end = time_end, max_time_gap = MAX_INTERACTION_GAP)    # function in "NetworkProperties_v082.R" -  # INTERACTIONS IN THIS FUNCTION ARE CALCULATED ACCORDING TO STROEYMEYT ET AL, SCIENCE 2018
        dead_by_REP <- metadata[which(metadata$IsAlive==FALSE & metadata$colony_id==colony_id),"antID"]                        # Remove interactions involving dead ants
        Interactions <- Interactions[!(Interactions$Tag1 %in% dead_by_REP | Interactions$Tag2 %in% dead_by_REP), ]
        
        # create time vars
        Interactions$time_hours   <- NA
        Interactions$time_of_day  <- NA
        # TIME_HOURS zero is the moment of exposed ants return
        # warning("interaction binning loop should be fixed as done later for the SpaceUse, see TIME_HOURS ")
        for (TIME_HOURS in Time_dictionary_PERIOD$time_hours) { ## increments by 3 hours for duration of experiment
          # TIME_HOURS <- -24                                                                                  #temporary filler if doing line by line testing
          TIME_OF_DAY <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "time_of_day"]
          
          # I had a rather simple design regarding timing and look only at a three hour block. Hence I deleted the long parts of linda or adriano (if needed look at their scripts)
          chunk_time_start <- time_start_GMT
          chunk_time_stop <- time_end_GMT
          
          #time windows
          From_TIME_HOURS <- as.numeric(chunk_time_start,3)
          To_TIME_HOURS <- as.numeric(chunk_time_stop,3)
          
          # assing time labels
          Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_hours"]   <- TIME_HOURS
          Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_of_day"]  <- TIME_OF_DAY
        }
  
        # Add metadata info
        Interactions <- cbind(data.frame(
          period = PERIOD,Interactions, REP_treat = REP_TREAT, colony=colony, treatment=treatment_code,
          treatment_detailed=metadata_colonies[which(metadata_colonies$colony_id==colony_id),"treatment"],
          
          stringsAsFactors = F
        ))
        
        Interactions <- Interactions[,c("Tag1","Tag2","Startframe","Stopframe","Starttime","Stoptime","Box","Xcoor1","Ycoor1","Angle1","Xcoor2","Ycoor2","Angle2","Direction","Detections","time_hours","time_of_day","colony","treatment","treatment_detailed","REP_treat","period","ant1.zones","ant2.zones","duration")]
       
        # remove extra -3h gap leftovers (few mins)
        Interactions <- Interactions[which(Interactions$time_hours%in%Time_dictionary_PERIOD$time_hours),] 
        
        # Interactions save (saved INSIDE the Network_analysis folder)
  
        write.table(Interactions, file = INTERACTIONS_FULL, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
        
        # split output into bins
          # 3-HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically (so in my case not absolutely necessary as I only look at roughly 3h)
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
        INTERACT_loop_end_time <- Sys.time()
        print(paste("Interactions 3h chunk took ", round(as.numeric(difftime(INTERACT_loop_end_time, INTERACT_loop_start_time, units = "secs")),1), " sec to complete"))
        } # End of RUN_INTERACT
        
        #### 6. SPACE USE ####
        # check if file already exists
        if (file.exists(SPACE_USE)) {
          print("file exists")
          spaceUse_done <- read.table(SPACE_USE, header = T, stringsAsFactors = F, sep = "")
          spaceUse_done <- spaceUse_done[which(spaceUse_done$colony==colony),"time_hours"]
        }else{
          spaceUse_done <- c() # placeholder
        }
          
        if (RUN_SPACEUSE) {
          print("Computing SpaceUse based on full time-window pre AND post exposure")
          SpaceUse_loop_start_time <- Sys.time()
            
          for (TIME_HOURS in Time_dictionary_PERIOD$time_hours) { ## increments by 3 hours for duration of experiment # TIME_HOURS zero is the moment of exposed ants return
            TIME_OF_DAY <- Time_dictionary_PERIOD[which(Time_dictionary_PERIOD$time_hours == TIME_HOURS), "time_of_day"]
  
              if (TIME_HOURS %in% spaceUse_done) {
                cat("\rDONE UP TO: TIME_HOURS", TIME_HOURS,"TIME_OF_DAY", TIME_OF_DAY,"PERIOD",PERIOD,"| SKIP >>")
                        }else{
                          print(paste("TIME_HOURS", TIME_HOURS,"TIME_OF_DAY", TIME_OF_DAY,"PERIOD",PERIOD,sep=" "))
                          # I had a rather simple design regarding timing and look only at a three hour block. Hence I deleted the long parts of linda or adriano (if needed look at their scripts)
                          chunk_time_start <- time_start_GMT
                          chunk_time_stop <- time_end_GMT
                          # format for myrmidon
                          chunk_time_start <- fmTimeCreate(chunk_time_start)
                          chunk_time_stop  <- fmTimeCreate(chunk_time_stop)
          
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
              time_hours = TIME_HOURS,
              time_of_day = TIME_OF_DAY,
              SpaceUsage,
              REP_treat = REP_TREAT,
              stringsAsFactors = F
            ))
  
              ## Space Use save (saved INSIDE the Network_analysis folder)
              if (file.exists(SPACE_USE)) { # append to the existing two tables 
                # pre_vs_post_treatment/individual_behavioural_data.txt
                write.table(SpaceUsage, file = SPACE_USE, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
                # pre_treatment/network_position_vs_time_outside.dat
                write.table(SpaceUsage[which(SpaceUsage$period=="pre"),], file = SPACE_USE_PRE, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
              } else { # if they do not exist yet write new tables. 
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
  
    # #### Run Netwworks ####
    #   # is set to FALSE and not run here so we ignore this bit for now. 
    # if (RUN_NETWORKS) { 
    #   # add status (large, small) info to NetworkProp_collective
    #   NetworkProp_collective <- dplyr::left_join(NetworkProp_collective, unique(metadata[c("treatment", "treatment_simple", "colony_id")]), by = "colony_id") # Apply left_join dplyr function
    # 
    #   # add status (large, small) info to NetworkProp_individual
    #   NetworkProp_individual <- dplyr::left_join(NetworkProp_individual, unique(metadata[c("time_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function
    # 
    #   ########################################
    #   ##### SAVE FILES IN FOLDER #############
    # 
    #   ## Network properties Collective save (saved INSIDE the Network_analysis folder)
    #   if (file.exists(NET_properties_collective)) {
    #     write.table(NetworkProp_collective, file = NET_properties_collective, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    #   } else {
    #     write.table(NetworkProp_collective, file = NET_properties_collective, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
    #   }
    # 
    #   ## Network properties Individual save (saved INSIDE the Network_analysis folder)
    #   if (file.exists(NET_properties_individual)) {
    #     write.table(NetworkProp_individual, file = NET_properties_individual, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    #   } else {
    #     write.table(NetworkProp_individual, file = NET_properties_individual, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
    #   }
    # } # end of RUN NETWORKS
      
      
      # cleaning
      rm(list = ls()[which(!ls() %in% to_keep)])
      gc()
      mallinfo::malloc.trim(0L)
      
  } # end of the loop -> next colony
  
  loop_end_time <- Sys.time()
  print(paste("loop took ", as.numeric(difftime(loop_end_time, loop_start_time, units = "mins")), " minutes to complete"))
  


  
  

#### To dos #### 

# move on to the next analysis scripts i.e. vital_main_analysis

# Dont forget to correct when going through the downstream pipeline:

# downstream: use 13_network_analysis.R and 14_summarise_interactions.R in the version linda sent around and include her fixes 
# In another script, 12_simulate_transmission.R (which is part of the 2018 Social Plasticity code and can be found here on Adriano's github), we discovered that the output of the Surv() function has changed the name for one of its columns, which is being referred to in lines 261 & 265: "rmean" (new) instead of "*rmean" (old); due to a library update perhaps? Please be aware of this since it does not give an error but just NAs for those specific cases.

# P.S. sorry, already found another bug (I didn't test very thoroughly on Friday afternoon...) -  in 13_network_analysis add "edge_weights" to "to_keep". 
#to_keep <- c(ls(),"to_keep","input_folder","network_files","options","option","summary_collective","summary_individual","outputfolder","network_file","queenid", "edge_weights")

# When Calculating the ant tasks it might be better to do it over the full 24h of the pre treatment network?! Check with the script 19_facet_net 

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

  
### 
# Note to myself keep updating the list of scripts used as part of this analysis pipeline.
# if not specified differently, the script version used in the pipeline is in the top git hub folder or in the folder with trophallacis_classifier. 
# previous version is either in the _old folder or in a folder referring to another persion e.g Adriano or Linda.
# At some point make a cleaner github version of all this...
# Further, what ever I do with the interactions, make sure to do the same things with the trophallaxis interactions
  




