rm(list = ls())
gc()
Sys.sleep(3)
mallinfo::malloc.trim(0L) # forced memory cleaning after every gc() | # install in terminal: > sudo apt-get install libtcmalloc-minimal4
                          # and in R: > install.packages("mallinfo", repos = "http://www.rforge.net/") | and run the following after gc(): > mallinfo::malloc.trim(0L)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ###   BASE ANALYSES ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#### READ ME ####
# This is Script is based on Adrianos, Tom & Nathalies Exp1_base_analysis script and got adjusted to the data structure and needs of Daniel 

#### To Do's ####

# work through the whole script line by line
# get the exact end time for each of the experiments
# UPDATE all parts of the script that try to access BODYLENGTH_FILE by changing the code in a way that it refers to size in pixels based on colonies instead of tracking systems
# introduce all the fixes suggested via email:
# make sure everthing can run on laptop and uni computer (starting if statement)
# get a version of the myrmidon files with Nathalies 2018 capsule definitions for normal interactions. 
# What ever I do with the interaction network, make sure to get the same things with the trophallaxis networks
# Make sure to change frame rate to 6 (last in the last updated exp1 base analysis script)



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

# 
# Add the input from linda: 
#   here is one general error and a warning that i noticed in the EXP1_base_analysis. Maybe it's helpful for you too.
# 
# Everytime you define IF_frames (once in NetworkProperties_v082.R and once in SpaceUse_v082.R) you have to again order the time column (this will throw an error otherwise like it did in the previous script):
# 
# 
# IdentifyFrames      <- fmQueryIdentifyFrames(e,start, end,showProgress = FALSE)
# IF_frames           <- IdentifyFrames$frames
# IF_frames           <- IF_frames[order(IF_frames$time), ]
# 
# 
# Also in SpaceUse_v082.R right after you extract the trajectories like this:
# 
# 
# positions <- extract_trajectories(e, start , end, maximumGap = fmHour(24*365),computeZones = T,showProgress = F,IF_frames=IF_frames)
# positions_summaries       <- positions$trajectories_summary
# positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
# positions_list            <- positions$trajectories
# warning("order the trajectories to ensure that the following calculations are correct" ) 
# 
# there is a warning that was added by Adriano, telling you to make sure the trajectories are sorted (this would otherwise seriously mess with all consequent calculations performed on them). So I added a little loop to make sure the time of the trajectories of each ant in positions_list is sorted. (This is just a safety measure - the script runs without it). I used a try() statement because sometimes not all ants from 1 to 1:length(e$ants) exist in positions_list and if they don't, "to_sort" is neither TRUE nor FALSE and that error would stop the whole thing. When using try() it still gives you the error in the terminal but the script keeps running.
# 
# 
# # L: addition to the script to check/sort trajectories because of the warning ###
# 
# for (ant in   1: length(e$ants)){
#   ant <- paste("ant", ant, sep="_")
#   # Check if the column is sorted in ascending order
#   is_sorted <- !is.unsorted(positions_list[[ant]]$time)
#   #print(paste(ant, "is_sorted was created and it is", is_sorted))
#   
#   try({
#     if(is_sorted==FALSE)
#     {warning(paste(REP_TREAT,"- trajectories for", ant, "were not sorted initially"))
#       positions_list[[ant]]$time <- sort(positions_list[[ant]]$time)
#     }
#   })
# }




#### LIBRARIES #### 
library(data.table)
library(lubridate)
library(pals)
library(dplyr)
# library(reshape)
library(adehabitatHR)
library(adehabitatLT)
library(changepoint)
library(e1071)
library(igraph)
library(gtools)
library(Rcpp)
library(survival)

#### starting parameters ####

USER <- "AEL-laptop"  # Replace with the desired USER option: Nath_office, 2A13_Office_Adriano, 2A13_Office_Daniel, AEL-laptop
HD <- "/gismo_hd2" # One of the harddrives with the vital extrapolated data: possible values > Nathalies hd "/DISK_B" ; Daniels hds >  "/gismo_hd5" or  "/gismo_hd2"

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
  assign("usr", usr, envir = .GlobalEnv)  # Assign to the global environment
  hd <- NULL
  hd <- HD
  if (!is.null(hd)) {print(hd)} else {print("define new hd if necessary")}
  assign("hd", hd, envir = .GlobalEnv)  # Assign to the global environment
}
setUserAndHD(USER, HD)

#### set directories #### 
WORKDIR <- paste("/media/",usr, hd, "/vital/fc2",sep="") 
DATADIR <- paste(WORKDIR, sep = "/") # working directory and directory where the tracking data is saved (was not the same for Adriano, but for me the same
SAVEDIR <- paste("/media/",usr, hd,"/vital/fc2/vital_experiment/summary_data",sep="") # where to save the interactions - I added a folder to the vital
INTDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/intermediary_analysis_steps",sep="") # remember to use the same folder structure as as for Science 2018
BEHDIR <- paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment/processed_data/individual_behaviour",sep="")
SCRIPTDIR <- paste("/home/",usr,"/Documents/vital_rscripts_git",sep="") # place where the needed r scripts are stored
BEH_FUNCTIONS <-  paste(SCRIPTDIR, "/Behavioural_Inference_DS",sep="")


#### load metadata ####
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # will add colony_metadata data frame to the environment so it can be accessed within this script

# BODYLENGTH_FILE <- paste(BEH_FUNCTIONS,"Mean_ant_length_per_TrackingSystem.txt", sep = "/")

metadata <- read.table(paste(DATADIR, "/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",") # contains return time for each colony and end of experiment time! 



#### Run only 
# when creating the meta data for each ant add an additional thing to refer to the ants status 
# metadata$status_ant <- NA
# metadata$status_ant <-ifelse(metadata$Expose==TRUE,"treated","untreated")




  
  

### source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR, "SpaceUse_v082.R", sep = "/")) # SPACE_USE
source(paste(SCRIPTDIR, "NetworkProperties_v082.R", sep = "/")) # NET_properties collective + individual
suppressMessages(source(paste(BEH_FUNCTIONS, "BEH_libraries.R", sep = "/")))
source(paste(BEH_FUNCTIONS, "interaction_detection.R", sep = "/"))
source(paste(BEH_FUNCTIONS, "trajectory_extraction.R", sep = "/"))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # required by merge_interactions.cpp
sourceCpp(paste(BEH_FUNCTIONS, "merge_interactions.cpp", sep = "/"))
# list.functions.in.file(paste(SCRIPTDIR,"NetworkProperties_v082.R",sep="/"), alphabetic = TRUE)

#### FUNCTIONS
# list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n - 1)
    c(res, add)
  } else {
    res
  }
}

# body_lenghts
get_body_lengths <- function(e, all_body_lengths) {
  ### check how many spaces there are
  space_list <- e$spaces
  body_lengths <- NULL
  for (space_ID in 1:length(space_list)) {
    body_lengths <- rbind(body_lengths, data.frame(
      space_ID = space_list[[space_ID]]$ID,
      tracking_system = space_list[[space_ID]]$name,
      body_length_px = all_body_lengths[which(all_body_lengths$TS == space_list[[space_ID]]$name), "mean_worker_length_px"],
      body_length_mm = all_body_lengths[which(all_body_lengths$TS == space_list[[space_ID]]$name), "mean_worker_length_mm"]
    ))
  }
  return(body_lengths)
}


#### PARAMETERS
all_body_lengths <- read.table(BODYLENGTH_FILE, header = T, stringsAsFactors = F, sep = ",") ### body length information
TimeWind <- 3600 ## in seconds (3600 is an hour)
N_DECIMALS <- 3
MAX_INTERACTION_GAP <- 10 # fmSecond(10) # the maximum gap in tracking before cutting the interactions in two different object.  Used in ComputeGraph function
MAX_DIST_GAP_MM <- 0.6
FRAME_RATE <- 8
interactions_of_interest <- list(c("head", "body")) ### if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))

Time_dictionary <- data.frame(time_hours = -36:35, time_of_day = rep(0:23, 3)) # time_hours= time since exposure (negatives pre-exposire, positive post exposure), time_of_day=time hour slot of the day
Time_dictionary <- Time_dictionary[which(Time_dictionary$time_hours <= 21 & Time_dictionary$time_hours >= -27), ]
Time_dictionary$period <- ifelse(Time_dictionary$time_hours<0, "pre", "post")

#OUTPUT_FOLDER <-  "Exp1_Results_2023" # paste0("Exp1_Results_", Sys.Date()) #OUT OF DATE?
#SAVE FILES IN THE LOCATION WHERE THESE ARE SAVED FOR THE SCIENCE 2018 PIPELINE
SPACE_USE     <-  file.path(BEHDIR,"pre_vs_post_treatment","individual_behavioural_data.txt")
SPACE_USE_PRE <-  file.path(BEHDIR,"pre_treatment","network_position_vs_time_outside.dat")




#### FLAGS
RUN_INTERACT     <- TRUE
RUN_SPACEUSE     <- FALSE
RUN_NETWORKS     <- FALSE
warning(paste("RUN_INTERACT is set to:",RUN_INTERACT,
              "\nRUN_SPACEUSE is set to:",RUN_SPACEUSE,
              "\nRUN_NETWORKS is set to:",RUN_NETWORKS,sep="\t"))

## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED!
#window_shift <- 60 * 15 # approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition

# list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
# select REP folders
files_list <- files_list[grep("REP", files_list)]

### initialise general output folder
### remove folder if already exists to make sure we don't mix things up
# if (file.exists(file.path(DATADIR, "NetworkAnalysis_outcomes"))){
#   unlink(file.path(DATADIR, "NetworkAnalysis_outcomes"),recursive=T)
# }

# ### create folder
# dir.create(file.path(SAVEDIR, OUTPUT_FOLDER), recursive = T) 
# ### define name of general output files
# NET_properties_collective <- file.path(SAVEDIR, OUTPUT_FOLDER, "NetworkProp_collective.txt") # (saved INSIDE the folder)
# NET_properties_individual <- file.path(SAVEDIR, OUTPUT_FOLDER, "NetworkProp_individual.txt") # (saved INSIDE the folder)

##### RUNNING TIME
loop_start_time <- Sys.time()


###create object that will contain time limits for annotations_all (contains all annotations for all behaviour so true time limit for period analysed)
# time_window_all        <- merge(aggregate(T_start_UNIX ~ PERIOD + REPLICATE, FUN=min, data=annotations_all),aggregate(T_stop_UNIX  ~ PERIOD + REPLICATE, FUN=max, data=annotations_all))
# names(time_window_all) <- c("PERIOD","REPLICATE","time_start","time_stop")
time_window_all        <- data.frame(
  PERIOD = "all"
  ,REPLICATE = metadata_info$REP_treat
  ,return_time = metadata_info$ReturnExposed_time #time_start
)
# 
time_window_all$return_time <- as.POSIXct(time_window_all$return_time, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
time_window_all$time_start <- time_window_all$return_time - 27*60*60 #(24h + 3h gap)
time_window_all$time_stop  <- time_window_all$return_time + 24*60*60
# 

#### define to_keep variables to keep clearing memory between runs
to_keep <- c(ls(), c("to_keep"))

#### OPEN REPLICATE
# replicate folder
for (REP.n in 1:length(files_list)) {
  # REP.n <- 13    #temp
  REP.folder <- files_list[REP.n]
  REP.files <- list.files(REP.folder, pattern = "CapsuleDef2018_q") #AntsCreated_AutoOriented_withMetaData_NS_NS_q_
  ##
  ### select only for the science2018 capsules, the name to grep for is in NOTION
  ##
  #REP.files <- REP.files[!grepl("CapsuleDef", REP.files)]
  REP.filefolder <- paste(REP.folder, REP.files, sep = "/")

  # replicate file
  for (REP.FILES in REP.filefolder) {
    # REP.FILES <-  REP.filefolder[1]   #temp

    ## some initialization
    Period_dataframe <- NULL # checking time correspondances
    # start fresh
    NetworkProp_collective <- data.frame()
    NetworkProp_collective_hour <- data.frame()
    NetworkProp_individual <- data.frame()
    NetworkProp_individual_hour <- data.frame()
    
    #Interactions_total <- data.frame()
    Interactions_REP_TREAT <- data.frame()
    #SpaceUsage <- data.frame()

    cat(paste("########################################\n",basename(REP.FILES)),sep="") ## }}
    # open experiment
    e <- fmExperimentOpen(REP.FILES)
    # e.Ants <- e$ants
    print(paste0("Processing ", basename(REP.FILES)))

    ### get body length
    body_lengths <- get_body_lengths(e, all_body_lengths)

    # base file info
    # Split the string into a vector based on "_"
    input_vector <- strsplit(basename(REP.FILES), "_")[[1]]
    # Extract the first 2 elements of the vector
    # Use grep to find the element containing "R"
    REP_TREAT <- input_vector[grep("R", input_vector[1:2])]
    # REP_TREAT <- sub("\\_.*", "", )
    # SIZE_TREAT <- substr(REP_TREAT,(nchar(REP_TREAT)+1)-2,nchar(REP_TREAT))
  
  #   
  #   REPS <- c(REPS, REP_TREAT) #
  # 
  # }}

    
    #####
    # end time is return time + 24h
    exp_end <- time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"return_time"] + 60*60*24 # 24h
    #exp_end <- fmQueryGetDataInformations(e)$end - window_shift
    
#     if (RUN_NETWORKS) {
#     COLONY_SIZE <- unique(metadata[which(metadata$REP_treat == REP_TREAT), "colony_size"])
# }
# 
#     ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
#     print(paste0("Compute 3-hours analysis"))
#     # TIME_HOURS zero is the moment of exposed ants return
#     for (TIME_HOURS in Time_dictionary$time_hours[seq(1, length(Time_dictionary$time_hours), 3)]) { ## increments by 3 hours for 48 hours
#       # HOUR <- seq(from=0, to=48, by=3)
#       # From  <- fmQueryGetDataInformations(e)$end - 51*TimeWind + (HOUR * TimeWind) - window_shift
#       # To    <- fmQueryGetDataInformations(e)$end - 48*TimeWind + (HOUR * TimeWind) - window_shift
#       #
#       From <- fmQueryGetDataInformations(e)$end + (TIME_HOURS - 24) * TimeWind - window_shift
#       To <- fmQueryGetDataInformations(e)$end + (TIME_HOURS - 21) * TimeWind - window_shift
#       #############################
# 
#       print(paste0("computing hour ", TIME_HOURS))
#       print(paste("Time window, from", From, "to", To))
#       time_start <- fmTimeCreate(offset = From) # time_stop minus 48 hours plus incremental time
#       time_stop <- fmTimeCreate(offset = To) # time_stop minus 45 hours plus incremental time
# 
#       # base file information
#       # PERIOD
#       TimeDiff <- difftime(exp_end, To, units = "hours")
#       if (TimeDiff < 24) {
#         PERIOD <- "post"
#       } else if (TimeDiff >= 27 & TimeDiff < 51) {
#         PERIOD <- "pre"
#       } else {
#         PERIOD <- "EXPOSURE_GAP"
#       }
#       Period_dt <- data.frame(From, To, PERIOD)
#       Period_dataframe <- rbind(Period_dataframe, Period_dt)
#       # TIME_OF_DAY
#       TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
#       # }
#       #
#       # table(Period_dataframe$PERIOD) # shall be equal!
# 
#       if (RUN_SPACEUSE) {
# 
#       # RUN FOR PRE AND POST (skip the 3h exposure gap)
#       if (!Period_dt$PERIOD == "EXPOSURE_GAP") {
# 
# 
#         # ############ NETWORK PROPERTIES: INDIVIDUAL & COLONY LEVEL ##########################
#         #
#         # # COMPUTE NETWORK
#         # Graph <- compute_G(e = e, start = time_start, end = time_stop, gap = MAX_INTERACTION_GAP)
#         # # COMPUTE NETWORK PROPERTIES (collective and individual)
#         # Network_summ_prop_hour <- NetProperties(graph = Graph)
#         #
#         # ### Collective
#         # # Add metadata info
#         # NetworkProp_collective_hour <- cbind(data.frame(
#         #   randy = REP.FILES, REP_treat = REP_TREAT, colony_size = COLONY_SIZE, period = PERIOD, time_hours = TIME_HOURS, time_of_day = TIME_OF_DAY, From, To,
#         #   Network_summ_prop_hour$summary_collective,
#         #   stringsAsFactors = F
#         # ))
#         # # stack
#         # NetworkProp_collective <- rbind(NetworkProp_collective, NetworkProp_collective_hour)
#         #
#         # ### Individual
#         # # Add metadata info
#         # NetworkProp_individual_hour <- cbind(data.frame(
#         #   randy = REP.FILES, REP_treat = REP_TREAT, colony_size = COLONY_SIZE, period = PERIOD, time_hours = TIME_HOURS, time_of_day = TIME_OF_DAY, From, To,
#         #   Network_summ_prop_hour$summary_individual,
#         #   stringsAsFactors = F
#         # ))
#         # # stack
#         # NetworkProp_individual <- rbind(NetworkProp_individual, NetworkProp_individual_hour)
# 
#       }
#         }
#       } # REP LOOP
    

    ### PERFORM THE FULL INTERACTION AND NETWORK ANALYSIS OUTSIDE OF THE HOURLY LOOP #########################
    
    # # # Extract the minimum "From" time and maximum "To" time for each PERIOD
    Period_windows <- data.frame(Period = c("pre","post")
                                 , From = c(time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"time_start"], #skip 3h gap
                                            time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"return_time"]
                                 ))
    Period_windows$To <-   Period_windows$From + 24*60*60
    
    
    for (PERIOD in c("pre","post")) {
    
      #segment the time dictionary
      Time_dictionary_PERIOD <- Time_dictionary[which(Time_dictionary$period==PERIOD),]
      #conform naming to science2018
      # colony code
      REP_NUM           <- substring(REP_TREAT, 2, nchar(REP_TREAT))
      colony            <- paste("colony", paste(rep(0,4-nchar(REP_NUM)),collapse=""), REP_NUM,sep="")
      # colony_status
      status_char       <- substr(REP_TREAT, nchar(REP_TREAT), nchar(REP_TREAT))
      colony_status     <- ifelse(status_char == "P", "pathogen", ifelse(status_char == "S", "control", NA))
      # treatment code
      period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
      size_char         <- substr(REP_TREAT, nchar(REP_TREAT)-1, nchar(REP_TREAT)-1)
      size_status       <- ifelse(size_char == "S", "small", ifelse(size_char == "B", "big", NA))
      treatment_code    <- paste(colony_status,size_status,sep=".")
      
      # Select the full PERIOD (24h)
      cat(paste("#######","period:", PERIOD, sep= " "))
      #time_start <- as.numeric(Period_windows[which(Period_windows$Period==PERIOD),"From"]) # time_stop minus 48 hours plus incremental time
      time_start <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"From"]) # time_stop minus 48 hours plus incremental time
      time_stop <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"To"]) # time_stop minus 45 hours plus incremental time
      ##TEMP TIME STOP
      # warning("using tiny time window for testing")
      # time_stop <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"From"] + 2*60)
      
      
      if (RUN_INTERACT) {
        INTERACT_loop_start_time <- Sys.time()
      # DESTINATION FOLDER
      INTERACTIONS_FULL   <-  file.path(INTDIR,"full_interaction_lists",period_code,"observed", 
                                 paste(colony,treatment_code,period_code,"interactions.txt",sep="_"))

    ############ GET INTERACTIONS ##########################
    # if(file.exists(INTERACTIONS_FULL)){  #TEMP
    #   warning(paste(REP_TREAT,"|",colony,treatment_code,period_code,"| already present, SKIP >>", sep=" "))
    # }  #TEMP
    #   
    # if(!file.exists(INTERACTIONS_FULL)){  #TEMP
    
      # INTERACTIONS IN THIS FUNCTION ARE CALCULATED ACCORDING TO STROEYMEYT ET AL, SCIENCE 2018
      Interactions <- compute_Interactions(e = e, start = time_start, end = time_stop, max_time_gap = MAX_INTERACTION_GAP)
      
      # Remove interactions involving dead ants
      dead_by_REP <- metadata[which(metadata$IsAlive==FALSE & metadata$REP_treat==REP_TREAT),"antID"]
      Interactions <- Interactions[!(Interactions$Tag1 %in% dead_by_REP | Interactions$Tag2 %in% dead_by_REP), ]
      
      # create time vars
      Interactions$time_hours   <- NA
      Interactions$time_of_day  <- NA
      # TIME_HOURS zero is the moment of exposed ants return
      warning("interaction binning loop should be fixed as done later for the SpaceUse, see TIME_HOURS ")
      for (TIME_HOURS in Time_dictionary$time_hours[seq(1, length(Time_dictionary$time_hours), 3)]) { ## increments by 3 hours for 48 hours
        
        # TIME_OF_DAY
        TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
        #time windows
        From_TIME_HOURS <- as.numeric((time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind),N_DECIMALS)
        To_TIME_HOURS <- as.numeric((time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 21) * TimeWind),N_DECIMALS)
        
        #assing time labels
        Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_hours"]   <- TIME_HOURS
        Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_of_day"]  <- TIME_OF_DAY
      }
      
      # Add metadata info
      Interactions <- cbind(data.frame(
        period = PERIOD,Interactions, REP_treat = REP_TREAT, colony=colony, treatment=treatment_code,
        stringsAsFactors = F
      ))
      
      # extras not present in Science files: (added at the end of the file output)
      #"REP_treat","period","ant1.zones","ant2.zones","duration"
      Interactions <- Interactions[,c("Tag1","Tag2","Startframe","Stopframe","Starttime","Stoptime","Box","Xcoor1","Ycoor1","Angle1","Xcoor2","Ycoor2","Angle2","Direction","Detections","time_hours","time_of_day","colony","treatment","REP_treat","period","ant1.zones","ant2.zones","duration")]
      # remove extra -3h gap leftovers (few mins)
      Interactions <- Interactions[which(Interactions$time_hours!=-3),] 
      ## Interactions save (saved INSIDE the Network_analysis folder)
      #if (file.exists(INTERACTIONS_FULL)) {
      #  write.table(Interactions, file = INTERACTIONS_FULL, append = T, col.names = F, row.names = F, quote = F, sep = ",")
      #} else {
        write.table(Interactions, file = INTERACTIONS_FULL, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
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
        print("Computing SpaceUse based on 24h time-window pre AND post exposure")
        warning("the Space Use Script can be substantially sped up by computing trajectories beforehand, then cutting the result by 3h chunks and performing operations on them.
                \n -This is because trajectories stitching is slow so may be better to perform it once. see inspiration from 8_process_trajectory_files.R
                \n -This will require saving the summary stats per ant in a safe place (possibly the trajectories summary as now but with clearing of the vars at each time cycle)
                \n -parallelisation of ant computation greatly improve speed. Before parallelisation of SpaceUse calculations: 2.1 sec per ant per 3h chunk (70h), after 0.21 sec with 10 cores. this includes only the summary stats computations (traj stitching excluded)")
        SpaceUse_loop_start_time <- Sys.time()
        # TIME_HOURS zero is the moment of exposed ants return
        for (TIME_HOURS in Time_dictionary_PERIOD$time_hours[seq(1, length(Time_dictionary_PERIOD$time_hours), 3)]) { ## increments by 3 hours for 48 hours
          
          # TIME_OF_DAY
          TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
            # has the time_hour been computed already? 
            if (TIME_HOURS %in% spaceUse_done) {
              cat("\rDONE UP TO: TIME_HOURS", TIME_HOURS,"TIME_OF_DAY", TIME_OF_DAY,"PERIOD",PERIOD,"| SKIP >>")
                      }else{
                        print(paste("TIME_HOURS", TIME_HOURS,"TIME_OF_DAY", TIME_OF_DAY,"PERIOD",PERIOD,sep=" "))
                        
          #time windows
          time_start_h <- fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind)) # time_stop minus 48 hours plus incremental time
          time_stop_h  <- fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 21) * TimeWind)) # time_stop minus 45 hours plus incremental time
          ##TEMP TIME STOP
          # warning("using tiny time window for testing")
          #time_stop_h <-  fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind + 40*60) ) 
          
         
          ############ SPACE USE ##########################
          # SPACE USAGE
          SpaceUsage <- SpaceUse(e = e, start = time_start_h, end = time_stop_h)
          
          SpaceUsage <- dplyr::left_join(SpaceUsage, metadata[which(metadata$REP_treat==REP_TREAT),c("status_ant", "antID","IsAlive")], by = "antID") # Apply left_join dplyr function
          colnames(SpaceUsage)[which(colnames(SpaceUsage)=="antID")] <- "tag"
          colnames(SpaceUsage)[which(colnames(SpaceUsage)=="status_ant")] <- "status"
          
          #remove dead ants
          SpaceUsage <- SpaceUsage[which(SpaceUsage$IsAlive==T),]
          SpaceUsage$IsAlive <- NULL

          # Add metadata info
          SpaceUsage <- cbind(data.frame(
            colony = colony,
            colony_size = unique(metadata[which(metadata$REP_treat==REP_TREAT),"colony_size"]),
            treatment = treatment_code,
            age = NA,
            period = PERIOD,
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
          } # if present already
        } # hour loop
        SpaceUse_loop_end_time <- Sys.time()
        print(paste("spaceUse 3h chunk took ", as.numeric(difftime(SpaceUse_loop_end_time, SpaceUse_loop_start_time, units = "secs")), " sec to complete"))
        }
      } # PERIOD
    }
    
  
    #### ADD EXTRA REPLICATE INFORMATIONS
  if (RUN_NETWORKS) {
    # add status (large, small) info to NetworkProp_collective
    NetworkProp_collective <- dplyr::left_join(NetworkProp_collective, unique(metadata[c("size_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function

    # add status (large, small) info to NetworkProp_individual
    NetworkProp_individual <- dplyr::left_join(NetworkProp_individual, unique(metadata[c("size_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function

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
}

    # start fresh
    NetworkProp_collective <- data.frame()
    NetworkProp_individual <- data.frame()
    Interactions           <- data.frame()
    SpaceUsage             <- data.frame()

    # cleaning
    rm(list = ls()[which(!ls() %in% to_keep)])
    gc()
    mallinfo::malloc.trim(0L)
  }

loop_end_time <- Sys.time()
print(paste("loop took ", as.numeric(difftime(loop_end_time, loop_start_time, units = "mins")), " minutes to complete"))
