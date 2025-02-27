rm(list = ls())
gc()
mallinfo::malloc.trim(0L)

# NOTE ON DIFFERENCES FROM SCIENCE 2018 PAPER:
# period is ("pre","post"), not ("after","before")
# treatment is ("Pathogen","Sham"), not ("pathogen","sham")
# status is linked to size of the colony, not ("treated","untreated") ant

# to do in MARCH 2023
# - proportion_time_active
# - average_bout_speed_pixpersec
# total_distance_travelled_pix

# ------------------------------------
# DONE:
# COLUMNS:
# colony    colony_size x
# treatment (two values: pathogen and control)  x
# tag    (only include treated workers) x
# age    x
# status    ( replace the content of the "status" column with either "large" or "small" (rather than "treated" or "untreated") ) x
# period    ( pre/post chunks corresponding to the same time of day should have the same value in column "time_of_day") x
# time_hours    x
# time_of_day x
#

# #reference
# individual_behavioural_data <- read.table("/home/cf19810/Documents/TEMP/individual_behavioural_data.txt",header=T,stringsAsFactors = F)
# Time_dictionary_N <- unique(individual_behavioural_data[c("time_hours","time_of_day","period")])

########################################################################################################
################### IMPORTANT THINGS TO MODIFY #########################################################

# UNIFORM Plot_Grooming_Pre-Post.R to the output of this script and the sourced functions
# AGGREGATE SPACE USE FUNCTION FOR tag_hex_ID (we want 1 row per ant!)

# ########### DEFINE ZONES PROPERLY - DEFINED INSIDE THE ANT TASKS FUNCTION, SHOULD DO SAME FOR NEST USE
# #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
# zones <- e$spaces[[1]]$zones #function to show the Zones present in the Space
# zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
# foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
# nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
# print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))

# change window_shift to actual date? As per Grooming (new exact data info to be gathered)
# THE EXACT DATE IS IN /home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data/Grooming_Classifier_CrossVal_Adriano2022_RETURN_EXP_TIME_ZULU.csv

##############################################################################
######
#######
##########
#############
##########
#######
######
##############################################################################
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"

##### LIBRARIES
library(data.table)
library(lubridate)
library(pals)
library(dplyr)
# library(reshape)

# starting params
USER <- "2A13_Office" # Nath_office 

if (USER == "2A13_Office") {
  usr <- "ll16598"
} else {
  usr <- "bzniks"
}
  
SAVEDIR <- "/media/ll16598/One Touch/SICC_NET/INTERACTIONS"

  WORKDIR <- paste("/media",usr,"One Touch/SICC_NET",sep="/")
  DATADIR <- paste(WORKDIR, "EXPERIMENT_DATA", sep = "/")
  SCRIPTDIR <- paste("/media",usr,"DISK4/EXP1_base_analysis/EXP1_analysis scripts",sep="/") # "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP1_analysis scripts"
  metadata <- read.table(paste(DATADIR, "/Metadata_Exp1_2021_2022-10-20.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
  BEH_FUNCTIONS <-  paste("/media",usr,"DISK4/Ants_behaviour_analysis/ScriptsR",sep="/")
  BODYLENGTH_FILE <- paste(BEH_FUNCTIONS,"Mean_ant_length_per_TrackingSystem.txt", sep = "/")
  metadata_info         <- read.csv(paste(DATADIR,"Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU_All_colonies.csv",sep = "/"), sep = ",")


# #### ACCESS FILES
# if (USER == "Adriano") {
#   WORKDIR <- "/media/cf19810/DISK4/ADRIANO"
#   DATADIR <- paste(WORKDIR, "EXPERIMENT_DATA", sep = "/")
#   SCRIPTDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP1_analysis scripts"
#   metadata <- read.table(paste(DATADIR, "/Metadata_Exp1_2021_2022-10-20.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
#   warning("METADATA MISSING FOR R9SS AND R9BS")
#   BEH_FUNCTIONS <- "/home/cf19810/Documents/scriptsR/Ants_behaviour_analysis/ScriptsR"
#   BODYLENGTH_FILE <- paste(BEH_FUNCTIONS,"Mean_ant_length_per_TrackingSystem.txt", sep = "/")
#   metadata_info         <- read.csv(paste(DATADIR,"Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU_All_colonies.csv",sep = "/"), sep = ",")
#   
# }

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
MAX_INTERACTION_GAP <- 10 # fmSecond(10) # the maximum gap in tracking before cutting the interactions in two different object.  Used in ComputeGraph function
MAX_DIST_GAP_MM <- 0.6
FRAME_RATE <- 8
interactions_of_interest <- list(c("head", "body")) ### if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))
Time_dictionary <- data.frame(time_hours = -36:35, time_of_day = rep(0:23, 3)) # time_hours= time since exposure (negatives pre-exposire, positive post exposure), time_of_day=time hour slot of the day
Time_dictionary <- Time_dictionary[which(Time_dictionary$time_hours <= 21 & Time_dictionary$time_hours >= -27), ]
OUTPUT_FOLDER <-  "Exp1_Results_2023" # paste0("Exp1_Results_", Sys.Date())

#### FLAGS
RUN_SPACEUSE <- FALSE
RUN_NETWORKS <- FALSE
warning(paste("RUN_SPACEUSE is set to",RUN_SPACEUSE,"\nRUN_NETWORKS is set to",RUN_NETWORKS,sep=" "))

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

### create folder
dir.create(file.path(SAVEDIR, OUTPUT_FOLDER), recursive = T) 
### define name of general output files
NET_properties_collective <- file.path(SAVEDIR, OUTPUT_FOLDER, "NetworkProp_collective.txt") # (saved INSIDE the folder)
NET_properties_individual <- file.path(SAVEDIR, OUTPUT_FOLDER, "NetworkProp_individual.txt") # (saved INSIDE the folder)
SPACE_USE <- file.path(SAVEDIR, OUTPUT_FOLDER, "Space_Usage.txt") # (saved INSIDE the folder)

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
  # REP.n <- 1    #temp
  REP.folder <- files_list[REP.n]
  REP.files <- list.files(REP.folder, pattern = "CapsuleDef2018_q") #AntsCreated_AutoOriented_withMetaData_NS_NS_q_
  ##
  ### select only for the science2018 capsules, the name to grep for is in NOTION
  ##
  #REP.files <- REP.files[!grepl("CapsuleDef", REP.files)]
  REP.filefolder <- paste(REP.folder, REP.files, sep = "/")

  # replicate file
  for (REP.FILES in REP.filefolder) {
    # REP.FILES <-  REP.filefolder[2]   #temp

    ## some initialization
    Period_dataframe <- NULL # checking time correspondances
    # start fresh
    NetworkProp_collective <- data.frame()
    NetworkProp_collective_hour <- data.frame()
    NetworkProp_individual <- data.frame()
    NetworkProp_individual_hour <- data.frame()
    
    #Interactions_total <- data.frame()
    Interactions_REP_TREAT <- data.frame()
    SpaceUsage <- data.frame()

    print(REP.FILES) ## }}
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

    INTERACTIONS <- file.path(SAVEDIR, OUTPUT_FOLDER, paste0("INTERACTIONS_",REP_TREAT,".txt")) # (saved INSIDE the folder)
    
    
    #####
    # end time is return time + 24h
    exp_end <- time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"return_time"] + 60*60*24 # 24h
    #exp_end <- fmQueryGetDataInformations(e)$end - window_shift
    
    if (RUN_NETWORKS) {
    COLONY_SIZE <- unique(metadata[which(metadata$REP_treat == REP_TREAT), "colony_size"])
}

    # ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
    # print(paste0("Compute 3-hours analysis"))
    # # TIME_HOURS zero is the moment of exposed ants return
    # for (TIME_HOURS in Time_dictionary$time_hours[seq(1, length(Time_dictionary$time_hours), 3)]) { ## increments by 3 hours for 48 hours
    #   # HOUR <- seq(from=0, to=48, by=3)
    #   # From  <- fmQueryGetDataInformations(e)$end - 51*TimeWind + (HOUR * TimeWind) - window_shift
    #   # To    <- fmQueryGetDataInformations(e)$end - 48*TimeWind + (HOUR * TimeWind) - window_shift
    #   #
    #   From <- fmQueryGetDataInformations(e)$end + (TIME_HOURS - 24) * TimeWind - window_shift
    #   To <- fmQueryGetDataInformations(e)$end + (TIME_HOURS - 21) * TimeWind - window_shift
    #   #############################
    # 
    #   print(paste0("computing hour ", TIME_HOURS))
    #   print(paste("Time window, from", From, "to", To))
    #   time_start <- fmTimeCreate(offset = From) # time_stop minus 48 hours plus incremental time
    #   time_stop <- fmTimeCreate(offset = To) # time_stop minus 45 hours plus incremental time
    # 
    #   # base file information
    #   # PERIOD
    #   TimeDiff <- difftime(exp_end, To, units = "hours")
    #   if (TimeDiff < 24) {
    #     PERIOD <- "post"
    #   } else if (TimeDiff >= 27 & TimeDiff < 51) {
    #     PERIOD <- "pre"
    #   } else {
    #     PERIOD <- "EXPOSURE_GAP"
    #   }
    #   Period_dt <- data.frame(From, To, PERIOD)
    #   Period_dataframe <- rbind(Period_dataframe, Period_dt)
    #   # TIME_OF_DAY
    #   TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
    #   # }
    #   #
    #   # table(Period_dataframe$PERIOD) # shall be equal!
    # 
    #   if (RUN_SPACEUSE) {
    #     
    #   # RUN FOR PRE AND POST (skip the 3h exposure gap)
    #   if (!Period_dt$PERIOD == "EXPOSURE_GAP") {
    # 
    # 
    #     # ############ NETWORK PROPERTIES: INDIVIDUAL & COLONY LEVEL ##########################
    #     # 
    #     # # COMPUTE NETWORK
    #     # Graph <- compute_G(e = e, start = time_start, end = time_stop, gap = MAX_INTERACTION_GAP)
    #     # # COMPUTE NETWORK PROPERTIES (collective and individual)
    #     # Network_summ_prop_hour <- NetProperties(graph = Graph)
    #     # 
    #     # ### Collective
    #     # # Add metadata info
    #     # NetworkProp_collective_hour <- cbind(data.frame(
    #     #   randy = REP.FILES, REP_treat = REP_TREAT, colony_size = COLONY_SIZE, period = PERIOD, time_hours = TIME_HOURS, time_of_day = TIME_OF_DAY, From, To,
    #     #   Network_summ_prop_hour$summary_collective,
    #     #   stringsAsFactors = F
    #     # ))
    #     # # stack
    #     # NetworkProp_collective <- rbind(NetworkProp_collective, NetworkProp_collective_hour)
    #     # 
    #     # ### Individual
    #     # # Add metadata info
    #     # NetworkProp_individual_hour <- cbind(data.frame(
    #     #   randy = REP.FILES, REP_treat = REP_TREAT, colony_size = COLONY_SIZE, period = PERIOD, time_hours = TIME_HOURS, time_of_day = TIME_OF_DAY, From, To,
    #     #   Network_summ_prop_hour$summary_individual,
    #     #   stringsAsFactors = F
    #     # ))
    #     # # stack
    #     # NetworkProp_individual <- rbind(NetworkProp_individual, NetworkProp_individual_hour)
    # 
    #     ############ SPACE USE: INDIVIDUAL ##########################
    #     # SPACE USAGE
    #     SpaceUsage_hour <- SpaceUse(e = e, start = time_start, end = time_stop)
    # 
    #     # Add metadata info
    #     SpaceUsage_hour <- cbind(data.frame(
    #       randy = REP.FILES, REP_treat = REP_TREAT, colony_size = COLONY_SIZE, period = PERIOD, time_hours = TIME_HOURS, time_of_day = TIME_OF_DAY, From, To,
    #       SpaceUsage_hour,
    #       stringsAsFactors = F
    #     ))
    # 
    #     # stack
    #     SpaceUsage <- rbind(SpaceUsage, SpaceUsage_hour)
    # 
    #     #####################################################################################################
    #   } 
    #     }
    #   } # REP LOOP
    

    
    ### PERFORM THE FULL INTERACTION AND NETWORK ANALYSIS OUTSIDE OF THE HOURLY LOOP
    
    ############ GET INTERACTIONS ##########################
    if(file.exists(INTERACTIONS)){
      print(paste("File",INTERACTIONS,"already present, skip", sep=" "))
    }
      
    if(!file.exists(INTERACTIONS)){
    
    # # # Extract the minimum "From" time and maximum "To" time for each PERIOD
    # Period_windows <- Period_dataframe %>%
    #   group_by(PERIOD) %>%
    #   summarise(From = min(c(From,To)),To = max(c(From,To)))
  
  Period_windows <- data.frame(Period = c("pre","post")
                               , From = c(time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"time_start"], #skip 3h gap
                                          time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"return_time"]
                               ))
  
  Period_windows$To <-   Period_windows$From + 24*60*60
  
  
    
    for (PERIOD in c("pre","post")) {
      # Select the full PERIOD (24h)
      print(paste("period:", PERIOD, sep= " "))
      time_start <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"From"]) # time_stop minus 48 hours plus incremental time
      time_stop <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"To"]) # time_stop minus 45 hours plus incremental time
      # # TEMP TIME STOP
      # warning("using tiny time window for testing")
      # time_stop <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"From"] + 1*60)
      # 
      # INTERACTIONS IN THIS FUNCTION ARE CALCULATED ACCORDING TO STROEYMEYT ET AL, SCIENCE 2018
      Interactions <- compute_Interactions(e = e, start = time_start, end = time_stop, max_time_gap = MAX_INTERACTION_GAP)
      
      ### Individual
      # Add metadata info
      Interactions_REP_TREAT <- cbind(data.frame(
         REP_treat = REP_TREAT, period = PERIOD,Interactions, randy = REP.FILES,
        stringsAsFactors = F
      ))
      # stack
      #Interactions_total <- rbind(Interactions_total, Interactions_REP_TREAT)
      
      ## Interactions save (saved INSIDE the Network_analysis folder)
      if (file.exists(INTERACTIONS)) {
        write.table(Interactions_REP_TREAT, file = INTERACTIONS, append = T, col.names = F, row.names = F, quote = T, sep = ",")
      } else {
        write.table(Interactions_REP_TREAT, file = INTERACTIONS, append = F, col.names = T, row.names = F, quote = T, sep = ",")
      }
      
    } # PERIOD
    } # only run if file does not exist
    
  
    #### ADD EXTRA REPLICATE INFORMATIONS
  if (RUN_NETWORKS) {
    # add status (large, small) info to NetworkProp_collective
    NetworkProp_collective <- dplyr::left_join(NetworkProp_collective, unique(metadata[c("size_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function

    # add status (large, small) info to NetworkProp_individual
    NetworkProp_individual <- dplyr::left_join(NetworkProp_individual, unique(metadata[c("size_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function
}

    # add task, exposure and status (large, small)  info to SpaceUsage
    if (RUN_SPACEUSE) {
    SpaceUsage <- dplyr::left_join(SpaceUsage, metadata[c("size_treat", "status", "treatment", "REP_treat", "antID", "Exposed", "AntTask")], by = c("REP_treat", "antID")) # Apply left_join dplyr function
}

    ########################################
    ##### SAVE FILES IN FOLDER #############

    # NOTE ON DIFFERNCES FROM SCIENCE 2018 PAPER:
    # period is ("pre","post"), not ("after","before")
    # treatment is ("Pathogen","Sham"), not ("pathogen","sham")
    # status is linked to size of the colony, not ("treated","untreated") ant

  

  if (RUN_NETWORKS) {    
    ## Network properties Collective save (saved INSIDE the Network_analysis folder)
    if (file.exists(NET_properties_collective)) {
      write.table(NetworkProp_collective, file = NET_properties_collective, append = T, col.names = F, row.names = F, quote = T, sep = ",")
    } else {
      write.table(NetworkProp_collective, file = NET_properties_collective, append = F, col.names = T, row.names = F, quote = T, sep = ",")
    }

    ## Network properties Individual save (saved INSIDE the Network_analysis folder)
    if (file.exists(NET_properties_individual)) {
      write.table(NetworkProp_individual, file = NET_properties_individual, append = T, col.names = F, row.names = F, quote = T, sep = ",")
    } else {
      write.table(NetworkProp_individual, file = NET_properties_individual, append = F, col.names = T, row.names = F, quote = T, sep = ",")
    }
}

    ## Space Use save (saved INSIDE the Network_analysis folder)
    if (RUN_SPACEUSE) {
    if (file.exists(SPACE_USE)) {
      write.table(SpaceUsage, file = SPACE_USE, append = T, col.names = F, row.names = F, quote = T, sep = ",")
    } else {
      write.table(SpaceUsage, file = SPACE_USE, append = F, col.names = T, row.names = F, quote = T, sep = ",")
    }
}


    # start fresh
    NetworkProp_collective <- data.frame()
    NetworkProp_individual <- data.frame()
    #Interactions_total     <- data.frame()
    SpaceUsage             <- data.frame()

    # cleaning
    rm(list = ls()[which(!ls() %in% to_keep)])
    gc()
    mallinfo::malloc.trim(0L)


    ######################################################
    ### PLOTTING (1 col..., should be all)
    # plot(AntTasks$delta_time_inside, col=as.factor(AntTasks$Exposed))

    ####################################################################################
    ### plot MEAN DELTA PER ANT SEPARATED BETWEEN NON EXPOSED AND EXPOSED, WITH COL. SIZE COMPARISON
    ###################################################################################
  }
}


loop_end_time <- Sys.time()
print(paste("loop took ", as.numeric(difftime(loop_end_time, loop_start_time, units = "mins")), " minutes to complete"))










################################################################################################################
########## THE INFORMATION ON METADATA IS USELESS AND MISLEADING, REMOVE ######################################
###############################################################################################################

# ########## GET EXPOSED ANTS from METADATA
# e.Ants <- e$ants
# metadata$Exposed <- "no"
# for (ant in e.Ants){
#   individual  <- ant$ID
#   #print(ant)
#     if (TRUE %in% ant$getValues("Exposed")[,"values"]) { metadata[individual,"Exposed"] <- "exposed" }
# }
#
# ########## GET QUEENS
# metadata$IsQueen <- "no"
# for (ant in e.Ants){
#   individual  <- ant$ID
#   #print(ant)
#   if (TRUE %in% ant$getValues("IsQueen")[,"values"]) { metadata[individual,"IsQueen"] <- "queen" }
# }
