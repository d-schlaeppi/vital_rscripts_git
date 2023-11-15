##########################################################################################
############## BEH MAIN Behaviours Analysis ##############################################
##########################################################################################

#### THIS VERSION IS FORT 0.8.1 COMPATIBLE ####

# Script created by Adriano Wanderlingh, Nathalie Stroeymeyt and Tom Richardson, with contributions by Enrico Gavagnign

#this should be the version of the script maintained for long term use.
#For previous versions of this script and the sourced one, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral


#clean start
rm(list=ls())
gc()

###parameter to set at start
USER                       <- "Nathalie"#"Adriano"
training_set               <- "Daniel" ##Adriano, Vasudha or both
BEH                        <- "T"
FRAME_RATE                 <- 6
FLUSH                      <- F###set to T if you want to delete the existing Machine Learning output folder to start fresh    
ITER                       <- 5
exclude_queen              <- T
queen_id                   <- 1
###TO DO: also edit the path to BODYLENGTH_FILE when defining folders from line 45 onwards

if (BEH =="G"){
  interactions_of_interest <- list(c("head","body"))
  ###if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))
}else if (BEH=="T"){
  interactions_of_interest <- list(c("head","head"))
}

###############################################################################
###### GLOSSARY ###############################################################
###############################################################################

# MAN         : manual interactions/trajectories deriving from the hand annotated data (from the file annotations)
# AUTO        : automatically extracted interactions/trajectories deriving from fmQueryComputeAntInteractions
# REPLICATE   : nest ID, either "R3SP" or "R9SP" (SP: small pathogen)
# PERIOD      : the treatment period, either "pre" or "post" pathogen exposure
# REP_PER     : each of the 4 blocks analised, crossing the REPLICATE and PERIOD

###############################################################################
###### LOAD LIBRARIES AND FUNCTIONS #####################################
###############################################################################

####### navigate to folder containing myrmidon file
WORKDIR <- "/media/bzniks/gismo_hd2/vital/fc2/"
DATADIR <- WORKDIR
SCRIPTDIR <- "~/Dropbox/SeniorLectureship_Bristol/Students_postdocs/Post-Docs/Daniel Schlaeppi/Trophallaxis_Classifier"
SAVEOUTPUT <- "/media/bzniks/gismo_hd6"
# BODYLENGTH_FILE <- paste(WORKDIR,"Data","/Mean_ant_length_per_TrackingSystem.txt",sep="/")


###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"BEH_Extract_movement_variables_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_Auto_Man_agreement_matrix_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_PCA_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_self_defined_functions.R",sep="/"))
source(paste(SCRIPTDIR,"interaction_detection.R",sep="/"))
suppressMessages(source(paste(SCRIPTDIR,"BEH_libraries.R",sep="/")))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(SCRIPTDIR,"add_angles.cpp",sep="/"))
sourceCpp(paste(SCRIPTDIR,"merge_interactions.cpp",sep="/"))

duplicate_annotations <- function(annotation_subset,time_limit,total_duplicated_events){
  ###time_limit will work well if those two events don't overlap
  ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
  ###first list events to duplicate (if any)
  to_duplicate  <- annotation_subset[which(  (as.numeric(annotation_subset$T_start_UNIX)<=time_limit & as.numeric(annotation_subset$T_stop_UNIX)>time_limit) ),]
  
  ###perform duplication if necessary
  if (nrow(to_duplicate)>=1){
    
    total_duplicated_events <- total_duplicated_events+nrow(to_duplicate)
    
    to_keep_as_is <- annotation_subset[which(!(as.numeric(annotation_subset$T_start_UNIX)<=time_limit & as.numeric(annotation_subset$T_stop_UNIX)>time_limit) ),]
    to_duplicate_before             <-  to_duplicate
    to_duplicate_before$T_stop_UNIX <- as.POSIXct(time_limit,  origin="1970-01-01", tz="GMT" )
    to_duplicate_before$duration    <- as.numeric(to_duplicate_before$T_stop_UNIX - to_duplicate_before$T_start_UNIX)
    
    to_duplicate_after              <-  to_duplicate
    to_duplicate_after$T_start_UNIX <- as.POSIXct(time_limit,  origin="1970-01-01", tz="GMT" )+1/FRAME_RATE
    to_duplicate_after$duration     <- as.numeric(to_duplicate_after$T_stop_UNIX - to_duplicate_after$T_start_UNIX)
    
    annotation_subset <- rbind(to_keep_as_is,to_duplicate_before,to_duplicate_after)
    annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
  }
  
  return(list(annotation_subset=annotation_subset,total_duplicated_events=total_duplicated_events))
}


split_by_number<- function(annotations_all,time_window_all,seed){
  set.seed(seed)
  ##############################################################################
  ######### SPLIT ANNOTATIONS DATASET INTO TEST AND TRAINING  ##################
  ##############################################################################
  ### Be careful about this - in the previous you had randomly allocated behaviours to test or training throughout the period
  ### This means true Hits were wrongly identified as misses when looking at automatic interactions, because the time span of automatic interaction detection was unchanged
  ### Instead you need to define contiguous periods of time that contain half the events, for each colony/period
  print("Splitting manual annotations into two  chunks...")
  ###first list nb of events for each behaviour, each replicate and each period
  nb_events <- aggregate ( ant1 ~ Behaviour + PERIOD + REPLICATE, FUN=length, data=annotations_all)
  names(nb_events)[which(names(nb_events)=="ant1")] <- "Nb"
  ###narrow down to behaviour of interest only
  nb_events <- nb_events[which(nb_events$Behaviour==BEH),]
  
  ###initialise new test and training objects
  annotations_training <- NULL
  annotations_test     <- NULL
  time_window_training <- NULL
  time_window_test     <- NULL
  
  ###then loop over nb_events
  
  
  total_duplicated_events <- 0
  for (i in 1:nrow(nb_events)){
    ###subset annotations_all to period/replicate of interest
    PERIOD    <- nb_events[i,"PERIOD"]
    REPLICATE <- nb_events[i,"REPLICATE"]
    annotation_subset <- annotations_all[which(annotations_all$Behaviour==BEH & annotations_all$PERIOD==PERIOD & annotations_all$REPLICATE==REPLICATE),]
    
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to half the events for this period  
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to 3/4 of the events for this period  
    annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
    
    ### now we split annotation_subset in two chunks before and after time_limit,
    ### randomly allocate each time chunk to test or training dataset,
    ### and store the time windows for those time chunks
    if (runif(1,0,1)<0.5){
      
      time_limit <- min (c(as.numeric(annotation_subset[floor(proportion_training_events*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+proportion_training_events*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor(proportion_training_events*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor(proportion_training_events*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
      if (length(time_limit)==0){
        time_limit <- min (c(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+0.5*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor(0.5*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
        
      }
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      
      
      
      annotations_training <- rbind( annotations_training , annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),]  )
      annotations_test     <- rbind( annotations_test     , annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),])
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                       time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      )
      )
      
      time_window_test     <- rbind( time_window_test     , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit + 1/FRAME_RATE,origin="1970-01-01", tz="GMT") - 1/FRAME_RATE, 
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
    }else{
      time_limit <- min (c(as.numeric(annotation_subset[floor((1-proportion_training_events)*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+(1-proportion_training_events)*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor((1-proportion_training_events)*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor((1-proportion_training_events)*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
      if (length(time_limit)==0){
        time_limit <- min (c(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+0.5*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor(0.5*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
        
      }
      
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      annotations_training <- rbind( annotations_training, annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),]  )
      annotations_test     <- rbind( annotations_test,annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),])
      
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit+1/FRAME_RATE,origin="1970-01-01", tz="GMT")- 1/FRAME_RATE,
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
      
      
      time_window_test     <- rbind( time_window_test , data.frame(PERIOD=PERIOD, 
                                                                   REPLICATE=REPLICATE,
                                                                   time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                   time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      ))
      
    }
    
  }
  
  ###finally, recalculate duration for all annotations_all objects and define new time columns express in seconds
  for (annotation_object in c("annotations_all","annotations_training","annotations_test")){
    annot <- get(annotation_object)
    #convert Zulu time to GMT
    annot$duration      <- as.numeric(annot$T_stop_UNIX - annot$T_start_UNIX)
    # #transform zulu time in GMT
    # annot$T_start_UNIX <- as.POSIXct(annot$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    # annot$T_stop_UNIX  <- as.POSIXct(annot$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    #assign time in sec to avoid issues on time management and matching
    annot$T_start_sec <- round(as.numeric(annot$T_start_UNIX),N_DECIMALS)
    annot$T_stop_sec <- round(as.numeric(annot$T_stop_UNIX),N_DECIMALS)
    
    assign(annotation_object,annot)
  }
  print(paste("Had to duplicate",total_duplicated_events,"events"))
  return(list(annotations_all=annotations_all,annotations_training=annotations_training,annotations_test=annotations_test,time_window_all=time_window_all,time_window_training=time_window_training,time_window_test=time_window_test))
}


test_training_split_by_time <- function(annotations_all,time_window_all,seed){
  set.seed(seed)
  ##############################################################################
  ######### SPLIT ANNOTATIONS DATASET INTO TEST AND TRAINING  ##################
  ##############################################################################
  ### Be careful about this - in the previous you had randomly allocated behaviours to test or training throughout the period
  ### This means true Hits were wrongly identified as misses when looking at automatic interactions, because the time span of automatic interaction detection was unchanged
  ### Instead you need to define contiguous periods of time that contain half the events, for each colony/period
  print("Splitting manual annotations into test and training chunks...")
  
  ###initialise new test and training objects
  annotations_training <- NULL
  annotations_test     <- NULL
  time_window_training <- NULL
  time_window_test     <- NULL
  
  ###then loop over time_window_all
  total_duplicated_events <- 0
  for (i in 1:nrow(time_window_all)){
    ###subset annotations_all to period/replicate of interest
    PERIOD    <- time_window_all[i,"PERIOD"]
    REPLICATE <- time_window_all[i,"REPLICATE"]
    annotation_subset <- annotations_all[which(annotations_all$Behaviour==BEH & annotations_all$PERIOD==PERIOD & annotations_all$REPLICATE==REPLICATE),]
    
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to half the events for this period  
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to 3/4 of the events for this period  
    annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
    
    ### now we split annotation_subset in two chunks before and after time_limit,
    ### randomly allocate each time chunk to test or training dataset,
    ### and store the time windows for those time chunks
    time_chunk <- (as.numeric(time_window_all[i,"time_stop"])-as.numeric(time_window_all[i,"time_start"]))*proportion_training_events
    time_chunk <- round(time_chunk*FRAME_RATE)/FRAME_RATE
    if (runif(1,0,1)<0.5){ ###training = first_chunk
      time_limit <- as.numeric(time_window_all[i,"time_start"])+time_chunk
      
      
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      
      annotations_training <- rbind( annotations_training , annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),]  )
      annotations_test     <- rbind( annotations_test     , annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),])
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                       time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      )
      )
      
      time_window_test     <- rbind( time_window_test     , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit + 1/FRAME_RATE,origin="1970-01-01", tz="GMT") - 1/FRAME_RATE, 
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
    }else{
      time_limit <- as.numeric(time_window_all[i,"time_stop"])-time_chunk
      
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      annotations_training <- rbind( annotations_training, annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),]  )
      annotations_test     <- rbind( annotations_test,annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),])
      
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit+1/FRAME_RATE,origin="1970-01-01", tz="GMT")- 1/FRAME_RATE,
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
      
      
      time_window_test     <- rbind( time_window_test , data.frame(PERIOD=PERIOD, 
                                                                   REPLICATE=REPLICATE,
                                                                   time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                   time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      ))
      
    }
    
  }
  
  ###finally, recalculate duration for all annotations_all objects and define new time columns express in seconds
  for (annotation_object in c("annotations_all","annotations_training","annotations_test")){
    annot <- get(annotation_object)
    #convert Zulu time to GMT
    annot$duration      <- as.numeric(annot$T_stop_UNIX - annot$T_start_UNIX)
    # #transform zulu time in GMT
    # annot$T_start_UNIX <- as.POSIXct(annot$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    # annot$T_stop_UNIX  <- as.POSIXct(annot$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    #assign time in sec to avoid issues on time management and matching
    annot$T_start_sec <- round(as.numeric(annot$T_start_UNIX),N_DECIMALS)
    annot$T_stop_sec <- round(as.numeric(annot$T_stop_UNIX),N_DECIMALS)
    
    assign(annotation_object,annot)
  }
  return(list(annotations_all=annotations_all,annotations_training=annotations_training,annotations_test=annotations_test,time_window_all=time_window_all,time_window_training=time_window_training,time_window_test=time_window_test))
}

###############################################################################
###### PARAMETERS #############################################################
###############################################################################
###body length information
# all_body_lengths <-read.table(BODYLENGTH_FILE,header=T,stringsAsFactors = F,sep=",")

# #plotting limits used for the coordinates plotting
# Xmin <- 2000
# Xmax <- 7900
# Ymin <- 0000
# Ymax <- 5500
# 
#general
N_DECIMALS                  <- 3 ## when assigning time in seconds, the number of decimals to preserve when rounding to match between interactions & collisions
FUZZY_MATCH                 <- TRUE  ## fuzzy matching between data frames in collisions detection

#trajectories cutting gap, relevant for fmQueryComputeAntTrajectories
max_gap                     <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)

###NATH_FLAG: don't just use max and min, give some extra 
###NATH_FLAG: these parameters should be adjusted for each box as ant length may depend on focus /  camera distance
# AntDistanceSmallerThan      <- 300 #for higher accuracy, recalculate it from: max(interaction_MANUAL$straightline_dist_px,na.rm = T)
# AntDistanceGreaterThan      <- 70 #for higher accuracy, recalculate it from: min(interaction_MANUAL$straightline_dist_px,na.rm = T)
# ANT_LENGHT_PX               <- 153 #useful for matcher::meanAntDisplacement mean and median value are similar
# minimumGap                  <- fmSecond(0.5) ## THIS OPTION DOES NOT WORK IN INTERACTIONS SO DISABLED! for a given pair of interacting ants, when interaction is interrupted by more than minimumGap, interaction will check whether ants have moved since - and if so, will create new interaction

DISAGREEMENT_THRESH <- 0.5
###Fixed parameter

###############################################################################
###READ WHOLE ANNOTATIONS DATASET #############################################
###############################################################################
print("Loading manual annotations...")
###READ ANNOTATIONS
trophy_annotations <- read.csv(paste0(DATADIR,"trophallaxis_annotation.csv"), header = TRUE, stringsAsFactors = F)
###edit file to remove accidental capital letters from colony_id
trophy_annotations$colony_id <- decapitalize(trophy_annotations$colony_id)
##exclude all annotations whichy have NA in either actor or partner
trophy_annotations <- trophy_annotations[which(!is.na(trophy_annotations$actor) & !is.na(trophy_annotations$partner)),]
trophy_annotations$Behaviour <- "T"
# head(trophy_annotations)
if (exclude_queen){
  trophy_annotations <- trophy_annotations[which(trophy_annotations$actor!=queen_id & trophy_annotations$partner!=queen_id),]
}

###READ METADATA
meta_data <- read.table(file=paste(WORKDIR,"trophy_metadata.txt",sep=""),header=T,stringsAsFactors = F)
###sometimes the time stop in the annotations went beyond the actual end of the pre-defined time window, so we need to truncate that time in the annotation dataframe 
trophy_annotations        <- merge(trophy_annotations, meta_data[c("colony_id","annotation_start","annotation_stop")],all.x=T)
trophy_annotations[which(trophy_annotations$T_stop>trophy_annotations$annotation_stop),"T_stop"] <- trophy_annotations[which(trophy_annotations$T_stop>trophy_annotations$annotation_stop),"annotation_stop"]
trophy_annotations[which(trophy_annotations$T_start<trophy_annotations$annotation_start),"T_start"] <- trophy_annotations[which(trophy_annotations$T_start<trophy_annotations$annotation_start),"annotation_start"]
trophy_annotations       <-trophy_annotations[, which(!names(trophy_annotations)%in%c("annotation_start","annotation_stop"))]
###create T_start_UNIX and 
trophy_annotations$T_start_UNIX  <- as.POSIXct(trophy_annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
trophy_annotations$T_stop_UNIX   <- as.POSIXct(trophy_annotations$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
trophy_annotations$T_start <- NULL
trophy_annotations$T_stop <- NULL
trophy_annotations$T_start_sec <-  round(as.numeric(trophy_annotations$T_start_UNIX),N_DECIMALS)
trophy_annotations$T_stop_sec <-  round(as.numeric(trophy_annotations$T_stop_UNIX),N_DECIMALS)

###modify some column names
trophy_annotations$PERIOD <- "after"
names(trophy_annotations)[which(names(trophy_annotations)=="colony_id")] <- "REPLICATE"
names(trophy_annotations)[which(names(trophy_annotations)=="actor")] <- "ant1"
names(trophy_annotations)[which(names(trophy_annotations)=="partner")] <- "ant2"

###add duration column
trophy_annotations$duration <- as.numeric(trophy_annotations$T_stop_UNIX - trophy_annotations$T_start_UNIX)

###DEFINE FOCAL LIST
###define focal list. This should inlcude the "actor" from metadata for all experimental colonies + all ants for the additional "trophy" colony
focal_list_trophy <-  with(meta_data[which(!is.na(meta_data$actor)),], paste(colony_id,actor,sep="_")) ###first, use meta_data to extract focal ant for eahc of the experimental colonies
##find myrmidon file for additional colony
myrmidon_files <- list.files(path=DATADIR,pattern="myrmidon")
myrmidon_files <- myrmidon_files[ grepl(   meta_data[which(is.na(meta_data$actor)),"colony_id"] ,   myrmidon_files)]
myrmidon_files <- myrmidon_files[grepl("final",myrmidon_files)]
myrmidon_file <- paste(DATADIR,myrmidon_files[!grepl("Capsule",myrmidon_files)]  ,sep="")
e <- fmExperimentOpen(myrmidon_file)
replicate_full_ant_list <- unlist(lapply(e$ants,function(x)x$ID))
###exclude queen
if (exclude_queen){
  replicate_full_ant_list <-replicate_full_ant_list[which(replicate_full_ant_list!=queen_id)]
}

focal_list_trophy <- c(focal_list_trophy,     paste( meta_data[which(is.na(meta_data$actor)),"colony_id"],replicate_full_ant_list,sep="_")   )

##CREATE TIME_WINDOW OBJECT
time_window_trophy        <- data.frame(PERIOD="after",REPLICATE=meta_data$colony_id,time_start=as.POSIXct(meta_data$annotation_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ),time_stop=as.POSIXct(meta_data$annotation_stop, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ))

if (training_set=="Daniel"){
  FOCAL             <- T
  
  ##STEP 1 - split annotations for whole-occurrence colony into a fitting test and a validation test. Validation test will be only reused in step 3
  proportion_training_events <- 0.75 ##PROPORTION OF ANNOTATED TIME WE WANT TO USE FOR TRAINING THE CLASSIFIER
  validation_split <- split_by_number(trophy_annotations[which(trophy_annotations$REPLICATE=="c29"),],time_window_trophy[which(time_window_trophy$REPLICATE=="c29"),],seed=1)
  
  annotations_fitting     <- validation_split[["annotations_training"]]
  annotations_validation  <- validation_split[["annotations_test"]]
  time_window_fitting    <- validation_split[["time_window_training"]]
  time_window_validation   <- validation_split[["time_window_test"]]
  
  ##STEP 2 - exclude annotations_validation and time_window_validation from the annotations and time_window objects, so we ignore them for the rest of this script
  annotations_all        <- rbind (  trophy_annotations[which(trophy_annotations$REPLICATE!="c29"),]   ,   annotations_fitting     )
  time_window_all        <- rbind (  time_window_trophy[which(time_window_trophy$REPLICATE!="c29"),]   ,   time_window_fitting     )
  all_training_ants      <- focal_list_trophy
}

# ###make sure annotations don't overrun 30 minute window
# for (i in 1:nrow(annotations_all)){
#   annotations_all[i,"T_stop_UNIX"] <- min(annotations_all[i,"T_stop_UNIX"],time_window_all[which(time_window_all$REPLICATE==annotations_all[i,"REPLICATE"]&time_window_all$PERIOD==annotations_all[i,"PERIOD"]),"time_stop"]-1/FRAME_RATE)
#   annotations_all[i,"T_start_UNIX"] <- max(annotations_all[i,"T_start_UNIX"],time_window_all[which(time_window_all$REPLICATE==annotations_all[i,"REPLICATE"]&time_window_all$PERIOD==annotations_all[i,"PERIOD"]),"time_start"]+1/FRAME_RATE)
# }
annotations_all$duration      <- as.numeric(annotations_all$T_stop_UNIX - annotations_all$T_start_UNIX)
annotations_all               <- annotations_all[which(annotations_all$duration>0),] 

###finally, subset annotations_all for the behaviour of interest
annotations_all <- annotations_all[which(annotations_all$Behaviour==BEH),]

###split into test and training sets
proportion_training_events <- 0.8 
random_split <- test_training_split_by_time(annotations_all,time_window_all,seed=2)
# random_split <- test_training_split_by_time(annotations_all,time_window_all)

annotations_all      <- random_split[["annotations_all"]]
annotations_training <- random_split[["annotations_training"]]
annotations_test     <- random_split[["annotations_test"]]
time_window_all      <- random_split[["time_window_all"]]
time_window_training <- random_split[["time_window_training"]]
time_window_test     <- random_split[["time_window_test"]]

#start fresh
Grooming_LDA_output             <- data.frame()
Grooming_LDA_eachRun            <- data.frame()

###initialise general output folder
##remove folder if already exists to amke sure we don't mix things up
if(FLUSH){
  if (file.exists(file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL"))){
    unlink(file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL"),recursive=T)
  }
}
###create folder
dir.create(file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL"),recursive = T)

###define name of general output table containing quality scores
output_name <- file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",paste("quality_scores_",ITER,".txt",sep=""))
if (file.exists(output_name)){
  Grooming_LDA_output <- read.table(output_name,stringsAsFactors = F,header=T)
  partial <- T
  classifier_LIST <- unique(Grooming_LDA_output$classifier)
}else{
  partial <- F
  classifier_LIST <- NULL
}



###############################################################################
###### OUTER PARAMETERS LOOP ##################################################
###############################################################################
#Varying capsule shapes
# CAPSULE_FILE_LIST           <- paste("CapsuleDef0",c(0:8),sep="") 
CAPSULE_FILE_LIST           <- paste("CapsuleDef0",c(3:8),sep="") 
MAX_INTERACTION_GAP_LIST    <- c(25,30,35)
DT_dist_THRESHOLD_BL_LIST   <- c(0,0.0005,0.001,0.0015,0.002,0.0025,0.003)
DT_frame_THRESHOLD_LIST     <- c(32,40,48)
trim_length_sec_LIST        <- c(0,3,4,6)
beta_LIST                   <- c(1)
# #####Arguments to loop over - to comment out when running the loop
# CAPSULE_FILE                <- "BODY-HEAD_17March22_auto_oriented_meanlength"
# DT_dist_THRESHOLD_BL           <- 0 # NOT HIGHER THAN 0.5 as it will cut a very large portion of data (most movements are on a very small scale) #tag length is 62 px approx (measured on full size pics in R9SP)
# MAX_INTERACTION_GAP         <- 10 ## in SECONDS, the maximum gap in tracking before cutting the trajectory or interactions in two different object
# DISAGREEMENT_THRESH         <- 0.4#Assign Hit or Miss after matrix difference according to threshold 
# trim_length_sec             <- 1  ####all automated interactions that last trim_length_sec or less will be removed from the automated detection to fit the LDA as they increase noise!
# DT_frame_THRESHOLD          <- 16  #arbitrary value #trajectories jumps/gaps thresholds to avoid getting skewed means summary values (see their use in params extraction scripts: BEH_Traj_from_Man_annotations_fort081.R and BEH_Parameters_Auto_fort081.R)
####

####initialise Loop_ID
Loop_ID       <- (ITER-1)*1000+1
Trunk_Loop_ID  <- (ITER-1)*1000+1


CAPSULE_FILE_LIST_operational         <- CAPSULE_FILE_LIST[ceiling(ITER / (length(MAX_INTERACTION_GAP_LIST)))]
MAX_INTERACTION_GAP_LIST_operational  <- MAX_INTERACTION_GAP_LIST[  ITER - length(MAX_INTERACTION_GAP_LIST)* floor((ITER-1) / (length(MAX_INTERACTION_GAP_LIST))  )   ]
# MAX_INTERACTION_GAP_LIST_operational  <- MAX_INTERACTION_GAP_LIST
# MAX_INTERACTION_GAP_LIST_operational  <- MAX_INTERACTION_GAP_LIST
DT_dist_THRESHOLD_BL_LIST_operational <- DT_dist_THRESHOLD_BL_LIST
for (CAPSULE_FILE in CAPSULE_FILE_LIST_operational) { #list of CAPUSLE FILES TO BE USED  ###NATH_FLAG: CAPSULE_FILE_LIST has not been defined
  if (
    (partial==F) 
    |
    (
      (partial==T)
      &
      (max(c(0,which( Grooming_LDA_output[nrow(Grooming_LDA_output),"CAPSULE_FILE"]==CAPSULE_FILE_LIST_operational)))<=which(CAPSULE_FILE==CAPSULE_FILE_LIST_operational))
      &
      (length(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE))< (length(MAX_INTERACTION_GAP_LIST_operational)*length(DT_dist_THRESHOLD_BL_LIST_operational)*length(DT_frame_THRESHOLD_LIST)*length(trim_length_sec_LIST)*length(classifier_LIST)*length(beta_LIST)))
    )
  ){
    
    
    
    #Sequentially vary the interaction gap-filling to check what effect this has on the agreement between the MANUAL & AUTOMATIC interactions
    for (MAX_INTERACTION_GAP in MAX_INTERACTION_GAP_LIST_operational) { #IN SECONDS (5,10, 15 never selected)
      trunk_loop_start_time <- Sys.time()
      print(paste("TRUNK LOOP ID:",Trunk_Loop_ID))
      
      if (
        (partial==F)
        |
        (
          (partial==T)
          &
          (max(c(0,which(Grooming_LDA_output[max(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE)),"MAX_INTERACTION_GAP"]==MAX_INTERACTION_GAP_LIST_operational)))<=which(MAX_INTERACTION_GAP==MAX_INTERACTION_GAP_LIST_operational))
          &
          (length(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP)) < (length(DT_dist_THRESHOLD_BL_LIST_operational)*length(DT_frame_THRESHOLD_LIST)*length(trim_length_sec_LIST)*length(classifier_LIST)*length(beta_LIST)))
        )
      ){
        
        ##################################################################################
        ##################### EXTRACT VARIABLES FOR MANUAL AND AUTOMATIC INTERACTIONS ####
        ##################################################################################
        ###for these particular parameters, evaluate how successful the interaction detetection parameters are at detecting candidate frames
        ###extract interactions and manual annotations for the whole dataset
        print(paste("Evaluating quality of interaction parameters for trunk loop ID ",Trunk_Loop_ID,"...",sep=""))
        all      <- extraction_loop("all",extract_movement_variables=F,all_body_lengths=meta_data,focal=FOCAL,focal_list=all_training_ants)
        
        ###Run auto_manual_agreement on all data to get an idea of how well the Loop performs in discovering the candidate interaction frames
        loop_interaction_detection_TPTNFPFN <- auto_manual_agreement (all[["summary_AUTO"]] , all[["summary_MANUAL"]], all[["list_IF_Frames"]],all[["list_replicate_full_ant_list"]],all[["list_replicate_focal_list"]]  )[["true_false_positive_negatives"]]
        
        ##THRESHOLD to exclude jitter in the individuals' movement (DISTANCE)
        for (DT_dist_THRESHOLD_BL in DT_dist_THRESHOLD_BL_LIST_operational){
          
          if (
            (partial==F)
            |
            (
              (partial==T)
              &
              (max(c(0,which(Grooming_LDA_output[max(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP)),"DT_dist_THRESHOLD_BL"]==DT_dist_THRESHOLD_BL_LIST_operational)))<=which(DT_dist_THRESHOLD_BL==DT_dist_THRESHOLD_BL_LIST_operational))
              &
              (length(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP&Grooming_LDA_output$DT_dist_THRESHOLD_BL==DT_dist_THRESHOLD_BL)) < (length(DT_frame_THRESHOLD_LIST)*length(trim_length_sec_LIST)*length(classifier_LIST)*length(beta_LIST)))
            )
          ){
            
            
            
            # Assign Hit based on threshold
            # maybe can be put somewhere better as it involves a later stage of the analysis
            for (DT_frame_THRESHOLD in DT_frame_THRESHOLD_LIST){
              if (
                (partial==F)
                |
                (
                  (partial==T)
                  &
                  (max(c(0,which(Grooming_LDA_output[max(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP&Grooming_LDA_output$DT_dist_THRESHOLD_BL==DT_dist_THRESHOLD_BL)),"DT_frame_THRESHOLD"]==DT_frame_THRESHOLD_LIST)))<=which(DT_frame_THRESHOLD==DT_frame_THRESHOLD_LIST))
                  &
                  (length(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP&Grooming_LDA_output$DT_dist_THRESHOLD_BL==DT_dist_THRESHOLD_BL&Grooming_LDA_output$DT_frame_THRESHOLD==DT_frame_THRESHOLD)) < (length(trim_length_sec_LIST)*length(classifier_LIST)*length(beta_LIST)))
                )
              ){
                
                ###for these particular parameters, calculate variables of interest for manual annotations and automatic interactions for each of training and test dataset
                print("Extracting movement variables for manually-annotated data and automatic interactions:")
                print("-training data...")
                training <- extraction_loop("training",all_body_lengths=meta_data,focal=FOCAL,focal_list=all_training_ants)
                print("-test data...")
                test     <- extraction_loop("test",all_body_lengths=meta_data,focal=FOCAL,focal_list=all_training_ants)
                
                ##############################################################################
                ######### MANUAL/AUTO AGREEMENT  ###################################
                ##############################################################################
                ###Run auto_manual_agreement on training data to define Hits and Misses
                training[["summary_AUTO"]] <- auto_manual_agreement (training[["summary_AUTO"]] , training[["summary_MANUAL"]], training[["list_IF_Frames"]],training[["list_replicate_full_ant_list"]],training[["list_replicate_focal_list"]] )[["summary_AUTO"]]
                
                ##############################################################################
                ######### FIT CLASSIFIERS WITH LDA/QDA/RF  ####################################
                ##############################################################################
                ###Note for Adriano:
                ###For the fit of the classifiers I believe a lot can be gained from trimming the shortest interactions from summary_AUTO, as these will introduce noise
                ####define to_keep variables to keep clearing memory between runs
                to_keep <- c(ls(),c("to_keep","trim_length_sec"))
                
                for (trim_length_sec in trim_length_sec_LIST){##This level of the loop is AFTER extracting movement variables for training and test to avoid unnecessary repetition of this slow step.
                  if(
                    (partial==F)
                    |
                    (
                      (partial==T)
                      &
                      (max(c(0,which(Grooming_LDA_output[max(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP&Grooming_LDA_output$DT_dist_THRESHOLD_BL==DT_dist_THRESHOLD_BL&Grooming_LDA_output$DT_frame_THRESHOLD==DT_frame_THRESHOLD)),"trim_length_sec"]==trim_length_sec_LIST)))<=which(DT_frame_THRESHOLD==DT_frame_THRESHOLD_LIST))
                      &
                      (length(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP&Grooming_LDA_output$DT_dist_THRESHOLD_BL==DT_dist_THRESHOLD_BL&Grooming_LDA_output$DT_frame_THRESHOLD==DT_frame_THRESHOLD)) < (length(trim_length_sec_LIST)*length(classifier_LIST)*length(beta_LIST)))
                    )
                  ){
                    
                    
                    
                    
                    print(paste("Performing LOOP ID: ",Loop_ID))
                    ###prepare output directories for Loop_ID
                    subDir <- paste0("Loop_ID_",Loop_ID)
                    if (file.exists(file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",subDir))){
                      unlink(file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",subDir),recursive=T)
                    }
                    dir.create(file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",subDir,"fits"),recursive = T)
                    
                    ###fit classifiers
                    print(paste("Perform Classification for Loop_ID ",Loop_ID,"...",sep=""))
                    classifiers <- fit_classifiers(
                      trim_short_interactions ( ######trim_short_interactions: function that cuts interactions shorter than a certain duration
                        training[["summary_AUTO"]] ######interaction table: summary_AUTO from training object 
                        ,trim_length_sec ### maximum duration of interactions that will be trimmed (set to 0 to trim nothing)
                        ,"duration_sec" ###name of the column that contains the duration of each interaction in seconds (could vary between object so specify it here)
                        
                      )
                    )
                    
                    ###############################################################################
                    ###### SAVING LOOP-RELEVANT OBJECTS        ####################################
                    ###############################################################################
                    #save the SIRUS rules object
                    dput(classifiers[["SirusRules"]], file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",subDir,"SirusRules.txt"))
                    
                    #save the BN_list object containing all information regarding the variables selected by RELIEF and the method to calculate them for new values
                    dput(classifiers[["selected_variables_BNobject_list"]], file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",subDir,"BN_object_list.dat"))
                    
                    ####################################################################################################################
                    ###then for each classifier method, perform prediction on training and test data and get quality scores    #########
                    ####################################################################################################################
                    print(paste("Predict training and test data for Loop_ID ",Loop_ID,"...",sep=""))
                    
                    for (class_idx in 1:length(classifiers[["fits"]])){
                      if (!is.null(classifiers[["fits"]][[class_idx]])){#if non null
                        classifier <- list(classifiers[["fits"]][[class_idx]]); names(classifier) <- names(classifiers[["fits"]])[class_idx]
                        
                        ##check match with manual 
                        for (beta in c(1)){
                          ###beta: relates to the calculation of the Fbeta score (https://en.wikipedia.org/wiki/F-score)
                          ###      This is how much we value recall (=sensitivity = TP / (TP+FN)) over precision (=positive predictive value = TP / (TP + FP))
                          ### if beta = 1 the score calculated is a F1 score
                          ### based on lab meeting discussion we felt that we valued precision (minimising FP) over sensitivity to decrease the noise / avoid detecting grooming when there are none
                          ### so we want a beta < 1
                          ### arbitrarily here, it is set at 0.5
                          
                          quality_scores_pre_classifier       <- quality_scores(loop_interaction_detection_TPTNFPFN,beta)
                          if (class_idx==1){
                            print(paste("Automatic interaction detection in trunk loop ID ",Trunk_Loop_ID," had a sensitivity of ",round(quality_scores_pre_classifier["sensitivity"],digits=2)," and a precision of ",round(quality_scores_pre_classifier["precision"],digits=2)," PRE CLASSIFIER.",sep=""))
                          }
                          
                          to_keep2 <- ls()
                          for (what in c("training","test")){
                            ###extract objects of interest
                            summary_AUTO   <- get(what)[["summary_AUTO"]]
                            summary_MANUAL <- get(what)[["summary_MANUAL"]]
                            list_IF_Frames <- get(what)[["list_IF_Frames"]]
                            list_replicate_full_ant_list <- get(what)[["list_replicate_full_ant_list"]]
                            list_replicate_focal_list <- get(what)[["list_replicate_focal_list"]]
                            
                            ###Perform prediction on training data
                            summary_AUTO["predicted_Hit"] <- predict_class  (summary_AUTO =    summary_AUTO
                                                                             ,BN_list      =  classifiers[["selected_variables_BNobject_list"]]
                                                                             ,classifier   = classifier
                            )
                            
                            
                            ####get fit between manual and predicted groomings, both overall and separately for each treatment
                            ###overall
                            
                            if (sum(summary_AUTO$predicted_Hit,na.rm=T)>0){
                              aut_man_agreement <- auto_manual_agreement (summary_AUTO[which(summary_AUTO$predicted_Hit==1),]
                                                                          , summary_MANUAL
                                                                          , list_IF_Frames
                                                                          , list_replicate_full_ant_list
                                                                          , list_replicate_focal_list
                              )
                              assign (paste("true_false_positive_negatives",what,sep="_"),aut_man_agreement[["true_false_positive_negatives"]])
                              assign (paste("quality_scores",what,sep="_"),round(quality_scores(aut_man_agreement[["true_false_positive_negatives"]],beta),digits=3 ))
                            }else{
                              assign (paste("quality_scores",what,sep="_"),round(c(CSI=0,Fbeta=0,precision=0,sensitivity=0),digits=3 ))
                              assign (paste("true_false_positive_negatives",what,sep="_"),data.frame(true_negatives=NA,true_positives=0,false_negatives=NA,false_positives=0))
                            }
                            
                            
                            # if(what=="test"){
                            print(paste("Automatic interaction detection in trunk loop ID ",Trunk_Loop_ID," has a sensitivity of ",round(get(paste("quality_scores",what,sep="_"))["sensitivity"],digits=2)," and a precision of ",round(get(paste("quality_scores",what,sep="_"))["precision"],digits=2)," POST CLASSIFIER ",class_idx," in the ",what," data.",sep=""))
                            # }
                          }
                          #create a row per each Loop run 
                          Grooming_LDA_eachRun  <- data.frame(Loop_ID=Loop_ID,
                                                              CAPSULE_FILE=CAPSULE_FILE,
                                                              DT_dist_THRESHOLD_BL=DT_dist_THRESHOLD_BL, 
                                                              MAX_INTERACTION_GAP=MAX_INTERACTION_GAP,
                                                              DISAGREEMENT_THRESH=DISAGREEMENT_THRESH,
                                                              DT_frame_THRESHOLD = DT_frame_THRESHOLD,
                                                              trim_length_sec=trim_length_sec,
                                                              beta=beta,
                                                              classifier=names(classifier),
                                                              CSI_pre_classifier = quality_scores_pre_classifier["CSI"],
                                                              Fbeta_pre_classifier  = quality_scores_pre_classifier["Fbeta"],
                                                              precision_pre_classifier = quality_scores_pre_classifier["precision"],
                                                              sensitivity_pre_classifier  = quality_scores_pre_classifier["sensitivity"],                                                  
                                                              CSI_training = quality_scores_training["CSI"],
                                                              Fbeta_training  = quality_scores_training["Fbeta"],
                                                              precision_training = quality_scores_training["precision"],
                                                              sensitivity_training  = quality_scores_training["sensitivity"],
                                                              CSI_test = quality_scores_test["CSI"],
                                                              Fbeta_test  = quality_scores_test["Fbeta"],
                                                              precision_test = quality_scores_test["precision"],
                                                              sensitivity_test  = quality_scores_test["sensitivity"]
                                                              ,
                                                              # 
                                                              # CSI_training_CrossTreatmentMean = quality_scores_training_CrossTreatmentMean["CSI"],
                                                              # Fbeta_training_CrossTreatmentMean  = quality_scores_training_CrossTreatmentMean["Fbeta"],
                                                              # precision_training_CrossTreatmentMean = quality_scores_training_CrossTreatmentMean["precision"],
                                                              # sensitivity_training_CrossTreatmentMean  = quality_scores_training_CrossTreatmentMean["sensitivity"],
                                                              # CSI_test_CrossTreatmentMean = quality_scores_test_CrossTreatmentMean["CSI"],
                                                              # Fbeta_test_CrossTreatmentMean  = quality_scores_test_CrossTreatmentMean["Fbeta"],
                                                              # precision_test_CrossTreatmentMean = quality_scores_test_CrossTreatmentMean["precision"],
                                                              # sensitivity_test_CrossTreatmentMean  = quality_scores_test_CrossTreatmentMean["sensitivity"],
                                                              # 
                                                              # CSI_training_CrossTreatmentsd = quality_scores_training_CrossTreatmentsd["CSI"],
                                                              # Fbeta_training_CrossTreatmentsd  = quality_scores_training_CrossTreatmentsd["Fbeta"],
                                                              # precision_training_CrossTreatmentsd = quality_scores_training_CrossTreatmentsd["precision"],
                                                              # sensitivity_training_CrossTreatmentsd  = quality_scores_training_CrossTreatmentsd["sensitivity"],
                                                              # CSI_test_CrossTreatmentsd = quality_scores_test_CrossTreatmentsd["CSI"],
                                                              # Fbeta_test_CrossTreatmentsd  = quality_scores_test_CrossTreatmentsd["Fbeta"],
                                                              # precision_test_CrossTreatmentsd = quality_scores_test_CrossTreatmentsd["precision"],
                                                              # sensitivity_test_CrossTreatmentsd  = quality_scores_test_CrossTreatmentsd["sensitivity"]
                                                              # ,
                                                              
                                                              stringsAsFactors = F,row.names = NULL)
                          
                          ###Add this line in general output table
                          if (partial ==F){
                            Grooming_LDA_output   <- rbind(Grooming_LDA_output,     Grooming_LDA_eachRun)
                            
                          }else{
                            if (
                              length(which(Grooming_LDA_output$CAPSULE_FILE==CAPSULE_FILE&Grooming_LDA_output$MAX_INTERACTION_GAP==MAX_INTERACTION_GAP&Grooming_LDA_output$DT_dist_THRESHOLD_BL==DT_dist_THRESHOLD_BL&Grooming_LDA_output$DT_frame_THRESHOLD==DT_frame_THRESHOLD&Grooming_LDA_output$trim_length_sec==trim_length_sec&Grooming_LDA_output$beta==beta&Grooming_LDA_output$classifier==names(classifier)))==0
                            ){
                              Grooming_LDA_output   <- rbind(Grooming_LDA_output,     Grooming_LDA_eachRun)
                              partial  <- F
                            }
                          }
                          
                          ###############################################################################
                          ###### SAVING CLASSIFIER-RELEVANT OBJECTS        ##############################
                          ###############################################################################
                          ###1. save new line in quality score object, then clear it
                          if (partial==F){
                            if (file.exists(output_name)){
                              write.table(Grooming_LDA_eachRun,file=output_name,append=T,col.names=F,row.names=F,quote=T)
                            }else{
                              write.table(Grooming_LDA_eachRun,file=output_name,append=F,col.names=T,row.names=F,quote=T)
                            }
                          }
                          rm(list=ls()[which(!ls()%in%to_keep2)])
                        }
                        
                        ###2. save fit
                        saveRDS(classifier[[names(classifier)]], file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",subDir,"fits",paste(names(classifier),".rds",sep="")))
                        
                      }#if non null
                    }#class_idx
                  }
                  ###add 1 to Loop_ID counter
                  Loop_ID <- Loop_ID +  1
                  rm(   list   =  ls()[which(!ls()%in%to_keep)]    )
                  gc()
                }#trim_length_sec
              }else{
                Loop_ID <- Loop_ID + length(trim_length_sec_LIST)
              }
            }#DT_frame_THRESHOLD
          }else{
            Loop_ID <- Loop_ID + length(DT_frame_THRESHOLD_LIST)*length(trim_length_sec_LIST)
          }
        }#DT_dist_THRESHOLD_BL
      }else{
        Loop_ID <- Loop_ID + length(DT_dist_THRESHOLD_BL_LIST_operational)*length(DT_frame_THRESHOLD_LIST)*length(trim_length_sec_LIST)
      }
      Trunk_Loop_ID <- Trunk_Loop_ID +1
      trunk_loop_end_time <- Sys.time()
      print (paste("Trunk loop took ",trunk_loop_end_time-trunk_loop_start_time," seconds to complete"))
    }#MAX_INTERACTION_GAP
  }else{
    Loop_ID <- Loop_ID + length(MAX_INTERACTION_GAP_LIST_operational)*length(DT_dist_THRESHOLD_BL_LIST_operational)*length(DT_frame_THRESHOLD_LIST)*length(trim_length_sec_LIST)
    Trunk_Loop_ID <- Trunk_Loop_ID + length(MAX_INTERACTION_GAP_LIST_operational)
  }
}#CAPSULE_FILE
rm(list=ls())
gc()
mallinfo::malloc.trim(0L)