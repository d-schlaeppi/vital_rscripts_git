# Task definition

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### TASK EXTRACTION ADJUSTED FOR DANIEL ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####
# based on Adriano, Nathalie, Tom and Linda Scripts. 

# THIS SCRIPT IS USED INSIDE: Extract_Metadata_v082.R
# In this script we access the zone definitions of the myrmidon files: i.e. nest and arena
# Then, space use is calculated for each ant and based on the amount of time spent outside (needs to be adjusted for each user) we define foragers and nurses
# The is written as a function that takes as argument a single myrmidon file referred to as "e" so it can be run in a loop in metadata extraction for each of your colonies
# The output is a Data frame that can be accessed for the Metadata extraction. 

# Note colony_metadata should be loaded. however that should already be the case as it was loaded in extract_metadata

AntTasks_DS <- function(e, file){
  print(paste0("Computing AntTasks based on 24h time-window before exposure for", file)) # in my case I always removed ants for the feeding treatment at 09:00 and thus ant task is calculated for the 24h windo before that
  #access zones: 1 = nest and 2 = arena  (Could be expanded to other EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
  zones <- e$spaces[[1]]$zones #function to show the Zones present in the Space
  zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
  foraging_zone <- zones_tab[which(grepl("arena",zones_tab$name)),"ID"] #fool-proofing - this way we are sure to always select the right zone
  nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"] #fool-proofing - this way we are sure to always select the right zone
  print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))
  #Get complete list of Ants
  AntID_list <- NULL
  for (ant in   1: length(e$ants)) {
      if(length(e$ants[[ant]]$identifications) == 0) {next}  # skip potential errors due to missing ants (i.e. a lost tag that was miss-identified as an ant and then removed during manual post processing results in a missing ant-id)
    AntID_list <- rbind(AntID_list,
                        data.frame (
                          antID = e$ants[[ant]]$ID,
                          tag_hex_ID = e$ants[[ant]]$identifications[[1]]$tagValue,
                          alive = e$ants[[ant]]$getValue("IsAlive",fmTimeNow()),
                          treated = e$ants[[ant]]$getValue("IsTreated",fmTimeNow())
                        ))
  }
  
  # extract treatment_start_time  #### At some point get true tracking stop time. for now set end time 3 hours later treatment start time
  colony_id <- sub(".*/final_(c\\d{2})\\.myrmidon", "\\1", file)
  time_treatment_start <- colony_metadata[which (grepl(paste0(colony_id) ,colony_metadata$colony_id)), "time_treatment_start" ] # get treatment start time 
  time_treatment_start <- as.POSIXct(time_treatment_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) # transform it into time R understands
  # time_exp_end <- time_treatment_start + 3*3600
  time_start_taskdef <- time_treatment_start -24*3600 # get the start time for the task definition 24 h before the treatment (after acclimatization period)
  time_stop_taskdef <- time_treatment_start -2*3600 # get the stop for time for task definition (just before workers get removed for treatment)
  # transform back into time objects for FortMyrmidon
  time_start <- fmTimeCreate(time_start_taskdef)
  time_end  <- fmTimeCreate(time_stop_taskdef)
  
  # get information on frames
  IdentifyFrames      <- fmQueryIdentifyFrames(e,start=time_start, end=time_end,showProgress = FALSE)
  IF_frames           <- IdentifyFrames$frames
  IF_frames$frame_num <- as.numeric(seq.int(nrow(IF_frames))) # Assign a frame to each time since start and use it as baseline for all matching and computation
  IF_frames$time_sec <- round(as.numeric(IF_frames$time),3) # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
  
  # extract trajectories in the selected time period (for me 22h) using Nathalie's function
  positions <- extract_trajectories(e = e,
                                    start = time_start,
                                    end = time_end,
                                    maximumGap = fmHour(24*365) ,
                                    IF_frames = IF_frames) 
  positions_summaries       <- positions$trajectories_summary
  positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
  positions_list            <- positions$trajectories
  
  # before going back and forth between positions_summaries and positions_list:
  nrow(positions_summaries)==length(positions_list)                             # 1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
  positions_summaries <- positions_summaries[order(positions_summaries$index),] # 2/ always make sure that positions_summaries is ordered correctly, using the index column
  # this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
  
  # calculate the zone usage for each ant (number of frames in nest or arena) and from that define nurse or foragers (and queen)
  for ( ant_index in 1:length(positions_list)) {
    positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]%in%foraging_zone))
    positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]%in%nest_zone))
  }
  positions_summaries$nb_frames_inside[which(is.na(positions_summaries$nb_frames_inside))] <-0
  positions_summaries$nb_frames_outside[which(is.na(positions_summaries$nb_frames_outside))] <-0
  positions_summaries$prop_time_outside <- positions_summaries$nb_frames_outside/(positions_summaries$nb_frames_outside+positions_summaries$nb_frames_inside)
  
  positions_summaries[which(positions_summaries$prop_time_outside<=0.01),"AntTask"] <- "nurse" # in this case the threshold is set so that ants that spend more than 1% of their time outside of the arena are considered foragers
  positions_summaries[which(positions_summaries$prop_time_outside>0.01),"AntTask"]  <- "forager"
  positions_summaries[which(positions_summaries$antID==1),"AntTask"]                <- "queen"
  
  merged_data <- merge(positions_summaries, AntID_list, by = "antID", all.x = TRUE)
  positions_summaries$tag_hex_ID <- merged_data$tag_hex_ID

  AntTasks <- data.frame(antID=positions_summaries[,"antID"],
                         tag_hex_ID=positions_summaries[,"tag_hex_ID"],
                         AntTask= positions_summaries[,"AntTask"],
                         prop_time_outside= positions_summaries[,"prop_time_outside"])
  
  print("AntTasks computed")
  
  # add missing ants and define them as nurses per definition
  missing_ants <- subset(AntID_list$antID, !(AntID_list$antID %in% AntTasks$antID))
  missing_ants_table <- data.frame()
  for (MISSING in missing_ants) {
    for (id in length(e$ants[[MISSING]]$identifications)) {
      #print ( e$ants[[MISSING]]$identifications[[id]]$tagValue )
      missing_ants_table <- rbind(missing_ants_table, data.frame(antID=MISSING, tag_hex_ID= e$ants[[MISSING]]$identifications[[id]]$tagValue ,AntTask="nurse", prop_time_outside= NA))
    }}
  AntTasks <- rbind(AntTasks, missing_ants_table)
  
  # add misssing ants  as NURSE by DEFAULT
  # AntTasks[which(is.na(AntTasks$AntTask)),"AntTask"] <- "nurse"
  # warning("Ants that died before the considered time window (pre treatment) will not be assigned a Task by the function. Currently, no task will default to Nurse")
  
  AntTasks$AntTask_num <- NA
  AntTasks[which(AntTasks$AntTask=="nurse"),"AntTask_num"] <- 1
  AntTasks[which(AntTasks$AntTask=="forager"),"AntTask_num"] <- 2
  AntTasks[which(AntTasks$AntTask=="queen"),"AntTask_num"] <- 0
  AntTasks <- AntTasks[order(AntTasks$antID),]
  AntTasks <- merge(AntTasks,AntID_list,all.x=T) 
  
  rm(list=ls()[which(!ls()%in%c("positions_summaries_list","positions_SUMS","e","foraging_zone","AntID_list","AntTasks"))]) #close experiment
  gc()
  mallinfo::malloc.trim(0L)
  
  return(AntTasks)  ##RETURN OUTPUT
  
} 




