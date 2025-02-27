################## GET ANT TASK ###################################
#exp<-e
AntTasks1 <- function(exp){
  print("Computing AntTasks based on 48h time-window before exposure")
  
  #required packages
  require(FortMyrmidon)
  require(mallinfo)
  
  #PARAMS for function
  #hour_chunk_start represents the start time (end time defined by TimeWindow)
  hour_chunk_start <- sort(seq(39, 51, by = 12),decreasing = TRUE ) # c(75,63,51,39)
  TimeWindow <- 12
  
  
  #define zones
  #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
  zones <- exp$spaces[[1]]$zones #function to show the Zones present in the Space
  zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
  foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))
  
  #Get complete list of Ants
  AntID_list <- NULL
  for (ant in   1: length(exp$ants)) {
    AntID_list <- c(AntID_list,exp$ants[[ant]]$ID)}
  
  positions_summaries_list <- list()
  loop_N <- 0
  
  for (HOUR_start in hour_chunk_start) {
    loop_N <- loop_N + 1
    print(paste0("Computing chunk ",loop_N," of 4"))
    
    ## get 2 12Hours window for the Task calculation
    ## calcualte the task BEFORE the EXPOSURE
    start <- fmQueryGetDataInformations(exp)$end - HOUR_start*3600
    #start <- fmQueryGetDataInformations(exp)$start + 33*3600 ####first time in tracking plus 21 hours, to skip acclimation time + 12 HOURS
    time_start <- fmTimeCreate(offset=start)
    #time_start <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end - 24*3600)####last time in tracking minus 24 hours
    stop <- fmQueryGetDataInformations(exp)$end - (HOUR_start-TimeWindow)*3600
    #stop  <- fmQueryGetDataInformations(exp)$start + 45*3600 ####pre-tracking period 
    time_stop   <- fmTimeCreate(offset=stop)
    #time_stop  <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end) ####last time in tracking
    ###QUERY 3: fmQueryComputeAntTrajectories()
    positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
    positions_summaries       <- positions$trajectories_summary
    positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
    positions_list            <- positions$trajectories
    
    #for each trajectory, check whether the ant was observed in the foraging zone (positions_list$zone=2) or not. 
    # if so the ant is a forager, if not the ant is a nurse
    positions_summaries$AntTask <- NA
    
    ##before going back and forth between positions_summaries and positions_list:
    ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
    nrow(positions_summaries)==length(positions_list)
    ####2/ always make sure that positions_summaries is ordered correctly, using the index column
    positions_summaries <- positions_summaries[order(positions_summaries$index),]
    ###this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
    for ( ant_index in 1:length(positions_list)) {
      positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]%in%foraging_zone))
      positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]%in%nest_zone))
      
      if ( foraging_zone %in% positions_list[[ant_index]][,"zone"]){
        positions_summaries[ant_index,"AntTask"] <- "forager"
      }else{
        positions_summaries[ant_index,"AntTask"] <- "nurse"
      }
    }
    #match antID and tagID (live tracking gives tagID). 
    IDs <- exp$identificationsAt(fmTimeNow()) #this skips dead ants
    IDs[sapply(IDs, is.null)] <- NA # assign NA to dead ants
    IDs <- data.frame(tag_hex_ID=unlist(IDs), antID=1:length(IDs),stringsAsFactors = F)
    positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
    #positions_summaries1 <- positions_summaries
    #positions_summaries_LOOP <- aggregate(cbind(positions_summaries$nb_frames_outside,positions_summaries$nb_frames_inside), by = list(positions_summaries$antID), FUN = sum);  colnames(positions_summaries_LOOP) [match("Group.1",colnames(positions_summaries_LOOP))] <- "antID"
    positions_summaries_LOOP <- aggregate(cbind(nb_frames_outside,nb_frames_inside) ~ antID + tag_hex_ID, FUN = sum, na.rm=T,na.action=na.pass,positions_summaries )
    #add names that help in merging later on
    colnames(positions_summaries_LOOP) [match("nb_frames_outside",colnames(positions_summaries_LOOP))] <- paste0("nb_frames_outside",loop_N)
    colnames(positions_summaries_LOOP) [match("nb_frames_inside",colnames(positions_summaries_LOOP))] <- paste0("nb_frames_inside",loop_N)
    
    
    #positions_summaries_list <- append(positions_summaries_list, positions_summaries_LOOP)
    positions_summaries_list[[loop_N]] <-  positions_summaries_LOOP
    
    
    rm(list=ls()[which(!ls()%in%c("positions_summaries_list","exp","foraging_zone","nest_zone","AntID_list","hour_chunk_start","TimeWindow","loop_N"))]) #close experiment
    gc()
    mallinfo::malloc.trim(0L)
    
  }
  
  #merge all data frames together
  
  # # non-good looking recursive merging
  # positions_summaries_mergA <- merge(positions_summaries_list[[1]][c("antID","nb_frames_outside1","nb_frames_inside1")], # , "tag_hex_ID"
  #                                    positions_summaries_list[[2]][c("antID","nb_frames_outside2","nb_frames_inside2")], # , "tag_hex_ID"
  #                                    all.x=T,all.y=T)
  # 
  # positions_summaries_mergB <-  merge(positions_summaries_list[[3]][c("antID","nb_frames_outside3","nb_frames_inside3")], # , "tag_hex_ID"
  #                                     positions_summaries_list[[4]][c("antID","nb_frames_outside4","nb_frames_inside4")], # , "tag_hex_ID"
  #                                     all.x=T,all.y=T)
  
  positions_summaries <- Reduce(function(x, y) merge(x, y, all=TRUE), positions_summaries_list)
  
  positions_summaries <- as.data.frame(sapply(positions_summaries,as.numeric))
  
  positions_summaries[is.na(positions_summaries)] <- 0
  
  #positions_summaries$prop_time_outside <- (positions_summaries$nb_frames_outside1+positions_summaries$nb_frames_outside2)/(positions_summaries$nb_frames_outside1+positions_summaries$nb_frames_outside2+positions_summaries$nb_frames_inside1+positions_summaries$nb_frames_inside2)
  
  #sum inside & outside
  positions_SUMS<- data.frame(antID=positions_summaries$antID,tag_hex_ID=positions_summaries$tag_hex_ID, outside=rowSums(positions_summaries[, grep("outside", colnames(positions_summaries))]),
                              inside=rowSums(positions_summaries[, grep("inside", colnames(positions_summaries))]))
  
  positions_SUMS$prop_time_outside <- positions_SUMS$outside/(positions_SUMS$outside+positions_SUMS$inside)
  
  positions_SUMS[which(positions_SUMS$prop_time_outside<=0.01),"AntTask"] <- "nurse"
  positions_SUMS[which(positions_SUMS$prop_time_outside>0.01),"AntTask"] <- "forager"
  
  AntTasks <- data.frame(antID=positions_SUMS[,"antID"],tag_hex_ID=positions_SUMS[,"tag_hex_ID"],AntTask= positions_SUMS[,"AntTask"])
  
  print("AntTasks computed")
  
  # #add missing ants
  missing_ants <- subset(AntID_list, !(AntID_list %in% AntTasks$antID))
  missing_ants_table <- data.frame()
  
  for (MISSING in missing_ants) {
    print(MISSING)
    for (id in length(exp$ants[[MISSING]]$identifications)) {
      #print ( exp$ants[[MISSING]]$identifications[[id]]$tagValue )
      try(
      missing_ants_table <- rbind(missing_ants_table, data.frame(antID=MISSING, tag_hex_ID= exp$ants[[MISSING]]$identifications[[id]]$tagValue ,AntTask=NA))
      )
    }}
  
  AntTasks <- rbind(AntTasks, missing_ants_table)
  # add misssing ants  as NURSE by DEFAULT
  # AntTasks[which(is.na(AntTasks$AntTask)),"AntTask"] <- "nurse"
  
  AntTasks$AntTask_num <- NA
  AntTasks[which(AntTasks$AntTask=="nurse"),"AntTask_num"] <- 1
  AntTasks[which(AntTasks$AntTask=="forager"),"AntTask_num"] <- 2
  AntTasks <- AntTasks[order(AntTasks$antID),]
  
  rm(list=ls()[which(!ls()%in%c("positions_summaries_list","positions_SUMS","exp","foraging_zone","AntID_list","AntTasks"))]) #close experiment
  gc()
  mallinfo::malloc.trim(0L)
  
  ##RETURN OUTPUT
  # warning("Ants that died before the considered time window (pre treatment) will not be assigned a Task by the function. Currently, no task will default to Nurse")
  return(AntTasks)
}

