# Metadata extraction file
rm(list=ls())





gc()
####################################################################################
#### THIS SCRIPT CONTAINS:
#### METADATA EXTRACTION FROM MYRMIDON FILES
####################################################################################

#### LIBRARIES
library(FortMyrmidon) ####R bindings
library(stringr)

#### FUNCTIONS
#list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}


#### ACCESS FILES
USER <- "Super"

if (USER=="Luke") {
  WORKDIR <- "/media/ll16598/SeagateDesktopDrive/SICC_DATA" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR)#"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  SCRIPTDIR <- "/home/ll16598/Documents/SICC_ANALYSIS"
}
if (USER=="Super") {
  WORKDIR <- "/media/cf19810/SeagateDesktopDrive/SICC_DATA" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR)#"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  SCRIPTDIR <- "/media/cf19810/SeagateDesktopDrive/SICC_ANALYSIS"
}

###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"ANT_TASKS.R",sep="/"))



###define name of general output table
Metadata_Exp1 <- file.path(DATADIR,"Metadata_Exp1_2023.txt") 


# #inferred$return_time <- NA
# inferred$end_time <- NA
# inferred$exposed <- "no"
# inferred$dead <- "no"
#inferred$time_since_start <- NA

### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
#select REP folders
files_list <- files_list[grep('V',files_list)]
files_list <- files_list[!grepl('NW',files_list)]

# replicate folder
for (REP.n in 1:length(files_list)) {
  
  # REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = "STROEYM_QUEEN.myr")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  #replicate file
  for (REP.FILES in REP.filefolder) {
    #get substring in variable until R9BS_
    # REP.FILES <- REP.filefolder[2]
    #REP_treat <- sub("\\_.*", "", basename(REP.FILES))
    REP_treat <- stringr::str_extract(basename(REP.FILES), "[^_]*_[^_]*" )
    
    #treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))
    
    print(REP.FILES) ##}}
    #open experiment
    e <- fmExperimentOpen(REP.FILES)
    exp.Ants  <- e$ants
    exp_end   <- fmQueryGetDataInformations(e)$end
    exp_start <- fmQueryGetDataInformations(e)$start
   # Metadata_Exp1 <- file.path(paste(REP.folder,"/Metadata_Exp1_2022_",REP_treat,".txt",sep="")) 
    ########### COMPUTE THE ANT TASKS (48h before exposure)
    print(paste0("Computing Ant Tasks and Zone Uses"))
    AntTasks <-AntTasks1(exp=e)
    
    ############# CREATE BASE FILE #########################
    metadata <- NULL
    
    for (ant in exp.Ants){
      for (id in ant$identifications){
        exposure_treat  = substr(REP_treat,(nchar(REP_treat)+1)-2,nchar(REP_treat))
        metadata <- rbind(metadata,data.frame(REP_treat = REP_treat,
                                              #exposure_treat  = substr(REP_treat,(nchar(REP_treat)+1)-2,nchar(REP_treat)),
                                              exposure_treat  = gsub("[^a-zA-Z]", "", exposure_treat),
                                              status           = NA,
                                              antID            = ant$ID,
                                              tagIDdecimal     = id$tagValue,
                                              identifStart     = capture.output(print(id$start)), 
                                              identifEnd       = capture.output(print(id$end)), 
                                              AntTask          = AntTasks[which(AntTasks$antID==ant$ID),"AntTask"],
                                              #AntTask_num          = AntTasks[which(AntTasks$antID==ant$ID),"AntTask_num"],
                                              stringsAsFactors = F))
        
      } }
    
    metadata[grepl("∞", metadata$identifEnd),"identifEnd"] <- NA
    metadata[grepl("∞", metadata$identifStart),"identifStart"] <- NA
    
    #transform time in GMT
    metadata$identifEnd <- as.POSIXct( metadata$identifEnd , format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    
    #colony-wide metadata
    metadata$colony_size <- length(unique(metadata$antID))
    
    metadata$status[metadata$exposure_treat=="S"] <- "Sham"
    metadata$status[metadata$exposure_treat=="P"] <- "Pathogen"

    
    #split in two cols output
    treat_colz <-     data.frame(do.call('rbind',strsplit(metadata$status,' ',fixed=TRUE)))
    metadata$status <- treat_colz$X1
    metadata$treatment <- treat_colz$X2
    
    #empty metadata columns
   # metadata$Comment     <- NA
    metadata$CHALL_BB    <- NA 
    metadata$CHALL_BLUE    <- NA 
    metadata$CHALL_KVL    <- NA 
    metadata$CHALL_YELLOW    <- NA 
    metadata$SIBB    <- NA 
    metadata$SIKVL    <- NA 
    metadata$SISH    <- NA 
    
    metadata$dead     <- NA 
    metadata$queen     <- NA 
   # metadata$surviv_time <- NA # time of survival: time of death or end of exp
    
    ############# METADATA KEYS AND VALUES #########################
    list_keys <- list()
    #assign metadata keys
    for (KEY in   1:length(e$metaDataKeys)) {
      key <- names(e$metaDataKeys[KEY])
      #defaultvalue <- unname(e$metaDataKeys[KEY][[1]])
      list_keys <- c(list_keys,key)
    }
    
    
    ## assign metadata values
    for (ant in exp.Ants){
      individual  <- ant$ID
      if (is.null(individual) == TRUE) {
        next
      }
      print(individual)
      #extract metadata info per key
      for (METADATA_KEY in list_keys) {
        #for more than one row, always last row will be picked (the relevant one if there is a timed change or the default one if there is not)
        for(ROW in 1:nrow(ant$getValues(METADATA_KEY))) {
          #assign metadata_key value when ID corresponds
          metadata[which(metadata$antID==ant$ID),METADATA_KEY] <- ant$getValues(METADATA_KEY)[ROW,"values"]
                    #if the ant died
          #if (METADATA_KEY=="dead") {
            #if (ant$getValues("dead")[ROW,"values"]==TRUE) {
              #metadata[which(metadata$antID==ant$ID),"surviv_time"] <- ant$getValues("IsAlive")[ROW,"times"]
              # if didn't die    
            #}else if (ant$getValues("IsAlive")[ROW,"values"]==TRUE) {
             # metadata[which(metadata$antID==ant$ID),"surviv_time"] <- as.POSIXct( exp_end,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
            #} # IsAlive check
        } # ROW
      } #METADATA_KEY
    }#ant
    #metadata$surviv_time <- as.POSIXct(metadata$surviv_time, origin="1970-01-01", tz="GMT")
    
    # N_exposed in the colony
    metadata$N_CHALL_BB<- sum(metadata$CHALL_BB)
    metadata$N_CHALL_BLUE<- sum(metadata$CHALL_BLUE)
    metadata$N_CHALL_KVL<- sum(metadata$CHALL_KVL)
    metadata$N_CHALL_YELLOW<- sum(metadata$CHALL_YELLOW)
    metadata$N_SIBB<- sum(metadata$SIBB)
    metadata$N_SIKVL<- sum(metadata$SIKVL)
    metadata$N_SISH<- sum(metadata$SISH)

    #start and end times
    metadata$ExpStart  <-  exp_start
    attr(metadata$ExpStart,"tzone") <- "GMT"
    metadata$ExpEnd    <-  exp_end
    attr(metadata$ExpEnd,"tzone") <- "GMT"
    
    
    # save
    if (file.exists(Metadata_Exp1)){
      write.table(metadata,file=Metadata_Exp1,append=T,col.names=F,row.names=F,quote=T,sep=",")
    }else{
      write.table(metadata,file=Metadata_Exp1,append=F,col.names=T,row.names=F,quote=T,sep=",")
    }
  } 

    #gc()
    #remove experiment and all
    #clean up
   #rm(list=ls()[which(!ls()%in%c("REP.folder","REP.files","REP.filefolder","files_list","Metadata_Exp1","SCRIPTDIR","WORKDIR","DATADIR","list.dirs.depth.n"))]) #close experiment
}  
    
  
  
  
 # REP by REP

