##################################################################
####### TABLES TO MATCH STROEYMEYT ET AL., 2018, SCIENCE #########
##################################################################
rm(list = ls())
### tables from metadata information
library(reshape2)
library(dplyr)
library(FortMyrmidon)
library(circular)
library(stringr)
library(dplyr)
FRAME_RATE <- 8
  
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


USER<-'FOCACCIA'
if (USER == "2A13_Office") {
  usr <- "ll16598"
} else {
  usr <- "cf19810"
}

SAVEDIR <- paste("/media/",usr,"/SeagateDesktopDrive/SICC_DATA/main_experiment/original_data",sep="")
WORKDIR <- paste("/media",usr,"SeagateDesktopDrive",sep="/")
DATADIR <- paste(WORKDIR, "SICC_DATA", sep = "/")
EXP_DATA<-DATADIR
metadata <- read.csv(paste(DATADIR, "/Metadata_AUG23.csv", sep = ""), header = T, sep=',')
metadata_unique<-distinct(metadata)
nrow(metadata_unique)
metadata_present<-metadata_unique
metadata_present$colony_size
n=0
for (i in 1:nrow(metadata_present)){
  if (metadata_present[i, "queen"]==TRUE){
    print('yes')
    n=n+1
  }
}
n
if (USER=="linda") {
  usr <- "lsartori"
  WORKDIR <- paste("/media", usr, EXTERNAL, "circadian_rhythm_2022_experiment/scripts/EXP1_base_analysis", sep="/")
  DATADIR <- paste(WORKDIR, "EXP_summary_data", sep="/") # paste(WORKDIR,"EXPERIMENT_DATA",sep="/")
  SAVEDIR <- paste("/media", usr, EXTERNAL, "circadian_rhythm_2022_experiment/Lasius-Bristol_pathogen_experiment/main_experiment/original_data", sep="/")
  EXP_DATA <- paste("/media", usr, EXTERNAL, "circadian_rhythm_2022_experiment/tracking", sep ="/")
  #SCRIPTDIR <- "/media/cf19810/DISK4/EXP1_base_analysis/EXP1_analysis scripts"
}

###### LOAD METADATA

# LS: different column names from Adriano: 
# Adriano -> Linda
# size_status -> time_status
# AntTask -> AntTask1perc

## LS: add pathogen_load ####
pathogen_load <- read.csv("/media/cf19810/SeagateDesktopDrive/SICC_pathogen_data/Adriano_qPCR_pathogen_load_MASTER_REPORT.csv")

#remove dead ants

metadata_present <- metadata_present[which(metadata_present$dead==FALSE),]
nrow(metadata_present)
#remove duplicates (the tag rotated ant duplicates have a different stop-start so they have to be eliminated before dups removal)
metadata_present$identifStart <- NULL
metadata_present$identifEnd <- NULL

# LS: before 4432 obs.
metadata_present <- metadata_present %>% distinct()
# LS: after 4408 obs.
# for ants which are alive but because of the loss of the tag they did not get assigned a Task,
# assign task "nurse" (total of 8 ants). NOT DOING SO causes issues with assortativity_nominal in 13_network_measures.
# Excluding them altogether causes issues with the 11_transmission_simulation as interactions will be in the interactions lists but not in the tag_list
# LS: NULL cases for me, skip
#metadata_present[which(is.na(metadata_present$AntTask)),"AntTask"] <- "nurse"

##### ADD EXTRA COLS TO METADATA
#conform naming to science2018
# colony code
#Adriano:
#metadata_present$colony <- NULL
#metadata_present$REP_NUM           <- substring(metadata_present$REP_treat, 2, nchar(metadata_present$REP_treat))
#metadata_present$N_CHAR            <-  4-nchar(metadata_present$REP_NUM)
#rep function is unhappy if the reps are 0...
#paste("colony", paste(rep(0,N_CHAR),collapse=""), metadata_present$REP_NUM,sep="")
#metadata_present[which(metadata_present$N_CHAR==0),"colony"] <- paste("colony", metadata_present[which(metadata_present$N_CHAR==0),"REP_NUM"],sep="")
#metadata_present[which(metadata_present$N_CHAR==1),"colony"] <- paste("colony", 0, metadata_present[which(metadata_present$N_CHAR!=0),"REP_NUM"],sep="")
head(metadata_present$REP_treat)
# colony code
metadata_present$colony <- sub(".*_", "", metadata_present$REP_treat)
head(metadata_present$colony)

# colony_status
metadata_present$status_char <- gsub("[0-9]", "", metadata_present$colony)
metadata_present$colony_status     <- ifelse(metadata_present$status_char  == "P", "pathogen", ifelse(metadata_present$status_char  == "S", "control", NA))

# treatment code
#period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
# Adriano:
# metadata_present$size_char         <- substr(metadata_present$REP_treat, nchar(metadata_present$REP_treat)-1, nchar(metadata_present$REP_treat)-1)
# metadata_present$size_status       <- ifelse(metadata_present$size_char  == "S", "small", ifelse(metadata_present$size_char  == "B", "big", NA))
# metadata_present$treatment_code    <- paste(metadata_present$colony_status ,metadata_present$size_status,sep=".")
# treatment code
#metadata_present$time_char        <- substring(metadata_present$time_treat, 1, 1)
#metadata_present$time_status      <- ifelse(metadata_present$time_char  == "D", "day", ifelse(metadata_present$time_char  == "N", "night", NA))
desired_length=3

num_zeros <- desired_length - nchar(metadata_present$colony)
# Pad the values with leading zeros to the desired length
metadata_present$treatment_code <- str_pad(metadata_present$colony, width = desired_length, pad = "0")
metadata_present$colony <-metadata_present$treatment_code 


#list subdirectories in parent folder EXPERIMENT_DATA
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}
# list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(EXP_DATA, n = 1)

#select REP folders
files_list <- files_list[grep('V',files_list)]
files_list <- files_list[!grepl('NW',files_list)]

#### OPEN REPLICATE
EXP_list <- NULL
# replicate folder
for (REP.n in 1:length(files_list)) {
  REP.folder      <- files_list[REP.n]
  #REP.files       <- list.files(REP.folder, pattern = ".myrmidon")
  #REP.files       <- REP.files[!grepl('NW_',REP.files)]
  REP.files       <- list.files(REP.folder, pattern = "STROEYM_QUEEN.myr")
  # REP.files       <- REP.files[!grepl('manual_orient',REP.files)]#remove man orient files
  
  for (variable in REP.files) {
    
    REP.filefolder  <- paste(REP.folder,variable,sep="/")
    tracking_data <- fmExperimentOpen(REP.filefolder)
    space <- tracking_data$spaces #save these to a list and then get the extrapolation to only run for files in this list
    tracking_data_file <- list.files(REP.folder, pattern =".000") #may not work for all1 or others
    dir_file <- dirname(REP.filefolder)
    for (s in space){
      space_name <- s$name
    }
    
    EXP_list <- rbind(EXP_list ,data.frame(space_name,variable,path_name=paste(REP.folder,variable,sep="/")))
    
  }
}##DONE
#select REP folders
#files_list <- files_list[grep("REP",files_list)]
# replicate folder
metadata_present$box<-NULL
for (REP.n in 1:length(files_list)) {
  #REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files <- list.files(REP.folder, pattern = "STROEYM_QUEEN.myr") #AntsCreated_AutoOriented_withMetaData_NS_NS_q_
  ##
  ### select only for the science2018 capsules, the name to grep for is in NOTION
  ##
  #REP.files <- REP.files[!grepl("CapsuleDef", REP.files)]
  REP.filefolder <- paste(REP.folder, REP.files, sep = "/")
  
  
  #replicate file
  #for (REP.FILES in REP.filefolder) {
  #get substring in variable until R9BS_
  # REP.FILES <- REP.filefolder[1]
  REP.FILES <- REP.filefolder  
  input_vector <- strsplit(basename(REP.FILES), "_")[[1]]
  box=input_vector[5]
  # Extract the first 2 elements of the vector
  # Use grep to find the element containing "R"
  REP_TREAT <- paste(input_vector[1:2], collapse = "_")
  REP_treat<-REP_TREAT
  #REP_treat <- sub("\\_.*", "", basename(REP.FILES))
  # base file info
  
  
  #REP_name <- unlist(strsplit(REP.FILES,split="/"))[length(unlist(strsplit(REP.FILES,split="/")))]
  #REP_treatment <- unlist(strsplit(REP_name,split="_"))[1]
  REP_ts <- input_vector[5]
  
  ##LL
  split_REP_TREAT <- strsplit(REP_TREAT, "_")[[1]]
  REP_NUM           <- split_REP_TREAT[2]
  status_char <- gsub("[^A-Za-z]", "", REP_NUM)
  REP_NUM<- gsub("\\D", "", REP_NUM)
  
  #colony            <- paste("colony", paste(rep(0,4-nchar(REP_NUM)),collapse=""), REP_NUM,sep="")
  # colony_status
  
  colony<-paste(REP_NUM, status_char, sep="")
  colony_status     <- ifelse(status_char == "P", "pathogen", ifelse(status_char == "S", "control", NA))
  # treatment code
  #period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
  desired_length=3
  num_zeros <- desired_length - nchar(colony)
  # Pad with zeros
  treatment_code <- paste0(str_pad("", width = num_zeros, pad = "0"), colony)
  for (i in 1:nrow(metadata_present)){
    if (metadata_present[i, "colony"]==treatment_code){
      print('yes')
      metadata_present[i, "box"]<-box
      metadata_present[i, "treatment"]<-colony_status
      
      metadata_present[i, "treatment_code"]<-paste0(str_pad("", width = 1, pad = "0"), status_char)
      print(metadata_present[i, "treatment_code"])
    }
  }
  
}
metadata_present[3000,]
##################################################################
#######################    info.txt    ###########################
##################################################################
metadata_colony <- metadata_present[,c("colony","treatment_code","colony_size","box","REP_treat")]
#metadata_colony <- metadata_present[,c("colony","treatment_code","REP_treat")]

metadata_colony <- unique(metadata_colony)

info <- data.frame(colony = metadata_colony$colony,
                   treatment = metadata_colony$treatment_code,
                   box = metadata_colony$box,
                   colony_size = metadata_colony$colony_size,
                   Ynestmin = NA,
                   Ynestmax = NA,
                   nb_larvae = NA,
                   nb_pupae = NA,
                   colony_age = NA,
                   REP_treat= metadata_colony$REP_treat)

write.table(info, file = file.path(SAVEDIR,"info.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

##### ADD EXTRA COLS TO PATHOHEN LOAD DATA
# LS: add pathogen to data ####
# colony code
pathogen_load$REP_NUM           <- substring(pathogen_load$Colony, 2, nchar(pathogen_load$Colony))
pathogen_load$N_CHAR            <-  4-nchar(pathogen_load$REP_NUM)
#rep function is unhappy if the reps are 0...
#paste("colony", paste(rep(0,N_CHAR),collapse=""), pathogen_load$REP_NUM,sep="")
pathogen_load[which(pathogen_load$N_CHAR==0),"colony"] <- paste("colony", pathogen_load[which(pathogen_load$N_CHAR==0),"REP_NUM"],sep="")
pathogen_load[which(pathogen_load$N_CHAR==1),"colony"] <- paste("colony", 0, pathogen_load[which(pathogen_load$N_CHAR!=0),"REP_NUM"],sep="")
# colony_status
pathogen_load$status_char       <- substr(pathogen_load$Colony, nchar(pathogen_load$Colony), nchar(pathogen_load$Colony))
pathogen_load$colony_status     <- ifelse(pathogen_load$status_char  == "P", "pathogen", ifelse(pathogen_load$status_char  == "S", "control", NA))
# treatment code
#period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
pathogen_load$size_char         <- substr(pathogen_load$Colony, nchar(pathogen_load$Colony)-1, nchar(pathogen_load$Colony)-1)
pathogen_load$size_status       <- ifelse(pathogen_load$size_char  == "S", "small", ifelse(pathogen_load$size_char  == "B", "big", NA))
pathogen_load$treatment_code    <- paste(pathogen_load$colony_status ,pathogen_load$size_status,sep=".")

##################################################################
##################     task_groups.txt    ########################
##################################################################

task_groups <- data.frame(colony = metadata_present$colony,
                          tag =  metadata_present$antID,
                          task_group = metadata_present$AntTask,
                          treatment = metadata_present$treatment_code,
                          REP_treat= metadata_present$REP_treat)
write.table(task_groups, file = file.path(SAVEDIR,"task_groups.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

##################################################################
##############     treated_worker_list.txt    ####################
##################################################################

# LL add table for BB and KVL treated
metadata_KVL_Y_treated <- metadata_present[which(metadata_present$CHALL_KVL == TRUE | metadata_present$CHALL_YELLOW == TRUE),]

treated_worker_list_KVLY <- data.frame(colony = metadata_KVL_Y_treated$colony,
                                  tag =  metadata_KVL_Y_treated$antID,
                                  survived_treatment = T,
                                  treatment = metadata_KVL_Y_treated$treatment_code,
                                  REP_treat= metadata_KVL_Y_treated$REP_treat,
                                  AntTask= metadata_KVL_Y_treated$AntTask)

write.table(treated_worker_list_KVLY, file = file.path(SAVEDIR,"treated_worker_list_KVLY.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
##################################################################
##############           seed files           ####################
##################################################################

# flag to determine if selecting 10% of the colony size (T) or the same N as the real exposed ants (F)
select_fixed_N <- F

# same as above : treated_workers.txt
# treated workers simulation is to compare effects with the observed post-exposure changes
write.table(treated_worker_list_KVLY, file = file.path(SAVEDIR,"seeds/treated_worker_list_KVLY.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

# LL now for Beauveria/blue

metadata_BBB_treated <- metadata_present[which(metadata_present$CHALL_BB == TRUE | metadata_present$CHALL_BLUE == TRUE),]

treated_worker_list_BBB <- data.frame(colony = metadata_BBB_treated$colony,
                                       tag =  metadata_BBB_treated$antID,
                                       survived_treatment = T,
                                       treatment = metadata_BBB_treated$treatment_code,
                                       REP_treat= metadata_BBB_treated$REP_treat,
                                       AntTask= metadata_BBB_treated$AntTask)

write.table(treated_worker_list_BBB, file = file.path(SAVEDIR,"treated_worker_list_BBB.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")



#SI LIST
# LL add table for BB and KVL treated
metadata_SIBB <- metadata_present[which(metadata_present$SIBB == TRUE),]

worker_list_SIBB <- data.frame(colony = metadata_SIBB$colony,
                                       tag =  metadata_SIBB$antID,
                                       survived_treatment = T,
                                       treatment = metadata_SIBB$treatment_code,
                                       REP_treat= metadata_SIBB$REP_treat,
                                       AntTask= metadata_SIBB$AntTask)

write.table(worker_list_SIBB, file = file.path(SAVEDIR,"worker_list_SIBB.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
metadata_SIKVL <- metadata_present[which(metadata_present$SIKVL == TRUE),]

worker_list_SIKVL <- data.frame(colony = metadata_SIKVL$colony,
                               tag =  metadata_SIKVL$antID,
                               survived_treatment = T,
                               treatment = metadata_SIKVL$treatment_code,
                               REP_treat= metadata_SIKVL$REP_treat,
                               AntTask= metadata_SIKVL$AntTask)

write.table(worker_list_SIKVL, file = file.path(SAVEDIR,"worker_list_SIKVL.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
metadata_SISH <- metadata_present[which(metadata_present$SISH == TRUE),]

worker_list_SISH <- data.frame(colony = metadata_SISH$colony,
                               tag =  metadata_SISH$antID,
                               survived_treatment = T,
                               treatment = metadata_SISH$treatment_code,
                               REP_treat= metadata_SISH$REP_treat,
                               AntTask= metadata_SISH$AntTask)

write.table(worker_list_SISH, file = file.path(SAVEDIR,"worker_list_SISH.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
# simulation of transmission
# simulations on pre- and post-networks 
#
# Pre networks:
# - compare how transmission is compared to observed post (seed: treated workers). simulation shows what is the expected in the post-networks and it can be compared with the actual post-net qPCR data
# - what would have happened if the exposed where the foragers/nurses/random_workers. This simulation is to study the effects in the pre-period (constitutive properties)


# For simulations, ants where randomly selected from the pool of nurses, foragers or all of the workers, matching the number of exposed nurses in each colony.

# select the N of EXPOSED ants that where alive (present in metadata_present)
# Calculate N_exposed_alive as the number of Exposed per colony = TRUE
metadata_present$Exposed_KVLY<-metadata_present$CHALL_KVL+metadata_present$CHALL_YELLOW
metadata_present$Exposed_BBB<-metadata_present$CHALL_BB+metadata_present$CHALL_BLUE

metadata_present$N_exposed_alive_KVLY <- ave(metadata_present$Exposed_KVLY, metadata_present$colony, FUN = function(x) sum(x))
metadata_present$N_exposed_alive_BBB <- ave(metadata_present$Exposed_BBB, metadata_present$colony, FUN = function(x) sum(x))

# Randomly select X rows per colony
#selected_rows <- dataset[ave(seq_len(nrow(dataset)), dataset$colony, FUN = function(x) sample(x, size = unique(x$X))), ]
metadata_present$queen

for (GROUP in c("nurse","forager","random_worker")) {
  # if (!(GROUP %in% c("nurse","forager"))) {
    metadata_GROUP <- metadata_present[which(metadata_present$queen!=TRUE),]
    metadata_GROUP <- metadata_present[which(metadata_present$SIBB!=TRUE),]
    metadata_GROUP <- metadata_present[which(metadata_present$SIKVL!=TRUE),]
    metadata_GROUP <- metadata_present[which(metadata_present$SISH!=TRUE),]
    metadata_GROUP <- metadata_present[which(metadata_present$CHALL_BB!=TRUE),]
    metadata_GROUP <- metadata_present[which(metadata_present$CHALL_KVL!=TRUE),]
    metadata_GROUP <- metadata_present[which(metadata_present$CHALL_BLUE!=TRUE),]
    metadata_GROUP <- metadata_present[which(metadata_present$CHALL_YELLOW!=TRUE),]
  # }else{
  #   metadata_GROUP <- metadata_present[which(metadata_present$AntTask1perc==GROUP),]
  # }
  if (select_fixed_N) {
    # nurses.txt
    # Group the dataset by "colony" and sample rows based on "status"
    metadata_GROUPs_seed <- metadata_GROUP %>%
      group_by(colony) %>%
      sample_n(ifelse(status == "Big", 18, 3)) %>% # replace FALSE means do not resample rows if threre are less than N to be sliced
      ungroup()
  }else{
    metadata_GROUPs_seed_KVLY <- metadata_GROUP %>%
      group_by(colony) %>%
      sample_n(first(N_exposed_alive_KVLY)) %>%
      ungroup()
    metadata_GROUP_remaining <- anti_join(metadata_GROUP, metadata_GROUPs_seed_KVLY)
    
    # Sample rows for BBB from the remaining rows
    metadata_GROUPs_seed_BBB <- metadata_GROUP_remaining %>%
      group_by(colony) %>%
      sample_n(first(N_exposed_alive_BBB)) %>%
      ungroup()
      }

  metadata_GROUPs_seed_KVLY <- metadata_GROUPs_seed_KVLY[,c("colony","antID","treatment_code","REP_treat","AntTask")]
  names(metadata_GROUPs_seed_KVLY)[which(names(metadata_GROUPs_seed_KVLY)=="treatment_code")] <- "treatment"
  names(metadata_GROUPs_seed_KVLY)[which(names(metadata_GROUPs_seed_KVLY)=="antID")] <- "tag"

  write.table(metadata_GROUPs_seed_KVLY, file = file.path(SAVEDIR,paste0("seeds/",GROUP,"s.txt")), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
  
  metadata_GROUPs_seed_BBB <- metadata_GROUPs_seed_BBB[,c("colony","antID","treatment_code","REP_treat","AntTask")]
  names(metadata_GROUPs_seed_BBB)[which(names(metadata_GROUPs_seed_BBB)=="treatment_code")] <- "treatment"
  names(metadata_GROUPs_seed_BBB)[which(names(metadata_GROUPs_seed_BBB)=="antID")] <- "tag"
  
  write.table(metadata_GROUPs_seed_BBB, file = file.path(SAVEDIR,paste0("seeds/",GROUP,"s.txt")), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
}



##################################################################
#####################    qPCR_file.txt    ########################
##################################################################

head(pathogen_load[,c("colony","treatment_code","antID","Exposed","IsAlive","MbruDNA","above_detection_threshold", "Colony")])

qPCR_file <- data.frame(colony = pathogen_load$colony,
                        treatment = pathogen_load$treatment_code,
                        tag =  pathogen_load$antID,
                        age = NA,
                        status = pathogen_load$Exposed,
                        alive_at_sampling_time = pathogen_load$IsAlive,
                        measured_load_ng_per_uL = pathogen_load$MbruDNA,
                        above_detection_threshold = pathogen_load$above_detection_threshold,
                        REP_treat= pathogen_load$Colony)


write.table(qPCR_file, file = file.path(SAVEDIR,"qPCR","qPCR_file.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

##################################################################
#####################    tag_files.txt    ########################
##################################################################

### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!


for (REP.n in 1:length(files_list)) {
  #REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files <- list.files(REP.folder, pattern = "STROEYM_QUEEN.myr") #AntsCreated_AutoOriented_withMetaData_NS_NS_q_
  ##
  ### select only for the science2018 capsules, the name to grep for is in NOTION
  ##
  #REP.files <- REP.files[!grepl("CapsuleDef", REP.files)]
  REP.filefolder <- paste(REP.folder, REP.files, sep = "/")
  

  #replicate file
  #for (REP.FILES in REP.filefolder) {
    #get substring in variable until R9BS_
    # REP.FILES <- REP.filefolder[1]
  REP.FILES <- REP.filefolder  
  input_vector <- strsplit(basename(REP.FILES), "_")[[1]]
  box=input_vector[5]
  # Extract the first 2 elements of the vector
  # Use grep to find the element containing "R"
  REP_TREAT <- paste(input_vector[1:2], collapse = "_")
  REP_treat<-REP_TREAT
  #REP_treat <- sub("\\_.*", "", basename(REP.FILES))
  # base file info
  
  
  #REP_name <- unlist(strsplit(REP.FILES,split="/"))[length(unlist(strsplit(REP.FILES,split="/")))]
  #REP_treatment <- unlist(strsplit(REP_name,split="_"))[1]
  REP_ts <- input_vector[5]
  
  ##LL
  split_REP_TREAT <- strsplit(REP_TREAT, "_")[[1]]
  REP_NUM           <- split_REP_TREAT[2]
  REP_NUM<- gsub("\\D", "", REP_NUM)
  
  #colony            <- paste("colony", paste(rep(0,4-nchar(REP_NUM)),collapse=""), REP_NUM,sep="")
  # colony_status
  status_char       <- substr(REP_TREAT, nchar(REP_TREAT), nchar(REP_TREAT))
  colony<-paste(REP_NUM, status_char, sep="")
  colony_status     <- ifelse(status_char == "P", "pathogen", ifelse(status_char == "S", "control", NA))
  # treatment code
  #period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
  desired_length=3
  num_zeros <- desired_length - nchar(colony)
  # Pad with zeros
  treatment_code <- paste0(str_pad("", width = num_zeros, pad = "0"), colony)
  for (i in 1:nrow(metadata_present)){
    if (metadata_present[i, "colony"]==colony){
      print('yes')
      metadata_present[i, "box"]<-box
      print(metadata_present[i, "box"])
    }
  }
  
    #conform naming to science2018
    # Adriano:
    # colony code
    # REP_NUM           <- substring(REP_treat, 2, nchar(REP_treat))
    # colony            <- paste("colony", paste(rep(0,4-nchar(REP_NUM)),collapse=""), REP_NUM,sep="")
    # colony_status
    # status_char       <- substr(REP_treat, nchar(REP_treat), nchar(REP_treat))
    # colony_status     <- ifelse(status_char == "P", "pathogen", ifelse(status_char == "S", "control", NA))
    # treatment code
    # size_char         <- substr(REP_treat, nchar(REP_treat)-1, nchar(REP_treat)-1)
    # size_status       <- ifelse(size_char == "S", "small", ifelse(size_char == "B", "big", NA))
    # treatment_code    <- paste(colony_status,size_status,sep=".")
    
    print(paste(REP_treat,colony,treatment_code,sep=" | ")) ##}}
    #open experiment
    e <- fmExperimentOpen(REP.FILES)
    exp.Ants  <- e$ants
    # L: might need to change the times here because my experiment did end before I ended the tracking but variables are not used anywhere else ####
    exp_end   <- fmQueryGetDataInformations(e)$end
    exp_start <- fmQueryGetDataInformations(e)$start

    tag_stats <- fmQueryComputeTagStatistics(e)
    
    ############# CREATE BASE FILE
    tag_file <- NULL
    
    for (ant in exp.Ants){
      #ant <- exp.Ants[[1]] #temp
      # # FIX
    if (is.null(ant)==T){
      print('null ant')
      next
    }

      if(all(ant$getValues("dead")["values"]$values)==FALSE){ # check if there is any false
      #warning("do check that ants are assigned correctly, final_status will fail as the ant is present 2 times in the metadata")
      for (id in ant$identifications){
        group = metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$REP_treat==REP_TREAT),"AntTask"]
        if (length(group)==0){
          print(id$tagValue )
          next
        }
        fs = ifelse(metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$REP_treat==REP_TREAT),"dead"]==F,"alive","dead")
        # id <- ant$identifications[[1]] #temp
        #if the ant died, skip
        ### THIS EXCLUDES THE ANTS WITH ROTATED TAGS! if (capture.output(ant$identifications[[1]]$end)=="+∞") { 
        tag_file <- rbind(tag_file,data.frame(tag =  ant$ID,
                               count = tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"count"], #count considers also the acclimation time, therefore it is always higher than the N of frames
                               # L: this is also different for me ####
                               # AW: last_det = 51*60*60 * 8, #fixed= sll alive, so it is the n of frames for 51 hours.   tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"lastSeen"], 
                               #last_det = round(as.numeric(difftime(tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"lastSeen"],exp_start,units="secs"))*FRAME_RATE), #LS: change made by Nathalie to have a general definition
                               rot = round(deg(id$antAngle),3), # assumed to be the the relative angle of the ant to the tag
                               displacement_distance = 0, #id$antPosition one of the two measures
                               displacement_angle = 0,
                               antenna_reach = NA,
                               trapezoid_length = NA,
                               type = "N", #what does it mean?
                               size = 0 ,
                               headwidth = 0,
                               death = 0,  #all alive, dead not included
                               age = 0,
                               final_status = ifelse(metadata_present[which(metadata_present$tagIDdecimal == id$tagValue& metadata_present$REP_treat==REP_TREAT),"dead"]==F,"alive","dead"),
                               group = metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$REP_treat==REP_TREAT),"AntTask"],
                               REP_treat = REP_TREAT,
                               stringsAsFactors = F))
        tag_file$rot <- round(tag_file$rot,3)
        # group = metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$REP_treat==REP_TREAT),"AntTask"]
        # metadata_present$REP_treat
        #identifEnd <- ifelse(capture.output(print(id$end))=="+∞",NA,ifelse(capture.output(print(id$end))))
        #fmQueryGetDataInformations(e)$end - 51*60*60
        
  
          # individual  <- ant$ID
          # #extract metadata info per key
          #   #for more than one row, always last row will be picked (the relevant one if there is a timed change or the default one if there is not)
          #   for(ROW in 1:2) {
          #     #assign metadata_key value when ID corresponds
          #     tag_file[which(tag_file$tag==ant$ID),"final_status"] <- ant$getValues("IsAlive")["IsAlive","values"]
          #     #if the ant died
          #     if (METADATA_KEY=="IsAlive") {
          #       if (ant$getValues("IsAlive")[ROW,"values"]==FALSE) {
          #         metadata[which(metadata$antID==ant$ID),"surviv_time"] <- ant$getValues("IsAlive")[ROW,"times"]
          #         # if didn't die    
          #       }else if (ant$getValues("IsAlive")[ROW,"values"]==TRUE) {
          #         metadata[which(metadata$antID==ant$ID),"surviv_time"] <- as.POSIXct( exp_end,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
          #       }} # IsAlive check
          #   } # ROW
        #}
        }
    }
    }#ant is dead
    # #check if the tag_file file exists
    # if(file.exists(file.path(SAVEDIR,"tag_files",paste(colony,treatment_code,".txt")))){
    #     print(paste0(REP_treat," already present in tag_file, skip"))
    #   } else {
    
    
    #remove dead (none should be there anyway)
   # tag_file <- tag_file[which(tag_file$final_status=="alive"),]
    # remove duplicated tag ant rows and keep single one (first row)
    tag_file <- tag_file[!duplicated(tag_file$tag),]
    # for ants which are alive but because of the loss of the tag they did not get assigned a Task,
    # assign task "nurse" (total of 8 ants). NOT DOING SO causes issues with assortativity_nominal in 13_network_measures.
    # Excluding them altogether causes issues with the 11_transmission_simulation as interactions will be in the interactions lists but not in the tag_list
    # tag_file[which(is.na(tag_file$group)),"group"] <- "nurse"

        write.table(tag_file, file = file.path(SAVEDIR,"tag_files",paste0(colony,"_",treatment_code,".txt")), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
        
        
        # }
  #} # REP folders
} # REP by REP



##################################################################
#################    time_aggregation_info    ####################
##################################################################

# needed by the simulations

# file format: colony020_pathogen_PostTreatment.txt

# time frames time_since_start time_hours date_GMT time_of_day_exact time_of_day
# 1415277843.35 800219 111.14 0 2014-11-06.12:44:03.34 12.73 12
# 1415288643.47 821820 114.14 3 2014-11-06.15:44:03.47 15.73 15
# 1415299443.11 843420 117.14 6 2014-11-06.18:44:03.10 18.73 18
# 1415310243.24 865021 120.14 9 2014-11-06.21:44:03.24 21.73 21
# 1415321043.37 886622 123.14 12 2014-11-07.00:44:03.36 0.73 0
# 1415331843.48 908223 126.14 15 2014-11-07.03:44:03.48 3.73 3
# 1415342643.11 929823 129.14 18 2014-11-07.06:44:03.10 6.73 6
# 1415353443.25 951424 132.14 21 2014-11-07.09:44:03.25 9.73 9

# keep cols: time, time_hours, time_of_day

source_folders <- c(
"/media/cf19810/SeagateDesktopDrive/SICC_DATA/main_experiment/intermediary_analysis_steps/full_interaction_lists/PostTreatment/observed",
"/media/cf19810/SeagateDesktopDrive/SICC_DATA/main_experiment/intermediary_analysis_steps/full_interaction_lists/PreTreatment/observed")

for (SOURCE in source_folders){
  # Get a list of all .txt Interaction files in SOURCE
  interaction_files <- list.files(SOURCE, pattern = "\\.txt$", full.names = TRUE)
  
  for (INT in interaction_files){
    # read each .txt Interaction file in SOURCE
    Interaction <- read.table(INT, header = TRUE, sep = "\t")
    
  # select, per each Interaction$time_hours, the min(Interaction$Starttime) and return a table of Interaction$Starttime, Interaction$time_hours and Interaction$time_of_day
  aggregated_table <- aggregate(Interaction$Starttime, by = list(Interaction$time_hours, Interaction$time_of_day), FUN = min)
  
  # rename Interaction$Starttime to Interaction$time
  colnames(aggregated_table) <- c("time_hours", "time_of_day", "time")
  aggregated_table <- aggregated_table[order(aggregated_table$time),]
  
  # Save the table in /media/cf19810/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/original_data/time_aggregation_info using the Interaction file name but trimming the chars "_interactions"
  # example name:colony09SS_control.small_PreTreatment.txt
  output_filename <- gsub("_interactions", "", basename(INT))
  write.table(aggregated_table, file = paste0("/media/cf19810/SeagateDesktopDrive/SICC_DATA/main_experiment/original_data/time_aggregation_info/", output_filename), sep = "\t", row.names = FALSE)
}
}




















