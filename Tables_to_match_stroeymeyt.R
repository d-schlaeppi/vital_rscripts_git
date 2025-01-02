#rm(list = ls())
rm(list = setdiff(ls(), "first_time_use_working_directory"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### TABLES TO MATCH STROEYMEYT ET AL., 2018, SCIENCE ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### Background information | Read me ####

# Script created by Adriano Wanderlingh, with contributions of  Nathalie Stroeymeyt
# with adaptations by Linda Sartoris and finally adjusted to the needs and data structure of Daniel Schläppi
# Should be compatible with Fort 0.8.1 and 0.8.3

# When running this you will create a couple of tables in the right format with the right names for variables so it should match 
# the data structure of Stroeymeyt 2018 which will help running the remaining analysis pipeline.

# Index
# 1. Prerequisites
# 2. To dos and Notes
# 3. Load data
# 4. Create tables
# 4.1 task_groups.txt 
# 4.2 treated_worker_list.txt
# 4.3 seed files
# 4.4 info.txt
# 4.5 bead_file.txt
# 4.6 tag_files.txt
# 4.7 time_aggregation_info

#### 1. Prerequisites ####

# tables from metadata information

pacman::p_load(tcltk,reshape2,dplyr,FortMyrmidon,circular)
FRAME_RATE <- 6
# setwd(tk_choose.dir(default = "~/", caption = "Select Working Directory"))
if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") { # direct it to where you have config_user_and_hd.R (typically the script folder or github folder)
  selected_dir <- tcltk::tk_choose.dir(default = "~/", caption = "Select Working Directory")
  if (is.null(selected_dir) || selected_dir == "") {
    cat("No directory selected. Exiting.\n")
    return()}
  setwd(selected_dir)
  first_time_use_working_directory <<- getwd()
  cat(crayon::blue(getwd()))
} else { setwd(first_time_use_working_directory)
  cat(crayon::blue(getwd())) }

source("config_user_and_hd.R")
SAVEDIR   <- paste("/media" , usr, hd,"vital/fc2/vital_experiment/main_experiment/original_data",sep="/")


# if (usr != "mac_gismo") {
#   DATADIR   <- paste("/media", usr, hd, "vital/fc2",sep="/") 
#   SCRIPTDIR <- paste("/home" , usr, "Documents/vital_rscripts_git",sep="/")
#   SAVEDIR   <- paste("/media" , usr, hd,"vital/fc2/vital_experiment/main_experiment/original_data",sep="/")
#   } else {
#   DATADIR   <- paste("/Volumes", hd, "vital/fc2", sep="/") 
#   SCRIPTDIR <- paste("/Users", usr, "Documents/GitHub/vital_rscripts_git",sep="/")
#   SAVEDIR   <- paste("/Volumes", hd, "vital/fc2/vital_experiment/main_experiment/original_data",sep="/")
#   }


#### 2. To dos and Notes ####
# AT SOME POINT IT MIGHT BE NECESSARY TO CHANGE TO THE FACET NETWORK TASK DEFINITION VARIABLE INSTEAD OF THE ONE BASED ON PROPORTION OF TIME IN OR OUR THE NEST
# HOWEVER, THIS SCRIPT NEEDS TO BE RUN BEFORE THE FACET NET ONE... So we might want to come back to at at some point. --> Facet Net was just run afterwards and its information was added in addition to the space use one
# note, seed files etc are all calculated based on the old task definition which was based on the percentage of time spent inside or outside the nest. --> maybe this step still needs to be adjusted.


# CHECK WHICH BELOW IS REFERRING TO BEAD LOAD AND ADJUST ADRIANOS PATHOGEN LOAD. 
# pathogen_load <- read.csv("/media/cf19810/DISK4/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data/Adriano_qPCR_pathogen_load_MASTER_REPORT.csv")
# bead_load <- read.csv("path to table containing the information of the beads for each individual: blue, yellow") ### UPDATE AS SOON AS POSSIBLE and also update the part below regarding the bead load. 
# Anything further down in the pipeline referring to the pathogen_load table needs to be run with metadata instead.


# Column names might be slightly different based on experimental design...
# AntTask -> AntTask1perc
# I did not have pathogen load but bead load which is included in the metadata table. 
# As I had only 1 Treatment some variable will be different to Adriano and Linda.
# However, I will need to adjust things once I get to further down in the pipeline to address the issue of bead colors. --> should be sorted in the bead data analysis. 

# Time aggregation info might also be needed for the trophallaxis network... check later... if running simulations along the trophallaxis network. 


#### 3. Load data ####
#  meta data
metadata_present <- read.table(paste(DATADIR, "/individual_metadata_vital.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # will add colony_metadata data frame to the environment so it can be accessed within this script (in my case containing bodylenght information)

#remove dead ants / only keep the ones that are alive
metadata_present <- metadata_present[which(metadata_present$IsAlive==TRUE),]
# ants with rotated tags (retaged with the same tags) appear as duplicates ->
# remove duplicates: however have different start and stop times as well as surv time. As those variables are not used they can be set to zero and then we can eliminate lines that are not distinct. 
nrow(metadata_present)
metadata_present <- metadata_present %>% distinct(colony_id, antID, .keep_all = TRUE) # In my case: before 4242 and after 4235.
nrow(metadata_present)
# for ants which are alive but because of the loss of the tag or another reasong did not get assigned a Task,
# assign task "nurse" (total of 5 ants). NOT DOING SO causes issues with assortativity_nominal in 13_network_measures.
# Excluding them altogether causes issues with the 11_transmission_simulation as interactions will be in the interactions lists but not in the tag_list
metadata_present[which(is.na(metadata_present$AntTask1perc)),"AntTask1perc"] <- "nurse"

# ADD EXTRA COLS TO METADATA and conform naming to Stroeymeyt science 2018 
# colony code
metadata_present$colony <- metadata_present$colony_id

# colony_status + treatment (Adriano and Linda had more treatments so they used different naming... I simplified but kept some names just in case...)
metadata_present$status_char       <- metadata_present$treatment_simple
metadata_present$colony_status     <- metadata_present$treatment_simple

#### 4. Create tables ####


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####              4.1 task_groups.txt                          ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

task_groups <- data.frame(colony = metadata_present$colony_id,
                          tag =  metadata_present$antID,
                          task_group_prop = metadata_present$AntTask1perc,
                          treatment = metadata_present$treatment_simple,
                          REP_treat= paste0(metadata_present$colony_id, "_" ,metadata_present$treatment_simple))
# write.table(task_groups, file = file.path(SAVEDIR,"task_groups.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

# added check to not overwrite the task_group file in case you are coming back from the facet net community detection which adds additional variables to this table.
file_path <- file.path(SAVEDIR, "task_groups.txt")
if (file.exists(file_path)) {
  cat(red("\n !!! WARNING: !!! \n")); cat(yellow("\n The file task_groups.txt already exists! \n"))
  user_response <- readline(prompt = paste(" Do you want to overwrite it? (y/n): "))
  if (tolower(user_response) == "y") {
    write.table(task_groups, file = file_path, append = FALSE, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    cat(green("The file has been overwritten.\n"))
  } else {cat(blue("The file was not overwritten.\n"))
    task_groups <- read.table(paste0(SAVEDIR,"/task_groups.txt"), header = TRUE) }
} else {
  write.table(task_groups, file = file_path, append = FALSE, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  cat(blue("The file has been written.\n"))}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####          4.2 treated_worker_list.txt                      ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

metadata_treated <- metadata_present[which(metadata_present$IsTreated==TRUE),]

treated_worker_list <- data.frame(colony = metadata_treated$colony_id,
                                  tag =  metadata_treated$antID,
                                  survived_treatment = T,
                                  treatment = metadata_treated$treatment_simple,
                                  REP_treat= paste0(metadata_treated$colony_id, "_" ,metadata_treated$treatment_simple),
                                  AntTask= metadata_treated$AntTask1perc)

write.table(treated_worker_list, file = file.path(SAVEDIR,"treated_worker_list.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####              4.3 seed files                               ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


# flag to determine if selecting 10% of the colony size (T) or the same N as the real exposed ants (F)
select_fixed_N <- F
# !!!! that is something to think about: Do it once with the same number of workers that were in the treatment groups and once with the same number of workers that I deem worth including! i.e. that did not drown themselves?

# same as above : treated_workers.txt
# treated workers simulation is to compare effects with the observed post-exposure changes
write.table(treated_worker_list, file = file.path(SAVEDIR,"seeds/treated_workers.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

# simulation of transmission
# simulations on pre- and post-networks
# Pre networks:
# - compare how transmission is compared to observed post (seed: treated workers). simulation shows what is the expected in the post-networks and it can be compared with the actual post-net qPCR data
# - what would have happened if the exposed where the foragers/nurses/random_workers. This simulation is to study the effects in the pre-period (constitutive properties)

# For simulations, ants where randomly selected from the pool of nurses, foragers or all of the workers, matching the number of exposed workers in each colony.
# select the N of EXPOSED ants that where alive (present in metadata_present)
# Calculate N_exposed_alive as the number of Exposed per colony = TRUE
metadata_present$N_exposed_alive <- ave(metadata_present$IsTreated, metadata_present$colony, FUN = function(x) sum(x))

# Randomly select X rows per colony 
# selected_rows <- dataset[ave(seq_len(nrow(dataset)), dataset$colony, FUN = function(x) sample(x, size = unique(x$X))), ]

for (GROUP in c("nurse","forager","random_worker")) {
  if (!(GROUP %in% c("nurse","forager"))) {
    metadata_GROUP <- metadata_present[which(metadata_present$AntTask1perc!="queen"),]
  }else{
    metadata_GROUP <- metadata_present[which(metadata_present$AntTask1perc==GROUP),]
  }
  if (select_fixed_N) { # currently always false so DS goes straight to the else statement
    # nurses.txt
    # Group the dataset by "colony" and sample rows based on "status"
    # metadata_GROUPs_seed <- metadata_GROUP %>%
    #   group_by(colony) %>%
    #   sample_n(ifelse(status == "Big", 18, 3)) %>% # replace FALSE means do not re-sample rows if there are less than N to be sliced
    #   ungroup()
  }else{
    set.seed(69) # for reproducibility of random sampling
    metadata_GROUPs_seed <- metadata_GROUP %>%
      group_by(colony) %>%
      sample_n(first(N_exposed_alive)) %>%
      ungroup()
  }
  metadata_GROUPs_seed <- metadata_GROUPs_seed[,c("colony_id","antID","treatment_simple","AntTask1perc")]
  names(metadata_GROUPs_seed)[which(names(metadata_GROUPs_seed)=="colony_id")] <- "colony"
  names(metadata_GROUPs_seed)[which(names(metadata_GROUPs_seed)=="treatment_simple")] <- "treatment"
  names(metadata_GROUPs_seed)[which(names(metadata_GROUPs_seed)=="antID")] <- "tag"
  metadata_GROUPs_seed <- metadata_GROUPs_seed %>%
    mutate(REP_treat = paste0(colony, "_", treatment))
  
  write.table(metadata_GROUPs_seed, file = file.path(SAVEDIR,paste0("seeds/",GROUP,"s.txt")), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
}



### If the first part of the facet net script has been run to get the task group allocation instead space use we create additional seed files based on this 

if ("task_group_FACETNET_0.5" %in% colnames(task_groups)) { # if that already exists run the following
# rename variables from task_groups to match with metadata_present
  task_groups_sub <- task_groups %>% rename(
      antID = tag,  # tag to antID
      colony_id = colony)   %>%  # colony to colony_id
      dplyr::select(colony_id, antID, task_group_FACETNET_0.5, Forager_score) # relevant variables
# add the facet net variables task_group_FACETNET_0.5 and forager score to metadata_present 
  metadata_present <- metadata_present %>% 
    left_join(task_groups_sub , by = c("colony_id", "antID"))  # colony_id and antID used as identifiers for the join
# rerun the selection of workers for seed files based on facet net task allocation
for (GROUP in c("nurse","forager","random_worker")) { # GROUP <-"nurse"
  if (!(GROUP %in% c("nurse","forager"))) {
    metadata_GROUP <- metadata_present[which(metadata_present$task_group_FACETNET_0.5!="queen"),]
  }else{
    metadata_GROUP <- metadata_present[which(metadata_present$task_group_FACETNET_0.5==GROUP),]
  }
  # if (select_fixed_N) { # currently always false (see above) }else{
  set.seed(69) # for reproducibility of random sampling
  metadata_GROUPs_seed <- metadata_GROUP %>%
      group_by(colony) %>%
      sample_n(first(N_exposed_alive)) %>%
      ungroup() # }
  metadata_GROUPs_seed <- metadata_GROUPs_seed[,c("colony_id","antID","treatment_simple","task_group_FACETNET_0.5")]
  names(metadata_GROUPs_seed)[which(names(metadata_GROUPs_seed)=="colony_id")] <- "colony"
  names(metadata_GROUPs_seed)[which(names(metadata_GROUPs_seed)=="treatment_simple")] <- "treatment"
  names(metadata_GROUPs_seed)[which(names(metadata_GROUPs_seed)=="antID")] <- "tag"
  metadata_GROUPs_seed <- metadata_GROUPs_seed %>%
    mutate(REP_treat = paste0(colony, "_", treatment))
  
  write.table(metadata_GROUPs_seed, file = file.path(SAVEDIR,paste0("seeds/",GROUP,"s_facetnet.txt")), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
}} else {  cat(yellow("Facet net task allocation not defined yet \nTo get this, run the first part of 19_Facetnet_community_detection_DS.R\nOr skip for now and come back later"))}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####                4.4 info.txt                               ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

metadata_colony <- metadata_present[,c("colony_id","treatment_simple","box","colony_size")]
metadata_colony <- unique(metadata_colony)
metadata_colony <- metadata_colony %>%
  mutate(REP_treat = paste0(colony_id, "_", treatment_simple))

info <- data.frame(colony = metadata_colony$colony_id,
                   treatment = metadata_colony$treatment_simple,
                   box = metadata_colony$box,
                   colony_size = metadata_colony$colony_size,
                   Ynestmin = NA,
                   Ynestmax = NA,
                   nb_larvae = NA,
                   nb_pupae = NA,
                   colony_age = NA,
                   REP_treat= metadata_colony$REP_treat)

write.table(info, file = file.path(SAVEDIR,"info.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####                4.5 bead_file.txt                          ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# used to be called qPCR_file.txt and was based on Adrianos pathogen load file as I have the bead data included in the metadata file we just use the metadata file instead: 
# head(pathogen_load[,c("colony","treatment_code","antID","Exposed","IsAlive","MbruDNA","above_detection_threshold", "Colony")])
# qPCR_file <- data.frame(colony = pathogen_load$colony,
#                         treatment = pathogen_load$treatment_code,
#                         tag =  pathogen_load$antID,
#                         age = NA,
#                         status = pathogen_load$Exposed,
#                         alive_at_sampling_time = pathogen_load$IsAlive,
#                         measured_load_ng_per_uL = pathogen_load$MbruDNA,
#                         above_detection_threshold = pathogen_load$above_detection_threshold,
#                         REP_treat= pathogen_load$Colony)
# write.table(qPCR_file, file = file.path(SAVEDIR,"qPCR","qPCR_file.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

bead_file <- metadata_present  %>% dplyr::select(-treatment)
names(bead_file)[which(names(bead_file)=="treatment_simple")] <- "treatment"
names(bead_file)[which(names(bead_file)=="tagIDdecimal")] <- "tag"
bead_file <- bead_file %>% mutate(REP_treat = paste0(colony, "_", treatment))
bead_file <- bead_file %>% mutate(age = NA)
bead_file <- bead_file %>% mutate(status = IsTreated)
#save
bead_file_path <- file.path(SAVEDIR, "beads", "bead_file.txt")
if (!file.exists(dirname(bead_file_path))) {dir.create(dirname(bead_file_path), recursive = TRUE)}
write.table(bead_file, file = bead_file_path, append = F, col.names = T, row.names = F, quote = F, sep = "\t")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####                 4.6 tag_files.txt                         ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

files_list <- list.files(DATADIR)
files_list <- list.files(DATADIR, pattern="CapsuleDef2018.myrmidon")     
files_list <- files_list[which(!grepl("c29",files_list))] # getting a list of each colony
# files_list <- files_list[which(grepl("c11",files_list))]

for (file in files_list) { # file <- files_list[12]
  colony_id <- unlist(strsplit(file,split="_"))[grepl("c",unlist(strsplit(file,split="_")))]
  colony <- colony_id
  # get base info for this colony
  treatment <- colony_metadata[which(colony_metadata$colony_id==colony_id),"treatment_simple"]
  tracking_system <-colony_metadata[which(colony_metadata$colony_id==colony_id),"tracking_system_main"]
  REP_treat <- paste(treatment, tracking_system, colony_id, sep="_") 
  colony_status <- treatment
  print(paste(REP_treat,colony,treatment,sep=" | "))
  e <- fmExperimentOpen(paste(DATADIR,file, sep = "/"))
  exp.Ants  <- e$ants
  exp_end   <- colony_metadata[which(colony_metadata$colony_id==colony_id),"time_treatment_end"]#fmQueryGetDataInformations(e)$end
  treatment_start <-  as.POSIXct(colony_metadata[which(colony_metadata$colony_id==colony_id),"time_treatment_start"], format="%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) #fmQueryGetDataInformations(e)$start
  exp_start <- fmTimeCreate(treatment_start-24*3600)
  rec_start <- as.POSIXct(fmQueryGetDataInformations(e)$start, format="%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )  # the time at which the recording started (variable for each treatment as I sometimes included the acclimatisation time and sometimes did not)
  tag_stats <- fmQueryComputeTagStatistics(e)
  
  # create tag file and populate it...
  tag_file <- NULL
  for (ant in exp.Ants){  # ant <- exp.Ants[[92]] | ant <- exp.Ants[[2]]
    if(all(ant$getValues("IsAlive")["values"]$values)){ #making just making sure there is no false, i.e. ant is alive
      for (id in ant$identifications){ # id <- ant$identifications[[1]]
        if (capture.output(ant$identifications[[1]]$start)=="-∞") {   #DS: Changed from $start)=="+∞" to avoid loosing queen of colony 11 ### THIS EXCLUDES one instance of THE ANTS WITH ROTATED TAGS (or dead ants which should already be excluded anyways) ### skip lines for ants that had the tag rotated which were deleted as duplicates above in metadata present. 
        tag_file <- rbind(tag_file, data.frame(tag =  ant$ID,
                                              count = tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"count"],
                                              # version 1:
                                              # last_det = 27*3600*FRAME_RATE, #not entirely clear: simple like adriano: last_det = 51*60*60 * 8, #fixed= sll alive, so it is the n of frames for 51 hours.tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"lastSeen"],
                                              # version 2 or 3: number of frames since official starting of tracking 24 h before treatment until last seen? But could also be number of frames since recording starts which would be experiment start instead of treatment start -12 
                                              last_det = round(as.numeric(difftime(tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"lastSeen"], treatment_start-24*3600, units="secs"))*FRAME_RATE),
                                              # last_det = round(as.numeric(difftime(tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"lastSeen"], rec_start, units="secs"))*FRAME_RATE),
                                              rot = round(deg(id$antAngle),3), # assumed to be the the relative angle of the ant to the tag
                                              displacement_distance = 0, #id$antPosition one of the two measures
                                              displacement_angle = 0,
                                              antenna_reach = NA,
                                              trapezoid_length = NA,
                                              type = "N", #what does it mean?
                                              size = 0,
                                              headwidth = 0,
                                              death = 0,  #all alive, dead not included
                                              age = 0,
                                              #final_status = ifelse(metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$colony_id==colony_id),"IsAlive"]==T,"alive","dead"),   ### needed fixing as it caused an error for retagged ants with two identifications: removing one of the duplicates in metadata_present means that there not all tag values are still there -> thus we refer to antID instead.
                                              final_status = ifelse(metadata_present[which(metadata_present$antID == id$targetAntID & metadata_present$colony_id==colony_id),"IsAlive"]==T,"alive","dead"),
                                              #group = metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$colony_id==colony_id),"AntTask1perc"],
                                              group =               metadata_present[which(metadata_present$antID == id$targetAntID & metadata_present$colony_id==colony_id),"AntTask1perc"],
                                              REP_treat = REP_treat,
                                              stringsAsFactors = F))
 
        }
      }
    }
   }
  tag_file <- tag_file[which(tag_file$final_status=="alive"),] # just making sure that there are indeed no dead ants
  tag_file <- tag_file[!duplicated(tag_file$tag),]  # remove duplicated tag ant rows and keep single one (first row)
  write.table(tag_file, file = file.path(SAVEDIR,"tag_files",paste0(colony,"_",treatment,".txt")), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####                 4.7 time_aggregation_info                 ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# needed by the simulations
# ev. I need to do the same for the trophallaxis interaction network. However, but probably not as I only had one 3h time block for the pre and the post treatment duration.

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
  paste(DATADIR, "vital_experiment/main_experiment/intermediary_analysis_steps/full_interaction_lists/PostTreatment/observed", sep = "/"),
  paste(DATADIR, "vital_experiment/main_experiment/intermediary_analysis_steps/full_interaction_lists/PreTreatment/observed", sep = "/"))     

for (SOURCE in source_folders){ # SOURCE <-  source_folders[1]
  interaction_files <- list.files(SOURCE, pattern = "\\.txt$", full.names = TRUE) # Get a list of all .txt Interaction files in SOURCE
  for (INT in interaction_files){ # INT  <-  interaction_files[1]
    Interaction <- read.table(INT, header = TRUE, sep = "\t") # read each .txt Interaction file in SOURCE
    
    # select, per each Interaction$time_hours, the min(Interaction$Starttime) and return a table of Interaction$Starttime, Interaction$time_hours and Interaction$time_of_day
    aggregated_table <- aggregate(Interaction$Starttime, by = list(Interaction$time_hours, Interaction$time_of_day), FUN = min)
    
    # rename Interaction$Starttime to Interaction$time
    colnames(aggregated_table) <- c("time_hours", "time_of_day", "time")
    aggregated_table <- aggregated_table[order(aggregated_table$time),]
    
    # save the file 
    output_filename <- gsub("_interactions", "", basename(INT))
    write.table(aggregated_table, file = paste(DATADIR,"vital_experiment/main_experiment/original_data/time_aggregation_info", output_filename, sep = "/"), sep = "\t", row.names = FALSE)
  }
}

###
#' You should now have finished postprocessing, vital_base_analysis and tables_to_match_stroeymeyt meaning you can move on to vital_main_analysis
#' Where you first run the facet_net_community detectino to get a different definition of foragers and nurses.

















