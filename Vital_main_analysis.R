# rm(list = ls())
rm(list = setdiff(ls(), "first_time_use_working_directory"))
gc()
Sys.sleep(3)
mallinfo::malloc.trim(0L)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### VITAL MAIN ANALYSIS     ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### 1. Background information | Read me | ToDos ####

# This script contains code to run separate analysis scripts from the stroeymeyt 2018 pipeline: 

# Before starting this scrip, first start with the pre processing of the data following Vital.R and Adrianos guides 
# Then extract the meta data for your colonies and individuals based on the Extract meta data script
# Next do the base analyses (interactions and space use) in Vital_base_analysis.R
# Run Tables to match Stroeymeyt et al 2018 to get a couple more things in the right format for the following pipeline.
# Run 19_facet net to get forager nurse distribution via community detection
# Add pathogen/spore/bead load & new forager nurse allocation including space use to meta data (script 20) 

# Todo:
# Look at distribution of foragers and nurses to see which of the two task allocation methods is better? 
# Percentage per colony, overlap between the two different definitions. 

# Notes: 
# Due to difference in the framerate in comparison to adrianos experiments we so far did not generate the grooming interactions and anything regarding grooming is so far neglected down the line.

### Index ###
# 1. Read Me
# 2. Prerequisites - Input to define by User
# 3.  Analysis programs vital
# 3.1 11_randomise_interactions_DS.R 
# 3.2 12_simulate_transmission_DS.R
# 3.3 13_network_analysis.R
# 3.4 14_summarise_interactions.R
# 3.5 19_Facetnet_community_detection.R
# 4. All available analysis programs 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####     2. Prerequisites - Input to define by User        ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Make sure the below is directing to you data folder, the code/script folder, and the c++ executables folder (source scipts)

### input parameters:
FRAME_RATE <- 6

# Define which interactions to look at:                                                                                                  
{
RUN_CLASSIC_INTERACTIONS           <- TRUE
RUN_GROOMING_INTERACTIONS          <- FALSE
RUN_TROPHALLACTIC_INTERACTIONS     <- FALSE

# Define what analysis step to run: 
RUN_11_randomise_interactions_DS.R        <- FALSE
RUN_12_simulate_transmission_DS.R         <- FALSE
RUN_13_network_analysis_DS.R              <- TRUE
RUN_14_summarise_interactions_DS.R        <- FALSE
RUN_19_Facetnet_community_detection_DS.R  <- FALSE
RUN_19_SUBSECTION_ONLY                    <- FALSE  # Community detection on observed data only for task allocation --> run only the first section of the script (subsetted) -- Usually this is set to false, unless you are opening this script for the first time and want to get task allocation via the facet net method as alternative to the space use method in base analysis
}

# setwd("/home/ael/Documents/vital_rscripts_git")
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

source("config_user_and_hd.R") # contains getUserOptions() that defines usr and hd and the clean() function


# additional functions 
choose_data_path <- function() { # only works on linux
  list(
    CLASSIC_INTERACTIONS = if (RUN_CLASSIC_INTERACTIONS) paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment", sep="/") else NULL,
    GROOMING_INTERACTIONS = if (RUN_GROOMING_INTERACTIONS) paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment_grooming", sep="/") else NULL,
    TROPHALLACTIC_INTERACTIONS = if (RUN_TROPHALLACTIC_INTERACTIONS) paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment_trophallaxis", sep="/") else NULL
  )
}


# get the data and script paths as well as additional info tables, libraries, parameters and functions loaded: 
code_path   <- paste(SCRIPTDIR, "source_scripts", sep="/")  # paste("/home/",usr,"/Documents/vital_rscripts_git/source_scripts",sep="") # place where the needed r scripts are stored
data_paths  <- choose_data_path()
info        <- read.table(paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment/original_data/info.txt",sep="/"), header=T,stringsAsFactors = F) 
treated     <- read.table(paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment/original_data/treated_worker_list.txt",sep="/"),header=T,stringsAsFactors = F)
task_groups <- read.table(paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment/original_data/task_groups.txt",sep="/"),header=T,stringsAsFactors = F) # might need to be checked if this how the task groups are defined at the moment... % or facet net?  
tag_list    <- paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment/original_data/tag_files/",sep="/")
seed_path   <- paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment/original_data/seeds/",sep="/")
splitpath   <- paste("/media", usr, hd, "vital/fc2/vital_experiment/main_experiment/original_data/time_aggregation_info/",sep="/")
split_list  <- paste(splitpath,list.files(splitpath),sep="") # needed for transmission simulation
high_threshold <- 0.0411 * 0.5/0.3 # Science paper: high threshold =0.0411 where 1 = load of treated   # In Adriano's experiment, spore concentration was the same but volume was 0.3 microliter instead of 0.5 microliter --> needs to be adjusted so that it fits the virus or bead data (CHECK WHAT LUKE DID HERE?), Probably not required for my data as I had beads and not spores. 
to_keep <- c(ls(),"to_keep", "loop_start_time")



### ### ### ### ### ### ### ### ### ### ### ### 
####     3.  Analysis programs vital       ####
### ### ### ### ### ### ### ### ### ### ### ###


### START THE LOOP             ###
###                            ###
###                            ###
###                            ###

if (TRUE) {

#### 3.1 11_randomise_interactions_DS.R ####
# create randomized interaction networks based on the observed interactions 
  if (RUN_11_randomise_interactions_DS.R){ # RUN_11_randomise_interactions_DS.R <- FALSE
    print("Running 11_randomise_interactions_DS.R")
    for (interaction_type in names(data_paths)) { # interaction_type <- "CLASSIC_INTERACTIONS" 
      data_path <- data_paths[[interaction_type]]
      if (!is.null(data_path)) {
        print(paste("Processing file for", interaction_type, "\U0001F91D"))
        source(paste(code_path,"/11_randomise_interactions_DS.R",sep=""))
        print(paste(interaction_type, "ALL DONE", "\U0001F973", sep = " "))
        clean()
      } else {print(paste("Skipping randomisations for", interaction_type, sep = " ")) }
    }
    print(paste("Randomisations", "ALL DONE", "\U0001F973", sep = " "))
  } else {
    print("Skipping 11_randomise_interactions_DS.R")
  }
  

  
  
  
#### 3.2 12_simulate_transmission_DS.R ####
# simulate transmission of an agent from a list of originally contaminated workers to the rest of the colony. 

  if (RUN_12_simulate_transmission_DS.R){
    print("running 12_simulate_transmission_DS.R")
    for (interaction_type in names(data_paths)) { # interaction_type <- "TROPHALLACTIC_INTERACTIONS" or # interaction_type <- "CLASSIC_INTERACTIONS"
      data_path <- data_paths[[interaction_type]]
      if (!is.null(data_path)) {
        print(paste("Processing files for", interaction_type, "\U0001F91D"))
        source(paste(code_path,"/12_simulate_transmission_DS.R",sep=""))
        print(paste(interaction_type, "ALL DONE",  "\U0001F973"))
        clean()
      } else {print(paste("Skipping transmission simulations for", interaction_type)) }
    }
    print(paste("Simulations", "ALL DONE", "\U0001F973", sep = " "))
  } else {
    print("skipping 12_simulate_transmission_DS.R")
  }
  
  
  
  


#### 3.3 13_network_analysis.R ####
#  conducts network analysis on interaction data
if (RUN_13_network_analysis_DS.R) {
    print("running 13_network_analysis_DS.R")
    for (interaction_type in names(data_paths)) { # interaction_type <- "CLASSIC_INTERACTIONS" or interaction_type <- "TROPHALLACTIC_INTERACTIONS"
      data_path <- data_paths[[interaction_type]] 
      if (!is.null(data_path)) {
        print(paste("Processing files for", interaction_type, "\U0001F91D"))
        source(paste(code_path, "/13_network_analysis_DS.R", sep=""))
        print(paste(interaction_type, "ALL DONE", "\U0001F973"))
        clean()
      } else {
        print(paste("Skipping", interaction_type))}
      }} else {print("skipping 13_network_analysis_DS.R")}
  


  
#### 3.4 14_summarise_interactions.R ####
  
if (RUN_14_summarise_interactions_DS.R){ 
  print("summarize interactions currently being prepared and tested for Daniel")
  for (interaction_type in names(data_paths)) { # interaction_type <- "CLASSIC_INTERACTIONS" or interaction_type <- "TROPHALLACTIC_INTERACTIONS"
    data_path <- data_paths[[interaction_type]] 
    if (!is.null(data_path)) {
      print(paste("Processing files for", interaction_type, "\U0001F91D"))
      source(paste(code_path,"/14_summarise_interactions_DS.R",sep=""))
      clean()
      print(paste(interaction_type, "ALL DONE", "\U0001F973"))
    } else {
      print(paste("Skipping", interaction_type))}
  }} else {print("skipping 14_summarise_interactions_DS.R")}


  
#### 3.5 19_Facetnet_community_detection.R ####
if (RUN_19_Facetnet_community_detection_DS.R){
    print(paste(Sys.time(), ":    running 19_Facetnet_community_detection.R"))
    loop_start_time <- Sys.time()
    source(paste(code_path,"/19_Facetnet_community_detection_DS.R",sep=""))
    clean()
    loop_end_time <- Sys.time()
    print(paste(Sys.time(), ":    ALL DONE", "\U0001F973"))
    print(paste("script took ", as.numeric(difftime(loop_end_time, loop_start_time, units = "mins")), " minutes to complete"))
} else { print("skipping 19_Facetnet_community_detection_DS.R")} 
  
  
  

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
} 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### END OF THE LOOP OVER ALL FUNCTIONS ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### NEXT: Move on the analysis of the data using script "Vital_stats_and_plots_R"








#### 4. All available analysis programs ####
# check adrianos github for all the scripts should they be needed 

# source(paste(code_path,"/1_trackconverter.R",sep=""))
# clean()
# source(paste(code_path,"/2_define_deaths.R",sep=""))
# clean()
# source(paste(code_path,"/3_apply_rotation_to_datfiles.R",sep=""))
# clean()
# source(paste(code_path,"/4_retagged_ant_modifications.R",sep=""))
# clean()
# source(paste(code_path,"/5_zoneconverter_nest.R",sep=""))
# clean()
# source(paste(code_path,"/6_time_investment.R",sep=""))
# clean()
# source(paste(code_path,"/7_trajectory.R",sep=""))
# clean()
# source(paste(code_path,"/8_process_trajectory_files.R",sep=""))
# clean()
# source(paste(code_path,"/9_interaction_detection.R",sep=""))
# clean()
# source(paste(code_path,"/10_process_interaction_files.R",sep=""))
# clean()

# TO USE AW:
# source(paste(code_path,"/11_randomise_interactions.R",sep=""))
# clean()
# source(paste(code_path,"/12_simulate_transmission.R",sep=""))
# clean()
# source(paste(code_path,"/13_network_analysis.R",sep=""))
# clean()
# source(paste(code_path,"/14_summarise_interactions.R",sep=""))
# clean()
# source(paste(code_path,"/15_heatmaps_individual.R",sep=""))
# clean()
# source(paste(code_path,"/16_heatmaps_groups.R",sep=""))
# clean()
# source(paste(code_path,"/17_brood_location.R",sep=""))
# clean()
# source(paste(code_path,"/18_process_heatmaps.R",sep=""))
# clean()
# source(paste(code_path,"/19_Facetnet_community_detection.R",sep=""))
# clean()