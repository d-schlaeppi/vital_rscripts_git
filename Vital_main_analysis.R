rm(list = ls())
gc()
Sys.sleep(3)
mallinfo::malloc.trim(0L)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### VITAL MAIN ANALYSIS     ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### 1. Background information | Read me | ToDos ####
# Before starting this scrip, first start with the pre processing of the data following Vital.R and Adrianos guides 
# Then extract the meta data for your colonies and individuals based on the Extract meta data script
# Next do the base analyses (interactions and space use) in Vital_base_analysis.R
# Run Tables to match Stroeymeyt et al 2018 to get a couple more things in the right format for the following pipeline.
# If not already done add the data on pathogen or spore load to each individual (e.g. by merging it with the meta data file) 
# and when doing so also adjust blanks
# Then continue here. 

# This script contains code to run separate analysis scripts from the stroeymeyt 2018 pipeline: 

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
RUN_CLASSIC_INTERACTIONS           <- TRUE
RUN_GROOMING_INTERACTIONS          <- FALSE
RUN_TROPHALLACTIC_INTERACTIONS     <- FALSE

# Define what analysis step to run: 
{RUN_11_randomise_interactions_DS.R     <- FALSE
RUN_12_simulate_transmission_DS.R      <- FALSE
RUN_13_network_analysis.R              <- FALSE
RUN_14_summarise_interactions.R        <- FALSE
RUN_19_Facetnet_community_detection.R  <- FALSE}

setwd(tk_choose.dir(default = "~/", caption = "Select Working Directory")) # direct it to where you have config_user_and_hd.R (typically your scriptfolder)
source("config_user_and_hd.R") # contains getUserOptions() that defines usr and hd and the clean() function

choose_data_path <- function() { # does not work on the mac. 
  list(
    CLASSIC_INTERACTIONS = if (RUN_CLASSIC_INTERACTIONS) paste0("/media/", usr, hd, "/vital/fc2/vital_experiment/main_experiment") else NULL,
    GROOMING_INTERACTIONS = if (RUN_GROOMING_INTERACTIONS) paste0("/media/", usr, hd, "/vital/fc2/vital_experiment/main_experiment_grooming") else NULL,
    TROPHALLACTIC_INTERACTIONS = if (RUN_TROPHALLACTIC_INTERACTIONS) paste0("/media/", usr, hd, "/vital/fc2/vital_experiment/main_experiment_trophallaxis") else NULL
  )
}


# Call the functions to get the data and script paths as well as additional libraries and functions loaded: 
code_path <- paste("/home/",usr,"/Documents/vital_rscripts_git/source_scripts",sep="") # place where the needed r scripts are stored
data_paths <- choose_data_path()
source(paste(code_path,"/libraries_DS.R",sep=""))
source(paste(code_path,"/functions_and_parameters_DS.R",sep=""))

# executables_path <- "~/executables" # not clear what this is and where it is used. ?
to_keep <- c(ls(),"to_keep")


### ### ### ### ### ### ### ### ### ### ### ### 
####     3.  Analysis programs vital       ####
### ### ### ### ### ### ### ### ### ### ### ###


#### 3.1 11_randomise_interactions_DS.R ####

# create randomized interaction networks based on the observed interactions 
if (RUN_11_randomise_interactions_DS.R){
  print("runnung 11_randomise_interactions_DS.R")
  for (interaction_type in names(data_paths)) {
    data_path <- data_paths[[interaction_type]]
    if (!is.null(data_path)) {
      print(paste("Processing file for", interaction_type, "\U0001F91D"))
      source(paste(code_path,"/11_randomise_interactions_DS.R",sep=""))
      clean()
    }
  print(paste("ALL DONE", "\U0001F973"))
  }} else {
  print("skipping 11_randomise_interactions_DS.R") 
}
  

#### 3.2 12_simulate_transmission_DS.R ####
# simulate transmission of an agent from a list of originally contaminated workers to the rest of the colony. 
for (interaction_type in names(data_paths)) {
  data_path <- data_paths[[interaction_type]]
  #data_path <- data_paths[["TROPHALLACTIC_INTERACTIONS"]]
  #data_path <- data_paths[["CLASSIC_INTERACTIONS"]]
  if (!is.null(data_path)) {
    print(paste("Processing files for", interaction_type, "\U0001F91D"))
    source(paste(code_path,"/12_simulate_transmission_DS.R",sep=""))
    clean()
  }
  print(paste("ALL DONE", "\U0001F973"))
}

#### 3.3 13_network_analysis.R ####
source(paste(code_path,"/13_network_analysis.R",sep=""))
clean()

#### 3.4 14_summarise_interactions.R ####
source(paste(code_path,"/14_summarise_interactions.R",sep=""))
clean()


#### 3.5 19_Facetnet_community_detection.R ####
if (RUN_19_Facetnet_community_detection.R){
  print("running 19_Facetnet_community_detection.R")
  source(paste(code_path,"/19_Facetnet_community_detection.R",sep=""))
  clean()
  print(paste("ALL DONE", "\U0001F973"))
    } else { print("skipping 19_Facetnet_community_detection.R")}




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