rm(list = ls())
gc()
Sys.sleep(3)
mallinfo::malloc.trim(0L)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### VITAL MAIN ANALYSIS     ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Background information | Read me ####

# Before starting this scrip, first start with the pre processing of the data following Vital.R and Adrianos guides 
# Then extract the meta data for your colonies based on the Extract meta data script
# Next do the base analyses (interactions and space use) in Vital_base_analysis.R 
# Run Tables to match Stroeymeyt et al 2018 to get a couple more things in the right format for the following pipeline. 
# Then continue here. 

# This script contains:

# Index


# 1. Input definition 




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1. INPUT TO DEFINE BY USER ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Please fill in the path to the data folder, the code folder, and the c++ executables folder, as in the example below

#### Todo ####
# create the functions to run the user input super easily.

#### Prerequisites ####

# input parameters
USER <- "2A13_Office_Daniel"  # Replace with the desired USER option: Nath_office, 2A13_Office_Adriano, 2A13_Office_Daniel, AEL-laptop
HD <- "/DISK_B"               # One of the harddrives with the vital extrapolated data: possible values > Nathalies hd "/DISK_B" ; Daniels hds >  "/gismo_hd5" or  "/gismo_hd2" | just make sure to put the name of your HD in here
FRAME_RATE <- 6

# Define which interactions to look at:                                                                                                  
RUN_CLASSIC_INTERACTIONS           <- TRUE
RUN_GROOMING_INTERACTIONS          <- FALSE
RUN_TROPHALLACTIC_INTERACTIONS     <- FALSE


# functions
clean <- function(){
  rm(list=ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv)%in%to_keep], envir = .GlobalEnv)
  no_print <- gc(verbose=F)
  Sys.sleep(1) }

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
  hd <- NULL
  hd <- HD
  if (!is.null(hd)) {print(hd)} else {print("define new hd if necessary")}
  assign("hd", hd, envir = .GlobalEnv)  # Assign to the global environment
  assign("usr", usr, envir = .GlobalEnv)  # Assign to the global environment
}

choose_data_path <- function() {
  list(
    CLASSIC_INTERACTIONS = if (RUN_CLASSIC_INTERACTIONS) paste0("/media/", usr, hd, "/vital/fc2/vital_experiment/main_experiment") else NULL,
    GROOMING_INTERACTIONS = if (RUN_GROOMING_INTERACTIONS) paste0("/media/", usr, hd, "/vital/fc2/vital_experiment/main_experiment_grooming") else NULL,
    TROPHALLACTIC_INTERACTIONS = if (RUN_TROPHALLACTIC_INTERACTIONS) paste0("/media/", usr, hd, "/vital/fc2/vital_experiment/main_experiment_trophallaxis") else NULL
  )
}


# Call the functions to get the data and script paths as well as additional libraries and functions loaded: 
setUserAndHD(USER, HD)
code_path <- paste("/home/",usr,"/Documents/vital_rscripts_git/source_scripts",sep="") # place where the needed r scripts are stored
data_paths <- choose_data_path()
source(paste(code_path,"/libraries_DS.R",sep=""))
source(paste(code_path,"/functions_and_parameters_DS.R",sep=""))

# executables_path <- "~/executables" # not clear what this is and where it is used. ?
to_keep <- c(ls(),"to_keep")





#### Used analysis programs vital ####

# create randomized interaction networks based on the observed interactions 

# for (interaction_type in names(data_paths)) {
#   data_path <- data_paths[[interaction_type]]
#   if (!is.null(data_path)) {
#     print(paste("Processing file for", interaction_type, "\U0001F91D"))
#     source(paste(code_path,"/11_randomise_interactions_DS.R",sep=""))
#     clean()
#   }
#  print(paste("ALL DONE", "\U0001F973"))
# }


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


source(paste(code_path,"/13_network_analysis.R",sep=""))
clean()

source(paste(code_path,"/14_summarise_interactions.R",sep=""))
clean()

source(paste(code_path,"/19_Facetnet_community_detection.R",sep=""))
clean()




#### All available analysis programs ####
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