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

USER <- "2A13_Office_Daniel"  # Replace with the desired USER option: Nath_office, 2A13_Office_Adriano, 2A13_Office_Daniel, AEL-laptop
HD <- "/DISK_B" # One of the harddrives with the vital extrapolated data: possible values > Nathalies hd "/DISK_B" ; Daniels hds >  "/gismo_hd5" or  "/gismo_hd2" | just make sure to put the name of your HD in here

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
setUserAndHD(USER, HD)


choose_data_path <- function() {
  # Ask user which data they would like to process
  choice <- as.character(readline(prompt = "Which data would you like to process?\nAll standard interactions (ala Stroeymeyt et al 2018)  (1)\nGrooming interactions (2)\ \nTrophallactic interactions (3)\n"))
  # Set data_path depending on user's choice
  if (choice == "1") {
    data_path <- paste0("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment")
    print("ALL INTERACTIONS")
  } else if (choice == "2") {
    data_path <- paste0("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment_grooming")
    print("GROOMING INTERACTIONS")
  } else if (choice == "3") {
    data_path <- paste0("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment_trophallaxis")
    print("TROPHALLACTIC INTERACTIONS")
  } else {
    stop("Invalid choice entered. Please choose either '1, '2' or '3'.")
  }
  return(data_path)
}

# Call the function to choose the data path
data_path <- choose_data_path()


code_path  <- "/home/cf19810/Ant_Tracking/Scripts/code_Social_Network_Plasticity_Exp_2018_AW/1_data_post_processing/source"
executables_path <- "~/executables"
FRAME_RATE <- 8 #AW


########################################################################
###   END INPUT ########################################################
########################################################################

source(paste(code_path,"/libraries.R",sep=""))
source(paste(code_path,"/functions_and_parameters.R",sep=""))

to_keep <- c(ls(),"to_keep")

#####Run analysis programs ####
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