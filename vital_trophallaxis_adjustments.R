rm(list = ls())
gc()
Sys.sleep(3)
mallinfo::malloc.trim(0L) 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### VITAL - Adjusting inferred trophy interactions  ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# This script contais the code to get the inferred trophallactic interactions list into the same format as the classic interaction lists (3 h chunks, pre/post treatment and storage in the same format)

library(FortMyrmidon)
library(dplyr)


#### 1. Starting parameters & directories ####
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

WORKDIR <- paste("/media/",usr, hd, "/vital/fc2",sep="") 
DATADIR <- paste(WORKDIR, sep = "/") # working directory and directory where the tracking data is saved (was not the same for Adriano, but for me the same
THROPHY_DIR <-paste("/media/",usr, hd, "/vital/fc2/vital_experiment/main_experiment_trophallaxis/",sep="") 
SCRIPTDIR <- paste("/home/",usr,"/Documents/vital_rscripts_git",sep="") 

source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # will add colony_metadata data frame to the environment so it can be accessed within this script (in my case containing bodylenght information)
metadata_colonies <- colony_metadata

files_list <- list.files(paste0(THROPHY_DIR,"inferred_trophallaxis_whole_experiment/"))
files_list <- files_list[which(!grepl("c29",files_list))]

for (file in files_list) {
  colony_id <- unlist(strsplit(file,split="_"))[grepl("c",unlist(strsplit(file,split="_")))]
  cat(paste("########################################\n",basename(file)),sep="")    
  treatment <- colony_metadata[which(colony_metadata$colony_id==colony_id),"treatment_simple"]
  Interactions <- read.table(paste0(THROPHY_DIR,"inferred_trophallaxis_whole_experiment/", file), header = T, sep = ",",  stringsAsFactors = F)
  Interactions <- Interactions %>%
    rename(Tag1 = ant1, Tag2 = ant2, Starttime = T_start_sec, Stoptime = T_stop_sec)   # renaming to match adrianos version for future scripts... 
  for (PERIOD in c("pre","post")) {
    period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
    cat(paste("#######","period:", PERIOD, sep= " "))
    
    #define start and stop times to split the data
    treatment_start <-  as.POSIXct(metadata_colonies[    which   (   metadata_colonies$colony_id == colony_id   )     ,   "time_treatment_start"     ], format="%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    treatment_end   <-  as.POSIXct(metadata_colonies[    which   (   metadata_colonies$colony_id == colony_id   )     ,   "time_treatment_end"     ], format="%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    if (PERIOD=="pre"){
      time_start_GMT <- treatment_start - 24*3600
      time_end_GMT <- treatment_end - 24*3600
      time_label <- "_TH-24_TD11" # just used to match some of Adriano's labeling (needs to be adjusted if longer time chunks are analysed i.e. multiple 3 h blocks in the time dictionary) 
    }else{
      time_start_GMT <- treatment_start
      time_end_GMT <- treatment_end
      time_label <- "_TH0_TD11"
    }
    
    # subset the interaction list with trophallactiv interactions inferred for the whole recording time
    Interactions$start <- as.POSIXct(Interactions$start, format="%Y-%m-%d %H:%M:%S", tz="GMT")
    Interactions$end <- as.POSIXct(Interactions$end, format="%Y-%m-%d %H:%M:%S", tz="GMT")
    Interaction_subset <- Interactions[which(Interactions$start >= time_start_GMT & Interactions$end <= time_end_GMT),]
    
    # save ~3h subsets to pre and post interaction folders (full and binned - in my case the same, but I just want to make sure to have them in both places)
    outputfolder_1 <- paste0(THROPHY_DIR, "intermediary_analysis_steps/full_interaction_lists/", period_code, "/observed")
    if (!file.exists(outputfolder_1)) {dir.create(outputfolder_1, recursive = TRUE)}
    full_file_name <- paste0("/",colony_id, "_", treatment,"_", period_code, "_trophallactic_interactions.txt") 
    write.table(Interaction_subset, file = paste0(outputfolder_1, full_file_name), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
    
    outputfolder_2 <- paste0(THROPHY_DIR, "intermediary_analysis_steps/binned_interaction_lists/", period_code, "/observed")
    if (!file.exists(outputfolder_2)) {dir.create(outputfolder_2, recursive = TRUE)}
    binned_file_name <- paste0("/",colony_id, "_", treatment,"_", period_code, time_label, "_trophallactic_interactions.txt")
    write.table(Interaction_subset, file = paste0(outputfolder_2, binned_file_name), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
  }
}

