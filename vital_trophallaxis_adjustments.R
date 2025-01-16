# rm(list = ls())
rm(list = setdiff(ls(), "first_time_use_working_directory"))
gc(); mallinfo::malloc.trim(0L)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### VITAL - Adjusting inferred trophy interactions  ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# This script contais the code to get the inferred trophallactic interactions list into the same format as the classic interaction lists (3 h chunks, pre/post treatment and storage in the same format)


#### 1. Starting parameters & directories ####

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
THROPHY_DIR <- paste(DATADIR, "/vital_experiment/main_experiment_trophallaxis/",sep="") 

source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # will add colony_metadata data frame to the environment so it can be accessed within this script (in my case containing bodylenght information)


files_list <- list.files(paste0(THROPHY_DIR,"inferred_trophallaxis_whole_experiment/"))
files_list <- files_list[which(!grepl("c29",files_list))]

if(TRUE){
cat(blue ("\n Processing : \n"))  
 for (file in files_list) { # file <- files_list[1]
  colony_id <- unlist(strsplit(file,split="_"))[grepl("c",unlist(strsplit(file,split="_")))]
  cat(white(paste0("\r ", basename(file))))    
  treatment <- colony_metadata[which(colony_metadata$colony_id==colony_id),"treatment_simple"]
  Interactions <- read.table(paste0(THROPHY_DIR,"inferred_trophallaxis_whole_experiment/", file), header = T, sep = ",",  stringsAsFactors = F)
  Interactions <- Interactions %>%
    rename(Tag1 = ant1, Tag2 = ant2, Starttime = T_start_sec, Stoptime = T_stop_sec)   # renaming to match adrianos version for future scripts... 
  for (PERIOD in c("pre","post")) { # PERIOD <- "pre"
    period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
    #cat(paste("#######","period:", PERIOD, sep= " "))
    
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
    
    # Adjust variables so fortmatting matches perfectly with that of classic interactions
    Interaction_subset <- Interaction_subset %>%
      dplyr::select(-PERIOD) %>% # Remove the PERIOD variable
      rename(
        Startframe = startframe,
        Endframe = endframe,
        Xcoor1 = ant1.mean.x,
        Ycoor1 = ant1.mean.y,
        Angle1 = ant1.mean.angle,
        Xcoor2 = ant2.mean.x,
        Ycoor2 = ant2.mean.y,
        Angle2 = ant2.mean.angle,
        Detections = detections,
        Box = space,
        colony = REPLICATE
      ) %>%
      # Add new variables
      mutate(
        period = PERIOD,
        time_hours = ifelse(PERIOD == "pre", -24, 0),
        time_of_day = 11,
        treatment = colony_metadata[which(colony_metadata$colony_id == colony_id), "treatment_simple"],
        treatment_detailed = colony_metadata[which(colony_metadata$colony_id == colony_id), "treatment"],
        REP_treat = paste(treatment, Box, colony_id, sep = "_"),
        Direction = types
      )
    
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
  cat(yellow("\n All done... \n"))
}






            


            