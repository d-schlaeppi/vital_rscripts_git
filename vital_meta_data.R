### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Vital meta data ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####
# Quick script to read in some of the essential data from the vital tracking experiment and safe it in a data frame to be accessed in other scripts
# to load the data frame in another script just run the following line 
# source(paste0(directory_scripts,"vital_meta_data.R"))
#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"


dat <- read.csv(paste0(directory_scripts,"fc2_overview_data.csv"), header = TRUE, stringsAsFactors = F)
meta_data <- NULL
for(i in 1:nrow(dat)) {
  # collect variables
  nr                    <- i
  experiment            <- dat[i, "experiment"] 
  colony_id             <- dat[i, "colony_id"]
  block                 <- dat[i, "block"]
  colony_nr             <- paste0("c", sprintf("%02d", i))
  treatment             <- dat[i, "treatment"]
  treatment_simple      <- ifelse(dat[i, "treatment"] == "cc", "control", "virus")
  food_position_1       <- dat[i, "food_position_1"]
  food_position_2       <- dat[i, "food_position_2"]
  tracking_system_main  <- dat[i, "tracking_system_main"]
  tracking_system_feeding <- dat[i, "tracking_system_feeding"]
  annotation_start      <- dat[i, "annotation_start_time"] 
  annotation_stop       <- dat[i, "annotation_stop_time"]
  mean_ant_lenght_px    <- dat[i, "mean_ant_length_px"] # mean lenght in pixcel of the ants based on a manually oriented colony recorded with the same tracking system & setup (assuming that worker size is consistent accros colonies)
  
  # combine variables to a data frame  
  meta_data <-  rbind(meta_data, data.frame            (nr,
                                                        experiment, 
                                                        colony_id,
                                                        block, 
                                                        colony_nr,
                                                        treatment,
                                                        treatment_simple,
                                                        food_position_1,
                                                        food_position_2,
                                                        tracking_system_main,
                                                        tracking_system_feeding,
                                                        annotation_start, 
                                                        annotation_stop,
                                                        mean_ant_lenght_px,
                                                        stringsAsFactors = F))
}

# eventually to be updated with the complete meta data as adriano did with ant individual level data
