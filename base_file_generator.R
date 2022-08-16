rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      #### Base File Generator  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####

# automatically generate the myrmidon base files for all tracking files

#
setwd("/home/gw20248/Documents/vital_rscripts_git/")
dat <- read.csv("fc2_overview_data.csv", header = TRUE, stringsAsFactors = F)

# set working directory
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(directory)
  
# get all folders in the directory and compile them as a list (only folders containing the tracking data)
data_list <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
data_list <- grep(data_list, pattern = '.myrmidon', invert = TRUE, value = TRUE)

# idea: create a data frame containing information for each colony? and merge it to the data from the masterfile? 

data_collection <- NULL 

for(i in 1:nrow(dat)) {
  # collect variables
  nr                    <- i
  colony_id             <- dat[i, "colony_id"]
  block                 <- dat[i, "block"]
  colony_nr             <- paste0("c", sprintf("%02d", i))
  treatment             <- dat[i, "treatment"]
  food_position_1       <- dat[i, "food_position_1"]
  food_position_2       <- dat[i, "food_position_2"]
  main_tracking_system  <- "NA"
  # combine variables to a data frame  
  data_collection <-  rbind(data_collection, data.frame(nr, 
                                                        colony_id,
                                                        block, 
                                                        colony_nr,
                                                        treatment,
                                                        food_position_1,
                                                        food_position_2,
                                                        main_tracking_system, 
                                                        stringsAsFactors = F))
  }


# starting to copy the ideas from enricos script?!
# https://github.com/EnricoGavagnin/auto_vs_manual_orientation_check/blob/main/auto_orientation_loop.R

for (tracking_data_file in data_list) {
  if(substr(tracking_data_file, nchar(tracking_data_file)-0, nchar(tracking_data_file)) !="0") {next}
  print(tracking_data_file)
  # base file name
  
  
  base_file_name <- paste0(substring(tracking_data_file,11,nchar(tracking_data_file)-8),'_base.myrmidon')
  
}

print(substr(tracking_data_file, nchar(tracking_data_file)-1, nchar(tracking_data_file)))


substr(tracking_data_file, nchar(tracking_data_file)-(nchar(tracking_data_file)), nchar(tracking_data_file))
substr(tracking_data_file, nchar(tracking_data_file)-0, nchar(tracking_data_file))



print(tracking_data_file)
nchar(tracking_data_file)

?substr

substr(dataset_name, 1, nchar(dataset_name)-9),'_metaID.myrmidon')


print(tracking_data_file)