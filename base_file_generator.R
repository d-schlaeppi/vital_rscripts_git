rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      #### Base File Generator  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####

# automatically generate the myrmidon base files for all tracking files

#### prerequisites ####

# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()

setwd("/home/gw20248/Documents/vital_rscripts_git/")
dat <- read.csv("fc2_overview_data.csv", header = TRUE, stringsAsFactors = F)

# set working directory
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(directory)
list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE)


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
  tracking_system_main  <- dat[i, "tracking_system_main"]
  tracking_system_feeding <- dat[i, "tracking_system_feeding"] 
  # combine variables to a data frame  
  data_collection <-  rbind(data_collection, data.frame(nr, 
                                                        colony_id,
                                                        block, 
                                                        colony_nr,
                                                        treatment,
                                                        food_position_1,
                                                        food_position_2,
                                                        tracking_system_main,
                                                        tracking_system_feeding,
                                                        stringsAsFactors = F))
  }


# get all folders in the directory and compile them as a list (only folders containing the tracking data)
data_list <- list.files(path=directory, pattern=NULL, all.files=FALSE, full.name=FALSE) # all files
data_list <- grep(data_list, pattern = '.myrmidon', invert = TRUE, value = TRUE)
data_list <- grep(data_list, pattern = '.txt', invert = TRUE, value = TRUE)

# first an empty .myrmidon file is created manually as a source file using fort studio: base_source.myrmidon this is then used to create the rest of the files... 

tracking_type <- c("main", "feeding")
for (i in 1:nrow(data_collection)) {
  for (type in tracking_type) {
    if (type == "main") {
      file_name <- paste0(data_collection[i,"colony_nr"], '_m_base.myrmidon') # file name will be of the following structure: c01_m_base.myrmidon (m = main tracking vs f = feeding)
      tracking_data <- fmExperimentOpen("base_source.myrmidon")
      s <- tracking_data$createSpace(data_collection[i,"tracking_system_main"])
      printf("Space '%s' has ID: %d\n",s$name,s$ID)
      tracking_data$name <- paste0("vital fc2 ",data_collection[i,"colony_nr"], " main")
      # assign the tracking data
      for (folder_name in data_list) {
        if (grepl(data_collection[i,"colony_nr"], folder_name, fixed = TRUE) & !grepl("feeding", folder_name, fixed = TRUE)) {
          tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name))
        } else {next}
      }
      # save the file base file with created ants 
      tracking_data$save(paste0(directory, data_collection[i,"colony_nr"], '_main.myrmidon'))
    } else {
      file_name <- paste0(data_collection[i,"colony_nr"], '_f_base.myrmidon') # file name will be of the following structure: c01_m_base.myrmidon (m = main tracking vs f = feeding)
      tracking_data <- fmExperimentOpen("base_source.myrmidon")
      s <- tracking_data$createSpace(data_collection[i,"tracking_system_feeding"])
      printf("Space '%s' has ID: %d\n",s$name,s$ID)
      tracking_data$name <- paste0("vital fc2 ",data_collection[i,"colony_nr"], " feeding")
      # assign the tracking data
      for (folder_name in data_list) {
        if (grepl(data_collection[i,"colony_nr"], folder_name, fixed = TRUE) & grepl("feeding", folder_name, fixed = TRUE)) {
          tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name))
        }
      }
      # save the file base file with created ants 
      tracking_data$save(paste0(directory, data_collection[i,"colony_nr"], '_feeding.myrmidon'))
    }
  }
}



#### once the issue with the tracking data valid method error is resolved The loop above should work to create all the base files..
# if it does not work the files can be created automatically, but the directories need to be added manually which would be annoying. 
# run the above code on the extrapolated data.
# if the code above runs, include the ant orientation code and the capsule code for a quick post processing, 
# then apply the ant pose cloner to synchronize feeding and main tracking files.
# (check if there is indeed no scaling required to make the capsule cloning work)
# then, run the ant orient_express!






# below is just some practice code that can be deleted once the aboce loop works. 


file_name <- paste0(data_collection[i,"colony_nr"], '_f_base.myrmidon') # file name will be of the following structure: c01_m_base.myrmidon (m = main tracking vs f = feeding)
tracking_data <- fmExperimentOpen("base_source.myrmidon")
s <- tracking_data$createSpace(data_collection[i,"tracking_system_feeding"])
printf("Space '%s' has ID: %d\n",s$name,s$ID)
tracking_data$name <- paste0("vital fc2 ",data_collection[i,"colony_nr"], " feeding")
tracking_data$save(paste0(directory, data_collection[i,"colony_nr"], '_feeding.myrmidon'))


# assign the tracking data
tracking_data$addTrackingDataDirectory(s$ID, paste0(directory, folder_name))
tracking_data$addTrackingDataDirectory(s$ID, x)

x <- '/home/gw20248/Documents/data_copy_for_trials/vital_fc2_esterhase_c01_feeding_DS.0000'
paste0(directory, folder_name)

for (folder_name in data_list) {
  if (grepl(data_collection[i,"colony_nr"], folder_name, fixed = TRUE) & grepl("feeding", folder_name, fixed = TRUE)) {
    tracking_data$addTrackingDataDirectory(s$ID, paste0(directory,folder_name))
  }
}
# save the file base file with created ants 
tracking_data$save(paste0(directory, data_collection[i,"colony_nr"], '_feeding.myrmidon'))










# create ants
tag_statistics <- fmQueryComputeTagStatistics(tracking_data)
for (y in 1:nrow(tag_statistics)) {  # loop over each tag
  if (tag_statistics[y,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T)) {
    a <- tracking_data$createAnt();
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[y,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
}

tag_statistics <- fmQueryComputeTagStatistics(tracking_data)
for (x in 1:nrow(tag_statistics)) {  # loop over each tag
  if (tag_statistics[x,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) {
    a <- tracking_data$createAnt();
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[x,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
}




# finding everything with colony nr
for (folder_name in data_list) {
  if (grepl("c02", folder_name, fixed = TRUE)) {print(folder_name)}
}

# finding everything with colony nr but excluding feeding
for (folder_name in data_list) {
  if (grepl("c02", folder_name, fixed = TRUE) & !grepl("feeding", folder_name, fixed = TRUE)) {print(folder_name)} 
}

# finding everything with colony nr and feeding
for (folder_name in data_list) {
  if (grepl("c02", folder_name, fixed = TRUE) & grepl("feeding", folder_name, fixed = TRUE)) {print(folder_name)}
}



?grepl
grepl("c27", "vital_fc2_trojan_c27_DS.0000", fixed = TRUE)

??addTrackingDataDirectory

###add tracking data directory
tddURI <- tracking_data$addTrackingDataDirectory(s$ID,paste(dir_data,tracking_data_file,sep=''))
tracking_data$save(paste(dir_data,auto_orient_file,sep='')) # file now exists

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
paste(substr(dataset_name, 1, nchar(dataset_name)-9),'_metaID.myrmidon')
print(tracking_data_file)













