# this script contains: 
# 1. ANT-GENERATOR to define ant indentifications in Fort

rm(list=ls())
gc()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1. ANT Generator  ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### INFORMATION ####

# create a list with all the tracking files for which a measurement is needed and run the loop to create the ants in the fort myrmidon files
# for more information on fort-myrmidon and fort-studio see: 
# https://formicidae-tracker.github.io/myrmidon/latest/index.html

#### PREREQUISITES ####

# load required libraries
library(FortMyrmidon) # R bindings

directory <-  '/media/gw20248/gismo_hd5/trophy_data/' 
setwd(directory)

# list your files  (for larger numbers directly access your files and extract file names from directory to create the list automatically instead of manually)
files <- list(
  paste(directory,"new_trophy_01_base.myrmidon",sep='') #, 
  # paste(directory,"vital_fc2_trojan_c27_DS_base.myrmidon",sep=''), and so on...
)

for (file in files) {
  tracking_data <- fmExperimentOpen(file) # files that need the ants created
  tag_statistics <- fmQueryComputeTagStatistics(tracking_data)   # extract the tag statistics to know how many times each tag was detected etz
  # create ants - decide on the cutoff: stronger cut off to reduce the chances of false positives if too many false ants are created | tag_statistics[,"count"] #check for the count numbers to see what range should be included (in this case the cutoff is reduced from 0.001 to 0,01)
  for ( i in 1:nrow(tag_statistics)) {  # loop over each tag
    if ( tag_statistics[i,"count"] >= 0.01*max(tag_statistics[,"count"],na.rm=T) ) { # antID only if the tag detection rate was more than 1/100 (adriano used 1/1000) of the best tag detection rate
      a <- tracking_data$createAnt(); # creates an antID, i.e. associates a decimal antID number to that particular tagID
      identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
      print(identification)
    }
  }
  tracking_data$save(paste0(substr(file, 1, nchar(file)-13), 'AntsCreated.myrmidon')) # Save to newly named myrmidon file: "previousname_AntsCreated"
}

