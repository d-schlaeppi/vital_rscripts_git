# this script contains
# 1. Automated ant creation for fort-studio

rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
     #### 1. AUTOMATED ANT CREATION FOR FORT STUDIO  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####

# for more information on fort-myrmidon and fort-studio see: 
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"

# When you create ants directly in fort-studio, antIDs will be created in the order in which you orient the ants
# This means that if two users were to orient the same experiment, but did it in a slightly different orders, they would end up with different antID lists where different antID would correspond to different tagID
# It also means that if for any reasons you wanted to replicate your results / reorient your ants from scratch, and wanted to compare results to your previous analysis, you may end up with a different antID list and this would make the comparison difficult
# One way to avoid that is to create the ants from R using the script below
# First, a base myrmidon file needs to be created with fill in first page (defining the space, author, tag size and comment). LABEL IT WITH "your_experiment_details" + _base.myrmidon (your_experiment_details_base.myrmidon)
# Then, run the following script in RStudio (remember to open it in right environment) to automatically create the ants in order (will be consistent over time)
# The script includes a  filter to remove tags that had very few detections compared to others (will be discarded as false detections / low-quality ants
# You will need to open the experiment as read/write as you will want to overwrite your "blank" myrmidon file with a file that contains antIDs

### Advantages of creating ants in R
# 1. Repeatably (as explained above)
# 2. Useful when you want to identify and display ants with specific properties halfway through the experiment, and you don't have the time to orient them manually (e.g. Adriano's experiment requiring to identify nurses)
### Disadvantages: 
# Does not allow to have several tags pointing to the same ant (retag). 
# BUT this can be altered manually in fort-studio afterwards by deleting the second ant and then adding the second tagID as identification to the first ant

#### prerequisites ####

### load necessary libraries
library(FortMyrmidon)   #### R bindings
library(R.utils)        #### contains printf()




#### things to modify !!! ####

### directory and file you will use 
# file_name <- 'vital_fc2_esterhase_c01_feeding_DS_base.myrmidon'  #insert your base file name here. NOTE: Files names are supposed to end with base.myrmidon for the paste functions to work
file_name <- 'vital_fc2_esterhase_c02_feeding_DS_base.myrmidon'
file_name <- 'vital_fc2_guillam_c12_feeding_DS_base.myrmidon'
file_name <- 'vital_fc2_guillam_c27_feeding_DS_base.myrmidon'




# notes for tomorrow: files to do for manual orientation: 1x esterhase feeding, 1 oder 2x guillam feeding, 2tes mal guillam after feeding 
# directory <- '/media/gw20248/gismo_hd2/vital/fc2/' # insert your directory here
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(directory)





#### code to run without further modification ####

### tracking file you would like to create ants for: 
tracking_data <- fmExperimentOpen(file_name)

### extract the tag statistics to know how many times each tag was detected:
tag_statistics <- fmQueryComputeTagStatistics(tracking_data)

### create ants
for ( i in 1:nrow(tag_statistics)) {  #####loop over each tag
  if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
}

### create ants with a stronger cut off to reduce the chances of false positives if too many false ants are created 
#tag_statistics[,"count"] #check for the count numbers to see what range should be included (in this case the cutoff is reduced from 0.001 to 0,01)
for ( i in 1:nrow(tag_statistics)) {  #####loop over each tag
  if ( tag_statistics[i,"count"] >= 0.01*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
 }


### Check: print identifications
ants <- tracking_data$ants
for (a in ants) {
  printf("Ant %s is identified by:\n", fmFormatAntID(a$ID))
  for (i in a$identifications){
    printf(" * %s\n", capture.output(i))
  }
}

### Save to newly named myrmidon file: "previousname_AntsCreated"
tracking_data$save(paste0(directory, substr(file_name, 1, nchar(file_name)-13), 'AntsCreated.myrmidon'))

# paste0 prints the strings given to it with no separator (if you need a separator use paste() and define sep = "")
# substr() can be used to subset a string - e.g. a subset of a string: substr(string, starting_position, last_character)  --> in our case we started with the first character and printed until the full length minus 13 characters to substract base.myrmidon

rm(list=c("tracking_data"))
gc()







