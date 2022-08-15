rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
     #### 1. ANT-ID-MATCHER  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####

# If ants are used in two separate tracking systems it is useful (or even essential) to make sure that based on the tagID the same AntID's are created in both files
# (e.g. tracking of a colony in one tracking system and then a few workers get sampled for treatment that is also tracked and then the workers get returned to their colony)
# when creating antID automatically based on tagID, antID are automatically created starting with 001. That means you will have e.g 200 in your main colony and you take a random sample of X workers from that
# in the new tracking file you also need to assign antID ant they will be automatically created from 001-0XX --> antID will not automatically be matched with the original antID 
# Because antID is read only and can not be overwritten we will create a new variable - meta_ID both for the main tracking and for the treatment tracking files
# For anything following this step we will need to call meta_ID instead of ant_ID if we refer to ants to avoid any confusion

# for more information on fort-myrmidon and fort-studio see: 
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"


#### prerequisites ####

# load libraries
library(FortMyrmidon) #R bindings

# set working directory
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(directory) # adjust the script the two files are not in the same folder

# select the files you want to work with
main_file_name <- "vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon"
secondary_file_name <- "vital_fc2_esterhase_c02_feeding_DS_AntsCreated.myrmidon"
no_ants <- "vital_fc2_esterhase_c02_feeding_DS_base.myrmidon"

# tracking file you would like to create ants for: 
main_tracking_data <- fmExperimentOpen(main_file_name)
secondary_tracking_data <-  fmExperimentOpen(secondary_file_name)
no_ant_data <-  fmExperimentOpen(no_ants)

# define a new key (metadata) that will contain the new ID

for main and secondary 
main_tracking_data$setMetaDataKey(key = "meta_ID", default_Value = 001)
main_tracking_data$setMetaDataKey(key = "queen", default_Value = FALSE)
main_tracking_data$setMetaDataKey(key = "name", default_Value = "NA")

main_tracking_data$setMetaDataKey(key = "queen", default_Value = FALSE)

then safe


# define/get the starting time of the experiment
# t <- fmTimeCreate(offset = 0) #SET TIME TO 1970 which is per definition way before the actual start of the experiment.
t <- fmTimeCreate(offset = fmQueryGetDataInformations(main_tracking_data)$start)
# other useful syntax for other things when dealing with tracking system time: 
# from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start) # experiment start time
# to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start + 12*3600  ) # start time plus the duration in seconds

# loop to create create/assign each ant of the main set ID and use the tag 



main_tracking_data$ants[[1]]$getValue(key="name", time = t)
main_tracking_data$ants[[2]]$getValue(key="name", time = t)

main_tracking_data$ants[[1]]$setValue(key="name", time = t)

??setValue
main_tracking_data$ants[[1]]$setValue(key="name", value = "vilya", time = t)
main_tracking_data$ants[[1]]$setValue(key="name", value = "vilya", time = )

main_tracking_data$ants[[1]]$getValue(key="name", time = t)

main_tracking_data$ants[[5]]$getValue(key="name", time = t)
main_tracking_data$ants[[6]]$getValue(key="name", time = t)


fmQueryGetDataInformations(main_tracking_data)$start


main_tracking_data$












# save meta data in new myrmidon file
main_tracking_data$save(paste0(directory, substr(main_file_name, 1, nchar(main_file_name)-9), '_meta.myrmidon'))















??setMetaDataKey





??getValue()


secondary_tracking_data$ants[[1]]$
secondary_tracking_data$ants[[1]]$ID <- 189

secondary_tracking_data$ants[[1]]$setValue(069, secondary_tag_statistics[1,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())

??setValue

no_ant_data$ants
str(no_ant_data$ants)
attributes(no_ant_data)

no_ant_data$ants <- secondary_tracking_data$ants







### playing around with lists

x <- list(1, "a", TRUE, 1+4i, list("a", 2))
x  
x[[5]][[2]] <- 3



tracking_data$ants[[1]]$ID
secondary_tracking_data$ants[[1]]$identifications[[1]]$tagValue

secondary_tracking_data$ants[[1]]$identifications[[1]]$targetAntID <- 003
secondary_tracking_data$ants[[1]]$identifications[[1]]$




secondary_tracking_data$ants[[1]]$identifications[[1]]

hello <- secondary_tracking_data$ants
str(hello)
hello[[1]]$identifications[[1]]$targetAntID <- 999







for ( i in 1:nrow(tag_statistics)) {  #####loop over each tag
  if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
}






  
# extract the tag statistics to know how many times each tag was detected:
main_tag_statistics <- fmQueryComputeTagStatistics(main_tracking_data)
secondary_tag_statistics <- fmQueryComputeTagStatistics(secondary_tracking_data)






# define name of output file
output_name <- file.path(paste0('/media/gw20248/gismo_hd2/vital/fc2/', "Mean_ant_length_colonies.txt")) 




##this changes metadata if tag IDs match for SI ants 02/08/22
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (b in SIBB$ants){
    if (a$identifications[[1]]$tagValue == b$identifications[[1]]$tagValue){
      a$setValue("SIBB",TRUE, t)}}
  for (k in SIKVL$ants){
    if (a$identifications[[1]]$tagValue == k$identifications[[1]]$tagValue){
      a$setValue("SIKVL",TRUE, t)}}
  for (s in SISHAM$ants){
    if (a$identifications[[1]]$tagValue == s$identifications[[1]]$tagValue){
      a$setValue("SISH",TRUE, t)}}
}



#for ( i in 1:nrow(number of tags)) { # for each ant in the secondary file
#  find the same tag id in main file, copy the antID assigned to that main-tag-id copy it and overwrite the previous one 

for (i in 1:nrow(secondary_tag_statistics)) 


for ( i in 1:nrow(tag_statistics)) {  #####loop over each tag
  if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
}


str(tag_statistics)
names(tag_statistics)
tracking_data$ant

ants

#### stuff from Luke ####
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)

dir_data <- "/media/ll16598/Seagate Desktop Drive/"
dir_SI <- "/home/ll16598/Documents/SICC_METADATA/SI_ID/13P/"
dir_chall <- "/home/ll16598/Documents/SICC_METADATA/Pathogen_ID/"

myrmidon_file1 <- paste(dir_data,"13P_IV31_110522_HDN/13P_IV31_110522_HDN_manual_orient_test.myrmidon",sep='')
Metadata_exp <- fmExperimentOpen(myrmidon_file1)

#FORPATHOGEN CHALLENGED COLONIES
CHALL_KVL_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
CHALL_KVL <- fmExperimentOpen(CHALL_KVL_)
CHALL_BB_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
CHALL_BB <- fmExperimentOpen(CHALL_BB_)

#PATHOGEN CHALLENGED IDS
BB_tag_statistics <- fmQueryComputeTagStatistics(CHALL_BB) #do for myrmidon file also
BB_IDs <- BB_tag_statistics$tagDecimalValue
for ( i in 1:nrow(BB_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- BB$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- CHALL_BB$addIdentification(a$ID,BB_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
BB_CHALL_ants <- CHALL_BB$ants #BEAUVERIA IDS

KVL_tag_statistics <- fmQueryComputeTagStatistics(KVL_CHALL) #do for myrmidon file also
KVL_IDs <- KVL_tag_statistics$tagDecimalValue
for ( i in 1:nrow(KVL_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- KVL_CHALL$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- KVL_CHALL$addIdentification(a$ID,KVL_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
KVL_CHALL_ants <- KVL_CHALL$ants

#FOR BEAD CHALLENGED COLONIES
BLUE_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
BLUE <- fmExperimentOpen(BLUE_)
YELLOW_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
YELLOW <- fmExperimentOpen(YELLOW_)

#BEAD CHALLENGE IDS
BLUE_tag_statistics <- fmQueryComputeTagStatistics(BLUE) #do for myrmidon file also
BLUE_IDs <- BLUE_tag_statistics$tagDecimalValue
for ( i in 1:nrow(BLUE_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- BLUE$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- BLUE$addIdentification(a$ID,BLUE_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
BLUE_ants <- BLUE$ants #BLUE ANT IDS

YELLOW_tag_statistics <- fmQueryComputeTagStatistics(YELLOW) #do for myrmidon file also
YELLOW_IDs <- YELLOW_tag_statistics$tagDecimalValue
for ( i in 1:nrow(YELLOW_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- YELLOW$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- YELLOW$addIdentification(a$ID,YELLOW_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
YELLOW_ants <- YELLOW$ants

#READ SI FILES
SIBB_ <- paste(dir_SI, "13P_SIBB.myrmidon",sep='') #Beauveria SI IDs
SIBB <- fmExperimentOpen(SIBB_)
SIKVL_ <- paste(dir_SI, "13P_SIKVL.myrmidon",sep='') #Metarhizium SI IDs
SIKVL <- fmExperimentOpen(SIKVL_)
SISH_ <- paste(dir_SI, "13P_SISH.myrmidon",sep='') #Sham SI IDs
SISH <- fmExperimentOpen(SISH_)



#CREATE SIBB ANTS
SIBB_tag_statistics <- fmQueryComputeTagStatistics(SIBB) #do for myrmidon file also
SIBB_IDs <- SIBB_tag_statistics$tagDecimalValue
for ( i in 1:nrow(SIBB_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- SIBB$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- SIBB$addIdentification(a$ID,SIBB_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
SIBB_ants <- SIBB$ants

#CREATE SIKVL ANTS
SIKVL_tag_statistics <- fmQueryComputeTagStatistics(SIKVL) #do for myrmidon file also
SIKVL_IDs <- SIKVL_tag_statistics$tagDecimalValue
for ( i in 1:nrow(SIKVL_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- SIKVL$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- SIKVL$addIdentification(a$ID,SIKVL_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
SIKVL_ants <- SIKVL$ants

#CREATE SISH ANTS
SISH_tag_statistics <- fmQueryComputeTagStatistics(SISH) #do for myrmidon file also
SISH_IDs <- SISH_tag_statistics$tagDecimalValue
for ( i in 1:nrow(SISH_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- SISH$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- SISH$addIdentification(a$ID,SISH_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
SISH_ants <- SISH$ants

#CREATE ANTS FOR MAIN EXP FILE
Metadata_ants <- Metadata_exp$ants
t <- fmTimeCreate(offset = 0)#SET TIME TO 1970

##this changes metadata if tag IDs match for SI ants 02/08/22
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (b in SIBB$ants){
    if (a$identifications[[1]]$tagValue == b$identifications[[1]]$tagValue){
      a$setValue("SIBB",TRUE, t)}}
  for (k in SIKVL$ants){
    if (a$identifications[[1]]$tagValue == k$identifications[[1]]$tagValue){
      a$setValue("SIKVL",TRUE, t)}}
  for (s in SISHAM$ants){
    if (a$identifications[[1]]$tagValue == s$identifications[[1]]$tagValue){
      a$setValue("SISH",TRUE, t)}}
}
# changes metadata for challenge ants
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (B in BB_CHALL$ants){
    if (a$identifications[[1]]$tagValue == B$identifications[[1]]$tagValue){
      a$setValue("BB_CHALL",TRUE, t)}}
  for (K in KVL_CHALL$ants){
    if (a$identifications[[1]]$tagValue == K$identifications[[1]]$tagValue){
      a$setValue("KVL_CHALL",TRUE, t)}}
}
#CONTROL BEAD ANT METADATA
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (N in BLUE$ants){
    if (a$identifications[[1]]$tagValue == N$identifications[[1]]$tagValue){
      a$setValue("BLUE",TRUE, t)}}
  for (Y in YELLOW$ants){
    if (a$identifications[[1]]$tagValue == Y$identifications[[1]]$tagValue){
      a$setValue("YELLOW",TRUE, t)}}
}

Metadata_exp$save(paste0(sub("\\..*", "", myrmidon_file1),"_withMetaData.myrmidon")) 









