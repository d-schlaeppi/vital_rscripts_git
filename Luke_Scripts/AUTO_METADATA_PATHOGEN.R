library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)
rm(list=ls())
gc()
dir_data <- "/media/ll16598/SeagateDesktopDrive/SICC_DATA/"
dir_SI <- "/media/ll16598/One Touch/SICC_METADATA/"
dir_chall <- "/media/ll16598/One Touch/SICC_METADATA/CHALL_ID/"

myrmidon_file1 <- paste("/media/ll16598/SeagateDesktopDrive/SICC_DATA/2P_IV51_080222_POL/LLEXPP2P_IV51_080222_KEY_DEATHS_AutoOriented.myrmidon",sep='')
Metadata_exp <- fmExperimentOpen(myrmidon_file1)

#FORPATHOGEN CHALLENGED COLONIES
CHALL_KVL_ <- paste(dir_chall, "2P_KVL.myrmidon",sep='') #Sham SI IDs
CHALL_KVL <- fmExperimentOpen(CHALL_KVL_)
CHALL_BB_ <- paste(dir_chall, "2P_BB.myrmidon",sep='') #Sham SI IDs
CHALL_BB <- fmExperimentOpen(CHALL_BB_)

#READ SI FILES
SIBB_ <- paste(dir_SI, "2P_SIBB.myrmidon",sep='') #Beauveria SI IDs
SIBB <- fmExperimentOpen(SIBB_)
SIKVL_ <- paste(dir_SI, "2P_SIKVL.myrmidon",sep='') #Metarhizium SI IDs
SIKVL <- fmExperimentOpen(SIKVL_)
SISH_ <- paste(dir_SI, "2P_SISH.myrmidon",sep='') #Sham SI IDs
SISH <- fmExperimentOpen(SISH_)



#PATHOGEN CHALLENGED IDS
BB_tag_statistics <- fmQueryComputeTagStatistics(CHALL_BB) #do for myrmidon file also
BB_IDs <- BB_tag_statistics$tagDecimalValue

for ( i in 1:nrow(BB_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- CHALL_BB$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- CHALL_BB$addIdentification(a$ID,BB_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}

###GET LIST
BB_CHALL_ants <- CHALL_BB$ants #BEAUVERIA IDS
BB_CHALL_LIST <- NULL
for (ant in BB_CHALL_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  BB_CHALL_LIST <- rbind(BB_CHALL_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(BB_CHALL_LIST)
write.table(BB_CHALL_LIST, file = paste0(dir_chall, "/BB_CHALL_list",basename(CHALL_BB_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)

###

KVL_tag_statistics <- fmQueryComputeTagStatistics(CHALL_KVL) #do for myrmidon file also
KVL_IDs <- KVL_tag_statistics$tagDecimalValue
for ( i in 1:nrow(KVL_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- CHALL_KVL$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- CHALL_KVL$addIdentification(a$ID,KVL_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
KVL_CHALL_ants <- CHALL_KVL$ants

###GET LIST
KVL_CHALL_LIST <- NULL
for (ant in KVL_CHALL_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  KVL_CHALL_LIST <- rbind(KVL_CHALL_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(KVL_CHALL_LIST)
write.table(KVL_CHALL_LIST, file = paste0(dir_chall, "/KVL_CHALL_list",basename(CHALL_KVL_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)


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
SIBB_LIST <- NULL
for (ant in SIBB_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  SIBB_LIST <- rbind(SIBB_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(SIBB_LIST)
write.table(SIBB_LIST, file = paste0(dir_SI, "/SIBB_list",basename(SIBB_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)
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
SIKVL_LIST <- NULL
for (ant in SIKVL_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  SIKVL_LIST <- rbind(SIKVL_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(SIKVL_LIST)
write.table(SIKVL_LIST, file = paste0(dir_SI, "/SIKVL_list",basename(SIKVL_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)

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
SISH_LIST <- NULL
for (ant in SISH_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  SISH_LIST <- rbind(SISH_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(SISH_LIST)

write.table(SISH_LIST, file = paste0(dir_SI, "/SISH_list",basename(SISH_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)
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
  for (s in SISH$ants){
    if (a$identifications[[1]]$tagValue == s$identifications[[1]]$tagValue){
      a$setValue("SISH",TRUE, t)}}
}
# changes metadata for challenge ants
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (B in CHALL_BB$ants){
    if (a$identifications[[1]]$tagValue == B$identifications[[1]]$tagValue){
      a$setValue("CHALL_BB",TRUE, t)}}
  for (K in CHALL_KVL$ants){
    if (a$identifications[[1]]$tagValue == K$identifications[[1]]$tagValue){
      a$setValue("CHALL_KVL",TRUE, t)}}
}

Metadata_exp$save(paste0(sub("\\..*", "", myrmidon_file1),"_withMetaData.myrmidon"))  
