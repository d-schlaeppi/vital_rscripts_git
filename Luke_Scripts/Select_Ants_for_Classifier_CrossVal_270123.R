
###########################################################################
## RANDOMLY SELECT NON-FOCAL INDIVIDUALS FOR CLASSIFIER CROSS-VALIDATION ##
###########################################################################

## The selected ants should be followed for: 
# 30 mins after exposure (between 2 and 32 mins post-isolation in 15 mins blocks as done with Vasudha)
# only for grooming

DATADIR <- "/media/ll16598/SeagateDesktopDrive/SICC_DATA/"

Metadata_Exp1 <- read.csv(paste(DATADIR,"Metadata_Exp1_2022.txt.csv",sep=""),header=T, sep=",")
#stringsAsFactors = F,
#select only non-exposed individuals
Metadata_Exp1_BB <- Metadata_Exp1[which(Metadata_Exp1$CHALL_BB==TRUE),]
Metadata_Exp1_BB <- Metadata_Exp1_BB[which(Metadata_Exp1_BB$dead==FALSE),]
Metadata_Exp1_KVL <- Metadata_Exp1[which(Metadata_Exp1$CHALL_KVL==TRUE),]
Metadata_Exp1_KVL <- Metadata_Exp1_KVL[which(Metadata_Exp1_KVL$dead==FALSE),]
Metadata_Exp1_BLUE <- Metadata_Exp1[which(Metadata_Exp1$CHALL_BLUE==TRUE),]
Metadata_Exp1_BLUE <- Metadata_Exp1_BLUE[which(Metadata_Exp1_BLUE$dead==FALSE),]
Metadata_Exp1_YELLOW <- Metadata_Exp1[which(Metadata_Exp1$CHALL_YELLOW==TRUE),]
Metadata_Exp1_YELLOW <- Metadata_Exp1_YELLOW[which(Metadata_Exp1_YELLOW$dead==FALSE),]
#exclude tag falls
Metadata_Exp1_BB <- Metadata_Exp1_BB[!grepl("TagFall", Metadata_Exp1$Comment),]

# RANDOMLY 1 ant per colony
# randomly choose only one row in each Replicate
Metadata_Exp1_BB$Chosen <- 0
Metadata_Exp1_BB[-tapply(-seq_along(Metadata_Exp1_BB$REP_treat),Metadata_Exp1_BB$REP_treat, sample, size=1),]$Chosen <- 1
Metadata_BB_Chosen <- Metadata_Exp1_BB[which(Metadata_Exp1_BB$Chosen==1),]
Metadata_Exp1_KVL$Chosen <- 0
Metadata_Exp1_KVL[-tapply(-seq_along(Metadata_Exp1_KVL$REP_treat),Metadata_Exp1_KVL$REP_treat, sample, size=1),]$Chosen <- 1
Metadata_KVL_Chosen <- Metadata_Exp1_KVL[which(Metadata_Exp1_KVL$Chosen==1),]
Metadata_Exp1_BLUE$Chosen <- 0
Metadata_Exp1_BLUE[-tapply(-seq_along(Metadata_Exp1_BLUE$REP_treat),Metadata_Exp1_BLUE$REP_treat, sample, size=1),]$Chosen <- 1
Metadata_BLUE_Chosen <- Metadata_Exp1_BLUE[which(Metadata_Exp1_BLUE$Chosen==1),]
Metadata_Exp1_YELLOW$Chosen <- 0
Metadata_Exp1_YELLOW[-tapply(-seq_along(Metadata_Exp1_YELLOW$REP_treat),Metadata_Exp1_YELLOW$REP_treat, sample, size=1),]$Chosen <- 1
Metadata_YELLOW_Chosen <- Metadata_Exp1_YELLOW[which(Metadata_Exp1_YELLOW$Chosen==1),]

# exclude 6S
Metadata_YELLOW_Chosen <- Metadata_YELLOW_Chosen[!Metadata_YELLOW_Chosen$REP_treat %in% c("EXP_6S"),]
Metadata_BLUE_Chosen <- Metadata_BLUE_Chosen[!Metadata_BLUE_Chosen$REP_treat %in% c("EXP_6S"),]

write.table(Metadata_BB_Chosen,file=paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_BB",Sys.Date(),".txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
write.table(Metadata_KVL_Chosen,file=paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_KVL",Sys.Date(),".txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
write.table(Metadata_BLUE_Chosen,file=paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_BLUE",Sys.Date(),".txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
write.table(Metadata_YELLOW_Chosen,file=paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_YELLOW",Sys.Date(),".txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")


