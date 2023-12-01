rm(list = ls())


getwd()
#directory

directory <- "/Volumes/DISK_B/vital/fc2" # directory on mac
directory <- "/media/ael/DISK_B/vital/fc2" #ael laptop
setwd(directory)

bead_data <- read.table("vital_bead_data_temp.txt", header = TRUE, sep = "\t")

head(bead_data)

unique(bead_data$flowjo_plate)
unique(subset(bead_data, flowjo_sampletype == "sample")$colony_id)
names(bead_data)
names(bead_data)[names(bead_data) == "TagID"] <- "tagID"

meta_data <- read.table("Metadata_vital_2023-07-04.txt", header = T, stringsAsFactors = F, sep = ",")
meta_data$status_ant <- NA
meta_data$status_ant <-ifelse(meta_data$IsTreated==TRUE,"treated","untreated")

merged_data <- merge(meta_data, bead_data, by = c("colony_id", "tagID"), all.x = TRUE)
names(merged_data)

new_filename<- file.path(directory, paste0("individual_metadata_vital.txt")) # define an output file path
write.table(merged_data, file=new_filename, append=F,col.names=T,row.names=F,quote=T,sep=",")
 
# maybe include the bead correction (subtracting blanks from actual samples already at this point?)