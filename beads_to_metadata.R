rm(list = ls())

### ### ### ### ### ### ### ### ### ### ### 
### 20_Metadata update     ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### 

#### Read Me ####
# Little script specific for Daniel to update the metadata (individuals) to add the information on beads, space use and the facet net definition of foragers and nurses 

#### Prerequisites ####
library(tcltk)
UPDATE_METADATA <- TRUE

# set directories
setwd(tk_choose.dir(default = "~/", caption = "Select Working Directory")) # direct it to where you have config_user_and_hd.R (typically your scriptfolder)
source("config_user_and_hd.R") # contains getUserOptions() that defines usr and hd and the clean() function
DATADIR <- paste("/media", usr, hd, "vital/fc2", sep="/")
SCRIPTDIR <- paste("/home",usr,"Documents/vital_rscripts_git/source_scripts", sep="/")
setwd(DATADIR)


#load and prep data
meta_data <- read.table("Metadata_vital_2023-07-04.txt", header = T, stringsAsFactors = F, sep = ",")
bead_data <- read.csv("vital_bead_data.csv", header = TRUE)
# sort(as.numeric(gsub("\\D", "", unique(bead_data$flowjo_plate), perl = TRUE))) # just to check that all 25 flowcytometry plates are included
names(bead_data)
names(bead_data)[names(bead_data) == "TagID"] <- "tagID"
head(bead_data)

# Correct the bead data numbers to correct for the blanks: 
mean_blue_blanks <- round(mean(bead_data$blue_count_NB[bead_data$flowjo_sampletype == "blank"], na.rm = TRUE))
mean_yellow_blanks <- round(mean(bead_data$yellow_count_YB[bead_data$flowjo_sampletype == "blank"], na.rm = TRUE))
bead_data$yellow_beads_cor <- ifelse(bead_data$flowjo_sampletype == "sample", bead_data$yellow_count_YB - mean_yellow_blanks, bead_data$yellow_count_YB)
bead_data$blue_beads_cor <- ifelse(bead_data$flowjo_sampletype == "sample", bead_data$blue_count_NB - mean_blue_blanks, bead_data$blue_count_NB)

file_list <- list.files(paste(DATADIR, "vital_experiment/main_experiment/original_data/tag_files", sep="/"))

#### Loop ####
# Run the loop to updata meta data with bead data and facet net caste definition
if(UPDATE_METADATA) {  
 meta_data$count <- "NA"
 meta_data$last_det <- "NA"
 meta_data$rot <- "NA"
 meta_data$AntTask_facet <- "NA"
 for (colony in file_list) { # colony <- file_list[1]
   colony_id <- substr(colony, 1, 3)
   print(colony_id)
   facet_data <- read.table(paste(DATADIR, "vital_experiment/main_experiment/original_data/tag_files", colony, sep="/"), header = TRUE)
   names(facet_data)[names(facet_data) == "group"] <- "AntTask_facet"
   facet_data$colony_id <- colony_id
   
   for (i in 1:nrow(facet_data)) { # i=2 # Add data from facet net script. 
     matching_row <- which(meta_data$colony_id == facet_data$colony_id[i] & meta_data$antID == facet_data$tag[i])
     if (length(matching_row) > 0) {       # Update the values in meta_data
       meta_data$count[matching_row] <- facet_data$count[i]
       meta_data$last_det[matching_row] <- facet_data$last_det[i]
       meta_data$rot[matching_row] <- facet_data$rot[i]
       meta_data$AntTask_facet[matching_row] <- facet_data$AntTask_facet[i]
     }
   }
 }
 merged_data <- merge(meta_data, bead_data, by = c("colony_id", "tagID"), all.x = TRUE)
 names(merged_data)[names(merged_data) == "comment.x"] <- "comments_meta_data"
 names(merged_data)[names(merged_data) == "comment.y"] <- "comments_flowcytometry"
 new_filename<- file.path(DATADIR, paste0("individual_metadata_vital_updated.csv")) # define an output file path
 write.table(merged_data, file=new_filename, append=F,col.names=T,row.names=F,quote=T,sep=",")
}



# for (colony in file_list) {# colony <- file_list[1]
#   colony_id <- substr(colony, 1,3)
#   print(colony_id)
#   facet_data <- read.table(paste(DATADIR, "vital_experiment/main_experiment/original_data/tag_files", colony, sep="/"), header = T)
#   names(facet_data)[names(facet_data) == "group"] <- "AntTask_facet"
#   facet_data$colony_id <- colony_id
#   merged_data <- merge(meta_data, facet_data, by.x = c("colony_id", "antID"), by.y = c("colony_id", "tag"))
#   # Add columns to meta_data
# 
#   update_rows <- match(paste(meta_data$colony_id, meta_data$antID), paste(merged_data$colony_id, merged_data$antID))
#   
#   meta_data$count[!is.na(update_rows)] <- merged_data$count[update_rows[!is.na(update_rows)]]
#   meta_data$last_det[!is.na(update_rows)] <- merged_data$last_det[update_rows[!is.na(update_rows)]]
#   meta_data$rot[!is.na(update_rows)] <- merged_data$rot[update_rows[!is.na(update_rows)]]
#   meta_data$AntTask_facet[!is.na(update_rows)] <- merged_data$AntTask_facet[update_rows[!is.na(update_rows)]]
# }