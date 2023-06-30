### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Vital feeder meta data ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####

# Quick script to read in some of the essential data from the vital tracking experiment and safe it in a data frame to be accessed in other scripts
# you will get a data frame containing a line with information for each ant that was part of the treatment in the vital feeding experiment (meaning it was sampled and then got access to one of two food sources)

# to load the data frame in another script just run the following line 
# source(paste0(directory,"/vital_metadata_feeders.R"))
#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

# fetching the feeding data

source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # get the colony level data colony_metadata df

dat <- read.csv(paste0(SCRIPTDIR,"/vital_treatment_feeding_annotation_2.csv"), header = TRUE, stringsAsFactors = F) #read data
head(dat)

# create a data frame with the missing times from the colonies that encountered an error in fort.myrmidon
colony <- c("c05", "c09", "c12", "c13", "c17", "c21")
time <- c("2022-03-16T09:18:34.196Z", "2022-03-23T08:46:06.150Z", "2022-03-26T09:47:19:881Z", "2022-03-30T09:13:39.464Z", "2022-04-06T09:14:40.527Z", "2022-04-13T09:36:14:331Z")
time_correction_df <- data.frame(colony, time)

# adjust the time format and calculate the duration of each feeding event 
for (i in 1:nrow(dat)) {
  if (is.na(dat$time_start[i])) { 
    dat$time_diff[i] <- 0 # for the ants that had NA set time feeding as 0 then move on.
    next
  }
  if (nchar(dat$time_start[i]) == 30) { # if time is in the format of fort calculate the duration and save it 
    dat$time_diff[i] <- abs(as.numeric(difftime(as.POSIXct(dat$time_start[i], format = "%Y-%m-%dT%H:%M:%S"), as.POSIXct(dat$time_stop[i], format = "%Y-%m-%dT%H:%M:%S"), units="secs")))
  } else { # if time is in the format of the mpv player reformat it so a format similar to the fort format then calculate the feeding duration 
    colony <- dat$colony[i]
    correction_factor <- time_correction_df$time[which(time_correction_df$colony == colony)]
    correction_time <- as.POSIXct(correction_factor, format = "%Y-%m-%dT%H:%M:%S")
    start_offset <- as.difftime(dat$time_start[i], format = "%H:%M:%OS")
    stop_offset  <- as.difftime(dat$time_stop[i], format = "%H:%M:%OS")
    dat$time_start[i] <- as.character(as.POSIXct(correction_time + start_offset, origin = "1970-01-01"))
    dat$time_stop[i] <-  as.character(as.POSIXct(correction_time + stop_offset , origin = "1970-01-01"))
    dat$time_diff[i] <- abs(as.numeric(time_diff_sec <- difftime(dat$time_start[i], dat$time_stop[i], units="secs")))
  }
}

# Group the data by ant ID and colony and calculate the total duration of feeding events per ant 
summed_dat <- dat %>%
  group_by(colony, position, focal_AntID, tagID) %>%
  summarise(
    total_duration = sum(ifelse(is.na(time_diff), 0, time_diff)), 
    excluder = mean(exclude))
summed_dat
summed_dat$focal_AntID <- as.character(summed_dat$focal_AntID)
summed_dat$tagID <- paste0("0x", summed_dat$tagID)
summed_dat$tagIDdecimal <- as.numeric(summed_dat$tagID)

# create the dataframe containing all the ants that count as treated (isTreated == TRUE)
metadata_feeders <- NULL

for (i in 1:nrow(summed_dat)) {
  colony_id         <- as.character(summed_dat[i, "colony"])
  position          <- summed_dat[i, "position"]
  AntID             <- summed_dat[i, "focal_AntID"]
  tagID             <- summed_dat[i, "tagID"]
  tagIDdecimal      <- summed_dat[i, "tagIDdecimal"]
  feeding_duration  <- summed_dat[i, "total_duration"]
  treatment         <- colony_metadata$treatment[colony_metadata$colony_nr == colony_id]
  excluder           <- summed_dat[i, "excluder"]
  if (treatment == "cc"){
    food  <- "control"
    ifelse(position == "p1", 
           beads <- "yellow", beads <- "blue")
  } else {
    if (position == "p1" ){
      food  <- "virus"
      ifelse(treatment == "vy", 
             beads <- "yellow", beads <- "blue")
    } else {
      food <- "control"
      ifelse(treatment == "vy", 
             beads <- "blue", beads <- "yellow")
    }
  }
  metadata_feeders <-  rbind(metadata_feeders, 
                             data.frame(colony_id,
                                        AntID,
                                        tagID,
                                        tagIDdecimal,
                                        position, 
                                        feeding_duration,
                                        treatment,
                                        food,
                                        beads,
                                        excluder, 
                                        stringsAsFactors = F ))
}





