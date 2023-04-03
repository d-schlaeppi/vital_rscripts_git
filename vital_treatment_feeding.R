#### VITAL Feeding analysis #### 

# Which ants did feed in the treatment of the vital tracking experiment and which ants did not. 

#### prerequisites ####
rm(list=ls())

# load necessary libraries
library(dplyr)

### read in table that contains a collection of information on the experiments and each colony, and create a data frame containing all the useful things for the myrmidon files
setwd("/home/gw20248/Documents/vital_rscripts_git/")
dat <- read.csv("vital_treatment_feeding_annotation.csv", header = TRUE, stringsAsFactors = F)
dat$time_start[422] <- NA
dat$time_stop[422] <- NA
dat$time_start[680] <- NA
dat$time_stop[680] <- NA

colony <- c("c05", "c09", "c12", "c13", "c17", "c21")
time <- c("2022-03-16T09:18:34.196Z", "2022-03-23T08:46:06.150Z", "2022-03-26T09:47:19:881Z", "2022-03-30T09:13:39.464Z", "2022-04-06T09:14:40.527Z", "2022-04-13T09:36:14:331Z")
time_correction_df <- data.frame(colony, time)


# step 1 adjust the time format and calculate the duration of each feeding event. 
for (i in 1:nrow(dat)) {
  if (is.na(dat$time_start[i])) {
    dat$time_diff[i] <- 0
    next
  }
  if (nchar(dat$time_start[i]) == 30) {
    #time1 <- format(as.POSIXct(df$time_start[1], format = "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d %H:%M:%S") # get starting time of feeding 
    #time2 <- format(as.POSIXct(df$time_stop[1], format = "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d %H:%M:%S") # get stop time 
    #time_diff_sec <- difftime(time1, time2, units="secs") #calculate time difference in seconds to get the duration of the feeding
    dat$time_diff[i] <- abs(as.numeric(difftime(as.POSIXct(df$time_start[i], format = "%Y-%m-%dT%H:%M:%S"), as.POSIXct(df$time_stop[i], format = "%Y-%m-%dT%H:%M:%S"), units="secs")))
  } else {
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
  group_by(colony, position, focal_AntID) %>%
  summarise(
    total_duration = sum(ifelse(is.na(time_diff), 0, time_diff)))
summed_dat

# step 3 for each colony get the proportion of ants above a certain treshold of feeding duration for both sides

# set the threshold for total duration
threshold <- 30

# calculate the proportion of ants that have total_duration >= threshold for each colony and position
prop_dat <- summed_dat %>%
  group_by(colony, position) %>%
  summarise(prop = mean(total_duration >= threshold))

# select the colonies that meet the criteria
selected_colonies <- prop_dat %>%
  group_by(colony) %>%
  filter(all(prop >= 0.75)) %>%
  filter(abs(diff(prop)) <= 0.25) %>%
  pull(colony)





# step 4 get an idea of how it looks overall

# step 5 come up with some criteria to select colonies based on poportion of feeders and even ness between the two positions

# step 6 correct the data by going over the comments in the raw data and exlcude/include critical cases


# Create sample dataframes
dat <- data.frame(colony = c("A", "B", "C", "D", "E"),
                  time_start = c("01:38:30.000", "03:15:20.000", "02:12:10.000", "04:00:00.000", "01:30:40.000"),
                  stringsAsFactors = FALSE)

time_correction_df <- data.frame(colony = c("A", "B", "C", "D", "E"),
                                 time_correction = c("2022-03-09 10:10:40", "2022-03-09 10:15:20", "2022-03-09 10:20:30", "2022-03-09 10:25:00", "2022-03-09 10:30:10"),
                                 stringsAsFactors = FALSE)

# Merge the two dataframes based on the "colony" variable
merged_df <- merge(dat, time_correction_df, by = "colony")

# Convert "time_start" to POSIXct format
merged_df$time_start <- as.POSIXct(merged_df$time_start, format = "%H:%M:%S.%OS")

# Add time correction to "time_start"
merged_df$time_start <- merged_df$time_start + as.POSIXct(merged_df$time_correction)

# Remove "time_correction" variable from the merged dataframe
merged_df$time_correction <- NULL

# View the corrected dataframe
merged_df



