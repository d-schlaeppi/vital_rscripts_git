#### VITAL Feeding analysis #### 

# Which ants did feed in the treatment of the vital tracking experiment and which ants did not. 

#### prerequisites ####
rm(list=ls())

# load necessary libraries
library(dplyr)
library(ggplot2)
library(ggdark)
library(plotly) # 3d graph
library(viridis) 
library(pscl) # zero inflated poisson distribution model
library(glmmTMB) # zero inflated poisson distribution model with random factors

### read in table that contains a collection of information on the experiments and each colony, and create a data frame containing all the useful things for the myrmidon files
setwd("/home/gw20248/Documents/vital_rscripts_git/") #Uni
#or 
setwd("/Users/gismo/Documents/GitHub/vital_rscripts_git/") #home mac

#dat <- read.csv("vital_treatment_feeding_annotation.csv", header = TRUE, stringsAsFactors = F)
dat <- read.csv("vital_treatment_feeding_annotation_2.csv", header = TRUE, stringsAsFactors = F)


# create a data frame with the missing times from the colonies with an error in fort
colony <- c("c05", "c09", "c12", "c13", "c17", "c21")
time <- c("2022-03-16T09:18:34.196Z", "2022-03-23T08:46:06.150Z", "2022-03-26T09:47:19:881Z", "2022-03-30T09:13:39.464Z", "2022-04-06T09:14:40.527Z", "2022-04-13T09:36:14:331Z")
time_correction_df <- data.frame(colony, time)

#### step 1 adjust the time format and calculate the duration of each feeding event ####

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

#### Group the data by ant ID and colony and calculate the total duration of feeding events per ant ####
summed_dat <- dat %>%
  group_by(colony, position, focal_AntID) %>%
  summarise(
    total_duration = sum(ifelse(is.na(time_diff), 0, time_diff)), 
    excluder = mean(exclude))
summed_dat
summed_dat$focal_AntID <- as.character(summed_dat$focal_AntID)



#### For each colony get the proportion of ants above a certain threshold of feeding duration, proportion of feeders and balance between position 1 & 2 ####
# old version without excluder
# manually excluding specific individuals that were not viable upon return after treatment and will thus be excluded
threshold_feeding <- 0 # threshold to define how long an ants were feeding to be included as feeder (the number of seconds an ant spends feeding to be classified as an actual feeder)
threshold_proportion <- 1 # define a threshold for the proportion of ants of that need to be classified as feeder for colonies to be classified as good enough for subsequent analyses (number between 0 and 1 with 1= 100% of ants were feeding and 0 = 0% of ants were feeding)
threshold_balance <- 0 # define a threshold for the balance between the proportion of feeders in position 1 and 2 for colonies to be classified as good enough for subsequent analyses (balance specified as the difference in the proportion of ants for the two feeding position 0 --> the proportion is the same for both colonies. 1 --> in one position 100% of ants were feeding while in the other it was 0%
# these are just examples and can be played around with 
feeders <- summed_dat %>%
  group_by(colony) %>%
  summarize(
    feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1"),
    feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2"),
    feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding & position == "p1"), "/", 0.5*n()), 
    feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding & position == "p2"), "/", 0.5*n()), 
    prop_feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1") / (0.5*n()), 
    prop_feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2") / (0.5*n())
  )
# Filter colonies based on criteria
selected_colonies <- feeders %>%
  filter(prop_feeders_p1 >= threshold_proportion & prop_feeders_p2 >= threshold_proportion & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance) %>%
  pull(colony)
# Print list of selected colonies
selected_colonies
length(selected_colonies)


#### repeat the above but exclude the ants that need to excluded because they are not very viable after treatment (mostly because they managed to drown themselves)
threshold_feeding <- 0 # threshold to define how long an ants were feeding to be included as feeder (the number of seconds an ant spends feeding to be classified as an actual feeder)
threshold_proportion <- 1 # define a threshold for the proportion of ants of that need to be classified as feeder for colonies to be classified as good enough for subsequent analyses (number between 0 and 1 with 1= 100% of ants were feeding and 0 = 0% of ants were feeding)
threshold_balance <- 0 # define a threshold for the balance between the proportion of feeders in position 1 and 2 for colonies to be classified as good enough for subsequent analyses (balance specified as the difference in the proportion of ants for the two feeding position 0 --> the proportion is the same for both colonies. 1 --> in one position 100% of ants were feeding while in the other it was 0%
# these are just examples and can be played around with 
feeders <- NULL
feeders <- summed_dat %>%
  group_by(colony) %>%
  summarize(
    feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1),
    feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1),
    non_feeders_p1 = sum((position == "p1" & total_duration < threshold_feeding) | (position == "p1" & excluder == 1)),
    non_feeders_p2 = sum((position == "p2" & total_duration < threshold_feeding) | (position == "p2" & excluder == 1)),
    feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1), "/", 0.5*n()),
    feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1), "/", 0.5*n()), 
    prop_feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1) / (0.5*n()),
    prop_feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1) / (0.5*n()) 
  )
# Filter colonies based on criteria
selected_colonies <- feeders %>%
  filter(prop_feeders_p1 >= threshold_proportion & prop_feeders_p2 >= threshold_proportion & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance) %>%
  pull(colony)
# Print list of selected colonies
selected_colonies
length(selected_colonies)





#### selection of thresholds - loop over multiple values to get the best colonies ####

# Define different multiple threshold values to find the right combination
threshold_feeding <- c(0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 75, 90, 120, 180) 
threshold_proportion <- c(1, 0.9, 0.8, 0.75, 0.7, 0.66, 0.6, 0.5, 0.4, 0.25, 0)
threshold_balance <- c(0, 0.1, 0.2, 0.25, 0.3, 0.33, 0.4, 0.5, 0.6, 0.75, 1)

threshold_feeding <- c(0, 15, 30, 45, 60) 
threshold_proportion <- c(1,  0.85, 0.75, 0.66, 0.5)
threshold_balance <- c(0, 0.15, 0.25, 0.33, 0.5)

# create  empty list to store selected colonies for each combination of thresholds
selected_colonies_list <- list()
# create empty data frame to store excluded colonies
excluded_colonies <- data.frame(t_feeding = numeric(),
                                t_proportion = numeric(),
                                t_balance = numeric(),
                                included_colonies = numeric(),
                                excluded_colonies = numeric())

#it takes a while to run this loop. The list with the selected colonies and the dataframe with excluded colonies have been saved and can be called with the following code
#selected_colonies_list <- readRDS("selected_colonies_list.rds")
#excluded_colonies <- readRDS("excluded_colonies.rds")

# Iterate over threshold values
for (i in seq_along(threshold_feeding)) {
  for (j in seq_along(threshold_proportion)) {
    for (k in seq_along(threshold_balance)) {
      # recalculate the feeder table
      feeders <- NULL
      feeders <- summed_dat %>%
        group_by(colony) %>%
        summarize(
          feeders_p1 = sum(total_duration >= threshold_feeding[i] & position == "p1" & excluder != 1),
          feeders_p2 = sum(total_duration >= threshold_feeding[i] & position == "p2" & excluder != 1),
          non_feeders_p1 = sum((position == "p1" & total_duration < threshold_feeding[i]) | (position == "p1" & excluder == 1)),
          non_feeders_p2 = sum((position == "p2" & total_duration < threshold_feeding[i]) | (position == "p2" & excluder == 1)),
          feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding[i] & position == "p1" & excluder != 1), "/", 0.5*n()),
          feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding[i] & position == "p2" & excluder != 1), "/", 0.5*n()), 
          prop_feeders_p1 = sum(total_duration >= threshold_feeding[i] & position == "p1" & excluder != 1) / (0.5*n()),
          prop_feeders_p2 = sum(total_duration >= threshold_feeding[i] & position == "p2" & excluder != 1) / (0.5*n()) 
        )
      # Filter colonies based on criteria
      selected_colonies <- feeders %>%
        filter(prop_feeders_p1 >= threshold_proportion[j] & prop_feeders_p2 >= threshold_proportion[j] & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance[k]) %>%
        pull(colony)
      # Add selected colonies to the list
      selected_colonies_list[[paste(threshold_feeding[i], threshold_proportion[j], threshold_balance[k], sep = "_")]] <- selected_colonies
      # Calculate number of excluded colonies
      excluded_colonies_count <- length(unique(feeders$colony)) - length(selected_colonies)
      included_colonies_count <- length(unique(feeders$colony)) - excluded_colonies_count 
      # Add row to excluded_colonies data frame
      excluded_colonies <- rbind(excluded_colonies,
                                 data.frame(t_feeding = threshold_feeding[i],
                                            t_proportion = threshold_proportion[j],
                                            t_balance = threshold_balance[k],
                                            included_colonies = included_colonies_count,
                                            excluded_colonies = excluded_colonies_count))
    }
  }
}

# if needed the latest version of this can be saved so it can be called again without re-running the loop
saveRDS(selected_colonies_list, "selected_colonies_list.rds")
saveRDS(excluded_colonies, "excluded_colonies.rds")
selected_colonies_list[1]




#### visualize the exclusion power of the different thresholds ####
plot_ly(excluded_colonies, x = ~t_feeding, y = ~t_proportion, z = ~t_balance, 
        color = ~excluded_colonies, 
        colors = viridis_pal()(100),
        type = "scatter3d", mode = "markers") %>%
  add_markers(size = 5) %>%
  layout(scene = list(xaxis = list(title = "Threshold feeding (seconds)", nticks = 5),
                      yaxis = list(title = "Threshold proportion", nticks = 5),
                      zaxis = list(title = "Threshold balance", nticks = 5),
                      xaxis_title_font = list(size = 12),
                      yaxis_title_font = list(size = 12),
                      zaxis_title_font = list(size = 12),
                      axis.title = element_text(size = 14),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12),
                      axis.text.z = element_text(size = 12)))






#### which are the best colonies ####
selected_colonies_vec <- unlist(selected_colonies_list) # Combine all selected colonies into a single vector 
colony_counts <- table(selected_colonies_vec) #Count the number of times each colony appears
sorted_counts <- sort(colony_counts, decreasing = TRUE)# Sort the counts in descending order
# Print the sorted counts
print(sorted_counts)
plot(colony_counts)
plot(sorted_counts, las=2)





# filter the rows based on the excluded_colonies value
filtered_excluded <- excluded_colonies %>% 
  filter(included_colonies >= 14)
# print the combinations of thresholds that satisfy the condition
print(filtered_excluded)
filtered_excluded <- excluded_colonies %>% 
  filter(t_feeding == 30  & excluded_colonies <= 16)
print(filtered_excluded)
#choose one of the threshold sets.
# or simply select the 14 best colonies!!! 




#### check the colonies selected with those thresholds or the X best colonies from the treshold selection ####
{ # run the loop with the defined parameters 
threshold_feeding <- 30
threshold_proportion <- 0.5
threshold_balance <- 0.34
feeders <- summed_dat %>%
  group_by(colony) %>%
  summarize(
    feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1),
    feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1),
    non_feeders_p1 = sum((position == "p1" & total_duration < threshold_feeding) | (position == "p1" & excluder == 1)),
    non_feeders_p2 = sum((position == "p2" & total_duration < threshold_feeding) | (position == "p2" & excluder == 1)),
    feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1), "/", 0.5*n()),
    feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1), "/", 0.5*n()), 
    prop_feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1) / (0.5*n()),
    prop_feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1) / (0.5*n()) 
  )
selected_colonies <- feeders %>%
  filter(prop_feeders_p1 >= threshold_proportion & prop_feeders_p2 >= threshold_proportion & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance) %>%
  pull(colony)
selected_colonies
length(selected_colonies)
colonies_to_analyse <- subset(feeders, colony %in% selected_colonies)
colonies_to_analyse
}




# select the best performing colonies. 
x <- 14 #define the numbers of colonies to be selected
sorted_counts
selected_colonies <- names(sorted_counts)[1:14]
feeders <- summed_dat %>%
  group_by(colony) %>%
  summarize(
    feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1"),
    feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2"),
    feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding & position == "p1"), "/", 0.5*n()), 
    feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding & position == "p2"), "/", 0.5*n()), 
    prop_feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1") / (0.5*n()), 
    prop_feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2") / (0.5*n())
  )
colonies_to_analyse <- subset(feeders, colony %in% selected_colonies)
colonies_to_analyse



#### Check if selected colonies are equally balanced between the two treatments and bead colors ####

# import treatments and other data: 
experiment_data <- read.csv("fc2_overview_data.csv", header = TRUE, stringsAsFactors = F)
data_collection <- NULL 
for(i in 1:nrow(experiment_data)) {
  # collect variables
  nr                    <- i
  colony_id             <- experiment_data[i, "colony_id"]
  block                 <- experiment_data[i, "block"]
  colony             <- paste0("c", sprintf("%02d", i))
  treatment             <- experiment_data[i, "treatment"]
  food_position_1       <- experiment_data[i, "food_position_1"]
  food_position_2       <- experiment_data[i, "food_position_2"]
  tracking_system_main  <- experiment_data[i, "tracking_system_main"]
  tracking_system_feeding <- experiment_data[i, "tracking_system_feeding"] 
  # combine variables to a data frame  
  data_collection <-  rbind(data_collection, data.frame(nr, 
                                                        colony_id,
                                                        block, 
                                                        colony,
                                                        treatment,
                                                        food_position_1,
                                                        food_position_2,
                                                        tracking_system_main,
                                                        tracking_system_feeding,
                                                        stringsAsFactors = F))
}


# Merge the two data frames based on "colony" variable
merged_data <- merge(colonies_to_analyse, data_collection, by = "colony")
# check it the selected colonies are balanced across treatments 
table(merged_data$treatment)
# this is not perfect, there is one control too few and one vb too many... 


### plot sorted counts from the colony selection above 
treatments <- data_collection$treatment[match(names(sorted_counts), data_collection$colony)] # Match treatments with colonies
colors <- ifelse(treatments == "cc", "red", ifelse(treatments == "vy", "yellow", "blue")) # Create color vector
barplot(sorted_counts, col=colors, las=2)











#### Analyses of the feeding behavior ####

# import the general information of each colony regarding position and bead colors
df <- read.csv("fc2_overview_data.csv", header = TRUE, stringsAsFactors = F)
data_collection <- NULL 
for(i in 1:nrow(df)) {
  # collect variables
  nr                    <- i
  colony_id             <- df[i, "colony_id"]
  block                 <- df[i, "block"]
  colony_nr             <- paste0("c", sprintf("%02d", i))
  colony                <- paste0("c", sprintf("%02d", i))
  treatment             <- df[i, "treatment"]
  food_position_1       <- df[i, "food_position_1"]
  food_position_2       <- df[i, "food_position_2"]
  tracking_system_main  <- df[i, "tracking_system_main"]
  tracking_system_feeding <- df[i, "tracking_system_feeding"] 
  # combine variables to a data frame  
  data_collection <-  rbind(data_collection, data.frame(nr, 
                                                        colony_id,
                                                        block, 
                                                        colony_nr,
                                                        colony,
                                                        treatment,
                                                        food_position_1,
                                                        food_position_2,
                                                        tracking_system_main,
                                                        tracking_system_feeding,
                                                        stringsAsFactors = F))
}

# create new data frame containing all information needed
dynamic <- NULL
for (i in 1:nrow(summed_dat)) {
  nr                <- i
  colony_id         <- as.character(summed_dat[i, "colony"])
  position          <- summed_dat[i, "position"]
  focal_AntID       <- summed_dat[i, "focal_AntID"]
  feeding_duration  <- summed_dat[i, "total_duration"]
  treatment         <- data_collection$treatment[data_collection$colony_nr == colony_id]
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
  dynamic <-  rbind(dynamic, data.frame(nr,   # combine variables to a data frame  
                                        colony_id,
                                        position, 
                                        focal_AntID,
                                        feeding_duration,
                                        treatment,
                                        food,
                                        beads,
                                        stringsAsFactors = F ))
}






#### Stats on feeding behavior ####

### food control or virus
# is there an effect of food quality (control or virus) on the duration the ants spend feeding
boxplot(dynamic$total_duration ~ dynamic$food)
hist(dynamic$total_duration, breaks = 40)
# data is zero inflated poisson (e.g. counting number of seconds the ants spend feeding similar to e.g. the number of days a patient

# Fit zero-inflated Poisson regression model first test
model <- zeroinfl(total_duration ~ food | 1, data = dynamic, dist = "poisson")
summary(model) # ok but we need random factors


# Fit mixed-effects zero-inflated Poisson regression model with random factors
model <- glmmTMB(total_duration ~ food + (1 | colony_id) + (1 | treatment)+ (1 | beads),
                 ziformula = ~1, data = dynamic, family = "truncated_poisson")
summary(model)
head(dynamic)
# same but with virus fraction only: 
dynamic_virus_only <- subset(dynamic, treatment != "cc")
model <- glmmTMB(total_duration ~ food + (1 | colony_id) + (1 | treatment)+ (1 | beads),
                 ziformula = ~1, data = dynamic_virus_only, family = "truncated_poisson")
summary(model)



# test model assumptions
library(DHARMa)
res = simulateResiduals(model)
plot(res, asFactor = T)
# violates assumption of normally distributed residuals --> find alternative solution

#Check if the number of zeros is similar for both food sources!
dynamic_virus_only %>%
  group_by(food) %>%
  summarize(count = sum(total_duration == 0))



ggplot(dynamic_virus_only, aes(x = food, y = total_duration, fill = food)) +
  geom_boxplot(aes(fill= food),colour ="white", width = 0.3, alpha = 0.6)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.5, color = "white") +
  labs(title = "Feeding Duration by Food Source",
       x = "Food type", y = "Duration [log(seconds+0.5)])", 
       fill = "Food type") +
  dark_theme_classic() +
  theme(legend.position = "none")

# Count the number of observations for each food group
dynamic_virus_only %>%
  group_by(food) %>%
  summarize(count = n())







#### BEADS ####

# Fit mixed-effects zero-inflated Poisson regression model with random factors
model <- glmmTMB(total_duration ~ beads + (1 | colony_id) + (1 | food),
                 ziformula = ~1, data = dynamic, family = "truncated_poisson")
summary(model)
ggplot(dynamic, aes(x = beads, y = total_duration, fill = beads)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5, color = "white") +
  labs(title = "Feeding Duration by bead Color",
       x = "Bead color", y = "Duration (seconds)", 
       fill = "Bead color") +
  dark_theme_classic() +
  theme(legend.position = "none")



# Fit mixed-effects zero-inflated Poisson regression model with random factors
model <- glmmTMB(total_duration ~ beads + food + (1 | colony_id),
                 ziformula = ~1, data = dynamic, family = "truncated_poisson")
summary(model)




# next run quick stats to show that there is no difference in the feeding duration between the two colonies. 
# next go over the manual annotation file and insert a correction variable to exclude ants which then later on died because they drowned themselves in food (either dead or because their behavior completely off)


dynamic$total_duration
ggplot(data = dynamic, aes(x = total_duration)) +
  geom_histogram(binwidth = 50) +
  ggtitle("Histogram of Total Duration")

ggplot(data = dynamic, aes(sample = total_duration)) +
  stat_qq() +
  ggtitle("Q-Q Plot of Total Duration")

library(fitdistrplus)
fit_gamma <- fitdist(dynamic$total_duration, "gamma")
fit_zigamma <- fitdist(dynamic$total_duration, "ziggamma")

# Compare AIC values
AIC(fit_gamma, fit_zigamma)









#### selection of FCII ants to be manually annotated for the trophy classifier ####

# Select ants with total_duration >= 30 and excluder == 0
selected_ants <- summed_dat %>% 
  filter(total_duration >= 30, excluder == 0) %>% 
  group_by(colony, position) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  select(colony, position, antID = focal_AntID)
table(selected_ants$colony)

# Save selected ants to a text file
write.table(selected_ants, file = "vital_ants_to_annotate.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)




