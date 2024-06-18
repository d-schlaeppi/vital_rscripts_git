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
library(entropy)

### read in table that contains a collection of information on the experiments and each colony, and create a data frame containing all the useful things for the myrmidon files
setwd("/home/gw20248/Documents/vital_rscripts_git/") #Uni
#or 
setwd("/Users/gismo/Documents/GitHub/vital_rscripts_git/") #home mac

#dat <- read.csv("vital_treatment_feeding_annotation.csv", header = TRUE, stringsAsFactors = F)
dat <- read.csv("vital_treatment_feeding_annotation_2.csv", header = TRUE, stringsAsFactors = F)
head(dat)
#dat <- read.csv("/Users/gismo/Documents/GitHub/vital_rscripts_git/vital_treatment_feeding_annotation_2.csv", header = TRUE, stringsAsFactors = F)


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
  group_by(colony, position, focal_AntID, tagID) %>%
  summarise(
    total_duration = sum(ifelse(is.na(time_diff), 0, time_diff)), 
    excluder = mean(exclude))
summed_dat
summed_dat$focal_AntID <- as.character(summed_dat$focal_AntID)
summed_dat$tagID <- paste0("0x", summed_dat$tagID)



#### For each colony get the proportion of ants above a certain threshold of feeding duration, proportion of feeders and balance between position 1 & 2 ####

#### just a test with non-exclusive thresholds ####
# manually excluding specific individuals that were not viable upon return after treatment and will thus be excluded
threshold_feeding <- 0 # threshold to define how long an ants were feeding to be included as feeder (the number of seconds an ant spends feeding to be classified as an actual feeder)
threshold_proportion <- 1 # define a threshold for the proportion of ants of that need to be classified as feeder for colonies to be classified as good enough for subsequent analyses (number between 0 and 1 with 1= 100% of ants were feeding and 0 = 0% of ants were feeding)
threshold_balance <- 0 # define a threshold for the balance between the proportion of feeders in position 1 and 2 for colonies to be classified as good enough for subsequent analyses (balance specified as the difference in the proportion of ants for the two feeding position 0 --> the proportion is the same for both colonies. 1 --> in one position 100% of ants were feeding while in the other it was 0%
feeders <- summed_dat %>%
  group_by(colony) %>%
  summarize(
    feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1"),
    feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2"),
    feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding & position == "p1"), "/", sum(position == "p1")),
    feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding & position == "p2"), "/", sum(position == "p2")),
    prop_feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1") / sum(position == "p1"),
    prop_feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2") / sum(position == "p2")
  )
summed_dat
 # Filter colonies based on criteria
selected_colonies <- feeders %>%
  filter(prop_feeders_p1 >= threshold_proportion & prop_feeders_p2 >= threshold_proportion & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance) %>%
  pull(colony)
# Print list of selected colonies
selected_colonies
length(selected_colonies)
#should list all colonies in the experiment because no ant is excluded and all ants are classified as feeders - works

#### same test but with excluder ####
# repeat the above but exclude the ants that need to excluded because they are not very viable after treatment (mostly because they managed to drown themselves)
feeders <- summed_dat %>%
  group_by(colony) %>%
  summarize(
    feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1),
    feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1),
    non_feeders_p1 = sum((position == "p1" & total_duration < threshold_feeding) | (position == "p1" & excluder == 1)),
    non_feeders_p2 = sum((position == "p2" & total_duration < threshold_feeding) | (position == "p2" & excluder == 1)),
    feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1), "/", sum(position == "p1")),
    feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1), "/", sum(position == "p2")),
    prop_feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1) / sum(position == "p1"),
    prop_feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1) / sum(position == "p2")
  )
selected_colonies <- feeders %>%
  filter(prop_feeders_p1 >= threshold_proportion & prop_feeders_p2 >= threshold_proportion & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance) %>%
  pull(colony)
selected_colonies
length(selected_colonies)
# only shows the colonies for which no ant got excluded - works





#### selection of thresholds - loop over multiple values to get the best colonies (with excluder) ####

# Define different multiple threshold values to find the right combination
threshold_feeding <- c(0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 75, 90, 120, 180) 
threshold_proportion <- c(1, 0.9, 0.8, 0.75, 0.7, 0.66, 0.6, 0.5, 0.4, 0.25, 0)
threshold_balance <- c(0, 0.1, 0.2, 0.25, 0.3, 0.33, 0.4, 0.5, 0.6, 0.75, 1)

threshold_feeding <- c(0, 10, 15, 20, 25, 30, 35, 40, 45, 60)
threshold_proportion <- c(1,  0.85, 0.75, 0.66, 0.5, 0)
threshold_balance <- c(0, 0.16, 0.26, 0.34, 0.51, 1)

# create  empty list to store selected colonies for each combination of thresholds
selected_colonies_list <- list()
# create empty data frame to store excluded colonies
excluded_colonies <- data.frame(t_feeding = numeric(),
                                t_proportion = numeric(),
                                t_balance = numeric(),
                                included_colonies = numeric(),
                                excluded_colonies = numeric())



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
          feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding[i] & position == "p1" & excluder != 1), "/", sum(position == "p1")),
          feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding[i] & position == "p2" & excluder != 1), "/", sum(position == "p2")),
          prop_feeders_p1 = sum(total_duration >= threshold_feeding[i] & position == "p1" & excluder != 1) / sum(position == "p1"),
          prop_feeders_p2 = sum(total_duration >= threshold_feeding[i] & position == "p2" & excluder != 1) / sum(position == "p2")
        )
      selected_colonies <- feeders %>%
        filter(prop_feeders_p1 >= threshold_proportion[j] & prop_feeders_p2 >= threshold_proportion[j] & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance[k]) %>%
        pull(colony)
      selected_colonies_list[[paste(threshold_feeding[i], threshold_proportion[j], threshold_balance[k], sep = "_")]] <- selected_colonies # Add selected colonies to the list
      excluded_colonies_count <- length(unique(feeders$colony)) - length(selected_colonies)       # Calculate number of excluded colonies
      included_colonies_count <- length(unique(feeders$colony)) - excluded_colonies_count 
      excluded_colonies <- rbind(excluded_colonies, # Add row to excluded_colonies data frame
                                 data.frame(t_feeding = threshold_feeding[i],
                                            t_proportion = threshold_proportion[j],
                                            t_balance = threshold_balance[k],
                                            included_colonies = included_colonies_count,
                                            excluded_colonies = excluded_colonies_count))
    }
  }
}

# if needed the latest version of this can be saved so it can be called again without re-running the loop
# saveRDS(selected_colonies_list, "selected_colonies_list.rds")
# saveRDS(excluded_colonies, "excluded_colonies.rds")
# it takes a while to run this loop with many cobinaitons of thresholds. IF the list with the selected colonies and the dataframe with excluded colonies have been saved and can be called with the following code
# selected_colonies_list <- readRDS("selected_colonies_list.rds")
# excluded_colonies <- readRDS("excluded_colonies.rds")

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


#### which are the best colonies based on thresholds ####
selected_colonies_vec <- unlist(selected_colonies_list) # Combine all selected colonies into a single vector 
colony_counts <- table(selected_colonies_vec) #Count the number of times each colony appears
sorted_counts <- sort(colony_counts, decreasing = TRUE)# Sort the counts in descending order
# Print the sorted counts
print(sorted_counts)
plot(colony_counts)
plot(sorted_counts, las=2)
selected_colonies_sorted_counts <- sorted_counts[1:14]

# Check which threshold combinations include more than 14 colonies or based on further filters 
filtered_excluded <- excluded_colonies %>% filter(included_colonies >= 14 & included_colonies <=18)
filtered_excluded <- excluded_colonies %>% 
  filter(t_feeding == 30  & included_colonies >= 14)
print(filtered_excluded)
#choose one of the threshold sets or  select the 14 best colonies based on treshold counts




#### check the colonies selected with a reasonable thresholds set ####
{ # run the loop with the defined parameters 
threshold_feeding <- 20
threshold_proportion <- 0.5
threshold_balance <- 0.26
feeders <- summed_dat %>%
  group_by(colony) %>%
  summarize(
    feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1),
    feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1),
    non_feeders_p1 = sum((position == "p1" & total_duration < threshold_feeding) | (position == "p1" & excluder == 1)),
    non_feeders_p2 = sum((position == "p2" & total_duration < threshold_feeding) | (position == "p2" & excluder == 1)),
    feeders_out_of_p1 = paste(sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1), "/", sum(position == "p1")),
    feeders_out_of_p2 = paste(sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1), "/", sum(position == "p2")),
    prop_feeders_p1 = sum(total_duration >= threshold_feeding & position == "p1" & excluder != 1) / sum(position == "p1"),
    prop_feeders_p2 = sum(total_duration >= threshold_feeding & position == "p2" & excluder != 1) / sum(position == "p2")
  )
selected_colonies <- feeders %>%
  filter(prop_feeders_p1 >= threshold_proportion & prop_feeders_p2 >= threshold_proportion & abs(prop_feeders_p1 - prop_feeders_p2) <= threshold_balance) %>%
  pull(colony)
selected_colonies
length(selected_colonies)
colonies_to_analyse <- subset(feeders, colony %in% selected_colonies)
colonies_to_analyse
}








#### Entropy - Alternative method of colony selection ####
# idea: using entropy as a measure of balance between the two feeding positions without using tresholds. 
# calculate a standardized measure of entropy and then plot the difference in entropy between the two colonies on one a

# step1 -  calculate log10 (+1sec) of the feeding. duration (subset to remove individuals to exclude)

# Exclude rows where excluder != 0 & log10 of total_duration (with +1 second)
df_filtered <- summed_dat %>% filter(excluder == 0)
df_filtered$log_duration <- log10(df_filtered$total_duration + 1)
# Calculate the mean of total_duration per colony
mean_duration <- df_filtered %>%
  group_by(colony) %>%
  summarise(mean_duration_log = mean(log_duration))
# Calculate the entropy for position p1
entropy_p1 <- df_filtered %>%
  group_by(colony, position) %>%
  filter(position == "p1") %>%
  summarise(entropy_p1 = entropy(log_duration)/entropy(rep(1, length(log_duration))))
# Calculate the entropy for position p2
entropy_p2 <- df_filtered %>%
  group_by(colony, position) %>%
  filter(position == "p2") %>%
  summarise(entropy_p2 = entropy(log_duration)/entropy(rep(1, length(log_duration))))
# Merge the results into a new data frame
new_df <- left_join(mean_duration, entropy_p1, by = "colony") %>%
  left_join(entropy_p2, by = "colony")
# Calculate the absolute difference between entropy_p1 and entropy_p2
new_df <- new_df %>%
  mutate(entropy_diff = abs(entropy_p1 - entropy_p2))
# Print the new data frame
print(new_df)


# select the 14 best colonies based on entropy differences between positions and mean(duration)
max_entropy_diff <- max(new_df$entropy_diff, na.rm = TRUE)
max_duration_log <- max(new_df$mean_duration_log, na.rm = TRUE)
new_df <- new_df %>%
  mutate(distance = (max_duration_log - mean_duration_log)/max_duration_log +  (1-(max_entropy_diff - entropy_diff)/max_entropy_diff))

# Sort the new_df based on the distance in ascending order & select the top 14 colonies from the sorted_df
sorted_df <- new_df %>% arrange(distance)
top_left_colonies <- head(sorted_df, 14)

# Plotting mean_duration_log against entropy_diff with colony labels and displaying distance values, with selected colonies highlighted
treshold_selected_colonies <- names(selected_colonies_sorted_counts)
ggplot(new_df, aes(x = entropy_diff, y = mean_duration_log, label = colony)) +
  geom_point() +
  geom_text(vjust = -0.5) +
  xlab("Entropy Difference") +
  ylab("Mean Duration (Log)") +
  ggtitle("Mean Duration Log vs. Entropy Difference") +
  theme_minimal() +
  geom_point(data = top_left_colonies, color = "blue", size = 5, shape = 16) +
  geom_point(data = new_df[new_df$colony %in% treshold_selected_colonies, ], color = "yellow", size = 3, shape = 16)

#### GPT approach ####
df_filtered <- summed_dat %>% filter(excluder == 0)
df_filtered$log_duration <- log10(df_filtered$total_duration + 1)
# Calculate the balance between Position 1 and Position 2 for each colony
balance <- df_filtered %>%
  group_by(colony) %>%
  summarise(balance = abs(mean(total_duration[position == "p1"]) - mean(total_duration[position == "p2"])))
# Calculate the number of ants with total_duration = 0 for each colony
num_zeros <- df_filtered %>%
  group_by(colony) %>%
  summarise(num_zeros = sum(total_duration == 0))
# Calculate the overall total_duration for each colony
overall_duration <- df_filtered %>%
  group_by(colony) %>%
  summarise(overall_duration = sum(total_duration))
# Merge the characteristics into a single data frame
characteristics <- left_join(left_join(balance, num_zeros, by = "colony"), overall_duration, by = "colony")
# Assign scores to each colony based on the characteristics
characteristics$score <- with(characteristics, balance / max(balance) + num_zeros / max(num_zeros) + overall_duration / max(overall_duration))
# Sort the colonies based on scores
sorted_colonies <- characteristics %>%
  arrange(desc(score))
# Select the top 14 colonies based on scores
selected_colonies_gpt <- characteristics %>%
  top_n(14, score)
# Print the selected colonies
print(selected_colonies_gpt)


#### Ranking approach based on total duration and position difference in tot_duration ####
# Exclude rows where excluder != 0 & log10 of total_duration (with +1 second)
df_filtered <- summed_dat %>% filter(excluder == 0)
df_filtered$log_duration <- log10(df_filtered$total_duration + 1)
# Calculate the mean of total_duration per colony
mean_duration <- df_filtered %>%
  group_by(colony) %>%
  summarise(mean_duration_log = mean(log_duration))
ranked_duration <- mean_duration %>% 
  mutate(rank_dur = rank(-mean_duration_log))
# Calculate the difference between mean log_duration for position p1 and p2
colony_difference <-df_filtered %>% 
  group_by(colony) %>% 
  summarize(difference = abs(mean(log_duration[position == "p1"]) - mean(log_duration[position == "p2"])))
ranked_difference <- colony_difference %>%
  mutate(rank_dif = rank(difference))
rank_df <- left_join(ranked_duration, ranked_difference, by = "colony")
rank_df <- rank_df %>% 
  mutate(summed_rank = rank_dur + rank_dif)
# Sort the dataframe based on the overall rank in ascending order and create new variable to mark the top
rank_df <- rank_df %>% arrange(rank_dur + rank_dif)
rank_df <- rank_df %>% mutate(top_colony = ifelse(row_number() <= 14, "Top 14 Colonies", "Other Colonies"))

example_threshold_selected_colonies <- colonies_to_analyse$colony # based on a reasonable set of three thresholds
gpt_selected_colonies <- selected_colonies_gpt$colony #really bad and thus not pursued
entropy_selected_colonies <- top_left_colonies$colony
threshold_count_selected_colonies <- names(selected_colonies_sorted_counts) # based on threshold counts
top_ranked_colonies_simple <- rank_df$colony[1:14]


ggplot(rank_df, aes(x = rank_dif, y = rank_dur, label = colony, color = top_colony)) +
  geom_point(size = 9) +
  geom_text(vjust = -1.5) +
  scale_x_reverse() +
  scale_y_reverse() +
  scale_color_viridis_d(option = "A", begin = 0.85, end = 0.2) + 
  xlab("Rank Duration Difference") +
  ylab("Rank Mean Duration (Log)") +
  ggtitle("Ranked Colonies based on Mean Duration Log and Duration Difference") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_point(data = rank_df[rank_df$colony %in% threshold_count_selected_colonies, ], color = "white", size = 7, shape = 16) + 
  geom_point(data = rank_df[rank_df$colony %in% example_threshold_selected_colonies, ], color = "black", size = 5, shape = 16) + 
  geom_point(data = rank_df[rank_df$colony %in% entropy_selected_colonies, ], color = "yellow", size = 3, shape = 16)


#### Ranking with new criteria taking into account the number of ants ####
# Exclude rows where excluder != 0 & log10 of total_duration (with +1 second)
df_filtered <- summed_dat %>% filter(excluder == 0)
df_filtered$log_duration <- log10(df_filtered$total_duration + 1)
# Calculate the mean of total_duration per colony
mean_duration <- df_filtered %>%
  group_by(colony) %>%
  summarise(mean_duration_log = mean(log_duration))
ranked_duration <- mean_duration %>% 
  mutate(rank_dur = rank(-mean_duration_log))
# Calculate the difference between mean log_duration for position p1 and p2
colony_difference <-df_filtered %>% 
  group_by(colony) %>% 
  summarize(difference = abs(mean(log_duration[position == "p1"]) * sum(position == "p1") - mean(log_duration[position == "p2"]) * sum(position == "p2")))
ranked_difference <- colony_difference %>%
  mutate(rank_dif = rank(difference))
rank_df <- left_join(ranked_duration, ranked_difference, by = "colony")
rank_df <- rank_df %>% 
  mutate(summed_rank = rank_dur + rank_dif)
# Sort the dataframe based on the overall rank in ascending order and create new variable to mark the top
rank_df <- rank_df %>% arrange(rank_dur + rank_dif)
rank_df <- rank_df %>% mutate(top_colony = ifelse(row_number() <= 14, "Top 14 Colonies", "Other Colonies"))

top_ranked_colonies_corrected <- rank_df$colony[1:14]

ggplot(rank_df, aes(x = rank_dif, y = rank_dur, label = colony, color = top_colony)) +
  geom_point(size = 9) +
  geom_text(vjust = -1.5) +
  scale_x_reverse() +
  scale_y_reverse() +
  scale_color_viridis_d(option = "A", begin = 0.85, end = 0.2) + 
  xlab("Rank Duration Difference") +
  ylab("Rank Mean Duration (Log)") +
  ggtitle("Ranked Colonies based on Mean Duration Log and Duration Difference") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_point(data = rank_df[rank_df$colony %in% threshold_count_selected_colonies, ], color = "white", size = 7, shape = 16) + 
  geom_point(data = rank_df[rank_df$colony %in% example_threshold_selected_colonies, ], color = "black", size = 5, shape = 16) + 
  geom_point(data = rank_df[rank_df$colony %in% entropy_selected_colonies, ], color = "yellow", size = 3, shape = 16)




# Create a list of vectors
vector <- c(example_threshold_selected_colonies,
  entropy_selected_colonies,
  threshold_count_selected_colonies,
  #gpt_selected_colonies,
  top_ranked_colonies_simple, 
  top_ranked_colonies_corrected
)
x <- table(vector)
sort(x, decreasing = TRUE)

### ### ###
colonies_to_analyse_certain  <- names(sort(x, decreasing = TRUE)[1:12])
#select between based on treatment and bead colors 
target_colonies <- c("c02", "c04", "c11")
# we need to select two of the target colonies so to have the best ones + to make sure we have an even balance between colours and treatments

### ### ###


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
data_selected_colonies  <- data_collection[data_collection$colony %in% colonies_to_analyse_certain, ]
# check it the selected colonies are balanced across treatments 
table(data_selected_colonies$treatment)
# there is one control missing and we can choose one of the two virus ones irrespective of color from the three target colonies
data_addon_colonies  <- data_collection[data_collection$colony %in% target_colonies, ]

# -> c11 is control and thus will be added as 14th colony to the selected colonies and we add C02 because it seems to be doing better on the threshold counts
colonies_to_analyse_final <- c(colonies_to_analyse_certain, "c02", "c11")



# Subset data for selected and non-selected colonies
selected_counts <- sorted_counts[names(sorted_counts) %in% colonies_to_analyse_final]
non_selected_counts <- sorted_counts[!(names(sorted_counts) %in% colonies_to_analyse_final)]

### plot sorted counts from the colony selection above
data_collection_selected_colonies <- subset(data_collection, colony %in% colonies_to_analyse_final)
treatments <- data_collection_selected_colonies$treatment[match(names(selected_counts), data_collection_selected_colonies$colony)] # Match treatments with colonies
num_colors <- 3
colors <- viridis(num_colors)[match(treatments, c("cc", "vy", "vb"))]
# Plot for selected colonies & non-selected colonies
par(mfrow = c(1, 2))
barplot(selected_counts, col = colors, ylim = c(0, max(sorted_counts)), las = 2, main = "Selected Colonies")
barplot(non_selected_counts, col = "gray", ylim = c(0, max(sorted_counts)), las = 2, main = "Non-Selected Colonies")


 












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
head(summed_dat)
selected_ants <- summed_dat %>% 
  filter(total_duration >= 30, excluder == 0) %>% 
  group_by(colony, position) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  select(colony, position, antID = focal_AntID, tagID)
table(selected_ants$colony)

# Save selected ants to a text file
write.table(selected_ants, file = "vital_ants_to_annotate.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)




