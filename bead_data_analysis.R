rm(list = ls())

### ### ### ### ### ### ### ### ### ### ### 
### Preliminary bead data analysis  ### ###
### ### ### ### ### ### ### ### ### ### ### 

#### Read Me ####
# First look at the flowcytometry bead data to get an idea on the food consumption and distribution through colonies

#### Todo ####

#### Prerequisites ####
library(tcltk)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(car)

#model testing
library(blmeco)
library(DHARMa)
plot_model_diagnostics <- function(model) {
  compareqqnorm(model)   # Compare QQ plot
  cat("Press Enter to continue...")
  readline()
  sim_res <- simulateResiduals(model)   # Simulate residuals using DHARMa
  plot(sim_res)
  readline()
  par(mfrow = c(1, 3))
  testUniformity(sim_res)  # Test for uniformity (KS test),
  text(0.3, 0.95, " Test for uniformity", cex = 2)
  text(0.43, 0.93,"(KS test - N.S.--> The two compared distributions are identical)", cex = 1)
  testOutliers(sim_res)  # Test for outliers
  testDispersion(sim_res) # Test for homoscedasticity
  #text(1.5, 20, " Test for homoscedasticity against heteroscedacity", cex = 1)
  #text(1.5, 19,"(N.S.--> residuals have constant variance = homoscedasticity)", cex = 1)
  layout(1)
} # simple finction to display some model plots - not sure what models it works for... should be ok for lmer

# set directories
set_directories <- function() {
  cat("Select an option:\n") # Prompt and read the user for input
  cat("1. This is Daniel working on his Mac\n")
  cat("2. This is a random user working on Linux\n")
  cat("3. Else\n")
  option <- as.integer(readLines(n = 1))
  
  if (option == 1) { #MAC
    DATADIR <- "/Volumes/DISK_B/vital/fc2/"
    SCRIPTDIR <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/"
    setwd(DATADIR)
    assign("DATADIR", DATADIR, envir = .GlobalEnv)
    assign("SCRIPTDIR", SCRIPTDIR, envir = .GlobalEnv)
    cat("Directories set for Daniel on Mac.\n")
    cat("DATADIR:", DATADIR, "\n")
    cat("SCRIPTDIR:", SCRIPTDIR, "\n")
  } else if (option == 2) { #LINUX - config_user_and_hd.txt must exist in the selected working directory
    working_dir <- tk_choose.dir(default = "~", caption = "Select Working Directory")
    setwd(working_dir)
    source("config_user_and_hd.R")
    DATADIR <- paste("/media", usr, hd, "vital/fc2", sep = "/")
    SCRIPTDIR <- paste("/home", usr, "Documents/vital_rscripts_git", sep = "/")
    setwd(DATADIR)
    assign("DATADIR", DATADIR, envir = .GlobalEnv)
    assign("SCRIPTDIR", SCRIPTDIR, envir = .GlobalEnv)
    cat("Directories set for random user on Linux.\n")
    cat("DATADIR:", DATADIR, "\n")
    cat("SCRIPTDIR:", SCRIPTDIR, "\n")
  } else if (option == 3) { #ELSE
    cat("Please set the directories manually.\n")
  } else {
    cat("Invalid option selected. Please run the script again and select a valid option.\n")
  }
}
set_directories()
1
#### Creating the different dataframes used for the Analyses ####

# read tables
bead_data <- read.csv("vital_bead_data.csv", header = TRUE) # flowcytometry bead data for all workers that were analysed
brood <- read.csv("vital_brood_bead_data.csv", header = TRUE) # flowcytometry bead data for all brood (pooled larva samples from each colony)
individual_metadata <- read.csv("individual_metadata_vital_updated.csv", header= TRUE) # metadata for all workers involved in the tracking experiment (technically already has the bead_data for workers but needs to be corrected due to negative values)
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # loads colony metadata to know background information such as treatments for each colony

### preparation of tables ###

### brood
# creating a new variable containing colony id
brood$colony_id <- NA 
brood$colony_id <- ifelse(brood$comment != "", brood$comment, NA)
# renaming variables so it can be added to the bead data (for proper correction of bead counts)
col_names <- c("flowjo_plate", "flowjo_row_well", "flowjo_column_strip", "flowjo_name" ,"yellow_count_YB", "blue_count_NB", "flowjo_sampletype", "comment", "colony_id")
names(brood) <- col_names

# add brood data frame to bead_data --> combined data frame for the calculation of corrected bead numbers (subtraction of the mean value of the all the blanks to minimize false positives)
missing_columns <- setdiff(names(bead_data), names(brood))
for (col in missing_columns) {
  brood[[col]] <- NA
}
combined_bead_data <- rbind(bead_data, brood)
combined_bead_data$flowjo_sampletype <- trimws(combined_bead_data$flowjo_sampletype) # Remove leading and trailing spaces from flowjo_sampletype
individual_metadata$flowjo_sampletype <- trimws(individual_metadata$flowjo_sampletype)


# Correct the bead data numbers to correct for the blanks: 
mean_blue_blanks <- round(mean(combined_bead_data$blue_count_NB[combined_bead_data$flowjo_sampletype == "blank"], na.rm = TRUE))
mean_yellow_blanks <- round(mean(combined_bead_data$yellow_count_YB[combined_bead_data$flowjo_sampletype == "blank"], na.rm = TRUE))
combined_bead_data$yellow_beads_cor <- ifelse(combined_bead_data$flowjo_sampletype %in% c("sample", "sample_brood"), 
                                              combined_bead_data$yellow_count_YB - mean_yellow_blanks, 
                                              combined_bead_data$yellow_count_YB)
combined_bead_data$blue_beads_cor <- ifelse(combined_bead_data$flowjo_sampletype %in% c("sample", "sample_brood"),
                                            combined_bead_data$blue_count_NB - mean_yellow_blanks,
                                            combined_bead_data$blue_count_NB)
# individual metadata already contains the corrected bead numbers as this was done in the scritp "add_beads_to_metadata-R"
# However, it is necessary to correct negative values to be 0 - could have been done in previously... but let's do it now
# function to replace negatives 
replace_negatives <- function(x) {
  x[x<0] <- 0
  return(x)
}
cols_to_adjust <- c("yellow_beads_cor", "blue_beads_cor")  
individual_metadata[cols_to_adjust] <- lapply(individual_metadata[cols_to_adjust], replace_negatives)
combined_bead_data[cols_to_adjust] <- lapply(combined_bead_data[cols_to_adjust], replace_negatives)

#### creating the virus bead variable ####
# disentangling bead color form food source for individual metadata as well as for brood

brood_data <- subset(combined_bead_data, combined_bead_data$flowjo_sampletype == "sample_brood")

# Merge colony_metadata and individual metadata based on colony_id
merged_data_individuals <- merge(individual_metadata[, !names(individual_metadata) %in% c("treatment", "treatment_simple")], 
                                 colony_metadata, by = "colony_id")
merged_data_brood <- merge(brood_data[, !names(brood_data) %in% c("treatment", "treatment_simple")], 
                           colony_metadata, by = "colony_id")

# As there were so many individuals it was not possible to analyse all colonies --> colonies were ranked based on manual feeding observations and the the top 7 colonies of the two treatments were selected for the bead analysis
# Hence the individuals data set is subsetted in the bead analysis
# all colonies were analysed for brood (pooled larva samples from each colony)
# identify colonies not analysed in flowcytometry ( have bead count zero for both food sources)
filtered_data <- updated_df_individuals %>%
  group_by(colony_id, treatment_simple) %>%
  summarise(
    total_beads_food1v = sum(food_1v, na.rm = TRUE),
    total_beads_food2c = sum(food_2c, na.rm = TRUE),
    .groups = 'drop'  # Ungroup after summarising
  )  %>%
  filter(total_beads_food1v != 0 & total_beads_food2c != 0)
colonies_to_analyse <- filtered_data$colony_id


# Create a function to extract the number of beads for the two food sources:  food_1v and food_2c based on treatment
calculate_beads <- function(treatment, yellow_beads, blue_beads) {
  if (is.na(treatment)) {
    food_1v <- NA
    food_2c <- NA
  } else {
    if (treatment == "cc") {  # Control colonies
      food_1v <- yellow_beads
      food_2c <- blue_beads
    } else {  # Virus treatment colonies
      if (substr(treatment, 2, 2) == "b") {  # Virus food source contains blue beads
        food_1v <- blue_beads
        food_2c <- yellow_beads
      } else {  # Virus food source contains yellow beads
        food_1v <- yellow_beads
        food_2c <- blue_beads
      }
    }
  }
  return(list(food_1v = food_1v, food_2c = food_2c))
}
process_df <- function(data){
  if (df == "individuals") {
    data <- data %>% filter(colony_id %in% colonies_to_analyse)
  }
  data$food_1v <- NA
  data$food_2c <- NA
  for (i in 1:nrow(data)) {  # Loop over each row and get beadnumber for food_1v; food_2c & beads_combined | i <- 1
    result <- calculate_beads(data$treatment[i], 
                              data$yellow_beads_cor[i], 
                              data$blue_beads_cor[i])
    data$food_1v[i] <- result$food_1v
    data$food_2c[i] <- result$food_2c
  }
  data$beads_combined <- data$food_1v + data$food_2c
  updated_df <- data   # Save the intermediate data frame with the new variables and beads_combined
  # Reshape the dataframe to have one column for bead_count and another for food_type
  long_df <- data %>%
    pivot_longer(cols = c(food_1v, food_2c),
                 names_to = "bead_source",
                 values_to = "bead_count")
  
  long_df$bead_color <- ifelse(long_df$treatment == "cc" & long_df$bead_source == "food_1v", "yellow", 
                               ifelse(long_df$treatment == "cc" & long_df$bead_source == "food_2c", "blue", 
                                      ifelse(substr(long_df$treatment, 2, 2) == "b" & long_df$bead_source == "food_1v", "blue", "yellow")))
  
  long_df$food_type <- ifelse(long_df$treatment == "cc", "control",
                              ifelse(substr(long_df$treatment, 2, 2) == "b" & long_df$bead_source == "food_1v", "virus", "control"))
  return(list(updated_df = updated_df, long_df = long_df))
}



dfs <- c("individuals", "brood")
for (df in dfs) {
  data <- get(paste0("merged_data_", df))
  results <- process_df(data)
  assign(paste0("updated_df_", df), results$updated_df)
  assign(paste0("long_df_", df), results$long_df)
}

# after the above there are 4 dataframes available: updated_df_individuals, long_df_individuals, updated_df_brood, long_df_brood


### ### ### ### ### ### ### ### ### ### ### 
#### Analyses #### 
### ### ### ### ### ### ### ### ### ### ###




#### BLUE vs YELLOW BEADS ####
# compare overall number of blue and yellow beads to see if the number of beads (proxy for consumed food) is different between the two colors
updated_df_individuals %>% 
  summarize(mean_blue = mean(blue_beads_cor, na.rm = TRUE),
            median_blue = median(blue_beads_cor, na.rm = TRUE),
            sd_blue = sd(blue_beads_cor, na.rm = TRUE),
            mean_yellow = mean(yellow_beads_cor, na.rm = TRUE),
            median_yellow = median(yellow_beads_cor, na.rm = TRUE),
            sd_yellow = sd(yellow_beads_cor, na.rm = TRUE))
wilcox.test(updated_df_individuals$yellow_beads_cor, updated_df_individuals$blue_beads_cor)
boxplot(log10(updated_df_individuals[updated_df_individuals$treatment_simple == "control", ]$yellow_beads_cor+0.5), log10(updated_df_individuals[updated_df_individuals$treatment_simple == "control", ]$blue_beads_cor+0.5), 
        main = "Beads number for each color in the controls",
        names = c("yellow", "blue"),
        xlab = "color", ylab = "count (log10)", 
        col = c("yellow", "blue"))
segments(1, 4, 2, 4, lwd = 2)
text(1.5, 4.2, "*", cex = 2)
### noteworthy --> the two beads are not equal there might be a systematic error in that overall there seem to be more blue beads - either because it was consumed more (food preference of the ants) or because it was slightly higher in the feeding solution or other explanations possible
# --> needs to be accounted for latter on by including bead color as a random factor

#### Colony level food consumption ####

# mean sum of beads per food source per colony for the two treatments
# Summarize and get a feel for data

#by colny
colony_sum <- updated_df_individuals %>%
  group_by(colony_id, treatment_simple) %>% 
  summarise(
    total_beads_f1 = sum(food_1v, na.rm = TRUE),
    total_beads_f2 = sum(food_2c, na.rm = TRUE),
    total_beads_combined = sum(beads_combined, na.rm = TRUE)
  )
colony_sum

#by treatment
summary_df <- colony_sum %>% 
  group_by(treatment_simple) %>%
  summarise(
    mean_f1 = mean(total_beads_f1, na.rm = TRUE),
    sd_f1 = sd(total_beads_f1, na.rm = TRUE),
    median_f1 = median(total_beads_f1, na.rm = TRUE),
    mean_f2 = mean(total_beads_f2, na.rm = TRUE),
    sd_f2 = sd(total_beads_f2, na.rm = TRUE),
    median_f2 = median(total_beads_f2, na.rm = TRUE),
    n = n()
  )
summary_df

# plot the data
colony_long <- colony_sum %>% pivot_longer(cols = c(total_beads_f1, total_beads_f2),
               names_to = "food",
               values_to = "bead_count")
boxplot(colony_long$bead_count ~ colony_long$treatment_simple * colony_long$food, 
        main = "Bead count by treatment and food source",
        xlab = "Treatment and Food Source",
        ylab = "Bead Count",
        names = c("Control - Food 1", "Control - Food 2", "Virus - Food 1", "Virus - Food 2"))

mod <- lmer(bead_count ~ treatment_simple * food + (1|colony_id), data = colony_long)
Anova(mod)
plot_model_diagnostics(mod)

# based on the above the bead count is slighly higher for the virus food sompared to the control food sources but the difference is not significant

### stats
# check if bead data from treated ants roughly correlates with the the data from the manual feeding observations




#### Brood ####
head(updated_df_brood)

# Check if the brood has beads
updated_df_brood %>%
  group_by(treatment_simple) %>%
  summarize(nr_f1_pos = sum(food_1v > 0, na.rm = TRUE),
            mean_f1 = mean(food_1v, na.rm = TRUE),
            nr_f2_pos = sum(food_2c > 0, na.rm = TRUE),
            mean_f2 = mean(food_2c, na.rm = TRUE),
            food_positive = sum(beads_combined > 0, na.rm = TRUE),
            all_food = mean(beads_combined, na.rm = TRUE))

# only in six colonies (3x control, 3x virus) there is brood with beads indicating that food was delivered to them via the trophallxis feeding chain
# the number of beads in all six cases is very low

# plot the bead count by treatment and food type
labels_df <- long_df_brood %>%
  group_by(treatment_simple, bead_source) %>%
  summarize(mean_bead_count = mean(bead_count, na.rm = TRUE)) %>%
  mutate(coordinates = 0)  # Add a new column with default height

ggplot(long_df_brood, aes(x = treatment_simple, y = bead_count, fill = bead_source)) +
  geom_boxplot() +
  geom_text(data = labels_df, aes(label = bead_source, y = coordinates - 0.7), 
            position = position_dodge(width = 0.75), vjust = 0) +
  labs(title = "Bead Count by Treatment and Food Type",
       x = "Treatment Simple",
       y = "Bead Count") +
  theme_minimal() +
  scale_fill_manual(values = c("food_1v" = "blue", "food_2c" = "green"))


#### Treated Workers ####
# get a feel for the data
updated_df_individuals %>% 
  filter(status_ant == "treated") %>%
  group_by(treatment_simple) %>%
  summarize(nr_treated_ants = n(),
            nr_f1_pos = sum(food_1v > 0, na.rm = TRUE),
            mean_f1 = mean(food_1v, na.rm = TRUE),
            nr_f2_pos = sum(food_2c > 0, na.rm = TRUE),
            mean_f2 = mean(food_2c, na.rm = TRUE),
            food_positive = sum(beads_combined > 0, na.rm = TRUE),
            all_food = mean(beads_combined, na.rm = TRUE))

# among treated workers the number of beads is probably higher for the virus food source but is it significant?

### plot
treated_workers_long <- subset(long_df_individuals, status_ant == "treated" & !is.na(bead_count))
ggplot(treated_workers_long, aes(x = treatment_simple, y = bead_count, fill = bead_source)) +
  geom_boxplot() +
  geom_text(data = labels_df, aes(label = bead_source, y = coordinates - 0.7), 
            position = position_dodge(width = 0.75), vjust = 0) +
  labs(title = "Bead Count by Treatment and Food Type",
       x = "Treatment Simple",
       y = "Bead Count") +
  theme_minimal() +
  scale_fill_manual(values = c("food_1v" = "blue", "food_2c" = "green"))

ggplot(treated_workers_long, aes(x = treatment_simple, y = log(bead_count+0.5), fill = bead_source)) +
  geom_boxplot() +
  geom_text(data = labels_df, aes(label = bead_source, y = coordinates - 1.5), 
            position = position_dodge(width = 0.75), vjust = 0) +
  labs(title = "Bead Count by Treatment and Food Type",
       x = "Treatment Simple",
       y = "log(Bead Count+0.5)") +
  theme_minimal() +
  scale_fill_manual(values = c("food_1v" = "blue", "food_2c" = "green"))

#virus only
treated_workers_long_virus_only <-subset(treated_workers_long, treatment_simple == "virus")
boxplot(log(treated_workers_long_virus_only$bead_count+0.5) ~ treated_workers_long_virus_only$bead_source)

# figure out the stats... 
# data not normally distributed and zero inflated and needs, random factors 
# generalized mixed model with poisson distribution and random factors for colony, bead color and individual but taking into account zero inflation... 









#### Nestmates (non-treated workers) ####

# feel for the data
updated_df_individuals %>% 
  filter(status_ant == "untreated") %>%
  group_by(treatment_simple) %>%
  summarize(nr_nestmates = n(),
            nr_f1_pos = sum(food_1v > 0, na.rm = TRUE),
            mean_f1 = mean(food_1v, na.rm = TRUE),
            nr_f2_pos = sum(food_2c > 0, na.rm = TRUE),
            mean_f2 = mean(food_2c, na.rm = TRUE),
            food_positive = sum(beads_combined > 0, na.rm = TRUE),
            all_food = mean(beads_combined, na.rm = TRUE))

### plot
nestmates_long <- subset(long_df_individuals, status_ant == "untreated" & !is.na(bead_count))

ggplot(nestmates_long, aes(x = treatment_simple, y = bead_count, fill = bead_source)) +
  geom_boxplot() +
  geom_text(data = labels_df, aes(label = bead_source, y = coordinates - 0.7), 
            position = position_dodge(width = 0.75), vjust = 0) +
  labs(title = "Bead Count by Treatment and Food Type",
       x = "Treatment Simple",
       y = "Bead Count") +
  theme_minimal() +
  scale_fill_manual(values = c("food_1v" = "blue", "food_2c" = "green"))

ggplot(nestmates_long, aes(x = treatment_simple, y = log(bead_count+0.5), fill = bead_source)) +
  geom_boxplot() +
  geom_text(data = labels_df, aes(label = bead_source, y = coordinates - 1.5), 
            position = position_dodge(width = 0.75), vjust = 0) +
  labs(title = "Bead Count by Treatment and Food Type",
       x = "Treatment Simple",
       y = "log(Bead Count+0.5)") +
  theme_minimal() +
  scale_fill_manual(values = c("food_1v" = "blue", "food_2c" = "green"))

### stats

# stats number of ants positive in the virus treatment: 
# Create a data frame with all the relevant information
data <- data.frame(
  treatment = rep(c("Control", "Virus"), each = 2),
  food_type = rep(c("Control", "Control", "Virus", "Control")),
  positive = c(365, 375, 465, 389),
  total = c(1015, 1015, 1023, 1023)
)
# Fit a logistic regression model
model <- glm(cbind(positive, total - positive) ~ treatment + food_type, data = data, family = binomial)
summary(model)
# ants are more significantly more likely to positive for the virus food source compared to the control food source
# indicates that virus food sources might be overproportinally shared with nestmates

# is the number of beads also higher?
# generalized mixed model with poisson distribution and random factors for colony, bead color and individual but taking into account zero inflation...
mod <- glmmTMB(bead_count ~ food_type + (1|colony_id) + (1|bead_color) , data = nestmates_long, ziformula = ~1)
summary(mod)

# nonparametric test to see if the number of beads is higher for the virus food source
wilcox.test(nestmates_long$bead_count[nestmates_long$treatment_simple == "virus" & nestmates_long$bead_source == "food_2c"], 
            nestmates_long$bead_count[nestmates_long$treatment_simple == "virus" & nestmates_long$bead_source == "food_1v"])



