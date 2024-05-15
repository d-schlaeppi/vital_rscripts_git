rm(list = ls())

### ### ### ### ### ### ### ### ### ### ### 
### Preliminary bead data analysis  ### ###
### ### ### ### ### ### ### ### ### ### ### 

#### Read Me ####
# First rough look at the beads data to get an idea on the food consumption and spread through colonies

#### Todo ####

# 1. Load the data: Colony meta data, brood bead data, and individual metadata
# Assign bead colours to virus vs control food sources
# Check for general differences in the colours (probably to include later as a random factor)
# Check for general differences in the values per food source 


#### Prerequisites ####
library(tcltk)
library(dplyr)
library(tidyr)

# set directories
setwd(tk_choose.dir(default = "~/", caption = "Select Working Directory")) # direct it to where you have config_user_and_hd.R (typically your scriptfolder)
source("config_user_and_hd.R") # contains getUserOptions() that defines usr and hd and the clean() function
DATADIR <- paste("/media", usr, hd, "vital/fc2", sep="/")
SCRIPTDIR <- paste("/home",usr,"Documents/vital_rscripts_git", sep="/")
setwd(DATADIR) #setwd("/media/gw20248/DISK_B/vital/fc2")

# read tables & prep tables
bead_data <- read.csv("vital_bead_data.csv", header = TRUE)
brood <- read.csv("vital_brood_bead_data.csv", header = TRUE)
individual_metadata <- read.csv("individual_metadata_vital_updated.csv", header= TRUE)
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # loads colony metadata

brood$colony_id <- NA
brood$colony_id <- ifelse(brood$comment != "", brood$comment, NA)
col_names <- c("flowjo_plate", "flowjo_row_well", "flowjo_column_strip", "flowjo_name" ,"yellow_count_YB", "blue_count_NB", "flowjo_sampletype", "comment", "colony_id")
col_names <- c("flowjo_plate", "flowjo_row_well", "flowjo_column_strip", "flowjo_name", "yellow_count_YB", "blue_count_NB", "flowjo_sampletype", "comment", "colony_id")
names(brood) <- col_names

# add brood data frame to bead_data
# Identify missing columns, Create missing columns in brood and fill them with NA values
missing_columns <- setdiff(names(bead_data), names(brood))
for (col in missing_columns) {
  brood[[col]] <- NA
}
combined_data <- rbind(bead_data, brood)

# Correct the bead data numbers to correct for the blanks: 
mean_blue_blanks <- round(mean(combined_data$blue_count_NB[combined_data$flowjo_sampletype == "blank"], na.rm = TRUE))
mean_yellow_blanks <- round(mean(combined_data$yellow_count_YB[combined_data$flowjo_sampletype == "blank"], na.rm = TRUE))
combined_data$yellow_beads_cor <- ifelse(combined_data$flowjo_sampletype == "sample", combined_data$yellow_count_YB - mean_yellow_blanks, combined_data$yellow_count_YB)
combined_data$blue_beads_cor <- ifelse(combined_data$flowjo_sampletype == "sample", combined_data$blue_count_NB - mean_blue_blanks, combined_data$blue_count_NB)

# correct negative values to be 0
# function to replace negatives 
replace_negatives <- function(x) {
  x[x<0] <- 0
  return(x)
}
# idea replace things in a loop... 
# to_correct <- c("combined_data$yellow_beads_cor", "combined_data$blue_beads_cor", "individual_metadata$yellow_beads_cor", "individual_metadata$blue_beads_cor")
# for (column in to_correct) {
#   df_name <- sub("\\$.*", "", column)
#   col_name <- sub(".*\\$", "", column)
#   df <- get(df_name)
#   replacement_data <- replace_negatives(df[[col_name]])
#   df[[col_name]] <- replacement_data
#   assign(df_name, df, envir = .GlobalEnv)
# }
# easier alternative... 

# Apply the function to selected columns in individual_metadata and comined_data
cols_to_adjust <- c("yellow_beads_cor", "blue_beads_cor")  
individual_metadata[cols_to_adjust] <- lapply(individual_metadata[cols_to_adjust], replace_negatives)
combined_data[cols_to_adjust] <- lapply(combined_data[cols_to_adjust], replace_negatives)

# compare overall number of blue and yellow beads 
individual_metadata %>% 
  summarize(mean_blue = mean(blue_beads_cor, na.rm = TRUE),
            median_blue = median(blue_beads_cor, na.rm = TRUE),
            sd_blue = sd(blue_beads_cor, na.rm = TRUE),
            mean_yellow = mean(yellow_beads_cor, na.rm = TRUE),
            median_yellow = median(yellow_beads_cor, na.rm = TRUE),
            sd_yellow = sd(yellow_beads_cor, na.rm = TRUE))
wilcox.test(individual_metadata$yellow_beads_cor, individual_metadata$blue_beads_cor)
boxplot(log10(individual_metadata[individual_metadata$treatment_simple == "control", ]$yellow_beads_cor+0.5), log10(individual_metadata[individual_metadata$treatment_simple == "control", ]$blue_beads_cor+0.5), 
        main = "Beads number for each color in the controls",
        names = c("yellow", "blue"),
        xlab = "color", ylab = "count (log10)", 
        col = c("yellow", "blue"))
segments(1, 4, 2, 4, lwd = 2)
text(1.5, 4.2, "*", cex = 2)
### noteworthy --> the two beads are not equal there might be a systematic error in that overall there seem to be more blue beads - either because it was consumed more (food preference of the ants) or because it was slightly higher in the feeding solution or other explanations possible
# --> needs to be accounted for latter on by including bead color as a random factor
# Might need some Brain power for me to figure out how to include this in the model? maybe for each ant there will be two lines in the data set one wit the food 1 and the variables type, colour, bead number

#### creating the virus bead variable ####
names(colony_metadata)
names(individual_metadata)


# Merge colony_metadata and individual_metadata based on colony_id
merged_data <- merge(individual_metadata[, !names(individual_metadata) %in% c("treatment", "treatment_simple")], 
                     colony_metadata, by = "colony_id")

# Create a function to calculate food_1v and food_2c based on treatment
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

# Initialize new columns
merged_data$food_1v <- NA
merged_data$food_2c <- NA

# Loop over each row and calculate food_1v and food_2c
for (i in 1:nrow(merged_data)) {
  result <- calculate_beads(merged_data$treatment[i], 
                            merged_data$yellow_beads_cor[i], 
                            merged_data$blue_beads_cor[i])
  merged_data$food_1v[i] <- result$food_1v
  merged_data$food_2c[i] <- result$food_2c
}

# Calculate beads_combined
merged_data$beads_combined <- merged_data$food_1v + merged_data$food_2c
head(merged_data)
names(merged_data)

# making a "long" form of the dataframe
# Reshape the dataframe to have one column for bead_count and another for food_type
long_df <- merged_data %>%
  pivot_longer(cols = c(food_1v, food_2c),
               names_to = "food_type",
               values_to = "bead_count")
long_df$bead_color <- ifelse(long_df$treatment == "cc" & long_df$food_type == "food_1v", "yellow", 
                             ifelse(long_df$treatment == "cc" & long_df$food_type == "food_2c", "blue", 
                                    ifelse(substr(long_df$treatment, 2, 2) == "b" & long_df$food_type == "food_1v", "blue", "yellow")))
long_df$food_type <- ifelse(long_df$treatment == "cc", "control",
                            ifelse(substr(long_df$treatment, 2, 2) == "b" & long_df$food_type == "food_1v", "virus", "control"))

head(virus)
virus <- subset(long_df, treatment_simple == "virus")
boxplot(log10(bead_count+1) ~ food_type, data = virus)
names(virus)

#'Comparisons to make: 
#'Overall food (nr of beads) when comparing control with treatment colonies.
#'Flow of virus vs control food source in within the treatment colonies: overall, just the treated ones, just the non-treated nestmates
#'Check if and how much of the beads reached the 
#'