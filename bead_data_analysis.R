rm(list = setdiff(ls(), "first_time_use_working_directory_bead_data"))

### ### ### ### ### ### ### ### ### ### ### 
### Preliminary bead data analysis  ### ###
### ### ### ### ### ### ### ### ### ### ### 

#### 1. Read Me ####
# First look at the flow cytometry bead data to get an idea on the food consumption and distribution in ant colonies
# Note: Controls had one for their food source defined as the "virus corresponding food source" in a ~ paired fashion rather than randomly... ?

# Preliminary bead data analysis INDEX

#' 1. Read Me
#' 2. Prerequisites & Dataprep
#' 2.1 Compiling of Data frames
#' 2.2 Draw up the bead variable in its final form
#' 3. INITIAL BEAD ANALYSES
#' 3.1 BLUE vs YELLOW BEADS
#' 3.2 Colony level food consumption
#' 3.3 Feeding Duration (annotations) vs. nuber of Beads
#' 3.4 Feeding duration vs. Food source
#' 3.5 Beads in Brood
#' 3.6 Beads in  Workers - Individual level analyse
#' 3.6.1 Treated Workers
#' 3.6.2 Nestmates (non-treated workers)


#_________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Prerequisites ####
if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") { # direct it to where you have config_user_and_hd.R (typically the script folder or github folder)
  standard <- "/media/ael/gismo_hd2/vital/vital_rscripts_git" # if you are always working from the same directory just put its name here and it will save you some clicking.  
  selected_dir <- if  (dir.exists(standard)) {standard} else {tcltk::tk_choose.dir(default = "~/", caption = "Select Working Directory")}
  if (is.null(selected_dir) || selected_dir == "") {
    cat("No directory selected. Exiting.\n")
    return()}
  setwd(selected_dir)
  first_time_use_working_directory <<- getwd()
  cat(crayon::blue(getwd()))
} else { setwd(first_time_use_working_directory)
  cat(crayon::blue(getwd())) }


source("config_user_and_hd.R") # contains getUserOptions() that defines usr and hd and the clean() function
source(paste0(SCRIPTDIR,"/vital_fc1/func_test_norm.R"))

# # should now also work on windows and if not quickly define inputs manually:
# DATADIR <- "D:/DISK_B/vital/fc2"
# SCRIPTDIR <- "D:/DISK_B/vital/vital_rscripts_git"

overdispersion_test <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson_chisq <- sum(rp^2)
  prat <- Pearson_chisq / rdf
  pval <- pchisq(Pearson_chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson_chisq, ratio = prat, rdf = rdf, p = pval)
}




#_________________________________________________________________________________________________________________________________________________________________________________________
#### 2.1 Creating the different dataframes used for the Analyses ####
# read tables
bead_data <- read.csv("vital_bead_data.csv", header = TRUE) # flowcytometry bead data for all workers that were analysed
brood <- read.csv("vital_brood_bead_data.csv", header = TRUE) # flowcytometry bead data for all brood (pooled larva samples from each colony)
individual_metadata <- read.csv("individual_metadata_vital_updated.csv", header= TRUE) # metadata for all workers involved in the tracking experiment (technically already has the bead_data for workers but needs to be corrected due to negative values)
source(paste0(SCRIPTDIR,"/vital_meta_data.R")) # loads colony metadata to know background information such as treatments for each colony


#_________________________________________________________________________________________________________________________________________________________________________________________
#### Check colony metadata to see which colonies to analyse are affected by Fluorescent bead contamination ####
# First iteration of food sources was good 
# Second iteration of food sources: the bead color was switched for the two control foods
# and for the virus foods: the color is also switched and the one that is supposed to be yellow there is a contamination 
# Step 1 - see how many and which colonies are affected 
# Step 2 - Correct the dataframes 
# Step 2.5 - move fix to colony metadata script so that later on the right color can be used for other analyses
# Step 2.6 - Consider also implementing the fix for the variable "bead_colour" of treated_ants in the individual_metadata beforehand.
# Step 3 - rerun the bead analysis with the right colors. 
# (Step 4) - maybe clean up script by updating the source meta file to the correct colors so the fix is no longer required here 

# When was the switch of food sources?
# 23.03.2022 --> new food was used starting with colony 9+ 
# after that the colours have been switched and for the ones that were supposed to be virus blue there is a lot of contamination with the other color meaning we probably do not used them
# we are left with 7 control colonies and 5 virus colonies that we want to analyse.

colnames(colony_metadata)[which(names(colony_metadata) == "treatment")] <- "treatment_old"
# make an exclusion variable for colonies for which the beads were not analysed or for colonies that received second generation food where there was bead contamination in virus_yellow
colony_metadata$exclude_colony_for_beadanalyses <- ifelse(colony_metadata$block %in% 1:2,
                                                          ifelse(colony_metadata$fluorescence_measured == "no", "yes", "no"),
                                                          ifelse(colony_metadata$fluorescence_measured == "no","yes",
                                                                 ifelse(colony_metadata$fluorescence_measured == "yes" & colony_metadata$treatment_old == "vy", "yes","no")))

# create corrected treatment column
switch_letters <- function(x) {
  x <- gsub("b", "X", x)  # replace 'b' with 'X' (temporary placeholder)
  x <- gsub("y", "b", x)  # replace 'y' with 'b'
  x <- gsub("X", "y", x)  # replace 'X' (originally 'b') with 'y'
  return(x)}

colony_metadata$treatment <- ifelse(colony_metadata$block %in% 1:2,
                                    colony_metadata$treatment_old,
                                    switch_letters(colony_metadata$treatment_old))


colnames(individual_metadata)[which(names(individual_metadata) == "bead_colour")] <- "bead_colour_old"
individual_metadata$block <- ceiling(as.numeric(gsub("[^0-9]", "", individual_metadata$colony_id)) / 4) # transform c01 into a numeric value only, divide it by 4 and round up to next integer. 
individual_metadata$bead_colour <- ifelse(individual_metadata$block %in% 1:2,
                                      individual_metadata$bead_colour_old,
                                      ifelse(individual_metadata$bead_colour_old == "blue", "yellow", 
                                             ifelse(individual_metadata$bead_colour_old == "yellow", "blue", NA)))



### ### preparation of tables ### ###
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


# individual metadata already contains the corrected bead numbers as this was done in the script "add_beads_to_metadata-R"
# However, it is still necessary to correct negative values to be 0 - could have been done in before obviously, but I forgot so let's just do it now!
# function to replace negatives 
replace_negatives <- function(x) { x[x<0] <- 0
  return(x)}
cols_to_adjust <- c("yellow_beads_cor", "blue_beads_cor")  
individual_metadata[cols_to_adjust] <- lapply(individual_metadata[cols_to_adjust], replace_negatives)
combined_bead_data[cols_to_adjust] <- lapply(combined_bead_data[cols_to_adjust], replace_negatives)




#_________________________________________________________________________________________________________________________________________________________________________________________

#### 2.2 creating the virus bead variable ####

# disentangling bead color form food source for individual metadata as well as for brood
brood_data <- subset(combined_bead_data, combined_bead_data$flowjo_sampletype == "sample_brood")

# Merge colony_metadata and individual metadata based on colony_id
merged_data_individuals <- merge(individual_metadata[, !names(individual_metadata) %in% c("treatment", "treatment_simple")], 
                                 colony_metadata, by = "colony_id")
merged_data_brood <- merge(brood_data[, !names(brood_data) %in% c("treatment", "treatment_simple")], 
                           colony_metadata, by = "colony_id")

### Identify which colonies to analyse
# based on whether fluorescence has been measured and whether it was not fed with the contaminated virus food source. 

subset_colonies <- colony_metadata[colony_metadata$exclude_colony_for_beadanalyses == "no", ]
colonies_to_analyse <- subset_colonies$colony_id

### create random assignment of control food sources to virus and pseudo virus:
set.seed(4) # set.seed(4) or 8 
treatment_c_random <- subset_colonies$treatment
c_rows <- grep("^c", subset_colonies$treatment)
new_treatments <- c(rep("cb", length(subset_colonies$treatment[subset_colonies$treatment == "cb"])), rep("cy", length(subset_colonies$treatment[subset_colonies$treatment == "cy"])))
new_treatments <- sample(new_treatments)
treatment_c_random[c_rows] <- new_treatments
subset_colonies$treatment_c_random <- treatment_c_random
# table(subset_colonies$treatment_c_random)
# table(subset_colonies$treatment)

### over ride treatment if you want to use random assignment of the pseudo virus food in the controls.
USE_RANDOM_CONTROL_FOOD_SOURCE_ASSIGNMENT <- TRUE
if (USE_RANDOM_CONTROL_FOOD_SOURCE_ASSIGNMENT == TRUE) {
  # function to update merged_data_individuals with the new treatment values
  update_control_treatment <- function(df, subset_colonies) {
    df %>% left_join(subset_colonies %>% 
           dplyr::select(colony_id, treatment_c_random), by = "colony_id") %>%
           mutate(treatment = ifelse(!is.na(treatment_c_random), treatment_c_random, treatment)) %>%
           dplyr::select(-treatment_c_random)
      }
  merged_data_individuals <- update_control_treatment(merged_data_individuals, subset_colonies)
  merged_data_brood <- update_control_treatment(merged_data_brood, subset_colonies)
  colony_metadata <- update_control_treatment(colony_metadata, subset_colonies)
}




# Create a function to extract the number of beads for the two food sources:  food_1v and food_2c based on treatment
calculate_beads <- function(treatment, yellow_beads, blue_beads) {
  if (is.na(treatment)) {
    food_1v <- NA
    food_2c <- NA
  } else {
    if (substr(treatment, 2, 2) == "b") {  # Virus food source (or pseudo virus food source) contains blue beads
      food_1v <- blue_beads
      food_2c <- yellow_beads
    } else {  # Virus food source contains yellow beads
      food_1v <- yellow_beads
      food_2c <- blue_beads
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
  for (i in 1:nrow(data)) {  # Loop over each row and get bead number for food_1v; food_2c & beads_combined | i <- 1
    result <- calculate_beads(data$treatment[i], 
                              data$yellow_beads_cor[i], 
                              data$blue_beads_cor[i])
    data$food_1v[i] <- result$food_1v
    data$food_2c[i] <- result$food_2c
  }
  data$beads_combined <- data$food_1v + data$food_2c
  updated_df <- data   # Save the intermediate data frame with the new variables and beads_combined
  
  # Reshape the dataframe to have one column for bead_count and another for food_type
  long_df <- data %>% pivot_longer(cols = c(food_1v, food_2c),
                 names_to = "bead_source",
                 values_to = "bead_count")
  long_df$bead_color <- ifelse(
    (substr(long_df$treatment, 2, 2) == "y" & long_df$bead_source == "food_2c")|(substr(long_df$treatment, 2, 2) == "b" & long_df$bead_source == "food_1v"), "blue", "yellow")
  
  long_df$food_type <- ifelse( substr(long_df$treatment, 1,1) == "c", "control",
                              ifelse(substr(long_df$treatment, 1, 1) == "v" & long_df$bead_source == "food_1v", "virus", "control"))
  return(list(updated_df = updated_df, long_df = long_df))
}


dfs <- c("individuals", "brood")
for (df in dfs) { # df <- "individuals"
  data <- get(paste0("merged_data_", df))
  results <- process_df(data)
  assign(paste0("updated_df_", df), results$updated_df)
  assign(paste0("long_df_", df), results$long_df)
}
# after the above there are 4 dataframes available: updated_df_individuals, long_df_individuals, updated_df_brood, long_df_brood



# we clear out some of the variables not used here just to make the dataframes a bit nicer to look at
my_dataframes <- c("updated_df_individuals", "long_df_individuals", "updated_df_brood", "long_df_brood")
not_used_here <-c("tagIDdecimal", "identifStart", "identifEnd", "treatment_time", "exp_stop_time", "exp_start_time",  
                  "sampling_time", "glass_beads", "meta_ID", "tag_reoriented", "N_treated",  
                  "count", "last_det", "rot", "ScanTime", "comments_flowcytometry", "colony_nr",  
                  "freezer_container", "freezer_container_position", "PosID", "flowjo_plate", "flowjo_row_well",  
                  "flowjo_column_strip", "flowjo_name", "experiment", "tracking_system_main",  
                  "tracking_system_feeding", "time_treatment_start", "time_treatment_end",  
                  "annotation_start", "annotation_stop", "mean_ant_lenght_px", "mean_ant_lenght_mm", "nr")
for (df_name in my_dataframes) {
  df <- get(df_name)
  df <- df[, setdiff(names(df), not_used_here)]
  assign(df_name, df)}


# create new variable saying for is virus positive: Yes for the virus food source in the virus treatment and yes in one randomly selected  of the two food sources in the controls.
updated_df_individuals$is_virus_positive <- ifelse(updated_df_individuals$food_1v > 0, "yes", "no") 
updated_df_individuals$is_food_positive  <- ifelse(updated_df_individuals$beads_combined > 0, "yes", "no") 








#_________________________________________________________________________________________________________________________________________________________________________________________
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3. INITIAL BEAD ANALYSES #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### 3.1 BLUE vs YELLOW BEADS ####
# compare overall number of blue and yellow beads to see if the number of beads (proxy for consumed food) is different between the two colors
# using controls only as there might to avoid a potential bias 
updated_df_individuals %>%
  filter(treatment_simple == "control") %>% 
  summarize(mean_blue = mean(blue_beads_cor, na.rm = TRUE),
            median_blue = median(blue_beads_cor, na.rm = TRUE),
            sd_blue = sd(blue_beads_cor, na.rm = TRUE),
            mean_yellow = mean(yellow_beads_cor, na.rm = TRUE),
            median_yellow = median(yellow_beads_cor, na.rm = TRUE),
            sd_yellow = sd(yellow_beads_cor, na.rm = TRUE))


wilcox.test(updated_df_individuals[updated_df_individuals$treatment_simple == "control", ]$yellow_beads_cor,updated_df_individuals[updated_df_individuals$treatment_simple == "control", ]$blue_beads_cor)
boxplot(log10(updated_df_individuals[updated_df_individuals$treatment_simple == "control", ]$yellow_beads_cor+0.5), log10(updated_df_individuals[updated_df_individuals$treatment_simple == "control", ]$blue_beads_cor+0.5), 
        main = "Beads number for each color in the controls",
        names = c("yellow", "blue"),
        xlab = "color", ylab = "count (log10)", 
        col = c("#FFF2AE", "skyblue"))
segments(1, 4, 2, 4, lwd = 2)
text(1.5, 4.2, "p=0.55 (wilcox)", cex = 1)

### noteworthy --> the two beads are equal - NICE, there appears to be no systematic error or random preference for one of the bead colors.
# --> still, will be including bead color as a random factor when modeling

#### number of ants and proportion of each colony positive for food ####
summary_prop <- updated_df_individuals %>% 
  group_by(colony_id) %>% 
  summarize(
    nr_positive_for_food = sum(is_food_positive == "yes", na.rm = TRUE),
    proportion_food_positive = nr_positive_for_food/first(colony_size), 
    treatment = first(treatment_simple)) %>% as.data.frame()
path_int_results <- paste0(DATADIR, "/vital_experiment/main_experiment/intermediary_analysis_steps/")
filename <- "col_summary_prop_food_positive.csv"
write.csv(summary_prop, paste0(path_int_results,filename), row.names = FALSE)

boxplot(summary_prop$proportion_food_positive ~ summary_prop$treatment)
shapiro.test(summary_prop$proportion_food_positive) # normally distributed
leveneTest(proportion_food_positive ~ treatment, data = summary_prop) # non sig. --> variances are equal (assumption holds)
t.test(proportion_food_positive ~ treatment, data = summary_prop, var.equal = TRUE) # Mean higher in virus than compared to controls but not statistically significant ,


#### 3.1.5 Identifiy treated ants and get information on what they were feeding on ####  

# capping any feeding event longer than 5 minutes. 
cap_threshold <- 6 * 60  # capping threshold in seconds

# identifying what food sources treated ants were feeding on: virus treatment - virus or control , control treatment: pseudo-virus or control
# had_access_to_food_X
updated_df_individuals <- updated_df_individuals %>%
  mutate(had_access_to_food_X = case_when(
    IsTreated == TRUE & treatment_simple == "virus" & food_source == "virus" ~ "food_1v",  
    IsTreated == TRUE & treatment_simple == "virus" & food_source == "control" ~ "food_2c",  
    IsTreated == TRUE & treatment_simple == "control" & treatment == "cb" & bead_colour == "blue" ~ "food_1v",  
    IsTreated == TRUE & treatment_simple == "control" & treatment == "cb" & bead_colour == "yellow" ~ "food_2c",  
    IsTreated == TRUE & treatment_simple == "control" & treatment == "cy" & bead_colour == "blue" ~ "food_2c",  
    IsTreated == TRUE & treatment_simple == "control" & treatment == "cy" & bead_colour == "yellow" ~ "food_1v",  
    TRUE ~ NA_character_  ), # Assign NA if no conditions are met
    feeding_duration_capped = ifelse(
      feeding_duration > cap_threshold,
      cap_threshold,
      feeding_duration)  # Capping feeding_duration
    )  



# all_treated_ants <- updated_df_individuals %>% 
#   filter(IsTreated == TRUE)
updated_df_individuals %>%
  group_by(colony_id) %>%
  summarize(num_treated_ants = sum(IsTreated == TRUE),
            num_treated_ants_nonexclusion = sum(IsTreated == TRUE & feeding_exclusion == 0),
            num_treated_excluded = sum(IsTreated == TRUE & feeding_exclusion == 1),
            num_treated_ants_dead = sum(IsTreated == TRUE & IsAlive == FALSE)) %>% as.data.frame()


# get ants which are treatment, that were not classified as "to exclude" during feeding annotations, and that survived to the end.
treated_ants <- updated_df_individuals %>% 
  filter(IsTreated == TRUE & feeding_exclusion == 0 & IsAlive == TRUE)

# fix retagged ants by merging double lines into one and deleting the second one: extract duplicate rows, adjust them, add them back to original data.
duplicate_rows     <- treated_ants %>% group_by(colony_id, antID) %>% filter(n() > 1) %>% ungroup()
non_duplicate_rows <- treated_ants %>%  group_by(colony_id, antID) %>% filter(n() == 1) %>% ungroup()
cleaned_duplicate_rows <- duplicate_rows %>% arrange(colony_id, antID) %>%
  mutate(across(everything(), ~ coalesce(.[1], .[2]))) %>%  # fill NAs from the second row to the first
  distinct(colony_id, antID, .keep_all = TRUE) %>% ungroup() # keep only one row per unique combination
treated_ants <- bind_rows(cleaned_duplicate_rows, non_duplicate_rows) %>% as.data.frame()




#____________________________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3.2 Colony level food consumption ####
# mean total sum of beads per food source per colony for the two treatments

colors <- c("#CCEBC5","#CCEBC5", "#FBB4AE","#CCEBC5")

for (who in c("untreated_only","treated_only", "everyone")){ # who = "everyone"     who = "untreated_only"     who = "treated_only"
  if (who=="everyone"){
    DF <- updated_df_individuals
  }else if (who=="untreated_only"){
    DF <- updated_df_individuals[which(updated_df_individuals$IsTreated==F),]
  }else{
    DF <- treated_ants
  }
  cat(red("###","\n") , blue("###",who, "###", "\n"))
  colony_sum <- DF %>%
    group_by(colony_id, treatment) %>% 
    summarise(
      total_beads_f1 = sum(food_1v, na.rm = TRUE),
      total_beads_f2 = sum(food_2c, na.rm = TRUE),
      total_beads_combined = sum(beads_combined, na.rm = TRUE)
    ) %>% as.data.frame()
  colony_sum$food1_colour <- substr(colony_sum$treatment,2,2)
  colony_sum$treatment_simple <- ifelse(substr(colony_sum$treatment,1,1) == "c", "control","virus")
  
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
  boxplot(colony_long$bead_count ~ colony_long$food*colony_long$treatment_simple, 
          main = paste("Bead count by treatment and food source",who,sep=" - "),
          xlab = "Treatment and Food Source",
          ylab = "Bead Count", col = colors
  );abline(v = 2.5, col = "black") 
  mod <- lmer(sqrt(bead_count) ~ treatment_simple * food + (1|colony_id) + (1|food1_colour ), data = colony_long)
  print(Anova(mod))
  # plot_model_diagnostics(mod)
  
  if(Anova(mod)["treatment_simple:food","Pr(>Chisq)"]<0.05){
    contrast_mat <- rbind("Control_F1-Control_F2" = c(0, 0, -1, 0),
                          "Control_F1-Virus_F1" = c(0, -1, 0, 0),
                          "Control_F1-Virus_F2" = c(0, -1, -1, -1),
                          "Control_F2-Virus_F1" = c(0, -1, 1, 0),
                          "Control_F2-Virus_F2" = c(0, -1, 0, -1),
                          "Virus_F1-Virus_F2" = c(0, 0, -1, -1))
    posthocs <- summary(glht(mod,linfct=contrast_mat),test=adjusted("fdr"))
    print(posthocs)
    emms <- emmeans(mod, ~ treatment_simple * food)
    pairwise_results <- contrast(emms, method = "pairwise", adjust = "fdr")
    print(summary(pairwise_results))
    cld_results <- cld(emms, Letters = letters, adjust = "fdr")
    print(cld_results)
  }
}




#' no difference detectable when looking at everyone or at treated only?
#' But when looking the summed up food that was shared with nestmates it there seems it appears that maybe more virus food reaches the nestmates!
#' See if such a pattern can be detected at individual level.







#_____________________________________________________________________________________________________________________________________________________________________________________
#### 3.3 Feeding Duration (annotations) vs. nuber of Beads  ####

# check if bead data reflects manual feeding observations:
# for treated individuals correlate bead count with individual feeding time.

  
# create a summary table for colony feeding amounts (based on color and food type) as well as feeding durations 
feeding_summary_colonies <- treated_ants %>% 
  group_by(colony_id) %>% 
  summarize(duration         = sum(feeding_duration, na.rm = TRUE),
            duration_capped  = sum(feeding_duration_capped, na.rm = TRUE),
            duration_food_1v_cap = sum(feeding_duration_capped[had_access_to_food_X == "food_1v"], na.rm = TRUE), 
            duration_food_2c_cap = sum(feeding_duration_capped[had_access_to_food_X == "food_2c"], na.rm = TRUE)) %>% 
  left_join(updated_df_individuals %>% 
      group_by(colony_id) %>% 
      summarize(amount_food_1v = sum(food_1v, na.rm = TRUE), 
                amount_food_2c = sum(food_2c, na.rm = TRUE),
                total_amount_of_food = amount_food_1v + amount_food_2c), by = "colony_id") %>% as.data.frame()

# merge with feeding colony bead data
feeding_data_merged <- merge(feeding_summary_colonies, colony_sum, by = "colony_id")
plots <- list()
variables <- c("duration", "duration_capped")
for (var in variables) { # var <-  "colony_total_feeding_duration"
  correlation <- cor(feeding_data_merged[[var]], feeding_data_merged$total_beads_combined, use = "complete.obs")
  #print(paste("Correlation coefficient:", correlation))
  p <- ggplot(feeding_data_merged, aes(x = !!sym(var), y = total_beads_combined)) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    ggtitle(paste("Scatter Plot and Regression Line (Correlation: ", round(correlation, 2), ")", sep = "")) +
    xlab(paste0("Total feeding ", var)) +
    ylab("Total Beads Combined") +
    theme_minimal()
  plots[[length(plots) + 1]] <- p}

food_type <-  c("food_1v", "food_2c")
for (var in food_type) {
  correlation <- cor(feeding_data_merged[[paste0("duration_", var, "_cap")]], feeding_data_merged[[paste0("amount_", var)]], use = "complete.obs")
  print(paste("Correlation coefficient for", var, ":", correlation))
  p <- ggplot(feeding_data_merged, aes(x = !!sym(paste0("duration_", var, "_cap")), y = !!sym(paste0("amount_", var)))) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    ggtitle(paste("Scatter Plot and Regression Line (Correlation: ", round(correlation, 2), ")", sep = "")) +
    xlab(paste("Feeding duration on", var)) +
    ylab("Amount of food") +
    theme_minimal()
  plots[[length(plots) + 1]] <- p}
grid.arrange(grobs = plots, ncol = 2)


plots_2 <- list()  
treaty <- c("control", "virus")
for (treat in treaty) { # treat <-  "control"
  feeding_data_merged_sub <- feeding_data_merged %>% filter(treatment_simple == treat)
  for (var in food_type) {
    correlation <- cor(feeding_data_merged_sub[[paste0("duration_", var, "_cap")]], feeding_data_merged_sub[[paste0("amount_", var)]], use = "complete.obs")
    print(paste("Correlation coefficient for", var, ":", correlation))
    p <- ggplot(feeding_data_merged_sub, aes(x = !!sym(paste0("duration_", var, "_cap")), y = !!sym(paste0("amount_", var)))) +
      geom_point() +
      geom_smooth(method = "lm", col = "blue") +
      ggtitle(paste("Treatment:", treat, "; Food:", var, "-- cor uptake vs duration (Cor:", round(correlation, 2), ")", sep = " ")) +
      xlab(paste("Feeding duration on", var)) +
      ylab("Amount of food") +
      theme_minimal()
    plots_2[[length(plots_2) + 1]] <- p}}
grid.arrange(grobs = plots_2, ncol = 2)





# feeding duraton maybe not the beist predictor but look at individual level first: 

# Looking at individual level feeders only and food source they fed on only:
cor_df <- NULL
for (i in 1:nrow(treated_ants)) { # i = 1
  colony_id        <- treated_ants$colony_id[i]
  treatment        <- treated_ants$treatment[i]
  treatment_simple <- treated_ants$treatment_simple[i]
  tagID            <- treated_ants$tagID[i]
  feeding_duration <- treated_ants$feeding_duration[i]
  feeding_duration_capped <- treated_ants$feeding_duration_capped[i]
  food_source      <- treated_ants$food_source[i]
  color            <- treated_ants$bead_colour[i]
  nr_beads         <- ifelse(color == "yellow", treated_ants$yellow_beads_cor[i], treated_ants$blue_beads_cor[i])
  nr_other_beads   <- ifelse(color == "yellow", treated_ants$blue_beads_cor[i], treated_ants$yellow_beads_cor[i])
  right_color_of_beads_dominant <- ifelse(nr_beads-nr_other_beads >= 0, "yes", "no")
  cor_df <-  rbind(cor_df,data.frame(colony_id,               treatment,
                                     treatment_simple,
                                     tagID,                   feeding_duration,
                                     feeding_duration_capped, food_source, 
                                     nr_beads,                nr_other_beads, 
                                     right_color_of_beads_dominant,
                                     stringsAsFactors = F ))
}

# correlate feeding duration and amount of the corresponding beads. 
cor(cor_df$feeding_duration, cor_df$nr_beads, use = "complete.obs")
cor.test(cor_df$feeding_duration, cor_df$nr_beads, method = "pearson") # still significant despite the fact that this is after the food is shared with nestmates so the manual feeding annotations were not completely off...
plot(cor_df$feeding_duration, cor_df$nr_beads,
     xlab = "Feeding Duration", ylab = "Number of Beads",
     main = "Scatter plot of Feeding Duration vs Number of Beads")
abline(lm(cor_df$nr_beads ~ cor_df$feeding_duration), col = "blue")


#### 3.4 Feeding duration vs. Food source ####
# compare feeding duration among treated ants 
boxplot(cor_df$feeding_duration ~ cor_df$food_source)
boxplot(cor_df$feeding_duration_capped ~ cor_df$food_source)
cor_df %>% 
  group_by(treatment_simple, food_source) %>% 
  summarize(
    mean_duration = mean(feeding_duration_capped, na.rm = TRUE),
    sd_duration   = sd(feeding_duration_capped, na.rm = TRUE)
  )

cor_df %>%
  filter(!is.na(food_source)) %>%
  group_by(treatment_simple, food_source) %>%
  summarize(
    mean_duration = mean(feeding_duration_capped, na.rm = TRUE),
    sd_duration   = sd(feeding_duration_capped, na.rm = TRUE),
    .groups = "drop") %>% as.data.frame()

colors <- c("control" = "#CCEBC5", "virus" = "#FBB4AE")
ggplot(treated_ants %>% filter(!is.na(food_source)), aes(x = treatment_simple, y = feeding_duration_capped, fill = food_source)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  labs(title = "Feeding Duration by Food Source and Treatment",
       x = "Treatment",
       y = "Feeding Duration") +
  theme_classic()

# decide what model is right: 
# mod <- lmer(feeding_duration_capped ~ food_source + (1|treatment_simple)  + (1|colony_id), data = treated_ants, family = "poisson") # data not normally distributed
# mod <- glmer(feeding_duration_capped ~ food_source + (1|treatment_simple)  + (1|colony_id), data = treated_ants, family = "poisson") # quite ok but data is overdispersed which is not ok for poisson?!
mod <- glmmTMB(feeding_duration_capped ~ food_source + (1|treatment_simple)  + (1|colony_id), data = treated_ants, family = nbinom2) # ration is between 0-1 so should mean that this model is fine
mod <- glmmTMB(feeding_duration_capped ~ food_source + (1|colony_id), data = treated_ants, family = nbinom2) # variance of treatment_simple is extremely small (1.272e-07), which suggests that this random effect might not be contributing much to the model and was thus removed to improve the model
#resids <- residuals(object = mod)
#hist(resids)
#summary(mod)
Anova(mod)
#compareqqnorm(mod); par(mfrow = c(1,1))
#test_norm(mod)
overdispersion_test(mod)





#_________________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3.5 Beads in Brood ####

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
# but in the time window looked at it did not happen very much.
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










#_________________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3.6 Beads in  Workers - Individual level analyses ####

long_df_individuals <- as.data.frame(long_df_individuals)
long_df_individuals$IsPositive <- as.numeric(long_df_individuals$bead_count>0)

#### 3.6.1 Treated Workers ####
# get a feel for the data
updated_df_individuals %>% 
  filter(status_ant == "treated") %>%
  group_by(treatment_simple) %>%
  summarize(nr_nestmates = n(),
            nr_f1_pos = sum(food_1v > 0, na.rm = TRUE),
            mean_f1 = mean(food_1v, na.rm = TRUE),
            prop_f1_pos = nr_f1_pos / nr_nestmates,
            nr_f2_pos = sum(food_2c > 0, na.rm = TRUE),
            mean_f2 = mean(food_2c, na.rm = TRUE),
            prop_f2_pos = nr_f2_pos / nr_nestmates,
            food_positive = sum(beads_combined > 0, na.rm = TRUE),
            all_food = mean(beads_combined, na.rm = TRUE)) %>% as.data.frame()



#treated_workers_long <- subset(long_df_individuals, status_ant == "treated" & !is.na(bead_count))
treated_workers_long <- long_df_individuals %>%
  filter(paste(colony_id, antID) %in% paste(treated_ants$colony_id, treated_ants$antID))

# among treated workers the number of beads might be higher for one of the food sources but is it significant?
# At least looking at is positive it appears that in the virus treatment more treated ants are positive for the control food.

aggregate(IsPositive ~ bead_source + treatment_simple, FUN=mean, data=treated_workers_long)


### plot
ggplot(treated_workers_long, aes(x = treatment_simple, y = log(bead_count+0.5), fill = bead_source)) +
  geom_boxplot() +
  geom_text(data = labels_df, aes(label = bead_source, y = coordinates - 1.5), 
            position = position_dodge(width = 0.75), vjust = 0) +
  labs(title = "Treated Workers - Bead Count by Treatment and Food Type",
       x = "Treatment Simple",
       y = "log(Bead Count+0.5)") +
  theme_minimal() +
  scale_fill_manual(values = c("food_1v" = "#FBB4AE", "food_2c" = "#CCEBC5"))

# figure out the stats... 
# data not normally distributed and zero inflated and needs, random factors 
# generalized mixed model with poisson distribution and random factors for colony, bead color and individual but taking into account zero inflation... 

### binomial modal to test for difference if positive or not... figure out how decide what "is positive" would mean in the controls with two control food sources... 


#### Nestmates (non-treated workers) ####
# feel for the data
scales <- c("treatment", "treatment_simple")

for (scale in scales) {
  to_print <- updated_df_individuals %>% 
    filter(status_ant == "untreated") %>%
    group_by(across(all_of(scale))) %>%
    summarize(
      nr_nestmates = n(),
      nr_f1_pos = sum(food_1v > 0, na.rm = TRUE),
      mean_f1 = mean(food_1v, na.rm = TRUE),
      prop_f1_pos = nr_f1_pos / nr_nestmates,
      nr_f2_pos = sum(food_2c > 0, na.rm = TRUE),
      mean_f2 = mean(food_2c, na.rm = TRUE),
      prop_f2_pos = nr_f2_pos / nr_nestmates,
      food_positive = sum(beads_combined > 0, na.rm = TRUE),
      all_food = mean(beads_combined, na.rm = TRUE)
    ) %>% 
    as.data.frame()
print(to_print)
}



### plot
nestmates_long <- subset(long_df_individuals, status_ant == "untreated" & !is.na(bead_count))
nestmates_long <- as.data.frame(nestmates_long)
nestmates_long$food1_colour <- substr(nestmates_long$treatment,2,2)

ggplot(nestmates_long, aes(x = treatment_simple, y = log(bead_count+0.5), fill = bead_source)) +
  geom_boxplot() +
  geom_text(data = labels_df, aes(label = bead_source, y = coordinates - 1.5), 
            position = position_dodge(width = 0.75), vjust = 0) +
  labs(title = "Bead Count by Treatment and Food Type",
       x = "Treatment Simple",
       y = "log(Bead Count+0.5)") +
  theme_minimal() +
  scale_fill_manual(values = c("food_1v" = "#FBB4AE", "food_2c" = "#CCEBC5"))

###_____________________________________________________________________________________________________________________________________________________________________________________
###
#### Binomial Analysis ####

subsets <- list(
  "nestmates" = long_df_individuals[long_df_individuals$IsTreated == FALSE,],
  "treated_only" = treated_workers_long,
  "all_individuals" = long_df_individuals
)
cat(red("WARNING: The following LOOPS NEEDS TO BE CHECKED... There is some issue with the stats"), "\n")
for (subset_name in names(subsets)) {
  cat("\n", green("_____________________________________________________________________________________________"), "\n")
  cat(blue("Subset: ", subset_name), "\n")
  current_data <- subsets[[subset_name]]
  current_data$food1_colour <- substr(current_data$treatment,2,2)
  mean_data <- aggregate(IsPositive ~ treatment_simple + bead_source, 
                         FUN = mean, 
                         data = current_data)
  std_error_data <- aggregate(IsPositive ~ treatment_simple + bead_source, 
                              FUN = std.error,  
                              data = current_data)
  names(std_error_data)[names(std_error_data) == "IsPositive"] <- "std_error"

  binomial_data_mean_std <- left_join(mean_data, std_error_data, by = c("treatment_simple", "bead_source"))
  
  p <- ggplot(binomial_data_mean_std, aes(x = treatment_simple, y = IsPositive, fill = bead_source)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_errorbar(aes(ymin = IsPositive - std_error, ymax = IsPositive + std_error), 
                  position = position_dodge(0.9), width = 0.25) +
    labs(title = paste("Mean IsPositive by Treatment and Bead Source (", subset_name, ")", sep = ""),
         x = "Treatment",
         y = "Mean IsPositive") +
    scale_fill_manual(values = c("food_2c" = "#CCEBC5", "food_1v" = "#FBB4AE")) +
    theme_classic()
  
  print(p)

  model_binomial <- glmer(IsPositive ~ treatment_simple * bead_source + 
                             (1|colony_id) + (1|colony_id/antID), # + (1|food1_colour) 
                          family = binomial, data = current_data)
  
  # Run ANOVA on the model
  anova_result <- Anova(model_binomial)
  print("ANOVA Results:")
  print(anova_result)
  
  #if (anova_result["treatment_simple:bead_source", "Pr(>Chisq)"] < 0.05) {
  if (TRUE) {
    contrast_mat <- rbind(
      "Control_F1 minus Control_F2" = c(0, 0, -1, 0),
      "Control_F1 minus Virus_F1" = c(0, -1, 0, 0),
      "Control_F1 minus Virus_F2" = c(0, -1, -1, -1),
      "Control_F2 minus Virus_F1" = c(0, -1, 1, 0),
      "Control_F2 minus Virus_F2" = c(0, -1, 0, -1),
      "Virus_F1 minus Virus_F2" = c(0, 0, -1, -1)
    )
    print(summary(glht(model_binomial, linfct = contrast_mat), test = adjusted("BH")))
  }
}

# binomial just within the treatments:

treaty <- c("control", "virus")
for(treat in treaty ){ # treat <- treaty[1]
  cat(red("\n TREATMENT: ", toupper(treat), "\n"))
  for (subset_name in names(subsets)) { # subset_name <- names(subsets)[1]
   cat("\n", green("_____________________________________________________________________________________________"), "\n")
   cat(blue("Subset: ", subset_name), "\n")
   current_data <- subsets[[subset_name]]
   current_data <- current_data %>% filter(treatment_simple == treat)
   current_data$food1_colour <- substr(current_data$treatment,2,2)
   mean_data <- aggregate(IsPositive ~ treatment_simple + bead_source, 
                         FUN = mean, 
                         data = current_data)
   std_error_data <- aggregate(IsPositive ~ treatment_simple + bead_source, 
                              FUN = std.error,  
                              data = current_data)
   names(std_error_data)[names(std_error_data) == "IsPositive"] <- "std_error"
  
   binomial_data_mean_std <- left_join(mean_data, std_error_data, by = c("treatment_simple", "bead_source"))
  
   p <- ggplot(binomial_data_mean_std, aes(x = treatment_simple, y = IsPositive, fill = bead_source)) +
     geom_bar(stat = "identity", position = position_dodge(), color = "black") +
     geom_errorbar(aes(ymin = IsPositive - std_error, ymax = IsPositive + std_error), 
                   position = position_dodge(0.9), width = 0.25) +
     labs(title = paste("Mean IsPositive by Treatment and Bead Source (", subset_name, ")", sep = ""),
          x = "Treatment",
          y = "Mean IsPositive") +
     scale_fill_manual(values = c("food_2c" = "#CCEBC5", "food_1v" = "#FBB4AE")) +
     theme_classic()
    print(p)
  
   model_binomial <- glmer(IsPositive ~  bead_source + 
                            (1|colony_id) + (1|colony_id/antID), # + (1|food1_colour) 
                          family = binomial, data = current_data)
   anova_result <- Anova(model_binomial)
   print("ANOVA Results:")
   print(anova_result)
 }
}

# It seems that the virus food reaches fewer individuals than the control food, but how to get at this correctly?!




###_____________________________________________________________________________________________________________________________________________________________________________________
###
#### Asymmetry Coefficient ####
#' Assigning pseudo virus in the controls might not be the best approach. 
#' Depending on the random assignment the results change --> thus alternative way of analyzing the data might be better. 
#' Alternative way of analyzing the data using an asymmetry coefficient that is blind to the type of food.

### Colony level
#' For each colony check if the proportion of ants positive for one of the two food sources is even or uneven. 
#' Absolute values - |p1-p2|/|p1+p2| to get the assymetry coefficient

summary_prop <- updated_df_individuals %>% 
  group_by(colony_id) %>% 
  summarize(
    treatment = first(treatment_simple),
    nr_positive_for_food = sum(is_food_positive == "yes", na.rm = TRUE),
    proportion_food_positive = nr_positive_for_food/first(colony_size),
    nr_individuals = n(),
    colony_size = first(colony_size),
    proportion_positive_f1v = sum(food_1v > 0, na.rm = TRUE) / sum(!is.na(food_1v)),
    proportion_positive_f2c = sum(food_2c > 0, na.rm = TRUE) / sum(!is.na(food_2c)),
    p1 = proportion_positive_f1v,
    p2 = proportion_positive_f2c,
    asy_coef = abs(p1 - p2)/(p1+p2),
    asy_coef_log_directional = log((p1 + 0.001) / (p2 + 0.001)) 
    ) %>% 
  dplyr::select(-p1,-p2) %>%  as.data.frame()




coefficients <- c("asy_coef", "asy_coef_log_directional")
for (coef in coefficients) { # coef <- coefficients[1]
  cat(green("\n", "__________________________________________________________________________________________________", "\n"))
  cat(green("\n", coef, "\n"))
  print( #plot
    ggplot(summary_prop, aes(x = treatment, y = !!sym(coef), fill = treatment)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(title = paste("Asymmetry Coefficient by Treatment:", coef), 
           x = "Treatment", y = coef) +
      theme_minimal() +
      scale_fill_manual(values = c("control" = "lightblue", "virus" = "pink"))
  )
  # test noramlity
  shapiro_result <- shapiro.test(summary_prop[[coef]])
  print(shapiro_result)
  
  # Choose the appropriate test
  if (shapiro_result$p.value < 0.05) { # data is NOT normally distributed →  Wilcox
    test_result <- wilcox.test(reformulate("treatment", coef), data = summary_prop)
    cat("\nUsing Wilcoxon Test (non-parametric) for", coef, ":\n")
  } else { # data is normal → t-test
    test_result <- t.test(reformulate("treatment", coef), data = summary_prop, var.equal = TRUE)
    cat("\nUsing t-test (parametric) for", coef, ":\n")}
  print(test_result)
}


# ev linear model instead... 
# lm_model <- lm(coef ~ treatment, data = summary_prop)
# summary(lm_model)


### Individual level analysis of assymetry based on food amounts? could be interesting. 










###_____________________________________________________________________________________________________________________________________________________________________________________
###
#### Update metadata tables to get acces to this info in other scripts ####
###
###

### Individual Metadata
individual_metadata_orig <- read.csv("individual_metadata_vital_updated.csv", header= TRUE) # metadata for all workers involved in the tracking experiment (technically already has the bead_data for workers but needs to be corrected due to negative values)

#remove duplicates...
individual_metadata_orig <- distinct(individual_metadata_orig)

individual_metadata_orig <- individual_metadata_orig %>%
  distinct(colony_id, antID, tagID, .keep_all = TRUE) %>%
  group_by(colony_id, antID, tagID) %>%
  mutate(
    identifStart = first(identifStart, na.rm = TRUE),
    identifEnd = first(identifEnd, na.rm = TRUE)
  ) %>%
  ungroup() %>%  as.data.frame()

# individual_metadata_orig %>%
#   group_by(colony_id, antID, tagID) %>%
#   filter(n() > 1) %>%  as.data.frame()

updated_df_individuals_duplicates <- updated_df_individuals %>%
  group_by(colony_id, antID, tagID) %>%
  filter(n() > 1) %>%
  slice_head(n = 1) %>% 
  ungroup() %>%  as.data.frame()

updated_df_individuals <- bind_rows(
  updated_df_individuals_duplicates, 
  updated_df_individuals %>% 
    group_by(colony_id, antID, tagID) %>%
    filter(n() == 1) %>%
    ungroup()
) %>% as.data.frame()



### making sure all the correct variables are added to the updated individual_metadata table. 
updated_vars <- updated_df_individuals %>%
  dplyr::select(colony_id, tagID, antID, food_1v, food_2c, beads_combined,
         is_virus_positive, is_food_positive)

individual_metadata_orig <- individual_metadata_orig %>% 
  mutate(bead_colour_old = bead_colour, 
         block = ceiling(as.numeric(gsub("[^0-9]", "", colony_id)) / 4), # transform c01 into a numeric value only, divide it by 4 and round up to next integer. 
         treatment = colony_metadata$treatment[match(colony_id, colony_metadata$colony_id)], 
         bead_colour = ifelse(block %in% 1:2,
                                            bead_colour_old,
                                            ifelse(bead_colour_old == "blue", "yellow", 
                                                   ifelse(bead_colour_old == "yellow", "blue", NA))),
         food_position_1 = case_when(treatment == "cb" ~ "control_blue", treatment == "cy" ~ "control_yellow",treatment == "vb" ~ "virus_blue",treatment == "vy" ~ "virus_yellow"),
         food_position_2 = case_when(treatment == "cb" ~ "control_yellow",treatment == "cy" ~ "control_blue",treatment == "vb" ~ "control_yellow",treatment == "vy" ~ "control_blue"),
         fluorescence_measured = colony_metadata$fluorescence_measured[match(colony_id, colony_metadata$colony_id)],
         exclude_colony_for_beadanalyses = colony_metadata$exclude_colony_for_beadanalyses[match(colony_id, colony_metadata$colony_id)], 
         had_access_to_food_X = case_when(
           IsTreated == TRUE & treatment_simple == "virus" & food_source == "virus" ~ "food_1v",  
           IsTreated == TRUE & treatment_simple == "virus" & food_source == "control" ~ "food_2c",  
           IsTreated == TRUE & treatment_simple == "control" & treatment == "cb" & bead_colour == "blue" ~ "food_1v",  
           IsTreated == TRUE & treatment_simple == "control" & treatment == "cb" & bead_colour == "yellow" ~ "food_2c",  
           IsTreated == TRUE & treatment_simple == "control" & treatment == "cy" & bead_colour == "blue" ~ "food_2c",  
           IsTreated == TRUE & treatment_simple == "control" & treatment == "cy" & bead_colour == "yellow" ~ "food_1v",  
           TRUE ~ NA_character_  ), # Assign NA if no conditions are met
         feeding_duration_capped = ifelse(
           feeding_duration > cap_threshold,
           cap_threshold,
           feeding_duration)
  ) %>% left_join(updated_vars, by = c("colony_id", "tagID", "antID")) %>% as.data.frame()


write.csv(individual_metadata_orig, "individual_metadata_vital_with_bead_data_with_random_control_food_assignment.csv", row.names = FALSE)







### Colony Metadata:
#' we could opt to update the colony metadata with the new treatment variable, especially if we are using random assignment of food sources in the controls... but we might not do so either way.






### Treated Workers 
# replacing old treated ants with new extended version of treated ants


# get subset of indiviudal metadata table which already contains all the information by filtering for treated. 
treated_ants_all <- individual_metadata_orig %>% 
  filter(IsTreated == TRUE)

# clean duplicates | fix retagged ants by merging double lines into one and deleting the second one: extract duplicate rows, adjust them, add them back to original data.
duplicate_rows     <- treated_ants_all %>% group_by(colony_id, antID) %>% filter(n() > 1) %>% ungroup() 
cleaned_duplicates <- duplicate_rows %>% # Merge duplicate rows by filling NAs from the other row
  group_by(antID) %>%
  summarise(across(everything(), ~ coalesce(.x[1], .x[2])), .groups = "drop")
# add cleaned lines back to the other lines:
non_duplicate_rows <- treated_ants_all %>%  group_by(colony_id, antID) %>% filter(n() == 1) %>% ungroup()
treated_ants_all <- bind_rows(cleaned_duplicate_rows, non_duplicate_rows) %>% as.data.frame()

# # count what is there for each colony and whether it makes sence. 
# treated_ants_all %>% 
#   group_by(colony_id) %>%
#   summarize(num_treated_ants = n(),
#             num_treated_ants_nonexclusion = sum(feeding_exclusion == 0),
#             num_treated_feeding_excluded = sum(feeding_exclusion == 1),
#             num_treated_ants_dead = sum(IsAlive == FALSE), 
#             num_treated_ants_final = sum(IsAlive == TRUE & feeding_exclusion == 0)) %>% as.data.frame()


# add some variables to match stroeymeyt and adriano synthax
treated_ants_final <- treated_ants_all %>% 
  mutate(
    colony = colony_id,
    tag = antID,
    survived_treatment = "NA",
    REP_treat = paste0(colony_id, "_" , treatment_simple),
    final_inclusion = ifelse(feeding_exclusion == 0 & IsAlive == TRUE, TRUE, FALSE))



treated_file_orig <- paste(DATADIR, "vital_experiment/main_experiment/original_data/treated_worker_list.txt", sep = "/")
# Replace old version of treated file from tables to match stroeymeyt
# overwrite_confirmation <- readline(prompt = "Do you want to replace the existing file (and save the old one with a new name)? (yes/no): ")
{overwrite_confirmation <- menu(c("Yes", "No"), title = "Do you want to replace the existing file (and save the old one with a new name)?")
if (tolower(overwrite_confirmation) == 1) {
  if (file.exists(treated_file_orig)) {  # rename old version of the file 
    current_date <- format(Sys.Date(), "%Y-%m-%d")
    treated_file_backup <- paste(DATADIR, 
                                 "/vital_experiment/main_experiment/original_data/treated_worker_list_old_", 
                                 current_date, 
                                 ".txt", 
                                 sep = "")
    file.rename(treated_file_orig, treated_file_backup)
    cat("Existing file renamed to:", treated_file_backup, "\n")}
  write.table(treated_ants_final, treated_file_orig, row.names = FALSE, col.names = TRUE, sep = "\t")   # safe new file under old name
  cat("New treated table successfully saved\n")
} else {cat("File not saved.\n")}}







