### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### Trophy Validation ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#' Script contains the code used to check for false positive trophallaxis events in the automatically classified interaction lists
#' Written by DS
#' In brief: 150 trophallaxis events were randomly sampled from the post treatment interaction lists, 
#' Those interactions were then manually checked in with the fort video footage. 
#' As short interactions appeared to have more false positives and optimal treshold for subsetting the data is calculated.

# note, when running this the first time I did not set a seed so the random selection cannot be repeated... make better next time... 

#### 1. Prerequisites ####
rm(list = setdiff(ls(), "first_time_use_working_directory"))

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

source("config_user_and_hd.R")
library(purrr)

path <- paste(DATADIR,"vital_experiment/main_experiment_trophallaxis/intermediary_analysis_steps/full_interaction_lists/PostTreatment/observed", sep="/")
files <-  list.files(path, full.names = TRUE)
all_interactions <- NULL

#### 2. Sample interactions ####

# get list wiht all trophallactic interactions
for (trophy_file in files) { # trophy_file <- trophy_files[1] 
  trophy_data <- read.table(trophy_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  all_interactions <- rbind(all_interactions, trophy_data)}

# sample 150 random events
set.seed(69)
trophy_data_to_check_manually <- all_interactions[sample(nrow(all_interactions), 150), ]

# edit and safe data for manual annotations
library(purrr)
trophy_data_to_check_manually <- trophy_data_to_check_manually %>% arrange(colony, start) %>%
  mutate(fort_time = map(start, ~ {
    time_fm <- fmTimeCreate(offset = as.numeric(as.POSIXct(.x, format="%Y-%m-%d %H:%M:%S", tz="GMT")))
    capture.output(time_fm)  # Capture the output as a string
  })) %>% 
  dplyr::select(colony, fort_time,Tag1, Tag2, start, end, pair, T_start_UNIX, T_stop_UNIX, Starttime, Stoptime, duration) %>% as.data.frame()

trophy_data_to_check_manually$fort_time <- unlist(trophy_data_to_check_manually$fort_time)
write.table(trophy_data_to_check_manually, file=paste(DATADIR, "vital_experiment/main_experiment_trophallaxis/intermediary_analysis_steps/trophy_data_to_check_manually.txt",sep="/"),append=F,quote=F,row.names=F,col.names=T)

#### Manual Annotation ####
# done manually and saved in SCRIPTDIR


#### Evaluation ####
validation_file <- paste(SCRIPTDIR, "trophy_validation_data.csv", sep = "/")
valium <- read.csv(validation_file)

tresholds <- seq(0, 10, by = 0.2)
results <- data.frame(treshold = numeric(0), 
                      prop_true_positives_over_t = numeric(0), 
                      prop_false_positives_over_t = numeric(0))
# Loop over the thresholds and calculate the metrics
for (treshold in tresholds) {
  nr_interactions_below_t <- sum(valium$duration < treshold)
  nr_interactions_above_t <- sum(valium$duration >= treshold)
  
  # True positives and false positives above threshold
  prop_true_positives_above_t <- sum(valium$interaction_trophy == "t" & valium$duration > treshold) / nr_interactions_above_t
  prop_false_positives_above_t <- sum(valium$interaction_trophy == "f" & valium$duration > treshold) / nr_interactions_above_t
  
  # Store the results in the dataframe
  results <- rbind(results, data.frame(treshold = treshold,
                                       prop_true_positives_over_t = prop_true_positives_above_t,
                                       prop_false_positives_over_t = prop_false_positives_above_t))
}

ggplot(results, aes(x = treshold, y = 1, fill = prop_true_positives_over_t)) +
  geom_tile() +
  scale_fill_viridis(name = "Proportion of True Positives", option = "C") + 
  labs(x = "Threshold", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Effect of Threshold on Proportion of True Positives")

results[which.max(results$prop_true_positives_over_t),]
optimal_t <- results$treshold[which.max(results$prop_true_positives_over_t)]
valium_sub <-valium %>% filter(duration >= optimal_t)

nr_interactions_tot         <- nrow(valium)
nr_false_positives_tot      <- sum(valium$interaction_trophy == "f")
prop_true_positives_tot     <- sum(valium$interaction_trophy == "t") / nrow(valium)
nr_interactions_sub_t       <- nrow(valium_sub)
nr_false_positives_new       <- sum(valium_sub$interaction_trophy == "f")
nr_interactions_cut         <- nr_interactions_tot - nrow(valium_sub)
nr_false_positives_cut      <- nr_false_positives_tot - nr_false_positives_new
prop_true_positives_new     <- sum(valium_sub$interaction_trophy == "t") / nrow(valium_sub)
prop_false_positives_new    <- sum(valium_sub$interaction_trophy == "f") / nrow(valium_sub)


# Print results
print_results <- TRUE
if (print_results) {
cat("Number of interactions                   =", nr_interactions_tot, "\n")
cat("Number of false positives                =", nr_false_positives_tot, "\n")
cat("Proportion of true positives             =", prop_true_positives_tot, "\n")
cat("Threshold to filter data [seconds]       =", optimal_t, "\n")
cat("Number of interactions post filtering    =", nr_interactions_sub_t, "\n")
cat("Number of false positives post filtering =", nr_false_positives, "\n")
cat("Number of interactions removed           =", nr_interactions_cut, "\n")
cat("Number of false positives removed        =", nr_false_positives_cut, "\n")
cat("Proportion of true positives new         =", prop_true_positives_new, "\n")
cat("Number of false positives new            =", nr_false_positives_new, "\n")
}

hist(valium$duration)
hist(valium$duration[which(valium$interaction_trophy=="f")])

# Create histogram with two colors based on 'interaction' (t or f)
ggplot(valium, aes(x = duration, fill = interaction_trophy)) +
  geom_histogram(binwidth = 8, position = "stack", alpha = 0.7) + 
  scale_fill_manual(values = c("t" = "blue", "f" = "red")) +  # Customize colors
  labs(x = "Duration (seconds)", y = "Frequency", fill = "Interaction Type") +
  ggtitle("Histogram of Duration with Interaction Types") +
  theme_minimal()




# nr_sub_treshold <- sum(valium_sub$duration < optimal_t) 
# nr_false_positives_sub_t <- sum(valium$interaction_trophy == "f" & valium$duration < optimal_t)
# prop_false_positives_sub_t_among_tot_false_positives <- nr_false_positives_sub_t/nr_false_positives_tot
# prop_false_positives_sub_t <- nr_false_positives_sub_t / sum(valium$duration < optimal_t)
# nr_false_positives_over_t <- sum(valium$interaction_trophy == "f" & valium$duration > optimal_t)
# prop_false_positives_over_t <- nr_false_positives_over_t / sum(valium$duration > optimal_t)
# prop_true_positives_over_t    <- sum(valium$interaction_trophy == "t" & valium$duration > optimal_t) / sum(valium$duration > optimal_t)
# prop_false_positives_over_t   <- sum(valium$interaction_trophy == "f" & valium$duration > optimal_t) / sum(valium$duration > optimal_t)

