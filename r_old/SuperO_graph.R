####     Super O     ####

#' Read me: 
#' This is a preliminary script experimenting with potential ways of visualizing social insect colony health in terms of functionality
#' The idea is that colony functionality gets defined by a colonies need to uphold tasks required for survival

#' Notes: On the visualizations and the modelling: 
#' Check out ApisRam and whether I can use it to get a better version of my graphs?

rm(list = setdiff(ls(), "first_time_use_working_directory"))

# install.packages("polyclip", type = "binary")
pacman::p_load(dplyr, ggplot2, ggalluvial,pheatmap, reshape2, Rtsne, tibble, tidyverse, ggforce, fmsb, viridis, patchwork)
set.seed(69)


#' Note to myself: Creat the plot of how a healthy colony looks like... 
#' Then create a plot of what a unhealthy stressed colony woudl look like.
#' Then subset it into the different tasks: 
#' 1. E.g. Stressor that reduces each bees ability by 10 % in every task
#' 2. Stressor that reduces one specific task
#' 3. Event that decreases work force... shift... and reduction of spare work force that jump into gap. 
#' Order it based on age specific... 
#' Make seasonal plots (show the occurence of seasonal patterns such as gynes and drones in summer)
#' Show different colony sizes but that the overall pattern looks the same!
#' Use viridis style colors! 
#' Use some sort of aggregation by dominant tasks to create patterns... 
#' Use a dominant color to represent stressor impact (black or red?)
#' With the pie chart size could also be used for something like overall capacity of each worker...
#' Currently worker profiles are rather similar... e.g. there is no real reserve like fraction that does mostly do nothing... a bit like the drones or foraging specialists.
#' Worker cohorts with different profiles: Soldier bees, Royal patrilines
#' Create other optimal profiles for winter - dominated by thermoregulation


#________________________________________________________________________________________________________
#### 01. Creating a nice example colony ####

# Define task categories
tasks <- c(
  "reproduction_queen", "reproduction_gynes_drones", "brood_care_nursing", 
  "brood_thermoregulation", "nest_construction_repair", "nest_thermoregulation",
  "social_immunity", "foraging", "food_processing_storage", 
  "guarding_patrolling_defense", "active_defense", 
  "patrolling_task_allocation", "communication_behaviors", 
  "inactive_reserve", "self_maintenance"
)

# Age-based task contribution profile for workers (in %)
age_profile <- list(
  "1" = c(0, 0, 0.25, 0.1, 0.05, 0.02, 0.01, 0, 0.02, 0, 0, 0.01, 0, 0.5, 0.04),
  "2" = c(0, 0, 0.28, 0.1, 0.06, 0.02, 0.02, 0, 0.03, 0, 0, 0.01, 0, 0.4, 0.08),
  "3" = c(0, 0, 0.25, 0.08, 0.06, 0.03, 0.03, 0.02, 0.05, 0.01, 0.01, 0.02, 0.01, 0.3, 0.1),
  "4" = c(0, 0, 0.15, 0.05, 0.05, 0.04, 0.04, 0.1, 0.07, 0.02, 0.01, 0.03, 0.01, 0.25, 0.07),
  "5" = c(0, 0, 0.1, 0.04, 0.04, 0.04, 0.04, 0.2, 0.08, 0.03, 0.02, 0.03, 0.02, 0.2, 0.06),
  "6" = c(0, 0, 0.08, 0.04, 0.03, 0.05, 0.04, 0.3, 0.08, 0.03, 0.03, 0.03, 0.03, 0.15, 0.06),
  "7" = c(0, 0, 0.05, 0.04, 0.03, 0.05, 0.03, 0.4, 0.07, 0.03, 0.04, 0.04, 0.04, 0.1, 0.06)
)

# Normalize each age profile (sums to 1)
age_profile <- lapply(age_profile, function(x) x / sum(x))

# Function to create a worker
create_worker <- function(id, age) {
  profile <- age_profile[[as.character(age)]]
  task_profile <- rnorm(length(profile), mean = profile, sd = 0.02)
  task_profile <- pmax(task_profile, 0)  # prevent negative values
  task_profile <- task_profile / sum(task_profile)  # normalize
  
  # Ensure reproduction tasks are explicitly 0
  task_profile[1:2] <- 0
  
  out <- c(ID = id, caste = "worker", age = age, incapability = 0)
  names(task_profile) <- tasks
  as.list(c(out, task_profile))
}

# Create 10,000 workers
num_workers <- 9649
worker_ages <- sample(1:7, num_workers, replace = TRUE)
worker_list <- lapply(1:num_workers, function(i) create_worker(i, worker_ages[i]))

# Create queen
queen_row <- as.list(c(
  ID = num_workers + 1, caste = "queen", age = NA, incapability = 0,
  setNames(c(1, rep(0, length(tasks) - 1)), tasks)
))


# Function to create simplified non-worker (gyne/drone)
create_nonworker_individual <- function(id, caste) {
  task_profile <- setNames(rep(0, length(tasks)), tasks)
  # Primary reproduction task for gynes and drones
  task_profile["reproduction_gynes_drones"] <- runif(1, 0.6, 1.0)
  # Minor other activity
  task_profile["self_maintenance"] <- runif(1, 0.2, 0.4)
  task_profile["nest_thermoregulation"] <- runif(1, 0.0, 0.05)
  task_profile["communication_behaviors"] <- runif(1, 0.0, 0.05)
  # Optionally normalize
  total <- sum(task_profile)
  if (total > 0) {
    task_profile <- task_profile / total
  }
  out <- c(ID = id, caste = caste, age = NA, incapability = 0)
  as.list(c(out, task_profile))
}


# Create drones and gynes
num_gynes <- 50
num_drones <- 300
gynes_list <- lapply(1:num_gynes, function(i) create_nonworker_individual(num_workers + 1 + i, "gyne"))
drones_list <- lapply(1:num_drones, function(i) create_nonworker_individual(num_workers + 1 + num_gynes + i, "drone"))

# Flatten all bee lists into a data frame
all_bees <- lapply(c(worker_list, list(queen_row), gynes_list, drones_list), function(x) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  # Transpose if necessary (in case x is a named vector instead of list)
  if (nrow(df) != 1) df <- as.data.frame(t(df), stringsAsFactors = FALSE)
  return(df)
})

# Bind rows together
colony_df <- bind_rows(all_bees)

# Convert columns to correct types
colony_df$ID <- as.integer(colony_df$ID)
colony_df$age <- as.integer(colony_df$age)
colony_df$incapability <- as.numeric(colony_df$incapability)
colony_df$caste <- as.factor(unlist(colony_df$caste))  # Fix for factor conversion

# Ensure task columns are numeric
for (task in tasks) {
  colony_df[[task]] <- as.numeric(colony_df[[task]])
}


# Optional check
str(colony_df)
table(colony_df$caste)


#________________________________________________________________________________________________________
#### 2. Alluvium graph #### 
# First, gather tasks into long format

long_bees <- colony_df %>%
  dplyr::select(ID, caste, all_of(tasks)) %>%
  pivot_longer(
    cols = all_of(tasks),
    names_to = "Task",
    values_to = "Effort"
  )


sampled_ids <- sample(unique(long_bees$ID), 500)

long_subsample <- long_bees %>%
  filter(ID %in% sampled_ids)

# Now plot using axis1 and axis2
ggplot(long_subsample, aes(axis1 = ID, axis2 = Task, y = Effort)) +
  geom_alluvium(aes(fill = Task), width = 0, alpha = 0.7) +
  geom_stratum(width = 0, fill = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Task Allocation Across Bees", x = "Bees and Tasks", y = "Effort")


#________________________________________________________________________________________________________
#### 3. Pie Graph ####

# --- sample a subset of bees for the grid to represent the colony ---
queen_df <- colony_df %>% filter(caste == "queen")
# sample workers gynes and drones proportionally in relation to what is in the colony. 
subset_size <- 1600
n_gynes <- round(subset_size * (sum(colony_df$caste == "gyne") / nrow(colony_df)))
n_drones <- round(subset_size * (sum(colony_df$caste == "drone") / nrow(colony_df)))
n_workers <- subset_size - 1 - n_gynes - n_drones  # subtract queen
gynes_sample <- colony_df %>% filter(caste == "gyne") %>% sample_n(n_gynes)
drones_sample <- colony_df %>% filter(caste == "drone") %>% sample_n(n_drones)
workers_sample <- colony_df %>% filter(caste == "worker") %>% sample_n(n_workers)
bee_sample <- bind_rows(queen_df, gynes_sample, drones_sample, workers_sample)


# --- transform to long format ---
bee_data <- bee_sample %>%
  dplyr::select(ID, all_of(tasks)) %>%
  pivot_longer(
    cols = all_of(tasks),
    names_to = "task",
    values_to = "proportion"
  ) %>%
  arrange(ID, task)


# --- compute pie-slices ---
bee_data <- bee_data %>%
  group_by(ID) %>%
  mutate(
    total = sum(proportion),
    proportion = ifelse(total == 0, 0, proportion / total),  # normalize in case of 0
    start = cumsum(lag(proportion, default = 0)) * 2 * pi,
    end = cumsum(proportion) * 2 * pi,
    r0 = 0,
    r = 0.51
  ) %>%
  ungroup()

# # --- put bees in a grid layout ---
# bee_data <- bee_data %>%
#   mutate(
#     index = dense_rank(ID),
#     x = (index - 1) %% 40,
#     y = (index - 1) %/% 40
#   )


# --- Prepare wide-format matrix for PCA ---
bee_wide <- bee_sample %>%
  dplyr::select(ID, all_of(tasks)) %>%
  column_to_rownames("ID") %>%
  as.matrix()

# Normalize task efforts to proportions (row-wise)
bee_wide <- t(apply(bee_wide, 1, function(x) if (sum(x) == 0) x else x / sum(x)))

# --- PCA ---
bee_pca <- prcomp(bee_wide, center = TRUE, scale. = FALSE)

# Extract PC scores
pca_coords <- as.data.frame(bee_pca$x[, 1:2])  # take first 2 PCs
pca_coords$ID <- as.integer(rownames(pca_coords))

# Rank bees by PC1 and PC2
pca_coords <- pca_coords %>%
  arrange(PC1, PC2) %>%
  mutate(index = row_number())

grid_size <- ceiling(sqrt(nrow(pca_coords)))  # e.g. 20 for 400 bees

pca_coords <- pca_coords %>%
  mutate(
    x = (index - 1) %% grid_size,
    y = (index - 1) %/% grid_size
  )



# # Scale to grid
# pca_coords <- pca_coords %>%
#   mutate(
#     x = scales::rescale(PC1, to = c(0, 39)),
#     y = scales::rescale(PC2, to = c(0, 39))
#   )
# 
# # --- Merge PCA layout back into long-format bee_data ---
# bee_data_clustered <- bee_data %>%
#   left_join(pca_coords, by = "ID")
bee_data_clustered <- bee_data %>%
  left_join(pca_coords, by = "ID")


# --- plot ---

ggplot(bee_data_clustered) +
  geom_arc_bar(aes(
    x0 = x, y0 = y, r0 = r0, r = r,
    start = start, end = end, fill = task
  ),
  color = "white", size = 0.1
  ) +
  scale_fill_viridis_d(option = "D") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Bee Task Allocation: Clustered Grid of Pies")



#________________________________________________________________________________________________________
#### 4. Radar / Spider chart ####

# Step 1: Sum each task across the colony
colony_total <- colony_df %>%
  dplyr::select(all_of(tasks)) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  unlist()

colony_profile <- colony_total / sum(colony_total)

# Top and bottom rows define axis scaling
radar_df <- rbind(
  rep(1, length(tasks)),   # Max value
  rep(0, length(tasks)),   # Min value
  colony_profile           # Actual data
)
colnames(radar_df) <- tasks
radar_df <- as.data.frame(radar_df)

radarchart(radar_df,
           pcol = "#1f77b4",                           # line color
           pfcol = scales::alpha("#1f77b4", 0.4),      # fill color with alpha
           plwd = 2,
           cglcol = "grey", cglty = 1,
           axislabcol = "black", caxislabels = NULL,
           title = "Colony-Level Task Allocation")

scaled_profile <- scales::rescale(colony_profile, to = c(0, 1))
radar_df <- rbind(
  rep(1, length(tasks)),
  rep(0, length(tasks)),
  scaled_profile
) %>% as.data.frame()

radarchart(radar_df,
           pcol = "#1f77b4",                           # line color
           pfcol = scales::alpha("#1f77b4", 0.4),      # fill color with alpha
           plwd = 2,
           cglcol = "grey", cglty = 1,
           axislabcol = "black", caxislabels = NULL,
           title = "Colony-Level Task Allocation")

### spider graph based on optional functioning ###
# --- Create optimal profile (this should sum to 1) ---
optimal_profile <- c(
  reproduction_queen = 0.05,
  reproduction_gynes_drones = 0.1,
  brood_care_nursing = 0.2,
  brood_thermoregulation = 0.05,
  nest_construction_repair = 0.05,
  nest_thermoregulation = 0.05,
  social_immunity = 0.05,
  foraging = 0.2,
  food_processing_storage = 0.05,
  guarding_patrolling_defense = 0.05,
  active_defense = 0.025,
  patrolling_task_allocation = 0.025,
  communication_behaviors = 0.025,
  inactive_reserve = 0.05,
  self_maintenance = 0.05
)

# Ensure order matches the task vector
optimal_profile <- optimal_profile[tasks]

# --- Compute actual profile from colony_df ---
actual_profile <- colony_df %>%
  select(all_of(tasks)) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  unlist()

# Normalize to sum to 1
actual_profile <- actual_profile / sum(actual_profile)
actual_profile <- actual_profile[tasks]  # Ensure same order

# --- Combine for radar chart (fmsb expects top & bottom row for scale) ---
radar_df <- rbind(
  rep(1, length(tasks)),   # max scale
  rep(0, length(tasks)),   # min scale
  optimal_profile,
  actual_profile
)
colnames(radar_df) <- tasks
radar_df <- as.data.frame(radar_df)
rownames(radar_df) <- c("max", "min", "optimal", "actual")

# Your radar_df: rows max, min, optimal, actual
# Extract max, min, optimal, actual rows for transformation
max_vals <- radar_df["max", ]
min_vals <- radar_df["min", ]
optimal_vals <- radar_df["optimal", ]
actual_vals <- radar_df["actual", ]

# Function to rescale so that optimal maps to 0.5 on 0-1 scale
# If actual < optimal, scale between 0 and 0.5
# If actual > optimal, scale between 0.5 and 1
rescale_deviation <- function(actual, optimal, min_val, max_val) {
  # Convert data frames to numeric vectors if needed
  actual <- as.numeric(actual)
  optimal <- as.numeric(optimal)
  min_val <- as.numeric(min_val)
  max_val <- as.numeric(max_val)
  
  rescaled <- numeric(length(actual))
  for (i in seq_along(actual)) {
    if (actual[i] <= optimal[i]) {
      # scale from min to optimal => 0 to 0.5
      rescaled[i] <- 0.5 * scales::rescale(actual[i], to = c(0, 0.5), from = c(min_val[i], optimal[i]))
    } else {
      # scale from optimal to max => 0.5 to 1
      rescaled[i] <- 0.5 + 0.5 * scales::rescale(actual[i], to = c(0.5, 1), from = c(optimal[i], max_val[i]))
    }
  }
  rescaled
}

# Apply transformation to optimal and actual
optimal_scaled <- rescale_deviation(optimal_vals, optimal_vals, min_vals, max_vals)
actual_scaled <- rescale_deviation(actual_vals, optimal_vals, min_vals, max_vals)

# Build new radar data frame with scaled values
radar_scaled <- rbind(
  max = rep(1, length(tasks)),  # max is 1 (top of scale)
  min = rep(0, length(tasks)),  # min is 0
  optimal = optimal_scaled,
  actual = actual_scaled
)

colnames(radar_scaled) <- tasks
radar_scaled <- as.data.frame(radar_scaled)

# Plot radar chart with scaled data
# radarchart(
#   radar_scaled[-c(1,2), ],  # exclude first two rows
#   pcol = c("#1f77b4", "#ff7f0e"),
#   pfcol = c(scales::alpha("#1f77b4", 0.3), scales::alpha("#ff7f0e", 0.4)),
#   plwd = 2,
#   plty = c(1, 1),
#   cglcol = "grey",
#   cglty = 1,
#   axislabcol = "black",
#   vlcex = 0.7,
#   title = "Colony Task Allocation: Deviations from Optimal (0.5)"
# )

# Plot radar chart
radarchart(radar_scaled[c(1,2,4), ],  # plotting only max, min, actual
           pcol = c("red"),          # line color for actual data
           pfcol = scales::alpha("red", 0.3),  # fill color for actual data
           plwd = 2,                 # line width for actual data
           cglcol = "grey",          # grid line color
           cglty = 1,                # grid line type
           axislabcol = "grey",      # axis label color
           caxislabels = seq(0,1,0.25), # axis labels
           vlcex = 0.8)              # variable labels size

radarchart(radar_scaled,
           pcol = c("red", "blue"),   # polygon border colors
           pfcol = scales::alpha(c("red", "blue"), 0.3),  # fill colors with transparency
           plwd = 2,                          # polygon line width
           plty = 1,                          # line type (solid)
           cglcol = "grey",                   # grid color
           cglty = 1,                        # grid line type
           axislabcol = "black",              # axis label color
           vlcex = 0.8                       # variable label size
)


# Add optimal points (optional)
points(radarchart::radar.polygon(radar_scaled["optimal", ], col = NA), col = "blue", pch = 16)


legend(
  "topright", 
  legend = c("Optimal (centered at 0.5)", "Actual"), 
  col = c("#1f77b4", "#ff7f0e"), 
  lty = 1, 
  lwd = 2,
  bty = "n"
)




#________________________________________________________________________________________________________
#### 5. Task specific colony heat maps ####

plot_task_heatmap <- function(task_name, data = bee_data_clustered) {
  data %>%
    filter(task == task_name) %>%
    ggplot(aes(x = x, y = y, fill = proportion)) +
    geom_tile(color = NA) +
    scale_fill_gradient(
      low = "white", 
      high = "darkblue", 
      #name = paste("Proportion\nof", task_name)
    ) +coord_fixed() +
    theme_void() +
    #theme(legend.position = "bottom") +
    theme(
      legend.position = "none",  # this line removes the legend
      plot.title = element_text(hjust = 0.5, size = 10)  # optional: center title
    ) +
    labs(title = paste(task_name))
}

# Example usage:

plot_task_heatmap("nest_construction_repair")


## Get unique tasks from your data
task_list <- unique(bee_data_clustered$task)
# Create all plots
task_plots <- lapply(task_list, function(t) plot_task_heatmap(t))
# Combine plots in a grid
# Choose number of columns, e.g., 4 for 4xN grid
plot_grid <- wrap_plots(task_plots, ncol = 4)
# Display the combined plot
plot_grid




#________________________________________________________________________________________________________
#### 6. Upregulation in Feeding requirements ####

colony_df_mild_stress <- colony_df %>%
  rowwise() %>%
  mutate(
    total_effort = sum(c_across(all_of(tasks)), na.rm = TRUE),
    
    # Set increase factor (e.g., 20% more foraging)
    foraging_new = foraging * 1.1,
    
    # Adjust inactivity to compensate
    delta = foraging_new - foraging,
    inactive_reserve_new = max(inactive_reserve - delta, 0),
    
    # Update foraging and inactivity
    foraging = foraging_new,
    inactive_reserve = inactive_reserve_new,
    
    # Recalculate total and renormalize to original effort
    new_total = sum(c_across(all_of(tasks)), na.rm = TRUE),
    across(all_of(tasks), ~ . * total_effort / new_total)
  ) %>%
  ungroup() %>%
  dplyr::select(-total_effort, -new_total)




colony_df_severe_stress <- colony_df %>%
  rowwise() %>%
  mutate(
    total_effort = sum(c_across(all_of(tasks)), na.rm = TRUE),
    # Larger foraging increase (e.g., double)
    foraging_new = foraging * 2,
    delta = foraging_new - foraging,
    # Try to compensate from inactivity first
    inactive_comp = min(inactive_reserve, delta),
    remaining_delta = delta - inactive_comp,
    # Reduce other tasks proportionally (except foraging and inactive)
    other_tasks_total = sum(c_across(setdiff(tasks, c("foraging", "inactive_reserve")))),
    # If there's remaining delta, reduce other tasks
    scaling_factor = ifelse(remaining_delta > 0 & other_tasks_total > 0,
                            (other_tasks_total - remaining_delta) / other_tasks_total,
                            1),
    # Apply changes
    foraging = foraging_new,
    inactive_reserve = inactive_reserve - inactive_comp,
    across(setdiff(tasks, c("foraging", "inactive_reserve")),
           ~ . * scaling_factor),
    # Renormalize
    new_total = sum(c_across(all_of(tasks)), na.rm = TRUE),
    across(all_of(tasks), ~ . * total_effort / new_total)
  ) %>%
  ungroup() %>%
  dplyr::select(-total_effort, -new_total, -delta, -inactive_comp, -remaining_delta, -other_tasks_total, -scaling_factor)



#________________________________________________________________________________________________________
#### 7. Task Heatmap comparison  ####
# --- Sample bee IDs only once ---
queen_df <- colony_df %>% filter(caste == "queen")
subset_size <- 1600
n_gynes <- round(subset_size * (sum(colony_df$caste == "gyne") / nrow(colony_df)))
n_drones <- round(subset_size * (sum(colony_df$caste == "drone") / nrow(colony_df)))
n_workers <- subset_size - 1 - n_gynes - n_drones

sample_ids <- bind_rows(
  queen_df,
  colony_df %>% filter(caste == "gyne") %>% sample_n(n_gynes),
  colony_df %>% filter(caste == "drone") %>% sample_n(n_drones),
  colony_df %>% filter(caste == "worker") %>% sample_n(n_workers)
) %>% pull(ID)

baseline_sample <- colony_df %>% filter(ID %in% sample_ids)
mild_sample <- colony_df_mild_stress %>% filter(ID %in% sample_ids)
severe_sample <- colony_df_severe_stress %>% filter(ID %in% sample_ids)

bee_wide <- baseline_sample %>%
  select(ID, all_of(tasks)) %>%
  column_to_rownames("ID") %>%
  as.matrix()
bee_wide <- t(apply(bee_wide, 1, function(x) if (sum(x) == 0) x else x / sum(x)))
pca_coords <- as.data.frame(prcomp(bee_wide)$x[, 1:2]) %>%
  mutate(ID = as.integer(rownames(.))) %>%
  arrange(PC1, PC2) %>%
  mutate(index = row_number()) %>%
  mutate(
    x = (index - 1) %% 40,
    y = (index - 1) %/% 40
  )

bee_wide <- baseline_sample %>%
  select(ID, all_of(tasks)) %>%
  column_to_rownames("ID") %>%
  as.matrix()
bee_wide <- t(apply(bee_wide, 1, function(x) if (sum(x) == 0) x else x / sum(x)))
pca_coords <- as.data.frame(prcomp(bee_wide)$x[, 1:2]) %>%
  mutate(ID = as.integer(rownames(.))) %>%
  arrange(PC1, PC2) %>%
  mutate(index = row_number()) %>%
  mutate(
    x = (index - 1) %% 40,
    y = (index - 1) %/% 40
  )

prepare_bee_data <- function(sample_df, label) {
  long <- sample_df %>%
    select(ID, all_of(tasks)) %>%
    pivot_longer(cols = all_of(tasks), names_to = "task", values_to = "proportion") %>%
    group_by(ID) %>%
    mutate(
      total = sum(proportion),
      proportion = ifelse(total == 0, 0, proportion / total),
      start = cumsum(lag(proportion, default = 0)) * 2 * pi,
      end = cumsum(proportion) * 2 * pi,
      r0 = 0,
      r = 0.51
    ) %>%
    ungroup() %>%
    left_join(pca_coords, by = "ID") %>%
    mutate(condition = label)
  return(long)
}

bee_data_baseline <- prepare_bee_data(baseline_sample, "Baseline")
bee_data_mild <- prepare_bee_data(mild_sample, "Mild Stress")
bee_data_severe <- prepare_bee_data(severe_sample, "Severe Stress")

bee_data_all <- bind_rows(bee_data_baseline, bee_data_mild, bee_data_severe)

plot_task_heatmap <- function(task_name, data = bee_data_all) {
  data %>%
    filter(task == task_name) %>%
    ggplot(aes(x = x, y = y, fill = proportion)) +
    geom_tile(color = NA) +
    scale_fill_gradient(low = "white", high = "darkblue") +
    facet_wrap(~condition) +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    labs(title = task_name)
}

# Generate one plot per task across all 3 conditions
task_list <- unique(bee_data_all$task)
task_plots <- lapply(task_list, function(t) plot_task_heatmap(t))

# Combine into a grid
library(patchwork)
wrap_plots(task_plots, ncol = 3)  # 3 = one for each condition per row


#________________________________________________________________________________________________________
#### Delta Task Heatmaps ####

prepare_proportions <- function(sample_df, label) {
  sample_df %>%
    select(ID, all_of(tasks)) %>%
    pivot_longer(cols = all_of(tasks), names_to = "task", values_to = "proportion") %>%
    group_by(ID) %>%
    mutate(total = sum(proportion),
           proportion = ifelse(total == 0, 0, proportion / total)) %>%
    ungroup() %>%
    mutate(condition = label)
}

baseline_long <- prepare_proportions(baseline_sample, "Baseline")
mild_long <- prepare_proportions(mild_sample, "Mild")
severe_long <- prepare_proportions(severe_sample, "Severe")

# Join datasets by ID and task
delta_data <- baseline_long %>%
  rename(proportion_baseline = proportion) %>%
  select(ID, task, proportion_baseline) %>%
  left_join(mild_long %>% rename(proportion_mild = proportion) %>% select(ID, task, proportion_mild), by = c("ID", "task")) %>%
  left_join(severe_long %>% rename(proportion_severe = proportion) %>% select(ID, task, proportion_severe), by = c("ID", "task"))

# Add layout (same PCA coords as before)
delta_data <- delta_data %>% left_join(pca_coords, by = "ID")

# Compute differences
delta_data <- delta_data %>%
  mutate(
    delta_mild = proportion_mild - proportion_baseline,
    delta_severe = proportion_severe - proportion_baseline
  )

plot_task_delta <- function(task_name, delta_column = "delta_mild") {
  delta_data %>%
    filter(task == task_name) %>%
    ggplot(aes(x = x, y = y, fill = .data[[delta_column]])) +
    geom_tile(color = NA) +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue", midpoint = 0,
      limits = c(-max(abs(delta_data[[delta_column]]), na.rm = TRUE),
                 max(abs(delta_data[[delta_column]]), na.rm = TRUE))
    ) +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    labs(title = paste0(task_name, "\n(", gsub("delta_", "Δ ", delta_column), ")"))
}

task_list <- unique(delta_data$task)

# Plot deltas
plots_mild <- lapply(task_list, function(t) plot_task_delta(t, "delta_mild"))
plots_severe <- lapply(task_list, function(t) plot_task_delta(t, "delta_severe"))

# Combine in grids
grid_mild <- wrap_plots(plots_mild, ncol = 4)
grid_severe <- wrap_plots(plots_severe, ncol = 4)

# Display separately
grid_mild + plot_annotation(title = "Task Allocation Changes: Mild Stress vs Baseline")
grid_severe + plot_annotation(title = "Task Allocation Changes: Severe Stress vs Baseline")


####  deltas next to each other ####
delta_long <- delta_data %>%
  select(ID, task, x, y, delta_mild, delta_severe) %>%
  pivot_longer(
    cols = c(delta_mild, delta_severe),
    names_to = "condition",
    values_to = "delta"
  ) %>%
  mutate(condition = recode(condition,
                            delta_mild = "Mild Stress",
                            delta_severe = "Severe Stress"))


ggplot(delta_long, aes(x = x, y = y, fill = delta)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-max(abs(delta_long$delta), na.rm = TRUE),
               max(abs(delta_long$delta), na.rm = TRUE)),
    name = "Δ Proportion"
  ) +
  facet_grid(task ~ condition, switch = "y") +
  coord_fixed() +
  theme_void(base_size = 8) +
  theme(
    legend.position = "right",
    strip.text.y = element_text(angle = 0),
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  labs(
    title = "Task Allocation Changes under Stress (vs. Baseline)",
    subtitle = "Red = reduction; Blue = increase"
  )
