rm(list = setdiff(ls(), "first_time_use_working_directory"))
#rm(list = ls())
# CFG first look into qPCR data

# load required libraries and functions:
pacman::p_load(dplyr, ggplot2, readr, lme4, lmerTest, patchwork, car)

# Set working directory
if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") { # direct it to where you have config_user_and_hd.R (typically the script folder or github folder)
  #standard <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/CFG" # if you are always working from the same directory just put its name here and it will save you some clicking.  
  standard <- "/media/ael/gismo_hd2/vital/vital_rscripts_git/CFG" # if you are always working from the same directory just put its name here and it will save you some clicking.  
  selected_dir <- if  (dir.exists(standard)) {standard} else {tcltk::tk_choose.dir(default = "~/", caption = "Select Working Directory")}
  if (is.null(selected_dir) || selected_dir == "") {
    cat("No directory selected. Exiting.\n")
    return()}
  setwd(selected_dir)
  first_time_use_working_directory <<- getwd()
  cat(crayon::blue(getwd()))
} else { setwd(first_time_use_working_directory)
  cat(crayon::blue(getwd())) }

list.files()

### get all relevant data into one frame... improvised as currently everything is all over the place: treatment, colonies analysed, spore load, identity
cfg_frozen_ants <- as.data.frame(read_csv("cfg_frozen_ants.csv"))
cfg_qpcr_all    <- read_csv("CFG_qPCR_all.csv") %>% as.data.frame()
cfg_colony_metadata    <- read_csv("CFG_colony_metadata.csv") %>% as.data.frame()

names(cfg_qpcr_all)
names(cfg_frozen_ants)


cfg_merged <- cfg_frozen_ants %>%
  left_join(cfg_qpcr_all, by = c("ID" = "sample_id")) %>% 
  dplyr::select(-Comment, -`...9`) %>% as.data.frame()

analysed <- cfg_merged %>%
  group_by(Colony) %>%
  summarize(n_qpcr_samples = sum(!is.na(`spore_concentration[ng/ul]`))) %>% as.data.frame()

colonies_analysed <- cfg_merged %>%
  group_by(Colony) %>%
  summarize(n_qpcr_samples = sum(!is.na(`spore_concentration[ng/ul]`))) %>%
  filter(n_qpcr_samples > 0, Colony != "c09") %>%
  pull(Colony)


cfg_subset <- cfg_merged %>%
  filter(Colony %in% colonies_analysed)%>%
  rename(colony_id = Colony) %>% as.data.frame()

cfg_colony_metadata_sub <- cfg_colony_metadata %>% 
  dplyr::select(colony_id, treatment)%>% as.data.frame()

cfg <- cfg_subset %>%
  left_join(cfg_colony_metadata_sub, by = "colony_id") %>%
  rename(spore_concentration = `spore_concentration[ng/ul]`)%>% as.data.frame()

cfg <- cfg %>%
  mutate(spore_concentration = ifelse(is.na(spore_concentration), 0, spore_concentration),
         high_exposure = spore_concentration > 4.32e-5)

# Plots
p1 <- ggplot(cfg, aes(x = treatment, y = spore_concentration)) +
  stat_summary(fun = mean, geom = "bar", fill = "steelblue", color = "black", width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  ylab("Spore concentration [ng/µl]") +
  xlab("Treatment") +
  ggtitle("Mean ± SE") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p2 <- ggplot(cfg, aes(x = treatment, y = log10(spore_concentration))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.1) +
  ylab("Spore concentration [ng/µl]") +
  xlab("Treatment") +
  ggtitle("Boxplot (log-transformed)") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p3 <- ggplot(cfg, aes(x = treatment, y = log10(spore_concentration))) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.1) +
  ylab("Spore concentration [ng/µl]") +
  xlab("Treatment") +
  ggtitle("Violin Plot (log-transformed)") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

(p1 | p2 | p3) + 
  plot_annotation(title = "Spore Concentration by Treatment (Multiple Visualizations)")

# first rough model: 
mini <- 0.5*min(cfg$spore_concentration[cfg$spore_concentration > 0], na.rm = TRUE)
model <- lmer(log10(spore_concentration+mini) ~ treatment + (1 | colony_id), data = cfg)
summary(model)
anova(model)




#### high vs. low threshold ####
# Define threshold from the paper
threshold <- 4.32e-5 # based on nathalies input - ask for citation

# Add exposure category
cfg <- cfg %>%  mutate(exposure_level = if_else(spore_concentration > threshold, "high", "low"))
# Count and proportion of high/low per treatment
exposure_summary <- cfg %>%
  group_by(treatment, exposure_level) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(treatment) %>%
  mutate(proportion = count / sum(count))

ggplot(exposure_summary, aes(x = treatment, y = proportion, fill = exposure_level)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("low" = "#CCEBC5", "high" = "#FBB4AE")) +
  labs(
    title = "Proportion of Ants by Exposure Level",
    x = "Treatment",
    y = "Proportion of Ants",
    fill = "Exposure Level"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Create contingency table
exposure_table <- table(cfg$treatment, cfg$exposure_level)
# Chi-square test
chisq.test(exposure_table)
# If expected frequencies are low, use Fisher's exact test
fisher.test(exposure_table)
cfg$exposure_binary <- if_else(cfg$exposure_level == "high", 1, 0)
mod <- glm(exposure_binary ~ treatment, data = cfg, family = binomial)
summary(mod)


### binary variable
cfg <- cfg %>% mutate(high_exposure = spore_concentration > threshold)
colony_summary <- cfg %>%
  group_by(colony_id, treatment) %>%
  summarise(
    n_total = n(),
    n_high = sum(high_exposure),
    prop_high = n_high / n_total,
    .groups = "drop"
  )

ggplot(colony_summary, aes(x = treatment, y = prop_high, fill = treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("control" = "#CCEBC5", "flupy" = "#FBB4AE")) +
  labs(
    title = "Proportion of High Exposure Ants per Colony",
    x = "Treatment",
    y = "Proportion > Threshold"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")


mod <- glmer(high_exposure ~ treatment + (1 | colony_id), data = cfg, family = binomial)
summary(mod)
Anova(mod)


# transform ng/ug to actual spore count?
# test out a couple of different thresholds-


thresholds <- c(0.1 * 4.32e-5, 4.32e-5, 10 * 4.32e-5)

# Loop through each threshold
for (threshold in thresholds) { #threshold <- thresholds[1]
  message("Running analysis for threshold: ", threshold)
  # Add exposure category
  cfg <- cfg %>% mutate(exposure_level = if_else(spore_concentration > threshold, "high", "low"))
  
  # Count and proportion of high/low per treatment
  exposure_summary <- cfg %>%
    group_by(treatment, exposure_level) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(treatment) %>%
    mutate(proportion = count / sum(count))
  
  print(ggplot(exposure_summary, aes(x = treatment, y = proportion, fill = exposure_level)) +
          geom_col(position = "fill") +
          scale_y_continuous(labels = scales::percent_format()) +
          scale_fill_manual(values = c("low" = "#CCEBC5", "high" = "#FBB4AE")) +
          labs(
            title = paste("Proportion of Ants by Exposure Level (threshold =", format(threshold, scientific = TRUE), ")"),
            x = "Treatment",
            y = "Proportion of Ants",
            fill = "Exposure Level"
          ) +
          theme_bw() +
          theme(panel.grid = element_blank()))
  
  # Statistical tests
  exposure_table <- table(cfg$treatment, cfg$exposure_level)
  print(chisq.test(exposure_table))
  print(fisher.test(exposure_table))
  
  # Logistic regression
  cfg$exposure_binary <- if_else(cfg$exposure_level == "high", 1, 0)
  mod <- glm(exposure_binary ~ treatment, data = cfg, family = binomial)
  print(summary(mod))
  
  # Summary per colony
  cfg <- cfg %>% mutate(high_exposure = spore_concentration > threshold)
  colony_summary <- cfg %>%
    group_by(colony_id, treatment) %>%
    summarise(
      n_total = n(),
      n_high = sum(high_exposure),
      prop_high = n_high / n_total,
      .groups = "drop"
    )
  
  print(ggplot(colony_summary, aes(x = treatment, y = prop_high, fill = treatment)) +
          geom_boxplot(width = 0.6, outlier.shape = NA) +
          geom_jitter(width = 0.15, alpha = 0.7, color = "black") +
          scale_fill_manual(values = c("control" = "#CCEBC5", "flupy" = "#FBB4AE")) +
          labs(
            title = paste("Proportion of High Exposure Ants per Colony (threshold =", format(threshold, scientific = TRUE), ")"),
            x = "Treatment",
            y = "Proportion > Threshold"
          ) +
          theme_bw() +
          theme(panel.grid = element_blank(), legend.position = "none"))
  
  # Mixed-effects logistic regression
  mod2 <- glmer(high_exposure ~ treatment + (1 | colony_id), data = cfg, family = binomial)
  print(summary(mod2))
  print(Anova(mod2))
}
