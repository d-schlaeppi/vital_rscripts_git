rm(list = setdiff(ls(), "first_time_use_working_directory"))

# CFG first look into qPCR data

# load required libraries and functions:
pacman::p_load(dplyr, ggplot2, readr, lme4, lmerTest)

# Set working directory
if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") { # direct it to where you have config_user_and_hd.R (typically the script folder or github folder)
  standard <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/CFG" # if you are always working from the same directory just put its name here and it will save you some clicking.  
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
  select(-Comment, -`...9`) %>% as.data.frame()

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
  select(colony_id, treatment)%>% as.data.frame()

cfg <- cfg_subset %>%
  left_join(cfg_colony_metadata_sub, by = "colony_id") %>%
  rename(spore_concentration = `spore_concentration[ng/ul]`)%>% as.data.frame()

# Plots
ggplot(cfg, aes(x = treatment, y = spore_concentration)) +
  stat_summary(fun = mean, geom = "bar", fill = "steelblue", color = "black", width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  ylab("Spore concentration [ng/µl]") +
  xlab("Treatment") +
  theme_minimal()

ggplot(cfg, aes(x = treatment, y = log10(spore_concentration))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  ylab("Spore concentration [ng/µl]") +
  xlab("Treatment") +
  theme_minimal()

ggplot(cfg, aes(x = treatment, y = log10(spore_concentration))) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  ylab("Spore concentration [ng/µl]") +
  xlab("Treatment") +
  theme_minimal()

# first rought model: 


model <- lmer(log10(spore_concentration) ~ treatment + (1 | colony_id), data = cfg)
summary(model)
anova(model)




