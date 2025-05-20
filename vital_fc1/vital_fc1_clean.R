# rm(list = ls())
rm(list = setdiff(ls(), "first_time_use_working_directory"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ###  VITAL - FC1 Analyses   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####
#' This script was written by Nathalie Stroeymeyt and Daniel Schläppi
#' Analysis of the binary food choice experiment conducted prior to the vital tracking experiment.
#' Feeding events of 16 colonies were colonies were could forage on a virus and/or a control food source are being analysed.
#' The data table contains annotations of ant feeding durations (Annotations till 1h until after the second food source has been discovered)

# feeding annotations done visually with start and end being defined as the moment mouth parts connected or disconnected with the food source)

# 
#' INDEX
#' 1. Prerequisites
#' 2. Analysis of binary food choice experiment
#' 2.1 Data Preparation
#' 2.2 Mean duration of feeding events based on food type
#' 2.2.1 Plot
#' 2.3 Feeding duration dependent on feeding time
#' 2.4 Colony level analysis feeding events: Number of feedings and overall feeding duration
#' 2.5 Feeding rate



#________________________________________________________________________________________________________________________________________________
#### 1. Prerequisites ####


### directory

# NATHALIE DIRECTORY
# first_time_use_working_directory <- "~/Dropbox/SeniorLectureship_Bristol/Students_postdocs/Post-Docs/Daniel Schlaeppi/vital_rscripts_git-main/vital_fc1"
# first_time_use_working_directory <- "/media/ael/gismo_hd2/vital/vital_rscripts_git/vital_fc1"

if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") {
  standard <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/vital_fc1" # if you are always working from the same directory just put its name here and it will save you some clicking.  
  selected_dir <- if  (dir.exists(standard)) {standard} else {tcltk::tk_choose.dir(default = "~/", caption = "Select Working Directory")}
  if (is.null(selected_dir) || selected_dir == "") {
    cat("No directory selected. Exiting.\n")
    return()}
  setwd(selected_dir)
  first_time_use_working_directory <<- getwd()
  cat(crayon::blue(getwd()))
} else { setwd(first_time_use_working_directory)
  cat(crayon::blue(getwd())) }


### load data
dat_duration <- read.csv("vital_fc1_data_feedingdurations.csv", header = TRUE)
dat_duration <- dat_duration[!is.na(dat_duration$feeding_duration), ] # remove 6 NA lines (artifact from annotations)
dat_numbers <- read.csv("vital_fc1_data_numbers.csv", header = TRUE) # quick annotation with number of ants present on food sources in 5 min intervals, only some of the metadata in the table is used.


### libraries
# install.packages("pacman")
pacman::p_load(lubridate, plotrix, scales, car, lme4, Hmisc, 
               dplyr, blmeco, lmtest, lsmeans, lubridate,
               emmeans, multcompView, multcomp, viridis, crayon, 
               e1071, DHARMa, merTools, tidyr, pheatmap, grid,
               progress, tidyverse, ggplot2, glmmTMB, patchwork)


### functions
sem <- function(x) {sd(x,na.rm=T)/sqrt(length(na.omit(x)))} # standard error of means
source("func_test_norm.R")  # adds test_norm() to the environment (it takes your model as an argument)
check_overdispersion <- function(model) {
  cat("\033[34m", "🔍 Testing for Overdispersion in Poisson Model", "\033[39m\n\n")
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  p <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  cat("Dispersion ratio:", round(ratio, 2), "\n")
  cat("p-value:", signif(p, 4), "\n\n")
  if (ratio < 1.2) {
    cat("✅ Residuals show no signs of overdispersion\n")
    cat("   --> Model is fine. You can go ahead 👍\n")
  } else {
    cat("❌ Residuals indicate overdispersion\n")
    cat("   --> Consider using a Negative Binomial model or adjusting the model.\n")
  }
}



#________________________________________________________________________________________________________________________________________________
#### 2. - Analysis of feeding - First feeding session fully annotated

#________________________________________________________________________________________________________________________________________________
#### 2.1 Data preparation ####
# get duration in seconds: 
dat_duration <- dat_duration %>%
  mutate(feeding_duration_seconds = as.numeric(as.difftime(feeding_duration, format="%H:%M:%S", units="secs")))

# get correct position of virus and data from other numbers data frame
# create new variable food source saying for each of the feeding events whether it was no a virus of control food source
for (i in 1:nrow(dat_duration)) { # i <- 1
  index <- which(dat_numbers$colony_id == dat_duration$colony_id[i] &
                   dat_numbers$feeding_session == dat_duration$feeding_session[i]) 
  dat_duration$position_virus[i] <- dat_numbers$position_virus_corrected[index]
  dat_duration$treatment[i]      <- dat_numbers$treatment[index]
  if (dat_duration$feeding_side[i] == "r") { # get the side at which the feeding event was observed
    mapped_feeding_side <- "right"
  } else if (dat_duration$feeding_side[i] == "l") {
    mapped_feeding_side <- "left"
  } else {
    mapped_feeding_side <- NA  
  }
  if (mapped_feeding_side == dat_numbers$position_virus_corrected[index]) { #define whether feeding event occurred on a virus or control food source
    dat_duration$food_source[i] <- "virus"
  } else {
    dat_duration$food_source[i] <- "control"
  }
}
#subset to focus only on first feeding session as the second one did not get fully annotated.
dat_duration <- dat_duration %>% filter(feeding_session == 1)



#________________________________________________________________________________________________________________________________________________
#### 2.2 Mean duration of feeding events based on food type ####

# first look at data
dat_duration %>% 
  group_by(food_source) %>%
  summarize(
    mean_feeding_duration = mean(feeding_duration_seconds, na.rm = TRUE),
    sd_feeding_duration = sd(feeding_duration_seconds, na.rm = TRUE), 
    median_feeding_duration = median(feeding_duration_seconds, na.rm = TRUE)) %>% as.data.frame()
# Feeding duration might be longer for the virus food sources.

# capping of feeding events longer than 5 minutes (threshold so far arbitrary, could use some experimental validation) - idea is that after a certain time with head dipped in food, food uptake will no longer increase as technically ants can fell up in like 30 seconds...
cap_threshold <- 5 * 60  # capping threshold in seconds
dat_duration$feeding_duration_seconds_capped <- ifelse( 
  dat_duration$feeding_duration_seconds > cap_threshold,
  cap_threshold,
  dat_duration$feeding_duration_seconds)

Capping_of_feeding_duration <- FALSE
if(Capping_of_feeding_duration) {
  dat_duration$feeding_duration_seconds <- dat_duration$feeding_duration_seconds_capped
}

#
mod <- glmer(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block), data = dat_duration, family = "poisson")
summary(mod)
Anova(mod)
residuals_poisson <- residuals(mod, type = "pearson")
test_norm(mod) # not ok, but not needed for poisson
# Non-Normal Residuals: For Poisson models, residuals are not expected to be normally distributed. The Poisson distribution is skewed, especially for small means.
# Other model checks required here: Over dispersion
check_overdispersion(mod)
# ratio should be around 1 (otherwise the model underestimates variability) which it is not --> so the data is overdispersed
#' Solutions:
#' Negative Binomial Regression: Use a negative binomial model instead of a Poisson model, which includes an extra parameter to account for overdispersion.
#' Quasi-Poisson Model: Use a quasi-Poisson model, which adjusts the standard errors to account for overdispersion.
# --> negative binomial GLMM
mod_nb <- glmmTMB(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block),
                  data = dat_duration, family = nbinom2)
summary(mod_nb)
check_overdispersion(mod_nb) # ratio is close to 1 with a significant p value, indicating mild but significant overdispersion. mod_nb is designed to handle this so no action required unless overdispersion increases (larger than 2)





#________________________________________________________________________________________________________________________________________________
#### 2.2.1 Plotting feeding duration: 4 times the same but with different representations ####
if(TRUE){
  # p1: boxplot
  p1 <- ggplot(dat_duration, aes(x = food_source, y = feeding_duration_seconds_capped, fill = food_source)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.35, color = "black", alpha = 0.1) +
    labs(title = "Boxplot") +
    scale_fill_manual(values = c("virus" = "#D04C5B", "control" = "#A1D99B")) +
    geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(feeding_duration_seconds_capped, na.rm = TRUE), 
                     yend = 1.05 * max(feeding_duration_seconds_capped, na.rm = TRUE)),
                 color = "black", size = 0.25) +
    annotate("text", x = 1.5, y = 1.1 * max(dat_duration$feeding_duration_seconds_capped, na.rm = TRUE), 
             label = "*** (glmer nbinom2)", color = "black", size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # p2: violin plot
  p2 <- ggplot(dat_duration, aes(x = food_source, y = feeding_duration_seconds_capped, fill = food_source)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_jitter(width = 0.35, color = "black", alpha = 0.1) +
    labs(title = "Violin Plot") +
    scale_fill_manual(values = c("virus" = "#D04C5B", "control" = "#A1D99B")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Summary for barplots
  summary_data <- dat_duration %>%
    group_by(food_source) %>%
    summarize(
      mean_duration = mean(feeding_duration_seconds_capped, na.rm = TRUE),
      se_duration = sd(feeding_duration_seconds_capped, na.rm = TRUE) / sqrt(n()),
      sd_duration = sd(feeding_duration_seconds_capped, na.rm = TRUE)
    )
  
  # p3: barplot with SE
  p3 <- ggplot(summary_data, aes(x = food_source, y = mean_duration, fill = food_source)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_duration - se_duration, ymax = mean_duration + se_duration),
                  width = 0.2, color = "black") +
    labs(title = "Mean ± SE") +
    scale_fill_manual(values = c("virus" = "#D04C5B", "control" = "#A1D99B")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # p4: barplot with SD and raw data
  p4 <- ggplot(summary_data, aes(x = food_source, y = mean_duration, fill = food_source)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_duration - sd_duration, ymax = mean_duration + sd_duration),
                  width = 0.2, color = "black") +
    geom_jitter(data = dat_duration,
                aes(x = food_source, y = feeding_duration_seconds_capped),
                width = 0.15, color = "black", alpha = 0.1, inherit.aes = FALSE) +
    labs(title = "Mean ± SD + Raw Data") +
    scale_fill_manual(values = c("virus" = "#D04C5B", "control" = "#A1D99B")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Combine all four
  (p1 | p2) / (p3 | p4) + 
    plot_annotation(title = "Feeding Duration by Food Source (Different Visualizations)")
}







#________________________________________________________________________________________________________________________________________________
#### 2.3 Feeding duration dependent on feeding time ####
#' Ants are feeding longer on virus food. Considering that virus is often discovered first (see below), 
#' it is worth checking whether feeding time and with that freshness of food affects feeding duration
#' Checking for potential effects of food freshness or water content and such on feeding duration: Correlation
# transform time to seconds
dat_duration$feeding_start_seconds <- period_to_seconds(hms(dat_duration$feeding_start))
dat_duration$feeding_end_seconds <- dat_duration$feeding_start_seconds + dat_duration$feeding_duration_seconds
dat_duration$feeding_end_seconds_capped <- dat_duration$feeding_start_seconds + dat_duration$feeding_duration_seconds_capped

### Correlation
shapiro.test(dat_duration$feeding_duration_seconds)
shapiro.test(dat_duration$feeding_start_seconds)
# Pearson correlation (assumes normal distribution of both variables)
cor.test(dat_duration$feeding_duration_seconds, dat_duration$feeding_start_seconds, method = "spearman")
# significant correlation between feeding duration and feeding start time
# however, correlation coefficient (rho = 0.066) is very small --> indicating a very weak positive relationship
plot(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration,
     main = "Feeding Duration vs Feeding Start Time",
     xlab = "Feeding Start Time (seconds)",
     ylab = "Feeding Duration (seconds)",
     pch = 16, col = "#69b3a2")
abline(lm(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration), col = "red", lwd = 2)


mod <- lmer(log10(feeding_duration_seconds) ~ food_source * feeding_start_seconds + (1|colony_id), data = dat_duration)
Anova(mod)
summary(mod)
test_norm(mod) 


# plot feeding duration over time capped and non capped. 
mean_start_times <- tapply(dat_duration$feeding_start_seconds, dat_duration$food_source, mean, na.rm = TRUE) # mean feeding start times per food_source
duration_vars <- c("feeding_duration_seconds", "feeding_duration_seconds_capped")
plots <- list()
for (var in duration_vars) {
    p <- ggplot(dat_duration, aes(x = feeding_start_seconds, y = .data[[var]], color = food_source)) +
    geom_point(alpha = 0.5, size = 2) +
    scale_color_manual(values = c("virus" = "#D04C5B", "control" = "#A1D99B")) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_smooth(data = subset(dat_duration, food_source == "control"), method = "lm", se = FALSE, color = "#A1D99B") +
    geom_smooth(data = subset(dat_duration, food_source == "virus"), method = "lm", se = FALSE, color = "#D04C5B") +
    geom_vline(xintercept = mean_start_times["control"], linetype = "dashed", color = "#A1D99B", size = 1) +
    geom_vline(xintercept = mean_start_times["virus"], linetype = "dashed", color = "#D04C5B", size = 1) +
    labs(
      title = paste(gsub("_", " ", var), "vs Feeding Start Time"),
      x = "Feeding Start Time (seconds)",
      y = gsub("_", " ", var)
    ) +
    theme_minimal() +
    theme(legend.title = element_text(face = "bold"))
  plots[[var]] <- p
}
combined_plot <- plots[["feeding_duration_seconds"]] | plots[["feeding_duration_seconds_capped"]]
print(combined_plot)

#'Feeding duration is slightly influenced by feeding time
#'Food type as well as feeding time affects feeding duration. Might be a proxy for food uptake and could indicate preference...?







#________________________________________________________________________________________________________________________________________________
#### 2.4 Colony level analysis - Number of feeding events and overall feeding duration ####
#________________________________________________________________________________________________________________________________________________
#### 2.4.1 Number of feeding events and overall feeding duration (total annotation) ####


dat_summary <- dat_duration %>%
  group_by(food_source, colony_id) %>%
  summarise(
    num_feeding_events = n(),
    total_feeding_duration = sum(feeding_duration_seconds, na.rm = TRUE)
  ) %>% as.data.frame()

dat_summar_plot <- dat_summary %>%
  group_by(food_source) %>%
  summarise(
    mean_nr_feeding_events = mean(num_feeding_events, na.rm = TRUE),
    sd_nr_feeding_events = sd(num_feeding_events, na.rm = TRUE),
    se_nr_feeding_events = sd(num_feeding_events, na.rm = TRUE) / sqrt(n()),
    mean_tot_feeding_duration = mean(total_feeding_duration, na.rm = TRUE),
    sd_tot_feeding_duration = sd(total_feeding_duration, na.rm = TRUE),
    se_tot_feeding_duration = sd(total_feeding_duration, na.rm = TRUE) / sqrt(n())
  ) %>% as.data.frame()


# number of feeding events per colony per food source 
mod <- lmer(num_feeding_events ~ food_source + (1|colony_id) , data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow = c(1,1))
test_norm(mod)

ggplot(dat_summary, aes(x = food_source, y = num_feeding_events, fill = food_source)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, color = "black", alpha = 0.2) +
  labs(title = "Number of Feeding Events",
       x = "Food Source",
       y = "Number of Feeding Events") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(num_feeding_events, na.rm = TRUE), 
                   yend = 1.05 * max(num_feeding_events, na.rm = TRUE)),
               color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 1.1 * max(dat_summary$num_feeding_events, na.rm = TRUE), 
           label = "p = 0.046 (lmer)", color = "black", size = 3)




# summed up feeding duration per colony and food source 
mod <- lmer(log10(total_feeding_duration) ~ food_source + (1|colony_id) , data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow = c(1, 1))
test_norm(mod)

ggplot(dat_summary, aes(x = food_source, y = total_feeding_duration, fill = food_source)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, color = "black", alpha = 0.2) +
  labs(title = "Summed up Duration of Feeding Events per Colony",
       x = "Food Source",
       y = "Duration of Feeding Events in Seconds") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(total_feeding_duration, na.rm = TRUE), 
                   yend = 1.05 * max(total_feeding_duration, na.rm = TRUE)),
               color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 1.1 * max(dat_summary$total_feeding_duration, na.rm = TRUE), 
           label = "p = 0.0008 (lmer)", color = "black", size = 3)

#' Summary
#' Feeding duration a bit longer for virus food sources
#' More feeding events on virus food sources at colony level
#' --> Leading to more overall exploitation (amount of time feeding on each food source) at colony level!
#' This should also be reflected in the feeding rate






#________________________________________________________________________________________________________________________________________________
#### 2.4.2 Same but with 60m time window with shifted T0 ####

### Repeat of  colony level analysis once more but this time capping to 1h of annotations otherwise there is a bias towards the earlier discovered food
#' Note: Annotations were not done for the complete 2h of recordings per colony.
#' Instead, annotations were done until 1h after the second food was discovered.
#' So if for example one food was discovered min 10 and one at min 30 we would get 60 min of annotations until min 90 for the food discovered last
#' However, for the other food that was discovered earlier we would have 90 minutes of annotations.
#' This would create a bias if one food source was consistently discovered earlier/later, which now we know was the case.
#' --> by looking at only 60 minutes with shifted_discovery we take away some of this bias.
#' However there is still the natural bias towards the first discovered food that requires disentangling (see further below)
#' e.g. pheromone trails and effects of previous feedings (reduced hunger etc) on subsequent feedings


#' --> look at a number of feedings as well as the rate of feeding events during 1h post discovery of each food so rate per min in 60 minutes after a food has been discovered

# subset data to only include 60 minutes post discovery of for each food source
dat_duration$feeding_annotation_end_seconds <- period_to_seconds(hms(dat_duration$feeding_annotation_end))
summary_time_feeding <- dat_duration %>% group_by(colony_id, food_source) %>%
  summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
subsetted_data <- dat_duration %>% left_join(summary_time_feeding, by = c("colony_id", "food_source"))
# filter to only include feeding events starting within the 60 min (3600s) window since discovery (first actual feeding) of the food source of interest (different discovery times)
subsetted_data <- subsetted_data %>%
  mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
  mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>% 
  filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)


# Summary statistics by food source and colony
dat_summary <- subsetted_data %>%
  group_by(food_source, colony_id) %>%
  summarise(
    num_feeding_events = n(),
    total_feeding_duration = sum(feeding_duration_seconds, na.rm = TRUE)
  ) %>% as.data.frame()

# Summary for plotting
dat_summar_plot <- dat_summary %>%
  group_by(food_source) %>%
  summarise(
    mean_nr_feeding_events = mean(num_feeding_events, na.rm = TRUE),
    sd_nr_feeding_events = sd(num_feeding_events, na.rm = TRUE),
    se_nr_feeding_events = sd(num_feeding_events, na.rm = TRUE) / sqrt(n()),
    mean_tot_feeding_duration = mean(total_feeding_duration, na.rm = TRUE),
    sd_tot_feeding_duration = sd(total_feeding_duration, na.rm = TRUE),
    se_tot_feeding_duration = sd(total_feeding_duration, na.rm = TRUE) / sqrt(n())
  ) %>% as.data.frame()

# Model: number of feeding events per colony per food source
mod <- lmer(log(num_feeding_events) ~ food_source + (1 | colony_id), data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow = c(1, 1))
test_norm(mod) #is fine with log but else if not normally distributed --> we could try with a glmm with poisson as it is count data. 
# mod_glmm <- glmer(num_feeding_events ~ food_source + (1 | colony_id), 
#                   data = dat_summary, family = poisson)
# check_overdispersion(mod_glmm) # overdispersed so we would change to negative binomial
# #negative binomial instead
# mod_nb <- glmmTMB(num_feeding_events ~ food_source + (1 | colony_id), data = dat_summary, family = nbinom2)
# summary(mod_nb)
# car::Anova(mod_nb)
# sim_res_nb <- simulateResiduals(fittedModel = mod_nb)
# plot(sim_res_nb)


# Plot: number of feeding events
ggplot(dat_summary, aes(x = food_source, y = num_feeding_events, fill = food_source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, color = "black", alpha = 0.2) +
  labs(title = "Number of Feeding Events",
       x = "Food Source",
       y = "Number of Feeding Events") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(num_feeding_events, na.rm = TRUE), 
                   yend = 1.05 * max(num_feeding_events, na.rm = TRUE)),
               color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 1.1 * max(dat_summary$num_feeding_events, na.rm = TRUE), 
           label = "p = 0.08 (lmer)", color = "black", size = 3)

# Model: summed feeding duration
mod <- lmer(log10(total_feeding_duration) ~ food_source + (1 | colony_id), data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow = c(1, 1))
test_norm(mod)

# Plot: total feeding duration
ggplot(dat_summary, aes(x = food_source, y = total_feeding_duration, fill = food_source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, color = "black", alpha = 0.2) +
  labs(title = "Summed up Duration of Feeding Events per Colony",
       x = "Food Source",
       y = "Duration of Feeding Events in Seconds") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(total_feeding_duration, na.rm = TRUE), 
                   yend = 1.05 * max(total_feeding_duration, na.rm = TRUE)),
               color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 1.1 * max(dat_summary$total_feeding_duration, na.rm = TRUE), 
           label = "p = 0.02 (lmer)", color = "black", size = 3)

#' Looking at it this way, now the number of feeding events is just a trend and no longer significant 
#' --> Hmm... despite virus on average being discovered first.
#' Overall duration (exploitation based on number of events and duration) is still significant.





#________________________________________________________________________________________________________________________________________________
#### 2.5 Feeding rate (feeding events per min) ####
# first look at feeding rate simply calculated as average number of feeding events per minute, which per definition will just look like overall number of feeding events / 60...
dat_summary <- subsetted_data %>%
  group_by(food_source, colony_id) %>%
  summarise(
    num_feeding_events = n(),
    total_feeding_duration = sum(feeding_duration_seconds, na.rm = TRUE),
    num_feeding_events_per_min = n()/60
  ) %>% as.data.frame()
dat_summary %>%
  group_by(food_source) %>% 
  summarise(
    mean_feeding_rate = mean(num_feeding_events_per_min),
    sd_feeding_rate   = sd(num_feeding_events_per_min)) %>% as.data.frame()
boxplot(num_feeding_events_per_min ~food_source, data = dat_summary, main = "feeding rate")

# stats
mod <- lmer(log(num_feeding_events_per_min) ~ food_source + (1|colony_id), data = dat_summary) #sqrt transformation might also work
summary(mod)
print(Anova(mod, type=3))
test_norm(mod)
compareqqnorm(mod); par(mfrow=c(1,1))

# plot
ggplot(dat_summary, aes(x = food_source, y = num_feeding_events_per_min, fill = food_source))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, color = "black", alpha = 0.2) +
  labs(title = "Feeding rate", 
       x = "Food source", 
       y = "Feeding rate (feedings/minute") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(num_feeding_events_per_min, na.rm = TRUE), 
                   yend = 1.05 * max(num_feeding_events_per_min, na.rm = TRUE)),
               color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 1.1 * max(dat_summary$num_feeding_events_per_min, na.rm = TRUE), 
           label = "p = 0.11 (lmer)", color = "black", size = 3)



#________________________________________________________________________________________________________________________________________________
#### 3. First feeding #### 

#### 3.1 Time of first feeding ####
first_feeding <- dat_duration %>% group_by(colony_id, food_source) %>% 
  summarise(time_first_feeding = as.numeric(min(feeding_start_seconds, na.rm = TRUE))) %>% as.data.frame()

# mean time until first feeding per food source
mod <- lmer(log10(time_first_feeding) ~ food_source + (1|colony_id), data = first_feeding)
Anova(mod)
sim_res <- simulateResiduals(mod)
plot(sim_res)

# on average the virus food is feasted upon earlier than the control food source

#### 3.1 Number of first discoveries per food source ####
food_first_eaten_from <- first_feeding %>% group_by(colony_id) %>% 
  arrange(colony_id, time_first_feeding) %>%
  slice(1) %>% dplyr::select(colony_id, food_first_eaten_from = food_source, time_earliest_feeding = time_first_feeding) %>% as.data.frame()
n_virus_eaten_first <- sum(food_first_eaten_from$food_first_eaten_from == "virus")
n_max <- nrow(food_first_eaten_from)
binom.test(n_virus_eaten_first, n_max, p = 0.5)
prop.test(n_virus_eaten_first, n_max, alternative = "two.sided", p = 0.5)
t2 <- table(food_first_eaten_from$food_first_eaten_from)
barplot(t2, main = "", xlab = "food source", ylab = "number of first feedings", ylim = c(0,14), xaxt="n")
axis(1, at=c(0.7, 1.95), labels=c("control", "virus"))
segments(x0 = 0.7, y0 =12.5, x1 = 1.9)
text(x = 1.3, y = 13, label = "p = 0.08 (prop.test)", font = 2, cex = 1.1)


#' there might be a trend but it is just not significant
#' But still: it is crucial to disentangle effect of first discovery from a potential preference...
#' randomization might be an approach but because it is imbalanced between samples it is not ideal.
#' instead we can look at a difference in exploitation rate based on which food was discovered first in a two step approach:
#' Check if there is a difference in exploitation rate between the first discovered food source and the second.
#' If so, check if this difference is different between the treatments.
#' --> Create a graph in which we can see the difference in exploitation rates between the two treatments 
#' depending on which of the two food sources was discovered first, and if there is a difference analyse it properly.







#________________________________________________________________________________________________________________________________________________
#### 4. Exploitation of Food sources over time #### 
# plot a time curve showing the exploitation of the two food sources
# plot a the rate of exploitation over time 
# --> create a full grid with seconds, food_source and colonies and add in the the original data:


#### 4.1 Prep data ####
dat_duration_modified <- dat_duration %>%
  group_by(colony_id) %>%
  mutate(
    first_feeding_virus     = min(feeding_start_seconds[food_source == "virus"], na.rm = TRUE),
    first_feeding_control   = min(feeding_start_seconds[food_source == "control"], na.rm = TRUE),
    first_source_fed_on = ifelse(
      first_feeding_virus < first_feeding_control, "virus",
      ifelse(first_feeding_control < first_feeding_virus, "control", NA_character_))) %>% 
  ungroup() %>% as.data.frame()
# get colony metadata
colony_metadata <- dat_duration_modified %>%
  dplyr::select(colony_id, treatment, position_virus, block, 
                first_feeding_virus, first_feeding_control, first_source_fed_on) %>%
  distinct(colony_id, .keep_all = TRUE)


#### 4.2 Main Analysis ####
#' could look things from different perspectives:
#' start_experiment - would look at the entire annotations. Meaning from discovery of first food until 60 min after second food was discovered
#' discovery_first_food - looking at what happens in the first 60 minutes after the discovery of the first food.
#' shifted_discovery - looking at a 60 minute window post discovery for each food source separately.
#' All would have their merits. However, as we are trying to disentangle effects of first discovery from potential preference shifted_discovery is best 
#' The first two options would create a strong bias because the data shows that the virus food on average is discovered first
#' By looking at 1h for each food separately we at least overcome the issue the the first food discovered will "be over-represented"

start_time <- "shifted_discovery"
#for (start_time in c("shifted_discovery")) { # "start_experiment", "discovery_first_food", if reintroducing them make sure to keep shifted for last...# start_time <- "shifted_discovery"

  #### 4.2.1 subset and prepare data ####
  
    # if (start_time == "start_experiment") {subsetted_data <- dat_duration_modified}
  if (start_time == "shifted_discovery") {
    summary_time_first_feeding <- dat_duration_modified %>% group_by(colony_id, food_source) %>%
      summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
    subsetted_data <- dat_duration_modified %>% left_join(summary_time_first_feeding, by = c("colony_id", "food_source"))
    # filter to only include feeding events starting within the 60 min (3600s) window since discovery (first actual feeding) of the food source of interest (different discovery times)
    subsetted_data <- subsetted_data %>%
      mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
      mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>%
      filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)  }
  # if (start_time == "discovery_first_food") {
  #   summary_time_first_feeding <- dat_duration_modified %>% group_by(colony_id) %>%
  #     summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
  #   subsetted_data <- dat_duration_modified %>% left_join(summary_time_first_feeding, by = c("colony_id"))
  #   # filter to only include first hour since discovery of first food. 
  #   subsetted_data <- subsetted_data %>%
  #     mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
  #     mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>%
  #     filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)  }
  
  #create a variable indicating whether a feeding happens on the food source that was discovered first or second
  subsetted_data <- subsetted_data %>% 
    mutate(food_source_timed = if_else(first_source_fed_on == food_source, "first", "second"))
  
  
  #### 4.2.2 Expand data to get a continuous exploitation rate per food source (real or timed) ####
  # Analyse data based on food type and based on first discovery: 
  # loop over the two grouping variables: 
  
  for (var in c("food_source", "food_source_timed")) { #var = "food_source"
    grid <- expand_grid(
      colony_id = unique(subsetted_data$colony_id), 
      food_source = unique(subsetted_data[[var]]),
      # time = seq(0, max(subsetted_data$feeding_end_seconds, na.rm = TRUE))
      time = seq(0, 3600) # hard cut at 1h to avoid slow trailing off of feeding events that started within 1h but run beyond it.
    ) %>% rename(!!var := food_source) %>% as.data.frame()
    cumulative_explotation_over_time <- NULL
    
    for (colony in unique(grid$colony_id)) { 
      for (group_value in unique(grid[[var]])) {
        # Filter grid and feeding data for the current colony and variable value
        current_grid <- grid %>% filter(colony_id == colony, !!sym(var) == group_value)
        current_feedings <- subsetted_data %>% filter(colony_id == colony, !!sym(var) == group_value)
        if (nrow(current_feedings) > 0) {
          current_feedings <- current_feedings %>% arrange(feeding_start_seconds)
          current_grid <- current_grid %>%
            rowwise() %>%
            mutate(
              cumulated_exploitation_time = sum(
                pmin(current_feedings$feeding_duration_seconds, 
                     pmax(0, time - current_feedings$feeding_start_seconds))
              ),
              exploitation_rate = sum(
                time >= current_feedings$feeding_start_seconds &
                  time <= current_feedings$feeding_end_seconds
              )
            )
        } else {
          current_grid <- current_grid %>%
            mutate(cumulated_exploitation_time = 0,
                   exploitation_rate = 0)
        }
        cumulative_explotation_over_time <- bind_rows(cumulative_explotation_over_time, current_grid)
      }
    }
    cumulative_explotation_over_time <- cumulative_explotation_over_time %>% left_join(colony_metadata, by = "colony_id")
    assign(paste0("cumulative_explotation_", var), cumulative_explotation_over_time)
  }

  
  
  
  #### 4.2.3 plot mean exploitation over time per food source or food discovery time ####
  # aggregate data and calculate means of duration for plotting of exploitation rate per food source
  # for (dataset in c("cumulative_explotation_food_source", "cumulative_explotation_food_source_timed")) {
  # data <- get(dataset)
  
  ### based on food type
  mean_exploitation_time <- cumulative_explotation_food_source %>%
    group_by(time, food_source) %>%
    summarize(mean_time = mean(cumulated_exploitation_time, na.rm = TRUE),
              sd = sd(cumulated_exploitation_time, na.rm = TRUE), 
              sem = sem(cumulated_exploitation_time))
  title_lab_1 <- paste0("Food based - Cum Exploit (T0 = ", start_time, ")")
  print(ggplot(mean_exploitation_time, aes(x = time, y = mean_time, color = food_source)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = mean_time - sem, ymax = mean_time + sem, fill = food_source), alpha = 0.2) +
        labs(title = title_lab_1,
        x = "Time",
        y = "Mean Cumulated Exploitation Time (seconds)") +
        scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
        theme_bw() +
        theme(legend.title = element_blank())+
        theme(legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()))
  
  ### based on discovery time 
  mean_exploitation_time <- cumulative_explotation_food_source_timed %>%
    group_by(time, food_source_timed) %>%
    summarize(mean_time = mean(cumulated_exploitation_time, na.rm = TRUE),
              sd = sd(cumulated_exploitation_time, na.rm = TRUE), 
              sem = sem(cumulated_exploitation_time))
  title_lab_1 <- paste0("Time based - Cum Exploitat (T0 = ", start_time, ")")
  print(ggplot(mean_exploitation_time, aes(x = time, y = mean_time, color = food_source_timed)) +
          geom_line(size = 1) +
          geom_ribbon(aes(ymin = mean_time - sem, ymax = mean_time + sem, fill = food_source_timed), alpha = 0.2) +
          labs(title = title_lab_1,
               x = "Time",
               y = "Mean Cumulated Exploitation Time (seconds)") +
          scale_color_manual(values = c("first" = "#FBB4AE", "second" = "#CCEBC5")) +
          theme_bw() +
          theme(legend.title = element_blank())+
          theme(legend.title = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
  #### 4.2.4 Plot exploitation rate over time per food source or food discovery time ####

  ### based on food type
  mean_exploitation_rate <- cumulative_explotation_food_source %>%
    group_by(time, food_source) %>%
    summarize(mean_exploitation = mean(exploitation_rate, na.rm = TRUE),
              sd = sd(exploitation_rate, na.rm = TRUE))
  title_lab_2 <- paste0("Food-based Exploitation Rate (T0 = ", start_time, ")")
  print(ggplot(mean_exploitation_rate, aes(x = time, y = mean_exploitation, color = food_source)) +
      geom_line(size = 0.5) +
      labs(title = title_lab_2,
        x = "Time",
        y = "Mean Cumulated Exploitation Time (seconds)") +
      scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
      theme_bw() +
      theme(legend.title = element_blank())+
      theme(legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()))
  
  # Optional consider adding a bar or box plot with just the overall mean (see previous versions of script)
  
  ### based on discovery time   
  mean_exploitation_rate <- cumulative_explotation_food_source_timed %>%
    group_by(time, food_source_timed) %>%
    summarize(mean_exploitation = mean(exploitation_rate, na.rm = TRUE),
              sd = sd(exploitation_rate, na.rm = TRUE))
  title_lab_2 <- paste0("Time-based Exploitation Rate (T0 = ", start_time, ")")
  print(ggplot(mean_exploitation_rate, aes(x = time, y = mean_exploitation, color = food_source_timed)) +
          geom_line(size = 0.5) +
          labs(title = title_lab_2,
               x = "Time",
               y = "Mean Cumulated Exploitation Time (seconds)") +
          scale_color_manual(values = c("first" = "#FBB4AE", "second" = "#CCEBC5")) +
          theme_bw() +
          theme(legend.title = element_blank())+
          theme(legend.title = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
  
  
  
  
  #### 4.2.5 Stats - Mean exploitation rate per food source ####
  summary_dat_duration <- cumulative_explotation_food_source %>% 
    group_by(colony_id, food_source) %>% 
    summarise(mean_exploitation_rate = mean(exploitation_rate)) %>% 
    left_join(colony_metadata, by = "colony_id") %>% as.data.frame()
  #stats
  print( paste("Statistics - T0 = ",start_time))
  summary_dat_duration$first_source_fed_on_binary <- summary_dat_duration$food_source==summary_dat_duration$first_source_fed_on
  
  # model <- lmer(log10(mean_exploitation_rate) ~ food_source + (1|colony_id) + (1|first_source_fed_on), data=summary_dat_duration )
  # print(Anova(model))
  # test_norm(model) 
  model <- lmer((mean_exploitation_rate) ~ food_source*first_source_fed_on_binary + (1|colony_id) , data=summary_dat_duration)
  test_norm(model)
  print(Anova(model)) # interaction non-significant --> it will be left away
  model <- lmer((mean_exploitation_rate) ~ food_source+first_source_fed_on_binary + (1|colony_id) , data=summary_dat_duration )
  test_norm(model)
  print(Anova(model)) # food source is not significant while time of discovery is 
  # based on this we would say that it is only first discovery and not food source that drives the preference.


  #### 4.2.6 PLOT: Delta exploitation rate based on first discovery ####
  # pivot data to wide format so virus and control exploitation rates are side-by-side and calculate 

  data_delta <- cumulative_explotation_food_source %>%
    dplyr::select(colony_id, time, food_source, exploitation_rate, cumulated_exploitation_time, first_source_fed_on) %>%
    pivot_wider(
      names_from = food_source,
      values_from = c(exploitation_rate, cumulated_exploitation_time)
    ) %>%
    mutate(delta_exploitation_rate = if_else(first_source_fed_on == "virus",
             exploitation_rate_virus - exploitation_rate_control,
             exploitation_rate_control - exploitation_rate_virus),
      delta_cumulated_exploitation_time = if_else(first_source_fed_on == "virus",
        cumulated_exploitation_time_virus - cumulated_exploitation_time_control,
        cumulated_exploitation_time_control - cumulated_exploitation_time_virus),
      first_discovered = first_source_fed_on) %>% as.data.frame()
  
  
  

  ggplot(data_delta, aes(x = time, y = delta_exploitation_rate, color = first_discovered, fill = first_discovered)) +
    stat_summary(fun = mean, geom = "line", size = 0.5, alpha = 0.5) + # raw
    # stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.5, color = NA) + # raw
    # geom_smooth(method = "loess", se = TRUE, size = 0.8, span = 0.2) + # slow when using a lot of data. 
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, size = 0.8) + #faster than method loess
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    scale_color_manual(values = c("virus" = "#E41A1C", "control" = "#4DAF4A")) +
    scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
    labs(title = "Delta Exploitation Rate",
         x = "Time (s)",
         y = "Δ Exploitation Rate (First - Second)",
         color = "Discovered first",
         fill = "Discovered first") +
    theme_minimal(base_size = 13) +
    theme(
      panel.border = element_rect(fill = NA),
      axis.line = element_line(),
      legend.position = "right")

 #### Stats delta exploitation based on food source ####
  
  ### Colony mean
  summary_exploitation_delta <- data_delta %>% 
    group_by(colony_id, first_discovered) %>% 
    summarise(mean_delta_exploitation_rate = mean(delta_exploitation_rate)) %>% as.data.frame()
  
  summary_exploitation_delta %>% 
    group_by(first_discovered) %>% 
    summarise(mean_delta = mean(mean_delta_exploitation_rate, na.rm = TRUE),
              se_delta = sd(mean_delta_exploitation_rate, na.rm = TRUE) / sqrt(n())) %>% as.data.frame()
  
  model <- lm((mean_delta_exploitation_rate) ~ first_discovered  , data=summary_exploitation_delta )
  test_norm(model)
  print(Anova(model))
  
  # not significant... 
  
  ### no mean: 
  # cumulative_explotation_food_source$first_source_fed_on_binary <- cumulative_explotation_food_source$food_source==cumulative_explotation_food_source$first_source_fed_on
  # cumulative_explotation_over_time_first  <- cumulative_explotation_food_source[which(cumulative_explotation_food_source$first_source_fed_on_binary),c("colony_id","time","cumulated_exploitation_time","exploitation_rate","first_source_fed_on")]
  # names(cumulative_explotation_over_time_first)[names(cumulative_explotation_over_time_first)%in%c("cumulated_exploitation_time","exploitation_rate")] <- paste(  names(cumulative_explotation_over_time_first)[names(cumulative_explotation_over_time_first)%in%c("cumulated_exploitation_time","exploitation_rate")] ,"_first_source_discovered",sep="")
  # cumulative_explotation_over_time_second <- cumulative_explotation_food_source[which(!cumulative_explotation_food_source$first_source_fed_on_binary),c("colony_id","time","cumulated_exploitation_time","exploitation_rate","first_source_fed_on")]
  # names(cumulative_explotation_over_time_second)[names(cumulative_explotation_over_time_second)%in%c("cumulated_exploitation_time","exploitation_rate")] <- paste(  names(cumulative_explotation_over_time_second)[names(cumulative_explotation_over_time_second)%in%c("cumulated_exploitation_time","exploitation_rate")] ,"_second_source_discovered",sep="")
  # cumulative_explotation_over_time_differential <- merge(cumulative_explotation_over_time_first,cumulative_explotation_over_time_second)
  # cumulative_explotation_over_time_differential$exploitation_rate_differential <- cumulative_explotation_over_time_differential$exploitation_rate_first_source_discovered - cumulative_explotation_over_time_differential$exploitation_rate_second_source_discovered
  # cumulative_explotation_over_time_differential$cumulated_exploitation_time_differential <- cumulative_explotation_over_time_differential$cumulated_exploitation_time_first_source_discovered - cumulative_explotation_over_time_differential$cumulated_exploitation_time_second_source_discovered
  # aggregate(exploitation_rate_differential~first_source_fed_on,function(x)cbind(mean(x),std.error(x)),data=cumulative_explotation_over_time_differential)

  exploitation_delta %>%
    group_by(first_source_fed_on) %>%
    summarise(
      mean_exploitation_rate_diff = mean(delta_exploitation_rate, na.rm = TRUE),
      se_exploitation_rate_diff = sd(delta_exploitation_rate, na.rm = TRUE) / sqrt(n()),
      .groups = "drop") %>% as.data.frame()
  

  model <- lmer((delta_exploitation_rate) ~ first_source_fed_on + (1|colony_id)  , data=data_delta)
  kurtosis(residuals(model))
  skewness(residuals(model))
  print(Anova(model)) #same as colony summary? 
  
  model <- lmer((delta_cumulated_exploitation_time) ~ first_source_fed_on*time + (1|colony_id) , data=data_delta)
  kurtosis(residuals(model))
  skewness(residuals(model))
  summary(model)
  print(Anova(model))
  pvalue <- round(Anova(model)$`Pr(>Chisq)`,2)
  print(shapiro.test(residuals(model))) # for shifted t0 not significant... 
  test_norm(model)
  
  #### CONTINUE HERE - JUST ASK GPT WHAT THE BEST WAY OF DOING THE STATS IS.
  
# }























#________________________________________________________________________________________________________________________________________________
#### X. End - What is next #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
cat(yellow("Move on to the next script: bead_data_analysis.R"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
