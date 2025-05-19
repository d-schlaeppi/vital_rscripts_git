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
#' 



#________________________________________________________________________________________________________________________________________________
#### 1. Prerequisites ####


### directory

# NATHALIE DIRECTORY
# first_time_use_working_directory <- "~/Dropbox/SeniorLectureship_Bristol/Students_postdocs/Post-Docs/Daniel Schlaeppi/vital_rscripts_git-main/vital_fc1"
# first_time_use_working_directory <- "/media/ael/gismo_hd2/vital/vital_rscripts_git/vital_fc1"

if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") {
  standard <- "/media/ael/gismo_hd2/vital/vital_rscripts_git/vital_fc1" # if you are always working from the same directory just put its name here and it will save you some clicking.  
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
# remove 6 NA lines (artifact from annotations)
dat_duration <- dat_duration[!is.na(dat_duration$feeding_duration), ]
dat_numbers <- read.csv("vital_fc1_data_numbers.csv", header = TRUE) # quick annotation with number of ants present on food sources in 5 min intervalls, only some of the metadata in the table is used. 


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



#________________________________________________________________________________________________________________________________________________
#### 1. - Analysis of feeding - First feeding session fully annotated

#### 1.1 Data preparation ####
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
  if (mapped_feeding_side == dat_numbers$position_virus_corrected[index]) { #define whether feeding event occured on a virus or control food source
    dat_duration$food_source[i] <- "virus"
  } else {
    dat_duration$food_source[i] <- "control"
  }
}

dat_duration <- dat_duration %>%
  filter(feeding_session == 1)



#________________________________________________________________________________________________________________________________________________
#### 2. - Mean duration of feeding events based on food source  

dat_duration %>% 
  group_by(food_source) %>%
  summarize(
    mean_feeding_duration = mean(feeding_duration_seconds, na.rm = TRUE),
    sd_feeding_duration = sd(feeding_duration_seconds, na.rm = TRUE)
  )

# Feeding duration might be longer for the virus food sources.

# capping of feeding events longer than 5 minutes 
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
# Alternative checks suitable for this:
# Over dispersion
overdispersion <- sum(residuals_poisson^2) / df.residual(mod)
print(overdispersion) # should be around 1 (otherwise the model underestimates variability) which it is not --> so the data is overdispersed which needs to be addressed.
#' Possible Solutions:
#' Negative Binomial Regression: Use a negative binomial model instead of a Poisson model, which includes an extra parameter to account for overdispersion.
#' Quasi-Poisson Model: Use a quasi-Poisson model, which adjusts the standard errors to account for overdispersion.
#' Additional Random Effects: Adding more random effects or accounting for hierarchical structures in your data might help.
# Fit a negative binomial GLMM
mod_nb <- glmmTMB(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block),
                  data = dat_duration, family = nbinom2)
summary(mod_nb)
# # Fit a negative binomial GLMM using lme4 and MASS
# mod_nb_lme4 <- glmer.nb(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block),
#                         data = dat_duration_first)
# summary(mod_nb_lme4)



# Check for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
overdisp_fun(mod_nb) # close to 1 and thus the model handles over dispersion adequately



#________________________________________________________________________________________________________________________________________________
#### 2.1  Plotting feeding duration: 4 times the same but with different representations
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
#### 2.1  Feeding duration in dependence of feeding time
# are feeding events different depending on when they happen? Checking for potential effects of food freshness or water content and such on feeding duration
# transform time to seconds
dat_duration$feeding_start_seconds <- period_to_seconds(hms(dat_duration$feeding_start))
dat_duration$feeding_end_seconds <- dat_duration$feeding_start_seconds + dat_duration$feeding_duration_seconds
dat_duration$feeding_end_seconds_capped <- dat_duration$feeding_start_seconds + dat_duration$feeding_duration_seconds_capped

mod <- lmer(log10(feeding_duration_seconds) ~ food_source * feeding_start_seconds + (1|colony_id), data = dat_duration)
Anova(mod)
summary(mod)
test_norm(mod) 

par(mfrow = c(2, 2))
scatter.smooth(fitted(mod), resid(mod)); abline(h = 0, lty = 2)
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main = "normal QQ-plot, residuals") 
qqline(resid(mod))
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))
par(mfrow = c(1, 1))


# Plot of feeding durations over time

plot(feeding_duration_seconds_capped ~ feeding_start_seconds, data = dat_duration,
     main = "Feeding Duration vs Feeding Start Time",
     xlab = "Feeding Start Time (seconds)",
     ylab = "Feeding Duration (seconds)",
     pch = 16,
     col = rgb(t(col2rgb(c("virus" = "#D04C5B", "control" = "#A1D99B")[as.factor(dat_duration$food_source)]))/255, alpha = 0.5),
     cex = 0.8)
abline(lm(feeding_duration_seconds_capped ~ feeding_start_seconds, data = dat_duration), col = "red", lwd = 2)
abline(lm(feeding_duration_seconds_capped ~ feeding_start_seconds, data = dat_duration[dat_duration$food_source == "control", ]), col = "#A1D99B", lwd = 2)  # updated control color
abline(lm(feeding_duration_seconds_capped ~ feeding_start_seconds, data = dat_duration[dat_duration$food_source == "virus", ]), col = "#D04C5B", lwd = 2)  # updated virus color
abline(v = fixef(mod)["(Intercept)"], lty = 4, col = "#A1D99B", lwd = 3)  # updated control color for intercept
abline(v = fixef(mod)["(Intercept)"] + fixef(mod)["food_sourcevirus"], lty = 4, col = "#D04C5B", lwd = 3)  # updated virus color for intercept
legend("topright", legend = levels(as.factor(dat_duration$food_source)), pch = 16,
       col = rgb(t(col2rgb(c("virus" = "#D04C5B", "control" = "#A1D99B")))/255, alpha = 0.5), title = "Food Source")


# plot(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration,
#      main = "Feeding Duration vs Feeding Start Time",
#      xlab = "Feeding Start Time (seconds)",
#      ylab = "Feeding Duration (seconds)",
#      pch = 16,
#      col = rgb(t(col2rgb(viridis::viridis(2)[as.factor(dat_duration$food_source)]))/255, alpha = 0.5),
#      cex = 0.8)
# abline(lm(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration), col = "red", lwd = 2)
# abline(lm(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration[dat_duration$food_source == "control", ]), col = "blue", lwd = 2)
# abline(lm(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration[dat_duration$food_source == "virus", ]), col = "orange", lwd = 2)
# abline(v = fixef(mod)["(Intercept)"], lty = 4, col = "blue", lwd = 3)
# abline(v = fixef(mod)["(Intercept)"] + fixef(mod)["food_sourcevirus"], lty = 4, col = "orange", lwd = 3)
# legend("topright", legend = levels(as.factor(dat_duration$food_source)), pch = 16,
#        col = rgb(t(col2rgb(viridis::viridis(2)))/255, alpha = 0.5), title = "Food Source")






#________________________________________________________________________________________________________________________________________________
#### 2.2 Colony level analysis - Number of feeding events and overall feeding duration ####
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
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #not significant --> good 
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
           label = "p = 0.054 (lmer)", color = "black", size = 3)


# summed up feeding duration per colony and food source 
mod <- lmer(log10(total_feeding_duration) ~ food_source + (1|colony_id) , data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod)
test_norm(mod)
par(mfrow = c(1, 1))


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




