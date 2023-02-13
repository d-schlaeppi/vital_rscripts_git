#### Neonic Fungus Co-exposure SNF 2 ####

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###                                                                                             ###
###   ### ### ### ###         ###     ###            ### ### ###  ###     ###     #####         ###
###   ### ### ### ###         ###     ###            ### ### ###  ###     ###   ###   ###       ###
###   ###         ###         ###     ###            ###          ###     ###  ###     ###      ###
###   ###         ###         ###     ###            ###          ###     ###  ###              ###
###   ### ###     ###         ###     ###   ### ###  ### ###      ###     ###  ###              ###
###   ### ###     ###         ###     ###   ### ###  ### ###      ###     ###  ###   #####      ###
###   ###         ###         ###     ###            ###          ###     ###  ###   #####      ###
###   ###         ###         ###     ###            ###          ###     ###  ###     ###      ###
###   ###         ### ### ###  ###   ###             ###           ###   ###    ###   ###       ###
###   ###         ### ### ###    #####               ###             #####        #####         ###
###                                                                                             ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Background information | Read me ####

# This will be the script containing all the analyses from the Flupyradifurone Fungus (Metharhizium brunneum) trial 

# this script (will) contain all the code for the second SNF proposal on fungus and neonic coexposure of ants

# 1. prerequisites
# 2. Fluorescin Feeding trial
# 3. Neonic only survival analysis
# 4. Flupy Funguy combined exposure survival 

#### 1. prerequisites ####
rm(list=ls())
#load required libraries functions and code
install.packages("tydir")
library(survival)
library(coxme)
library(scales)
library(lme4) #lmer
library(reshape2) #melt
library(stringr) #str_replace()
library(dplyr)
library(ggpubr)
library(gld)
library(PairedData) #contains paired()
library(emmeans) #lsmeans to compare pairwise differences after the modeling using bonferroni corrections
library(multcompView) #for the CLD functions
library(lsmeans) #cld
library(multcomp) #cld 2
library(rlang)
library(ggplot2)
library(viridis)
library(broom)


# set general working directory and get data
directory <- "/Users/gismo/Documents/GitHub/snf_2/"
setwd(directory) # homeoffice mac github folder


#### 2. Fluorescin Feeding trial ####
files <- list.files()
print(files)
dat <- read.table("full_fluorescin_feeding_quantification.txt", header = TRUE)
head(dat)

# fluorescence corrected for ant debris -> subtract values from ant-only wells
# create new easy to work with data frame with fluorescence only

mean_negative_ant_week1_run1 <- mean(dat$fluorescence_week1[dat$sample_type == "negative_ant" & dat$run == 1])
mean_negative_ant_week2_run1 <- mean(dat$fluorescence_week2[dat$sample_type == "negative_ant" & dat$run == 1])
mean_negative_ant_week1_run2  <- mean(dat$fluorescence_week1[dat$sample_type == "negative_ant" & dat$run == 2])
data <- NULL
for(i in 1:nrow(dat)) {
  # collect variables
  run                   <- as.factor(dat[i, "run"])
  nr                    <- dat[i, "nr"]
  colony                <- dat[i, "colony"]
  petridish             <- dat[i, "petridish"]
  treatment             <- as.factor(dat[i, "treatment"])
  sample_type           <- dat[i, "sample_type"]
  concentration         <- dat[i, "conc"]
  fluorescence_week1    <- dat[i, "fluorescence_week1"]
  fluorescence_week1_corrected <- ifelse(dat[i, "run"] == 1 & dat[i, "fluorescence_week1"] >= mean_negative_ant_week1_run1,
                                         dat[i, "fluorescence_week1"] - mean_negative_ant_week1_run1,
                                         ifelse(dat[i, "run"] == 2 & dat[i, "fluorescence_week1"] >= mean_negative_ant_week1_run2,
                                                dat[i, "fluorescence_week1"] - mean_negative_ant_week1_run2,
                                                0))
  fluorescence_week2    <- dat[i, "fluorescence_week2"]
  fluorescence_week2_corrected <- ifelse(dat[i, "run"] == 1 & !is.na(dat[i, "fluorescence_week2"]) & dat[i, "fluorescence_week2"] >= mean_negative_ant_week2_run1,
                                         dat[i, "fluorescence_week2"] - mean_negative_ant_week2_run1,
                                                ifelse(is.na(dat[i, "fluorescence_week2"]),
                                                       NA,
                                                       0))
   data <- rbind(data, data.frame(run, nr, colony, petridish, treatment, sample_type, concentration, fluorescence_week1_corrected, fluorescence_week2_corrected
  ))
}


#### 2.1 Fluorescin decay ####
# Is week one significantly different from week 2 #
# To test if there is a decay of the fluorescin signal if the samples are in the freezer at -80°C for a week we compare week 1 values with week2 values (repeated measures on the same samples)
# modify data so we have a long table with fluorescence for week 1 and week 2 in one column
data_long <- melt(data, id.vars=c("run", "nr", "colony", "petridish", "treatment", "sample_type", "concentration"))
# renaming variables
names(data_long)[names(data_long) == "variable"] <- "week"
names(data_long)[names(data_long) == "value"] <- "fluorescence"
data_long$week <- str_replace(data_long$week, "fluorescence_week1_corrected", "week_1")
data_long$week <- str_replace(data_long$week, "fluorescence_week2_corrected", "week_2")
data_long$week <- as.factor(data_long$week)

# first look 
group_by(data_long, week, run) %>%
  summarise(
    count = n(),
    mean = mean(fluorescence, na.rm = TRUE),
    sd = sd(fluorescence, na.rm = TRUE)
  )

# separate standard curves from data points because they have totally different values and check how they look
# std curves only
data_stdcurves <- subset(data_long, dat$treatment == "std_curve", drop = TRUE)
data_stdcurves$treatment <- droplevels(data_stdcurves$treatment)
# samples only 
data_samples <- subset(data_long, dat$sample_type == "ant")
data_samples$treatment <- droplevels(data_samples$treatment)

# reorder the levels of the treatment variable
data_samples$treatment <- droplevels(data_samples$treatment)
data_samples$treatment <- relevel(data_samples$treatment, "mid")
data_samples$treatment <- relevel(data_samples$treatment, "low")
data_samples$treatment <- relevel(data_samples$treatment, "control")




### visualisations ### 
#without loosing the paired information:
# Subset data week 1 and week 2
before <- subset(data_samples,  week == "week_1", fluorescence, drop = TRUE)
after <- subset(data_samples,  week == "week_2", fluorescence, drop = TRUE)
# Plot paired data
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()
ggpaired(data_samples, x = "week", y = "fluorescence", xlab = "time", ylab = "fluorescence",
         color = "week", line.color = "gray", line.size = 0.4,
         palette = "jco")
#only the first round: 
data_samples_run1 <- data_samples[data_samples$run == 1, ]
ggpaired(data_samples_run1, x = "week", y = "fluorescence", xlab = "time", ylab = "fluorescence",
         color = "week", line.color = "gray", line.size = 0.4,
         palette = "jco")
# --> week two seems lower... statistical test: paired t test if data is normally distributed.



# compute the difference
d <- with(data_samples, fluorescence[week == "week_1"] - fluorescence[week == "week_2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value ckearly smaller than 0.05 and thus the nullhypothesis of normally distributed data has to be rejected
hist(d)
# statistical test: 
paired_t_test <- t.test(before, after, paired = TRUE)
paired_t_test # significant, but because the data is not normally distributed a non parametric test has to be used
# paired Wilcox signed rank test
wilcox.test(before, after, paired = TRUE) 
# reject null hypothesis that week 1 and week 2 are the same --> less fluorescence in week 2










#### 2.2 Food consumption ####
#is there a significant difference between treatments regarding food consumption #
# we no longer need week 2 --> create a subset of week1 only for both samples and standard curves
data_samples <- subset(data_samples, week == "week_1")
data_stdcurves <- subset(data_stdcurves, week == "week_1")

# preliminary data exploration: 
group_by(data_samples, treatment) %>%
  summarise(
    count = n(),
    mean = mean(fluorescence, na.rm = TRUE),
    sd = sd(fluorescence, na.rm = TRUE)
  )
boxplot(fluorescence ~ treatment, data = data_samples)
# looks like high is lower and mid might be slightly lower as well...
# to make sure it is not an artifact due to differences in the two runs we need to 
# correct the two runs based on their respecitve standard curves for the calculations of consumed food volumes
# Also make the comparison between the full data set only vs. just round 2. 

# Step 1 use standard curves to calculate micro liter of food instead of fluorescence
# Plot the standard curves
ggplot(data_stdcurves, aes(x = concentration, y = fluorescence, color = sample_type)) + 
  geom_line() + 
  xlab("Concentration") + 
  ylab("Fluorescence") + 
  ggtitle("Standard Curves") + 
  scale_color_discrete(name = "Sample Type")
# -> standard curves for the second run are flatter --> first and second run need different coefficients for food calculation
# remove no ant standard curve, then calculate a mean std curve for the first and the second round and derive their coefficients for food consumption

data_stdcurves <- subset(data_stdcurves, sample_type != "std_curve_no_ant")
data_stdcurves$sample_type <- droplevels(as.factor(data_stdcurves$sample_type))

# Plot the pooled standard curves for each of the two runs
data_mean <- data_stdcurves %>% 
  group_by(concentration, run) %>% 
  summarise(mean_fluorescence = mean(fluorescence))

ggplot(data_mean, aes(x = concentration, y = mean_fluorescence, color = as.factor(run))) + 
  geom_point(size = 3) + 
  geom_smooth(aes(group = run), method = "lm", se = FALSE, size = 0.5) + 
  ggtitle("Standard Curves with Linear Model Fit") + 
  xlab("Concentration [μl fluorescin / 100 μl sugar water]") + 
  ylab("Mean Fluorescence") + 
  scale_color_viridis_d(name = "Run")



# calculation of how much food the ants consumed based on the standard curves (forcing std curves through zero)
data_mean_lm_0 <- data_mean %>% 
  group_by(run) %>% 
  do(mod = lm(mean_fluorescence ~ 0 + concentration, data = .))
coef_df_0 <- data_mean_lm_0 %>% 
  do(tidy(.$mod)) %>% 
  bind_rows()

data_samples$consumed_volume_0 <- NA
for (i in 1:nrow(data_samples)) {
  if(data_samples[i, "run"] == 1) {
    Slope = coef_df_0$estimate[1]
  } else {Intercept = coef_df_0$estimate[2]}
  data_samples[i, "consumed_volume_0"] <- data_samples[i, "fluorescence"] / Slope
}


# calculation of how much food the ants consumed based on the standard curves (but NOT forcing std curves through zero) -> intercept != 0
data_mean_lm_1 <- data_mean %>% 
  group_by(run) %>% 
  do(mod = lm(mean_fluorescence ~ concentration, data = .))
coef_df_1 <- data_mean_lm_1 %>% 
  do(tidy(.$mod)) %>% 
  bind_rows()

data_samples$consumed_volume_1 <- NA
for (i in 1:nrow(data_samples)) {
  if(data_samples[i, "run"] == 1) {
    Intercept = coef_df_1$estimate[1]
    Slope = coef_df_1$estimate[2]
  } else {
    Intercept = coef_df_1$estimate[3]
    Slope = coef_df_1$estimate[4]}
  data_samples[i, "consumed_volume_1"] <- (data_samples[i, "fluorescence"] - Intercept) / Slope
}

#comparison
head(data_samples)
data_samples$consumed_volume_0
data_samples$consumed_volume_1
boxplot(data_samples$consumed_volume_0 ~ data_samples$treatment)
boxplot(data_samples$consumed_volume_1 ~ data_samples$treatment)

data_run_1 <- subset(data_samples, run == 1)
data_run_2 <- subset(data_samples, run == 2)
boxplot(data_samples$consumed_volume_0 ~ data_samples$treatment)
boxplot(data_run_1$consumed_volume_0~data_run_1$treatment)
boxplot(data_run_2$consumed_volume_0~data_run_2$treatment)

# nice graph based on 0 intercept model 











# volume = slope*fluorescence + intercept
data_samples$consumed_volume <- data_samples$fluorescence*slope + intercept
boxplot(data_samples$consumed_volume ~ data_samples$treatment, xlab = "treatment", ylab = "sugarwater consumption [ul]")
#abline(median(data_samples$consumed_volume[data_samples$treatment == "high"]), 0)


#### next step do the statistical test on this data ####

mod <- lmer(log(consumed_volume) ~ treatment + (1|colony) + (1|petridish), data = data_samples, REML = FALSE)
summary(mod)
Anova(mod)
lsmeans(mod, pairwise ~ treatment, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment, data = data_samples)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

hist(log(data_samples$consumed_volume))

#check model assumptions
compareqqnorm(mod)
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
plot(mod)
par(mfrow=c(1,1))
#homogenity of variance # plot(mod, 1)
leveneTest(residuals(mod) ~ data$caste[!is.na(data$virus_titre)]) #non-significant = ok 
boxplot(residuals(mod) ~ data$caste[!is.na(data$virus_titre)])
#normality  # plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> ok$
# so far the high treatment is appears significantly lower than the controls and the low ones. But the model is not right (residuals indicate that model assumptions violated)
# --> next transform fluorescence values to estimates of the volumes consumed by the ants. 






#### 3. Neonic only survival analysis ####
dat_neo <- read.table("neo_survival.txt", fill=TRUE, header = TRUE)
head(dat_neo)



#### 4. Flupy Funguy combined exposure survival ####









#### Neonic only survival analysis ####

