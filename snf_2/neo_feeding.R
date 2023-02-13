### ### ### ### ### ### ### ### ### ### ### ###
### Fluorescin Neonic feeding trial ### ### ###
### ### ### ### ### ### ### ### ### ### ### ###


#### perrequisites ####

rm(list=ls())
# Does neonic contamination affect the sugarwater intake of ants?

suppressPackageStartupMessages(library(lme4)) #lmer
suppressPackageStartupMessages(library(reshape2)) #melt
suppressPackageStartupMessages(library(stringr)) #str_replace()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gld))
suppressPackageStartupMessages(library(PairedData)) #contains paired()
suppressPackageStartupMessages(library(emmeans)) #lsmeans to compare pairwise differences after the modeling using bonferroni corrections
suppressPackageStartupMessages(library(multcompView)) #for the CLD functions
suppressPackageStartupMessages(library(lsmeans)) #cld
suppressPackageStartupMessages(library(multcomp)) #cld 2

setwd("/Users/gismo/Documents/GitHub/snf_2") # homeoffice mac github folder
dat <- read.table("fluorescin_feeding_quantification.txt", header = TRUE)

# fluoresecence corrected for ant debris -> substract values from ant-only wells 
# create new easy to work with data frame
data <- NULL
mean_negative_ant_week1 <- mean(dat$fluorescence_week1[dat$sample_type == "negative_ant"])
mean_negative_ant_week2 <- mean(dat$fluorescence_week2[dat$sample_type == "negative_ant"])
for(i in 1:nrow(dat)) {
  # collect variables
  nr                    <- dat[i, "nr"]
  colony                <- dat[i, "colony"]
  petridish             <- dat[i, "petridish"]
  treatment             <- as.factor(dat[i, "treatment"])
  sample_type           <- dat[i, "sample_type"]
  concentration         <- dat[i, "conc"]
  fluorescence_week1    <- dat[i, "fluorescence_week1"]
  fluorescence_week1_corrected <- if (dat[i, "fluorescence_week1"] >= mean_negative_ant_week1){dat[i, "fluorescence_week1"]-mean_negative_ant_week1} else {0}
  fluorescence_week2    <- dat[i, "fluorescence_week2"]
  fluorescence_week2_corrected <- if (dat[i, "fluorescence_week2"] >= mean_negative_ant_week2){dat[i, "fluorescence_week2"]-mean_negative_ant_week2} else {0}
  data <- rbind(data, data.frame(nr, colony, petridish, treatment, sample_type, concentration, fluorescence_week1_corrected, fluorescence_week2_corrected
  ))
}


#### Questions: is week one significantly differnt from week 2 ####

# modify data so we have a long table with fluorescence for week 1 and week 2 in one column
data_long <- melt(data, id.vars=c("nr", "colony", "petridish", "treatment", "sample_type", "concentration"))
# renaming variables
names(data_long)[names(data_long) == "variable"] <- "week"
names(data_long)[names(data_long) == "value"] <- "fluorescence"
data_long$week <- str_replace(data_long$week, "fluorescence_week1_corrected", "week_1")
data_long$week <- str_replace(data_long$week, "fluorescence_week2_corrected", "week_2")
data_long$week <- as.factor(data_long$week)

# first feeling 
group_by(data_long, week) %>%
  summarise(
    count = n(),
    mean = mean(fluorescence, na.rm = TRUE),
    sd = sd(fluorescence, na.rm = TRUE)
  )

# separate standard curves from data points because they have totally different values... 
# std curves only
data_stdcurves <- subset(data_long, dat$treatment == "std_curve", drop = TRUE)
data_stdcurves$treatment <- droplevels(data_stdcurves$treatment)
# samples only 
data_samples <- subset(data_long, dat$sample_type == "ant")
data_samples$treatment <- droplevels(data_samples$treatment)
data_samples$treatment <- relevel(data_samples$treatment, "low")
data_samples$treatment <- relevel(data_samples$treatment, "control")


# visualisation: 
ggboxplot(data_stdcurves, x = "week", y = "fluorescence", 
          color = "week", palette = c("#00AFBB", "#E7B800"),
          order = c("week_1", "week_2"),
          ylab = "fluorescence", xlab = "time")
ggboxplot(data_samples, x = "week", y = "fluorescence", 
          color = "week", palette = c("#00AFBB", "#E7B800"),
          order = c("week_1", "week_2"),
          ylab = "fluorescence", xlab = "time")

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

# --> week two seems lower... statistical test: paired t test if data is normally distributed.
# compute the difference
d <- with(data_samples, fluorescence[week == "week_1"] - fluorescence[week == "week_2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value ckearly smaller than 0.05 and thus the nullhypothesis of normally distributed data has to be rejected
hist(d)
# statistical test: 
paired_t_test <- t.test(before, after, paired = TRUE)
paired_t_test # significant, but because the data is not nomrally distributed a non parametric test has to be used
# paired wilcoxon signed rank test
wilcox.test(before, after, paired = TRUE) 
# reject null hypothesis that week 1 and week 2 are the same --> less fluorescence in week 2


#### is there a significant difference between treatments regarding food consumption ####
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
# --> high seems to be lower... so lets to the analysis

# Step 1 use standard curves to calculate micro liter of food instead of fluorescence
data_stdcurves$sample_type <- droplevels(data_stdcurves$sample_type)

plot(data_stdcurves$fluorescence ~ data_stdcurves$nr)
abline(0, 0)
plot(data_stdcurves$fluorescence ~ data_stdcurves$concentration)
abline(0, 0)

# remove control control std curve
data_stdcurves <- subset(data_stdcurves, sample_type != "std_curve_no_ant")
plot(data_stdcurves$concentration[data_stdcurves$sample_type != "std_curve_no_ant"] ~ data_stdcurves$fluorescence[data_stdcurves$sample_type != "std_curve_no_ant"],
     xlab = "fluorescence", ylab = "volume of food in ul")
abline(linear_model <- lm(concentration ~ fluorescence, data = data_stdcurves))

# get the coefficients to calculate the volume of food in micro liters
cf <- coef(linear_model)
intercept <- cf[1]
slope <- cf[2]
# calculate consumed food volumes: 
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



