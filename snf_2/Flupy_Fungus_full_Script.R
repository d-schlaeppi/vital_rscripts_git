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

# this script (will) contain all the code for the second SNF project on fungus and neonic coexposure of ants

# 1. prerequisites
# 2. Fluorescin Feeding trial
# 3. Neonic only survival analysis
# 4. Flupy Funguy combined exposure survival 

#### 1.  prerequisites ####
rm(list=ls())
#load required libraries functions and code
{
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
library(car) # Anova()
library(blmeco) #compareqqnorm()
library(survminer) # used in the analysis of the survival curves
}


#### 2.  Fluorescin Feeding trial ####

# set general working directory and get data
directory <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/snf_2/"
setwd(directory) # homeoffice mac github folder
source('printme_coxme.R') # used in the analysis of the survival curves

# load data
files <- list.files()
print(files)
dat <- read.table("full_fluorescin_feeding_quantification.txt", header = TRUE)
head(dat)




###
# fluorescence corrected for ant debris -> subtract values from ant-only wells
# create new easy to work with data frame with fluorescence only

mean_negative_ant_week1_run1 <- mean(dat$fluorescence_week1[dat$sample_type == "negative_ant" & dat$run == 1])
mean_negative_ant_week2_run1 <- mean(dat$fluorescence_week2[dat$sample_type == "negative_ant" & dat$run == 1])
mean_negative_ant_week1_run2  <- mean(dat$fluorescence_week1[dat$sample_type == "negative_ant" & dat$run == 2])

#### 2.05 create an easy to work with data frame with corrected fluorescence values ####
data <- NULL
for(i in 1:nrow(dat)) {
  # collect variables
  run                   <- as.factor(dat[i, "run"])
  nr                    <- dat[i, "nr"]
  colony                <- as.factor(dat[i, "colony"])
  petridish             <- as.factor(dat[i, "petridish"])
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
# Is week one significantly different from week 2 (only available for the first run) #
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


############################################ CONTINUE HERE TO DO THE ANALYSIS #######################





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
# correct the two runs based on their respective standard curves for the calculations of consumed food volumes
# Also make the comparison between the full data set only vs. just round 2. 

# Step 1 use standard curves to calculate micro liter of food instead of fluorescence
# Plot the standard curves
ggplot(data_stdcurves, aes(x = concentration, y = fluorescence, color = sample_type)) + 
  geom_line() + 
  xlab("Concentration") + 
  ylab("Fluorescence") + 
  ggtitle("Standard Curves") + 
  scale_color_viridis_d(name = "Run")

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


#### 2.3 Food consumption graph ####
# nice graph based on 0 intercept model 

ggplot(data_samples, aes(x = treatment, y = consumed_volume_0)) +
  geom_boxplot(fill = alpha(viridis(4)[4], 0.5), color = "black", notch = TRUE) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, color = viridis(1)[1]) +
  labs(x = "Treatment", y = "Sugarwater Consumption [ul]", title = "Sugarwater Consumption by Treatment Group")+
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))

# separately for each of the two runs: 

data_samples$trt_run <- interaction(data_samples$treatment, data_samples$run)
levels_order <- c("control.1", "control.2", "low.1", "low.2", "mid.1", "mid.2", "high.1", "high.2")
data_samples$trt_run <- factor(data_samples$trt_run, levels = levels_order)

ggplot(data_samples, aes(x = trt_run, y = consumed_volume_0, fill = factor(run))) +
  geom_boxplot(color = "black", notch = TRUE, alpha = 0.5) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, color = viridis(1)[1]) +
  labs(x = "Treatment and Run", y = "Sugarwater Consumption [ul]", title = "Sugarwater Consumption by Treatment and Run")+
  scale_fill_manual(values = viridis(2)) +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 



#### 2.4 Fluoprescin Stats ####
hist(data_samples$consumed_volume_0)
# data is clearly not normally distributed. Log transformations did not work to get to a distribution suitable for lmer. 
# As the data is negative exponential distribution we use a glmer negative binomial model with thetha = 1 which coresponds to a negative exp distribution

mod <- glmer.nb(consumed_volume_0 ~ treatment + (1|colony) + (1|petridish) + (1|run), data = data_samples)
summary(mod)
Anova(mod)
cld(emmeans(mod, pairwise ~ "treatment", adjust = "tukey"), Letters = letters)
# a trend, but no significant differences



#### 3.  Neonic only survival analysis ####

exp1_data <- read.csv('flupy_survival_test.csv')
exp1_data$concentration <- factor(exp1_data$concentration )

#### 3.1 Stats neonic only ####
null_model <- coxme ( Surv (time = survival, event = censor) ~ 1                 + ( 1 | petri_dish) , data = exp1_data)
full_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + ( 1 | petri_dish) , data = exp1_data)
anova(null_model   ,  full_model )

summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH"))
letters <- cld(summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH")))



#### 3.1 Plot neonic only ####
surviplot1 <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration , data=exp1_data) #basic plot
legend <- c('0     - a', '0.5  - ab','1     - ab','5     - a','10   - ab','50   - ab','100  - b','500  - c')

surv_plot <- ggsurvplot(surviplot1, data = exp1_data,
                        censor = FALSE,
                        legend.title = 'concentration [ppm]',
                        legend.labs = legend, 
                        legend = c(0.15,0.25),
                        xlab = 'Time (days)',
                        ylab = 'Proportion Surviving',
                        break.time.by=2,
                        xlim = c(0,21),
                        ggtheme = theme_bw(),
                        palette = viridis(8, begin = 0, end = 0.8, option = 5),
                        conf.int= TRUE,
                        conf.int.alpha = 0.1
)
# remove the grid lines
surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank())
surv_plot

aggregate(  censor ~ concentration, FUN=mean,data=exp1_data)




#### 4.  Flupy Fungus combined exposure effects on individual survival ####
# prepare data
#reading data
antdata2 <- read.csv('flupy_fungus_individual_survival_data.csv')
str(antdata2)
antdata2$fungus <- as.factor(antdata2$fungus)
antdata2$petri_dish <- as.factor(antdata2$petri_dish)
antdata2$block <- as.factor(antdata2$block)

# creating a colony column, random factor must be controlled
antdata2 <- within(antdata2,colony <- substr(petri_dish,2,2))
antdata2$colony <- as.factor(antdata2$colony)

# making concentration a proper number, NOT factor
antdata2$concentration <- as.numeric(as.character(antdata2$concentration))

#### 4.1 Stats FluFu Indi ####

# full interaction model, compared to a model without interaction
flupy_fungus_interaction_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration  + fungus + concentration:fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
flupy_fungus_model             <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration  + fungus                        + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
fungus_model                   <- coxme ( Surv (time = survival, event = censor) ~ 1 + fungus + ( 1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
null_model                     <- coxme ( Surv (time = survival, event = censor) ~ 1                 + ( 1 | petri_dish) + (1 | colony) + (1 | block) , data = antdata2)
anova(null_model, fungus_model, flupy_fungus_model, flupy_fungus_interaction_model)
# clear effect of fungus and no effect of (explicitly sub lethal) Flupy concentration
# trend for interactive effect




#### 4.2 Flu Fu Indi plots #####
surviplot <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration + fungus, data=antdata2)
surv_plot2 <- ggsurvplot(surviplot, data = antdata2, 
           linetype = "fungus", 
           color = "concentration", break.time.by=2,
           palette = viridis(3, begin = 0, end = 0.9, option = 1),
           xlab = 'Time (days)', ylab = 'Proportion Surviving',
           xlim = c(0, max(antdata2$survival)),
           censor = FALSE,
           panel.labs = list(fungus = c("FUNGUS", "SHAM")), 
           legend.title = 'FPF Concentration [ppm]', legend.labs = c(), legend = c(0.17,0.25),
           ggtheme = theme_bw(), 
           conf.int = FALSE, 
           conf.int.alpha = 0.1
           )
surv_plot2$plot <- surv_plot2$plot + theme(panel.grid = element_blank()) +
                                     scale_linetype_discrete(name = "Fungus", labels = c("M. brunneum", "Sham"))
surv_plot2

theme <- theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour = "black", fill="NA", linewidth=1),
                     panel.background = element_blank())

surv_plot3 <- ggsurvplot(surviplot, data = antdata2,
           legend.title = 'FPF Concentration [ppm]', legend.labs = c(), legend = c(0.15,0.15),
           xlab = 'Time (days)', ylab = 'Proportion Surviving',
           palette = viridis(3, begin = 0, end = 0.9, option = 5), 
           ggtheme = theme,
           xlim = c(0,20), break.time.by=2,
           facet.by = "fungus",
           panel.labs = list(fungus = c("M. brunneum", "Sham")), 
           short.panel.labs = TRUE,
           panel.labs.font = list(face = "bold.italic", size = 11),
           panel.labs.background = list(color = "black", linewidth = 2, fill = "white"),
           censor = FALSE,
           conf.int = TRUE,
           conf.int.alpha = 0.1,
)
surv_plot3
aggregate(  censor ~ fungus + concentration, FUN=mean,data=antdata2)


#### 4.3 FluFu Indi analysed for 14 days ####
# Analyse dataset only for 14 days. 
# create a new column "survival_14" and "censor_14 in antdata2
antdata2$survival_14 <- ifelse(antdata2$survival >= 14, 14, antdata2$survival)
antdata2$censor_14 <- ifelse(antdata2$survival <= 13, 1, 0)
#baseline model which is then updated to become the full interaction model
null_model        <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1                          + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
flupy_model       <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration          + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
fungus_model      <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
interaction_model <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
anova(null_model, flupy_model, fungus_model, interaction_model)

antdata2$concentration_factor <- as.factor(antdata2$concentration)
null_model        <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1                          + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
flupy_model_f       <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration_factor          + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
fungus_model_f      <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration_factor + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
interaction_model_f <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration_factor * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
anova(null_model, flupy_model_f, fungus_model_f, interaction_model_f)

summary(interaction_model_f)
summary(glht(interaction_model_f, test=adjusted("BH")))
pairs_model <- emmeans(interaction_model_factor, ~ concentration_factor*fungus)
pairwise_model <- pairs(pairs_model)
CLD_model <- cld(pairwise_model, alpha = 0.05, adjust = "tukey")

levels(antdata2$concentration_factor)
levels(antdata2$fungus)

# get the pairwise differences?
coefficients(interaction_model_f)
#             0M 5M 50m 0s 5S 50S
mat <-  matrix(c( 1, 0, 0, 0, 0, 
                  0, 1, 0, 0, 0,
                  0, 0, 1, 0, 0,
                  0, 0, 0, 1, 0,
                  0, 0, 0, 0, 1,
                 -1, 1, 0, 0, 0,
                 -1, 0, 1, 0, 0,
                 -1, 0, 0, 1, 0,
                 -1, 0, 0, 0, 1,
                  0,-1, 1, 0, 0,
                  0,-1, 0,-1, 0,
                  0,-1, 0, 0, 1,
                  0, 0,-1, 1, 0,
                  0, 0,-1, 0, 1,
                  0, 0, 0, -1, 1),
               byrow=TRUE,
               ncol=5,
               dimnames=list(c(
                 "0M - 5M",
                 "0M - 50M",
                 "0M - 0S",
                 "0M - 5S",
                 "0M - 50S",
                 "5M - 50M",
                 "5M - 0S",
                 "5M - 5S",
                 "5M - 50S",
                 "50M - 0S",
                 "50M - 5S",
                 "50M - 50S",
                 "0S - 5S",
                 "0S - 50S",
                 "5S - 50S"
  ))
)

differences <- summary(glht(interaction_model_f, linfct = mat), test = adjusted("BH"))
pvalues <- as.numeric(differences$test$pvalues)
coefs <- as.numeric(differences$test$coefficients)
hazard_ratios <- sprintf("%.2f", exp(coefs))
names(pvalues) <- names(differences$test$coefficients)
names(coefs) <- names(differences$test$coefficients)
names(hazard_ratios) <- names(differences$test$coefficients)
differences
hazard_ratios








surviplot_14 <- survfit(Surv (time = survival_14, event = censor_14) ~ 1 + concentration + fungus, data=antdata2)
surv_plot14 <- ggsurvplot(surviplot_14, data = antdata2,
                         legend.title = 'FPF Concentration [ppm]', legend.labs = c(), legend = c(0.15,0.15),
                         xlab = 'Time (days)', ylab = 'Proportion Surviving',
                         palette = viridis(3, begin = 0, end = 0.9, option = 5), 
                         ggtheme = theme,
                         xlim = c(0,15), break.time.by=2,
                         facet.by = "fungus",
                         panel.labs = list(fungus = c("M. brunneum", "Sham")), 
                         short.panel.labs = TRUE,
                         panel.labs.font = list(face = "bold.italic", size = 11),
                         panel.labs.background = list(color = "black", linewidth = 2, fill = "white"),
                         censor = FALSE,
                         conf.int = TRUE,
                         conf.int.alpha = 0.1,
)
surv_plot14







#### Cutoffs other than 14 ####

# repeat the same for different cut offs to see if when we can see the interaction effect. 
# Set the range of cutoff values
cutoff_values <- 7:20
# Create a list to store the anova results for each model
anova_results <- vector("list", length(cutoff_values))
# Loop over the cutoff values
for (i in seq_along(cutoff_values)) {
  cutoff <- cutoff_values[i]
  # Create the survival and censor columns for the current cutoff
  survival_col <- paste0("survival_", cutoff)
  censor_col <- paste0("censor_", cutoff)
  antdata2[[survival_col]] <- ifelse(antdata2$survival >= cutoff, cutoff, antdata2$survival)
  antdata2[[censor_col]] <- ifelse(antdata2$survival <= (cutoff - 1), 1, 0)
  
  # Fit the four Cox proportional hazards models with the current cutoff
  null_model        <- coxme(Surv(time = antdata2[[survival_col]], event = antdata2[[censor_col]]) ~ 1 + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
  flupy_model       <- coxme(Surv(time = antdata2[[survival_col]], event = antdata2[[censor_col]]) ~ 1 + concentration + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
  fungus_model      <- coxme(Surv(time = antdata2[[survival_col]], event = antdata2[[censor_col]]) ~ 1 + concentration + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
  interaction_model <- coxme(Surv(time = antdata2[[survival_col]], event = antdata2[[censor_col]]) ~ 1 + concentration * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
  
  # Run the anova and store the output in the anova_results list
  anova_results[[i]] <- anova(null_model, flupy_model, fungus_model, interaction_model)
}

# Print the anova results for each cutoff value
for (i in seq_along(cutoff_values)) {
  cutoff <- cutoff_values[i]
  cat("Anova results for cutoff", cutoff, "\n")
  print(cutoff)
  print(anova_results[[i]])
}

# IT SEEMS THAT IT IS BEST MOST IDEAL FOR THE RESULTS IF WE LOOK AT THE EFFECTS AFTER 14 DAYS












#### Follow this up with the colony level effects ####