rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###                                                                                 ###
###   ### ### ###  ###         ###     ###     #####     ###     ###  ### ### ###   ###
###   ### ### ###  ###         ###     ###   ###   ###   ###     ###  ### ### ###   ###
###   ###          ###         ###     ###  ###     ###  ###     ###  ###           ###
###   ###          ###         ###     ###  ###          ###     ###  ###           ###
###   ### ###      ###         ###     ###  ###          ###     ###  ### ### ###   ###
###   ### ###      ###         ###     ###  ###   #####  ###     ###  ### ### ###   ###
###   ###          ###         ###     ###  ###   #####  ###     ###          ###   ###
###   ###          ###         ###     ###  ###     ###  ###     ###          ###   ###
###   ###          ### ### ###  ###   ###    ###   ###    ###   ###   ### ### ###   ###
###   ###          ### ### ###    #####        #####        #####     ### ### ###   ###
###                                                                                 ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Read me ####
# Author: Daniel Schläppi
# FLUGUS is short for FLUpyradifurone and funGUS
# This script contains all the analysis for the manuscript:
# Synergistic effects of the next generation insecticide flupyradifurone with a fungal pathogen
# Run it together with the following data files:
# exp1_flupy_susceptibility_test.csv
# exp2_flugus_interaction_survival.csv
# supplexp_food_uptake.txt


# Indexing
# 1. prerequisites
# 2. FPF susceptibility test
#   2.1 Stats
#   2.2 Survival plot
# 3. FPF Fungus co-Exposure
#   3.1 Interaction Stats
#   3.2 Interaction survival plot 
# _______________________________________________________________
# 4. Supplementary material
#   4.1 Prepare food uptake data
#   4.2 Standard curve & calculation of food uptake (volume in μl)
#   4.3 Food consumption stats 
#   4.4 Food consumption graph 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.  prerequisites ####

# load required libraries and functions:
{
library(dplyr) #contains pipe operator
library(ggplot2)
library(broom) #contains tidy()
library(viridis)
library(lme4)
library(car) # Anova()
library(multcomp) # contains cld
library(emmeans) # contains emmeans()
library(coxme)
library(survival)
library(survminer) # used in the analysis of the survival curves incl ggsurvplot
library(scales)
library(clipr)
}

source("/Users/gismo/Documents/GitHub/vital_rscripts_git/copy_paste_contrast_matrix.R")


# working directory 
directory <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/flugus_manuscript"
setwd(directory)

# load data
exp1_data <- read.csv('exp1_flupy_susceptibility_test.csv') #Experiment 1: flupy only survival data
exp1_data$concentration <- factor(exp1_data$concentration)

flugus_data <- read.csv('exp2_flugus_interaction_survival.csv') #Experiment 2: flupy fungus co-exposure survival data
flugus_data$fungus <- factor(flugus_data$fungus, levels=c("S","M"))

supl_data <- read.table("supplexp_food_uptake.txt", header = TRUE) # Supplementary Experiment: food uptake data



#### 2. FPF susceptibility test ####
#### 2.1 Stats ####
null_model <- coxme ( Surv (time = survival, event = censor) ~ 1                 + ( 1 | petri_dish), data = exp1_data)
full_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + ( 1 | petri_dish), data = exp1_data)
anova(null_model   ,  full_model )
summary(full_model)

posthocs_Tukey_1 <- summary(glht(full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH"))
posthocs_Tukey_2 <- summary(glht(full_model,            linfct = contrast_matrix)      ,test=adjusted("BH"))

letters <- cld(posthocs_Tukey_1)
letters <- cld(posthocs_Tukey_2)

copy_paste_contrast_matrix()
contrast_matrix <- rbind(
  "1minus2"=c(-1,0,0,0,0,0,0),
  "1minus3"=c(0,-1,0,0,0,0,0),
  "1minus4"=c(0,0,-1,0,0,0,0),
  "1minus5"=c(0,0,0,-1,0,0,0),
  "1minus6"=c(0,0,0,0,-1,0,0),
  "1minus7"=c(0,0,0,0,0,-1,0),
  "1minus8"=c(0,0,0,0,0,0,-1),
  "2minus3"=c(1,-1,0,0,0,0,0),
  "2minus4"=c(1,0,-1,0,0,0,0),
  "2minus5"=c(1,0,0,-1,0,0,0),
  "2minus6"=c(1,0,0,0,-1,0,0),
  "2minus7"=c(1,0,0,0,0,-1,0),
  "2minus8"=c(1,0,0,0,0,0,-1),
  "3minus4"=c(0,1,-1,0,0,0,0),
  "3minus5"=c(0,1,0,-1,0,0,0),
  "3minus6"=c(0,1,0,0,-1,0,0),
  "3minus7"=c(0,1,0,0,0,-1,0),
  "3minus8"=c(0,1,0,0,0,0,-1),
  "4minus5"=c(0,0,1,-1,0,0,0),
  "4minus6"=c(0,0,1,0,-1,0,0),
  "4minus7"=c(0,0,1,0,0,-1,0),
  "4minus8"=c(0,0,1,0,0,0,-1),
  "5minus6"=c(0,0,0,1,-1,0,0),
  "5minus7"=c(0,0,0,1,0,-1,0),
  "5minus8"=c(0,0,0,1,0,0,-1),
  "6minus7"=c(0,0,0,0,1,-1,0),
  "6minus8"=c(0,0,0,0,1,0,-1),
  "7minus8"=c(0,0,0,0,0,1,-1)
) 

posthocs_Tukey_1 <- summary(glht(interaction_model_qual,linfct=contrast_matrix_Tukey),test=adjusted("BH"))
posthocs_Tukey_2 <- summary(glht(interaction_model_qual_bis,linfct=mcp(interac="Tukey")),test=adjusted("BH"))
cld(posthocs_Tukey_2)






#### 2.2 Survival plot ####
# Survival plot for the FPF susceptibility test (exported as 800*600)
surv_plot <- survfit(Surv(time = survival, event = censor) ~ 1 + concentration, data = exp1_data)
legend <- c('0     - a', '0.5  - ab','1     - ab','5     - a','10   - ab','50   - ab','100  - b','500  - c')
line_types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "solid", "twodash")

surv_plot <- ggsurvplot(surv_plot, data = exp1_data,
                        censor = FALSE,
                        legend.title = 'Concentration [ppm]',
                        legend.labs = legend, 
                        legend = c(0.15, 0.25),
                        xlab = 'Time (days)',
                        ylab = 'Proportion Surviving',
                        break.time.by = 2,
                        xlim = c(0, 21),
                        ggtheme = theme_bw(),
                        palette = viridis(8, begin = 0, end = 0.9, option = 5),
                        conf.int = TRUE,
                        conf.int.alpha = 0.1,
                        linetype = line_types
)

surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank()) +  
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.ticks = element_line(linewidth = 1),
    strip.text = element_text(size = 14), 
    legend.key.width = unit(2.5, "line")
  )

surv_plot




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### ###

#### 3. FPF Fungus co-Exposure #### 

### Older version 
# #### 3.1 Stats ####
# fungus_model      <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
# interaction_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
# anova(fungus_model, interaction_model)
# summary(interaction_model)
# 
# contrast_matrix <- rbind("slope_fungusS"=c(1,0,0),"slope_fungusM"=c(1,0,1))
# summary(glht(interaction_model,linfct=contrast_matrix),test=adjusted("BH"))

#### 3.1 Stats - alternative ####
flugus_data_qual <- flugus_data
flugus_data_qual$concentration <- as.factor(flugus_data_qual$concentration)
flugus_data_qual$interac <- interaction(flugus_data_qual$fungus,flugus_data_qual$concentration)
fungus_model_qual          <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data_qual)
interaction_model_qual     <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data_qual)
interaction_model_qual_bis <- coxme ( Surv (time = survival, event = censor) ~ 1 + interac + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data_qual)

anova(fungus_model_qual, interaction_model_qual)
summary(interaction_model_qual)

contrast_matrix <- rbind(
"Delta50"=c(0,0,1,0,1),
"Delta5"=c(0,0,1,1,0),
"DeltaControl"=c(0,0,1,0,0),
"Delta50 minus Delta5"=c(0,0,0,-1,1),
"Delta50 minus DeltaControl"=c(0,0,0,0,1),
"Delta5 minus DeltaControl"=c(0,0,0,1,0))


posthocs <- summary(glht(interaction_model_qual,linfct=contrast_matrix_qual),test=adjusted("BH"))
HR_control <- exp(posthocs$test$coefficients["DeltaControl"])#increase in mortality rate caused by the fungus in the no-FPF group
HR_control_lo <- exp(posthocs$test$coefficients["DeltaControl"] - posthocs$test$sigma["DeltaControl"])
HR_control_hi <- exp(posthocs$test$coefficients["DeltaControl"] + posthocs$test$sigma["DeltaControl"])

HR_5 <- exp(posthocs$test$coefficients["Delta5"])#increase in mortality rate caused by the fungus in the FPF 5 group
HR_5_lo <- exp(posthocs$test$coefficients["Delta5"] - posthocs$test$sigma["Delta5"])
HR_5_hi <- exp(posthocs$test$coefficients["Delta5"] + posthocs$test$sigma["Delta5"])

HR_50 <- exp(posthocs$test$coefficients["Delta50"])#increase in mortality rate caused by the fungus in the FPF 50 group
HR_50_lo <- exp(posthocs$test$coefficients["Delta50"] - posthocs$test$sigma["Delta50"])
HR_50_hi <- exp(posthocs$test$coefficients["Delta50"] + posthocs$test$sigma["Delta50"])


contrast_matrix_Tukey <-  rbind(
  "1 minus 2"=c(-1,0,0,0,0)
  ,"1 minus 3"=c(0,-1,0,0,0)
  ,"1 minus 4"=c(0,0,-1,0,0)
  ,"1 minus 5"=c(-1,0,-1,-1,0)
  ,"1 minus 6"=c(0,-1,-1,0,-1)
  ,"2 minus 3"=c(1,-1,0,0,0)
  ,"2 minus 4"=c(1,0,-1,0,0)
  ,"2 minus 5"=c(0,0,-1,-1,0)
  ,"2 minus 6"=c(1,-1,-1,0,-1)
  ,"3 minus 4"=c(0,1,-1,0,0)
  ,"3 minus 5"=c(-1,1,-1,-1,0)
  ,"3 minus 6"=c(0,0,-1,0,-1)
  ,"4 minus 5"=c(-1,0,0,-1,0)
  ,"4 minus 6"=c(0,-1,0,0,-1)
  ,"5 minus 6"=c(1,-1,0,1,-1)
)


posthocs_Tukey_1 <- summary(glht(interaction_model_qual,linfct=contrast_matrix_Tukey),test=adjusted("BH"))
posthocs_Tukey_2 <- summary(glht(interaction_model_qual_bis,linfct=mcp(interac="Tukey")),test=adjusted("BH"))
cld(posthocs_Tukey_2)

#### 3.2 Plot ####
# # Survival plot for the FPF Fungus interaction test (exported as 800*600)
surviplot <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration + fungus, data=flugus_data)
aggregate(censor ~ fungus + concentration, FUN=mean,data=flugus_data)
chosen_colors <- viridis(2, option = 5, begin = 0, end = 0.8)
show_col(chosen_colors, labels = TRUE, borders = NULL)
color_vector <- rep(chosen_colors, times = 3)
y_offset <- 0.05 
line_types <- rep(c("dotdash", "longdash", "solid"), each = 2)

surv_plot <- ggsurvplot(surviplot, data = flugus_data,
                        pval = FALSE,
                        linetype = line_types,
                        lwd = 1,
                        xlab = 'Time (days)', ylab = 'Proportion Surviving',
                        palette = color_vector,
                        ggtheme = theme_bw(),
                        xlim = c(0, 15), break.time.by = 2,
                        censor = FALSE,
                        conf.int = TRUE,
                        conf.int.alpha = 0.2,
                        legend = c(0.4, 0.25),
                        legend.title = "Treatments"
)

surv_plot$plot <- surv_plot$plot + 
  scale_y_continuous(limits = c(0.2, 1), breaks = seq(0.2, 1, by = 0.2)) +
  theme(panel.grid = element_blank(), legend.position = "none") +  
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(linewidth = 1),
    strip.text = element_text(size = 14)) +
  geom_segment(aes(x = 14.3, y = 0.78, xend = 14.3, yend = 0.65))+
  geom_segment(aes(x = 14.3, y = 0.52, xend = 14.3, yend = 0.325)) +
  geom_segment(aes(x = 0.4, y = 0.33-y_offset, xend = 1, yend = 0.33-y_offset), color = chosen_colors[1], linetype = line_types[1]) +
  geom_segment(aes(x = 0.4, y = 0.30-y_offset, xend = 1, yend = 0.30-y_offset), color = chosen_colors[1], linetype = line_types[3]) +
  geom_segment(aes(x = 0.4, y = 0.271-y_offset, xend = 1, yend = 0.271-y_offset), color = chosen_colors[1], linetype = line_types[5]) +
  geom_segment(aes(x = 1.5, y = 0.33-y_offset, xend = 2.1, yend = 0.33-y_offset), color = chosen_colors[2], linetype = line_types[1]) +
  geom_segment(aes(x = 1.5, y = 0.30-y_offset, xend = 2.1, yend = 0.30-y_offset), color = chosen_colors[2], linetype = line_types[3]) +
  geom_segment(aes(x = 1.5, y = 0.271-y_offset, xend = 2.1, yend = 0.271-y_offset), color = chosen_colors[2], linetype = line_types[5]) +
  annotate("text", x = 0.3, y=0.36-y_offset, label="sham", color = "black", size = 4, fontface = "plain", hjust = 0) +
  annotate("text", x = 1.3, y=0.36-y_offset, label="fungus", color = "black", size = 4, fontface = "plain", hjust = 0) +
  annotate("text", x = 2.5, y=0.33-y_offset, label="0 ppm", color = "black", size = 4, fontface = "plain", hjust = 0) +
  annotate("text", x = 2.5, y=0.30-y_offset, label="5 ppm", color = "black", size = 4, fontface = "plain", hjust = 0) +
  annotate("text", x = 2.5, y=0.27-y_offset, label="50 ppm", color = "black", size = 4, fontface = "plain", hjust = 0) +
  annotate("text", x=14.7, y=0.71, label="n.s", color = "black", size = 4.5, fontface = "bold") +
  annotate("text", x=14.7, y=0.42, label="p = \n 0.01", color = "black", size = 4.5, fontface = "bold")

print(surv_plot)


#### 4. Supplementary material ####
#### 4.1 Prepare food uptake data ####

# adjust variables and correct fluorescence by subtracting the mean values of blank ants. 
supl_data$treatment <- as.factor(supl_data$treatment)
mean_negative_ant <- mean(supl_data$fluorescence[supl_data$sample_type == "negative_ant"])
for (i in 1:nrow(supl_data)) {
  if (supl_data$fluorescence[i] >= mean_negative_ant) {
    supl_data$fluorescence[i] <- supl_data$fluorescence[i] - mean_negative_ant
  } else {
    supl_data$fluorescence[i] <- 0
  }
}

# for simlicity create two subsets, one containing only the standard curves and one containing the samples
# std curves
data_stdcurves <- subset(supl_data, supl_data$treatment == "std_curve", drop = TRUE)
data_stdcurves$treatment <- droplevels(data_stdcurves$treatment)

{
  data_samples <- subset(supl_data, supl_data$sample_type == "ant")
  data_samples$treatment <- droplevels(data_samples$treatment) 
  data_samples$treatment <- relevel(data_samples$treatment, "mid") # reorder the levels of the treatment variable
  data_samples$treatment <- relevel(data_samples$treatment, "low")
  data_samples$treatment <- relevel(data_samples$treatment, "control")
}


#### 4.2 Standard curve & calculation of food uptake (volume in μl) ####
# Plot the pooled standard curve
data_mean <- data_stdcurves %>% 
  group_by(concentration) %>% 
  summarise(mean_fluorescence = mean(fluorescence))
ggplot(data_mean, aes(x = concentration, y = mean_fluorescence)) + 
  geom_line() + 
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  ggtitle("Mean Standard Curve with Linear Model Fit") + 
  xlab("Concentration [μl fluorescin / 100 μl sugar water]") + 
  ylab("Mean Fluorescence")
# get coeficients from linear model (forced through zero) to estimate how much food the ants consumed (based on the standard curves)
std_model <- lm(mean_fluorescence ~ 0 + concentration, data = data_mean)
coef <- coefficients(std_model)
slope <-  coef[[1]]
# calculate estimated volume of food uptake and add it to the dataframe
data_samples$consumed_volume <- data_samples$fluorescence / slope

#### 4.3 Food consumption stats ####
group_by(data_samples, treatment) %>%
  summarise(
    count = n(),
    mean = mean(consumed_volume, na.rm = TRUE),
    sd = sd(consumed_volume, na.rm = TRUE),
    min = min(consumed_volume, na.rm = TRUE), 
    max = max(consumed_volume, na.rm = TRUE)
  )

mod <- glmer.nb(consumed_volume ~ treatment + (1|petridish), data = data_samples)
summary(mod)
Anova(mod)
letters <- c("a", "b", "c", "d", "e", "f")
emmeans(mod, pairwise ~ "treatment", adjust = "tukey")
cld(emmeans(mod, pairwise ~ "treatment", adjust = "tukey"), Letters = letters)
residuals <- residuals(mod)
plot(residuals)
qqnorm(residuals)
qqline(residuals, col = "red")


#### 4.4 Food consumption graph ####

# nice graph based on 0 intercept model exportet at 600x700
# Define the labels and their positions
labels_df <- data.frame(treatment = c("high", "mid", "control", "low"),
                        label = c("ab", "a", "ab", "b"))
# Plot with labels (exported as )
gg <- ggplot(data_samples, aes(x = treatment, y = consumed_volume)) +
  geom_boxplot(fill = alpha("grey", 0.5), color = "black", notch = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, color = viridis(1)[1]) +
  labs(x = "Treatments", y = "Honey water consumption [μL]", title = "") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 0.7))
label_positions <- ggplot_build(gg)$data[[1]]$x  # calculate position of x to put the letters right above the boxplots. 
gg <- gg + geom_text(data = labels_df, aes(x = label_positions, y = 0.7, label = label, fontface = "bold"),
                     hjust = 0.5, vjust = 0.5, size = 5)
print(gg)

