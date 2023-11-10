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
# This script contains all the analysis for the manuscript:
# Synergistic effects of the insecticide Flupyradifurone (FPF) and the entomopathogen Metarhizium brunneum in ants
# Run it together with the following data files:

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

# Supplementary Material:
# data
# 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.  prerequisites ####

# load required libraries and functions:
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

# working directory 
directory <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/flugus_manuscript/"
setwd(directory)

#### 2. FPF susceptibility test ####
# load data
exp1_data <- read.csv('flupy_susceptibility_test.csv') #flupy only survival data
exp1_data$concentration <- factor(exp1_data$concentration)

flugus_data <- read.csv('FPF_fungus_survival.csv') #flupy fungus co-exposure survival data
flugus_data$fungus <- factor(flugus_data$fungus, levels=c("S","M"))

#### 2.1 Stats ####
null_model <- coxme ( Surv (time = survival, event = censor) ~ 1                 + ( 1 | petri_dish), data = exp1_data)
full_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + ( 1 | petri_dish), data = exp1_data)
anova(null_model   ,  full_model )
summary(full_model)

summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH"))
letters <- cld(summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH")))

#### 2.2 Survival plot ####
# Survial plot for the FPF susceptibility test (exported as 800*600)
surv_plot <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration , data=exp1_data) #basic plot
legend <- c('0     - a', '0.5  - ab','1     - ab','5     - a','10   - ab','50   - ab','100  - b','500  - c')

surv_plot <- ggsurvplot(surv_plot, data = exp1_data,
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

surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank()) +  
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.ticks = element_line(linewidth = 1),
    strip.text = element_text(size = 14)
  )
surv_plot






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### ###

#### 3. FPF Fungus co-Exposure #### 
#### 3.1 Stats ####
fungus_model      <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
interaction_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
anova(fungus_model, interaction_model)
summary(interaction_model)

contrast_matrix <- rbind("slope_fungusS"=c(1,0,0),"slope_fungusM"=c(1,0,1))
summary(glht(interaction_model,linfct=contrast_matrix),test=adjusted("BH"))

# how to write this in the stats section?


#### 3.2 Plot ####
surviplot <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration + fungus, data=flugus_data)
aggregate(censor ~ fungus + concentration, FUN=mean,data=flugus_data)
chosen_colors <- viridis(3, option = "inferno", begin = 0.35, end = 0.9)
show_col(chosen_colors, labels = TRUE, borders = NULL)
color_vector <- rep(chosen_colors, each = 2)

y_offset <- 0.05 

{surv_plot <- ggsurvplot(surviplot, data = flugus_data,
                          pval = FALSE,
                          linetype = "fungus",
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
    geom_segment(aes(x = 0, y = 0.4-y_offset, xend = 0.6, yend = 0.4-y_offset)) +
    geom_segment(aes(x = 0, y = 0.3699-y_offset, xend = 0.6, yend = 0.3699-y_offset), linetype = "dashed") +
    geom_segment(aes(x = 0, y = 0.33-y_offset, xend = 0.6, yend = 0.33-y_offset), color = chosen_colors[1]) +
    geom_segment(aes(x = 0, y = 0.30-y_offset, xend = 0.6, yend = 0.30-y_offset), color = chosen_colors[2]) +
    geom_segment(aes(x = 0, y = 0.271-y_offset, xend = 0.6, yend = 0.271-y_offset), color = chosen_colors[3]) +
    annotate("text", x = 0.8, y=0.4-y_offset, label="sham", color = "black", size = 4, fontface = "plain", hjust = 0) +
    annotate("text", x = 0.8, y=0.37-y_offset, label="fungus", color = "black", size = 4, fontface = "plain", hjust = 0) +
    annotate("text", x = 0.8, y=0.33-y_offset, label="0 ppm", color = "black", size = 4, fontface = "plain", hjust = 0) +
    annotate("text", x = 0.8, y=0.30-y_offset, label="5 ppm", color = "black", size = 4, fontface = "plain", hjust = 0) +
    annotate("text", x = 0.8, y=0.27-y_offset, label="50 ppm", color = "black", size = 4, fontface = "plain", hjust = 0) +
    annotate("text", x=14.7, y=0.71, label="n.s", color = "black", size = 4.5, fontface = "bold") +
    annotate("text", x=14.7, y=0.42, label="p = \n 0.01", color = "black", size = 4.5, fontface = "bold")
  print(surv_plot)
}





#### 4. Supplementary material ####
#### 4.1 Prepare food uptake data ####
# data
data <- read.table("food_uptake.txt", header = TRUE)
# adjust variables and correct fluorescence by subtracting the mean values of blank ants. 
data$treatment <- as.factor(data$treatment)
mean_negative_ant <- mean(data$fluorescence[data$sample_type == "negative_ant"])
for (i in 1:nrow(data)) {
  if (data$fluorescence[i] >= mean_negative_ant) {
    data$fluorescence[i] <- data$fluorescence[i] - mean_negative_ant
  } else {
    data$fluorescence[i] <- 0
  }
}

# for simlicity create two subsets, one containing only the standard curves and one containing the samples
# std curves
data_stdcurves <- subset(data, data$treatment == "std_curve", drop = TRUE)
data_stdcurves$treatment <- droplevels(data_stdcurves$treatment)
# data
{
  data_samples <- subset(data, data$sample_type == "ant")
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