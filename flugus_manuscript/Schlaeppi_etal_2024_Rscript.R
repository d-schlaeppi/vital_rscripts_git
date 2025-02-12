rm(list = setdiff(ls(), "first_time_use_working_directory"))

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
pacman::p_load(dplyr, ggplot2, broom, viridis, lme4, car, multcomp, multcompView,
               emmeans, coxme, survival, survminer, scales, clipr, gridExtra, grid)

# Set working directory
if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") { # direct it to where you have config_user_and_hd.R (typically the script folder or github folder)
  standard <- "/Users/gismo/Documents/GitHub/vital_rscripts_git" # if you are always working from the same directory just put its name here and it will save you some clicking.  
  selected_dir <- if  (dir.exists(standard)) {standard} else {tcltk::tk_choose.dir(default = "~/", caption = "Select Working Directory")}
  if (is.null(selected_dir) || selected_dir == "") {
    cat("No directory selected. Exiting.\n")
    return()}
  setwd(selected_dir)
  first_time_use_working_directory <<- getwd()
  cat(crayon::blue(getwd()))
} else { setwd(first_time_use_working_directory)
  cat(crayon::blue(getwd())) }


# source("config_user_and_hd.R") # contains getUserOptions() that defines usr and hd and some functions
directory <- paste(first_time_use_working_directory, "flugus_manuscript", sep = "/")
setwd(directory)

# source("/Users/gismo/Documents/GitHub/vital_rscripts_git/copy_paste_contrast_matrix.R")

# load data sets
exp1_data <- read.csv('exp1_flupy_susceptibility_test.csv') #Experiment 1: flupy only survival data
exp1_data$concentration <- factor(exp1_data$concentration)

flugus_data <- read.csv('exp2_flugus_interaction_survival.csv') #Experiment 2: flupy fungus co-exposure survival data
flugus_data$fungus <- factor(flugus_data$fungus, levels=c("S","M"))
flugus_data$concentration <- as.factor(flugus_data$concentration)

supl_data <- read.table("supplexp_food_uptake.txt", header = TRUE) # Supplementary Experiment: food uptake data


#### 2. FPF susceptibility test ####
#### 2.1 Stats ####
null_model <- coxme ( Surv (time = survival, event = censor) ~ 1                 + ( 1 | petri_dish), data = exp1_data)
full_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + ( 1 | petri_dish), data = exp1_data)
anova(null_model   ,  full_model )
summary(full_model)

posthocs_Tukey <- summary(glht(full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH"))
letters <- cld(posthocs_Tukey)



#### 2.2 Survival plot ####
# Survival plot for the FPF susceptibility test (exported as 800*600)

y_coordinates <- exp1_data %>% 
  group_by(concentration) %>%
  summarize(prop_surviving = sum(censor==0) / n()) %>% as.data.frame()
plot_data <- exp1_data %>% mutate(survival = ifelse(survival == 20 & censor == 0, 20.5, survival))


legend       <- c('    0', ' 0.5','    1','    5','  10','  50','100','500')
significance <- c('    0   -  a', ' 0.5   -  ab','    1   -  ab','    5   -  a','  10   -  ab','  50   -  ab','100*  -  b','500*  -  c')

dash             <- rep("-", 8)
significance_2   <- c('a','ab','ab','a','ab','ab','b','c')
treatment_labels <- c('0','0.5','1','5','10','50','100*','500*')
line_types       <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "solid", "twodash")

{surv_plot <- survfit(Surv(time = survival, event = censor) ~ 1 + concentration, data = plot_data)
  surv_plot <- ggsurvplot(surv_plot, data = exp1_data,
                          censor = FALSE,
                          legend.title = 'Concentration [ppm]',
                          legend.labs = legend, 
                          legend = c(0.15, 0.25),
                          xlab = 'Time (days)',
                          ylab = 'Proportion Surviving',
                          # break.time.by = 2,
                          xlim = c(0, 24),
                          ylim = c(0, 1),
                          ggtheme = theme_bw(),
                          palette = viridis(8, begin = 0, end = 0.9, option = 5),
                          conf.int = TRUE,
                          conf.int.alpha = 0.1,
                          linetype = line_types, 
                          size = 0.8,
                          alpha = 0.8)
  surv_plot$plot <- surv_plot$plot + 
    theme(panel.grid = element_blank()) +  
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.ticks = element_line(linewidth = 0.8),
      strip.text = element_text(size = 10), 
      legend.key.width = unit(2, "line"),
      legend.position = "none") + 
    scale_x_continuous(breaks = seq(0, 20, by = 2)) + 
    scale_y_continuous(
      breaks = seq(0.25, 1, by = 0.15),
      limits = c(0.16, 1)) +
    
    geom_text(data = y_coordinates, 
              aes(x = 23, y = prop_surviving, label = dash), 
              size = 3.5, hjust = 0, vjust = 0.5, inherit.aes = FALSE) +
    geom_text(data = y_coordinates, 
              aes(x = 23.5 , y = prop_surviving, label = significance_2), 
              size = 3.5, hjust = 0, vjust = 0.5, inherit.aes = FALSE) +
    geom_text(data = y_coordinates[1:6, ],
              aes(x = 21.5, y = prop_surviving, label = treatment_labels[1:6]),
              size = 3.5, hjust = 0, vjust = 0.5, inherit.aes = FALSE) +
    geom_text(data = y_coordinates[7:8, ],
              aes(x = 21.5, y = prop_surviving, label = treatment_labels[7:8]),
              size = 3.5, hjust = 0, vjust = 0.5, inherit.aes = FALSE,
              fontface = "bold") +
    
    
    # geom_text(data = y_coordinates[1:6, ], 
    #           aes(x = 21.5, y = prop_surviving, label = significance[1:6]), 
    #           size = 3.5, hjust = 0, vjust = 0.5, inherit.aes = FALSE) +
    # geom_text(data = y_coordinates[7:8, ], 
    #           aes(x = 21.5, y = prop_surviving, label = significance[7:8]), 
    #           size = 3.5, hjust = 0, vjust = 0.5, inherit.aes = FALSE, 
    #           fontface = "bold") +
    annotate("text", x = 20.75, y = 0.97, label = "c [ppm]", size = 3.5, hjust = 0, fontface = "bold")
  surv_plot
}

#700*620








### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### ###

#### 3. FPF Fungus co-Exposure #### 

#### 3.1 Stats ####

flugus_data$interac <- interaction(flugus_data$fungus,flugus_data$concentration)

null_model            <- coxme ( Surv (time = survival, event = censor) ~ 1                              + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
concentration_model   <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration              + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
additive_model        <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + fungus     + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
interaction_model     <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration * fungus     + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)
interaction_model_bis <- coxme ( Surv (time = survival, event = censor) ~ 1 + interac                    + (1 | petri_dish) + (1 | colony) + (1 | block), data = flugus_data)

anova(null_model, concentration_model, additive_model, interaction_model)
anova(additive_model, interaction_model)
Anova(interaction_model)
summary(interaction_model)

contrast_matrix <- rbind(
"Delta50"=c(0,0,1,0,1),
"Delta5"=c(0,0,1,1,0),
"DeltaControl"=c(0,0,1,0,0),
"Delta50 minus Delta5"=c(0,0,0,-1,1),
"Delta50 minus DeltaControl"=c(0,0,0,0,1),
"Delta5 minus DeltaControl"=c(0,0,0,1,0))

posthocs <- summary(glht(interaction_model,linfct=contrast_matrix),test=adjusted("BH"))

# for plot
# delta mortality chance by fungus depending on flupy concentration

HR_control <- exp(posthocs$test$coefficients["DeltaControl"]) #increase in mortality rate caused by the fungus in the no-FPF group
HR_control_lo <- exp(posthocs$test$coefficients["DeltaControl"] - posthocs$test$sigma["DeltaControl"])
HR_control_hi <- exp(posthocs$test$coefficients["DeltaControl"] + posthocs$test$sigma["DeltaControl"])

HR_5 <- exp(posthocs$test$coefficients["Delta5"]) #increase in mortality rate caused by the fungus in the FPF 5 group
HR_5_lo <- exp(posthocs$test$coefficients["Delta5"] - posthocs$test$sigma["Delta5"])
HR_5_hi <- exp(posthocs$test$coefficients["Delta5"] + posthocs$test$sigma["Delta5"])

HR_50 <- exp(posthocs$test$coefficients["Delta50"]) #increase in mortality rate caused by the fungus in the FPF 50 group
HR_50_lo <- exp(posthocs$test$coefficients["Delta50"] - posthocs$test$sigma["Delta50"])
HR_50_hi <- exp(posthocs$test$coefficients["Delta50"] + posthocs$test$sigma["Delta50"])

hr_data <- data.frame(
  Treatment = c("0 ppm", "5 ppm", "50 ppm"),
  HR = c(HR_control, HR_5, HR_50),
  HR_lo = c(HR_control_lo, HR_5_lo, HR_50_lo),
  HR_hi = c(HR_control_hi, HR_5_hi, HR_50_hi))





#### 3.2 Plot ####
# # Survival plot for the FPF Fungus interaction test (exported as 800*600)
surviplot <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration + fungus, data=flugus_data)
# aggregate(censor ~ fungus + concentration, FUN=mean,data=flugus_data)
chosen_colors <- viridis(2, option = 5, begin = 0, end = 0.8)
# show_col(chosen_colors, labels = TRUE, borders = NULL)
color_vector <- rep(chosen_colors, times = 3)
y_offset <- 0.05 
line_types <- rep(c("dotdash", "longdash", "solid"), each = 2)

y_coordinates <- flugus_data %>% group_by(concentration, fungus) %>%
  summarize(prop_surviving = sum(censor==0) / n(),
            mean_censor = mean(censor))  %>% 
  arrange(desc(prop_surviving)) %>% as.data.frame()
labeling <- c(" 5","50"," 0"," 0"," 5","50")

surv_plot <- ggsurvplot(surviplot, data = flugus_data,
                        pval = FALSE,
                        linetype = line_types,
                        lwd = 0.8,
                        xlab = 'Time (days)', ylab = 'Proportion Surviving',
                        palette = color_vector,
                        ggtheme = theme_bw(),
                        xlim = c(0, 16),
                        censor = FALSE,
                        conf.int = TRUE,
                        conf.int.alpha = 0.2,
                        legend = c(0.4, 0.25),
                        legend.title = "Treatments")

surv_plot$plot <- surv_plot$plot + 
  scale_y_continuous(limits = c(0.2, 1), breaks = seq(0.2, 1, by = 0.2)) +
  theme(panel.grid = element_blank(), legend.position = "none") +  
  theme(plot.margin = margin(20, 10, 10, 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.ticks = element_line(linewidth = 1),
        strip.text = element_text(size = 14)) +
  geom_segment(aes(x = 14.8, y = 0.78, xend = 14.8, yend = 0.65))+
  geom_segment(aes(x = 14.8, y = 0.52, xend = 14.8, yend = 0.325)) +
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
  annotate("text", x=15.5, y=0.71, label="n.s       ", color = "black", size = 3, fontface = "bold") +
  annotate("text", x=15.5, y=0.42, label=" p = 0.01", color = "black", size = 3, fontface = "bold") +
  geom_rect(aes(xmin = 0.2, xmax = 3.8, ymin = 0.2, ymax = 0.38 - y_offset), 
            fill = NA, color = "black", lwd = 0.05) +
    scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  geom_text(data = y_coordinates, 
            aes(x = 14.2, y = prop_surviving, label = labeling), 
            size = 3, hjust = 0, vjust = 0.5, inherit.aes = FALSE)

print(surv_plot)

surv_ggplot <- surv_plot$plot







# Hazard plot for M. brunneum in dependence of pesticide concentration 
hazard_plot <- ggplot(hr_data, aes(x = Treatment, y = HR)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.8) + 
  geom_errorbar(aes(ymin = HR_lo, ymax = HR_hi), width = 0.1, lwd = 1.2) + # Error bars for confidence intervals
  labs(y = "Fungus Hazard Ratio (HR)",
       x = "Flupyradifurone concentration") +
  theme_bw() +
  theme(plot.margin = margin(20, 10, 10, 10),
        panel.grid = element_blank(), legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.ticks = element_line(linewidth = 1),
        strip.text = element_text(size = 14)) +
  annotate("text", x = 0.96, y = 5.35, label = "a", color = "black", size = 6, fontface = "bold", hjust = 0) +
  annotate("text", x = 1.96, y = 5.35, label = "b", color = "black", size = 6, fontface = "bold", hjust = 0) +
  annotate("text", x = 2.96, y = 5.35, label = "b", color = "black", size = 6, fontface = "bold", hjust = 0) +
  coord_cartesian(ylim=c(1,5.5))  # Setting y-axis limit to start from 1
print(hazard_plot)

# print the two plots together | Exported as 1200*500 clip
layout_mat <- rbind(c(1,1,2), c(1,1,2))
label_a <- textGrob("a", x = unit(0, "npc") + unit(5, "mm"), y = unit(1, "npc") - unit(5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 20, fontface = "bold"))
label_b <- textGrob("b", x = unit(0, "npc") + unit(5, "mm"), y = unit(1, "npc") - unit(5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(surv_ggplot, hazard_plot, layout_matrix = layout_mat)
grid.arrange(
  arrangeGrob(surv_ggplot, top = label_a),
  arrangeGrob(hazard_plot, top = label_b),
  layout_matrix = layout_mat
)





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

# for simplicity create two subsets, one containing only the standard curves and one containing the samples
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

# nice graph based on 0 intercept model exported at 600x700
labels_df <- data.frame(treatment = c("high", "mid", "control", "low"),
                        label = c("ab", "a", "ab", "b"))

gg <- ggplot(data_samples, aes(x = treatment, y = consumed_volume)) +
  geom_boxplot(fill = alpha("grey", 0.5), color = "black", notch = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, color = viridis(1)[1]) +
  labs(x = "Treatments", y = "Honey water consumption [μL]", title = "") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 0.7))
label_positions <- ggplot_build(gg)$data[[1]]$x
gg <- gg + geom_text(data = labels_df, aes(x = label_positions, y = 0.7, label = label, fontface = "bold"),
                     hjust = 0.5, vjust = 0.5, size = 5)
print(gg)

