rm(list=ls())

### Group Flugus experiment with Ellie ###


directory <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/group_flugus/"
setwd(directory)

# libraries
# library(blmeco) #contains compareqqnorm  (multiple qq boxplots)
library(dplyr) #contains pipe operator
library(broom) #contains tidy()
library(viridis)
library(lme4)
library(car) # Anova()
library(multcomp) # contains cld
library(emmeans) # contains emmeans()
library(coxme)
library(survminer) # used in the analysis of the survival curves incl ggsurvplot
library(ggplot2)
library(survival)
library(emmeans)

source('printme_coxme.R') # used in the analysis of the survival curves


exp1_data <- read.csv("pre_exposure_survival.csv")
exp1_data$concentration <- factor(exp1_data$concentration )

#### pre fungus survival ####
null_model <- coxme ( Surv (time = survival, event = censor) ~ 1                 + ( 1 | petri_dish) , data = exp1_data)
full_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + ( 1 | petri_dish) , data = exp1_data)
anova(null_model   ,  full_model )
summary(full_model)

summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH"))
letters <- cld(summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH")))
letters

#### survival plot ####
surv_plot <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration , data=exp1_data) #basic plot
legend <- c('0  - a','5  - a','50 - b')

surv_plot <- ggsurvplot(surv_plot, data = exp1_data,
                        censor = FALSE,
                        legend.title = 'concentration [ppm]',
                        legend.labs = legend, 
                        legend = c(0.25,0.25),
                        xlab = 'Time (days)',
                        ylab = 'Proportion Surviving',
                        break.time.by=2,
                        xlim = c(0,13),
                        ggtheme = theme_bw(),
                        palette = viridis(3, begin = 0, end = 0.8, option = 5),
                        conf.int= TRUE,
                        conf.int.alpha = 0.1
)
# remove the grid lines
surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank()) +  
  theme(# Adjust font
    axis.title = element_text(size = 16),          # Font size of axis titles
    axis.text = element_text(size = 14),           # Font size of axis labels
    legend.title = element_text(size = 16),        # Font size of legend title
    legend.text = element_text(size = 14),         # Font size of legend labels
    axis.ticks = element_line(linewidth = 1),           # Size of axis tick marks
    strip.text = element_text(size = 14)           # Font size of facet labels (if applicable)
  )
surv_plot # export eg. as 800*600 




#### post fungus exposure survival  ####

antdata2 <- read.csv("post_exposure_survival.csv")
antdata2$fungus <- as.factor(antdata2$fungus)
antdata2$petri_dish <- as.factor(antdata2$petri_dish)
antdata2$concentration <- as.numeric(as.character(antdata2$concentration))
antdata2$concentration_factor <- as.factor(antdata2$concentration)
head(antdata2)

#### stats ####

#baseline model which is then updated to become the full interaction model
null_model          <- coxme ( Surv (time = survival, event = censor) ~ 1                                 + (1 | petri_dish) , data = antdata2)
flupy_model_f       <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration_factor          + (1 | petri_dish)  , data = antdata2)
fungus_model_f      <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration_factor + fungus + (1 | petri_dish) , data = antdata2)
interaction_model_f <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration_factor * fungus + (1 | petri_dish) , data = antdata2)
anova(null_model, flupy_model_f, fungus_model_f, interaction_model_f)

anova(fungus_model_f, interaction_model_f)
Anova(interaction_model_f)
# summary(interaction_model_f)
contrast_matrix <- rbind (
  "fungus:5FPF-fungus:0FPF" = c(1,0,0,1,0)
  ,
  "fungus:50FPF-fungus:0FPF"=c(0,1,0,0,1)
  ,
  "sham:5FPF - sham:0FPF"=c(1,0,0,0,0)
  ,
  "sham:50FPF - sham:0FPF"=c(0,1,0,0,0)
)
summary(glht(  interaction_model_f, linfct =  contrast_matrix),test=adjusted("BH"))
fungus_model_n      <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + fungus + (1 | petri_dish) , data = antdata2)
interaction_model_n <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration * fungus + (1 | petri_dish) , data = antdata2)
anova(fungus_model_n, interaction_model_n)
summary(interaction_model_n)
Anova(interaction_model_n, type="III")
Anova(interaction_model_n)
contrast_matrix <- rbind("slope_fungusS"=c(1,0,0),"slope_fungusM"=c(1,0,1))
summary(glht(interaction_model_n,linfct=contrast_matrix),test=adjusted("BH"))
# pairwise comparison
em_means <- emmeans(interaction_model_f, ~ concentration_factor * fungus)
letters <- c("a", "b", "c", "d", "e", "f")
# Compute the CLD (Comparisons of Least Squares Means) using the letters
cld_results <- cld(em_means, Letters = letters)
letters <- cld_results$.group


#### graph ####

surv_mod <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration + fungus, data=antdata2)
# because it seems impossible to create the legend the plot is created two times once to make the graph look right and once just to get the legend
colors <- rep(viridis(3, begin = 0, end = 0.9, option = 5), length.out = 3)
color_palette <- c(colors[1],colors[1], colors[2], colors[2],colors[3], colors[3])
surv_plot <- ggsurvplot(surv_mod, data = antdata2,
                        linetype = "fungus",
                        xlab = 'Time (days)', ylab = 'Proportion Surviving',
                        palette = color_palette, 
                        ggtheme = theme_bw(),
                        xlim = c(0, 15), break.time.by = 2,
                        censor = FALSE,
                        conf.int = TRUE,
                        conf.int.alpha = 0.1, 
                        legend = "none")
surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank())+  
  theme(# Adjust font
    axis.title = element_text(size = 16),          # Font size of axis titles
    axis.text = element_text(size = 14),           # Font size of axis labels
    legend.title = element_text(size = 16),        # Font size of legend title
    legend.text = element_text(size = 14),         # Font size of legend labels
    axis.ticks = element_line(linewidth = 1),           # Size of axis tick marks
    strip.text = element_text(size = 14)           # Font size of facet labels (if applicable)
  )
surv_plot

# redo the same and just copy paste a screenshot with the desired elements on top of the plot created earlier. Sounds stupid? Yes, it is! But I have wasted too much time already. And it is annoying because the solution is for sure super easy and one line of code with just fix it... 
# Create new factor variables with desired labels 
antdata_updated <- transform(antdata2,
                             fungus_label = ifelse(fungus == "F", "fungus", "sham"),
                             concentration_label = ifelse(concentration == 0, "0 ppm",
                                                          ifelse(concentration == 5, "5 ppm", "50 ppm")))
letters <- gsub(" ", "", letters)
# Create the survival plot with the updated factor variables
surv_mod <- survfit(Surv(time = survival, event = censor) ~ 1 + concentration_label + fungus_label, data = antdata_updated)
color_palette <- viridis(3, begin = 0, end = 0.9, option = 5)
surv_plot <- ggsurvplot(surv_mod, data = antdata_updated,
                        linetype = "fungus_label",
                        xlab = 'Time (days)', ylab = 'Proportion Surviving',
                        palette = color_palette,
                        color = "concentration_label",  # Specify the variable for color mapping
                        ggtheme = theme_bw(),
                        xlim = c(0, 15), break.time.by = 2,
                        censor = FALSE,
                        conf.int = TRUE,
                        conf.int.alpha = 0.1,
                        legend = c(0.1, 0.25),
                        legend.title = "Treatments")  # Update legend title
surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank()) +
  theme(# Adjust font
    axis.title = element_text(size = 16),          # Font size of axis titles
    axis.text = element_text(size = 14),           # Font size of axis labels
    legend.title = element_text(size = 16),        # Font size of legend title
    legend.text = element_text(size = 14),         # Font size of legend labels
    axis.ticks = element_line(linewidth = 1),      # Size of axis tick marks
    strip.text = element_text(size = 14),          # Font size of facet labels (if applicable)
    legend.key = element_blank()                   # Remove legend key (linetype)
  ) +
  geom_text(data = data.frame(x = rep(14.5, length(letters)), y = c(0.55, 0.49, 0.37, 0.18, 0.06, 0.03), label = letters),
            aes(x = x, y = y, label = label),
            color = "black",
            size = 5,
            hjust = 0.5,
            vjust = 0,
            fontface = "bold",
            inherit.aes = FALSE)
surv_plot




