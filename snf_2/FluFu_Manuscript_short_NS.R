rm(list=ls())

# FluFu short #
# analysis for the manuscript:
# Synergistic effects of the insecticide Flupyradifurone and the entomopathogen Metarhizium brunneum in ants

#### 1.  prerequisites ####


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

#### Data ####

# set general working directory and get data
directory <- "/Users/gismo/Documents/GitHub/vital_rscripts_git/snf_2/"
# directory <- "C:/Users/bzniks/Dropbox/Papers/5_Bristol/4_Pesticides_Individual_Immunity/data_and_stats"
setwd(directory) # homeoffice mac github folder


source('printme_coxme.R') # used in the analysis of the survival curves

# source("C:/Users/bzniks/Dropbox/Post-doc_Lausanne/source_programs/printme_coxme.R")
# load data
dat <- read.table("full_fluorescin_feeding_quantification.txt", header = TRUE)

# prepare data frame:
# fluorescence corrected for ant debris -> subtract values from ant-only wells (negative controls) (is different for run 1 ant 2)
mean_negative_ant_week1_run1 <- mean(dat$fluorescence_week1[dat$sample_type == "negative_ant" & dat$run == 1])
mean_negative_ant_week2_run1 <- mean(dat$fluorescence_week2[dat$sample_type == "negative_ant" & dat$run == 1])
mean_negative_ant_week1_run2  <- mean(dat$fluorescence_week1[dat$sample_type == "negative_ant" & dat$run == 2])

# data frame with corrected fluorescence values
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
  fluorescence <- ifelse(dat[i, "run"] == 1 & dat[i, "fluorescence_week1"] >= mean_negative_ant_week1_run1,
                                         dat[i, "fluorescence_week1"] - mean_negative_ant_week1_run1,
                                         ifelse(dat[i, "run"] == 2 & dat[i, "fluorescence_week1"] >= mean_negative_ant_week1_run2,
                                                dat[i, "fluorescence_week1"] - mean_negative_ant_week1_run2,
                                                0))
  data <- rbind(data, data.frame(run, nr, colony, petridish, treatment, sample_type, concentration, fluorescence
  ))
}

# subset the dataset to week 1 only (in the first week we did a second measurement of the samples after the plates have been in the freezer for a week --> showing clear decay of the signal)
# here only week 1 measurements are relevant
# also, for simplicity we split the data into samples and std. curves 

# std curves only data 
data_stdcurves <- subset(data, dat$treatment == "std_curve", drop = TRUE)
data_stdcurves$treatment <- droplevels(data_stdcurves$treatment)

# samples only data
data_samples <- subset(data, dat$sample_type == "ant")
data_samples$treatment <- droplevels(data_samples$treatment) 
{
data_samples$treatment <- relevel(data_samples$treatment, "mid")# reorder the levels of the treatment variable
data_samples$treatment <- relevel(data_samples$treatment, "low")
data_samples$treatment <- relevel(data_samples$treatment, "control")
}

#### preliminary data exploration ####

head(data_samples)

# preliminary data exploration: 
group_by(data_samples, treatment) %>%
  summarise(
    count = n(),
    mean = mean(fluorescence, na.rm = TRUE),
    sd = sd(fluorescence, na.rm = TRUE)
  )

boxplot(fluorescence ~ treatment, data = data_samples)
# looks like high is lower and mid might be slightly lower as well...
# calculate actual consumed volume (as we had two runs we well get the standard curves for both separately)



#### Calculate meal size in μl ####

# Step 1 use standard curves to calculate micro liter of food instead of fluorescence
# Plot raw standard curves (we had 4 std. curves for the first and two for the second run)
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
  theme_bw() +
  geom_point(size = 3) + 
  geom_smooth(aes(group = run), method = "lm", se = FALSE, linewidth = 0.5) + 
  ggtitle("Standard Curves with Linear Model Fit") + 
  xlab("Concentration [μl fluorescin / 100 μl sugar water]") + 
  ylab("Mean Fluorescence") + 
  scale_color_viridis_d(name = "Run")

# get coeficients from linear model to calculate of how much food the ants consumed based on the standard curves (forcing std curves through zero)
data_mean_lm_0 <- data_mean %>% 
  group_by(run) %>% 
  do(mod = lm(mean_fluorescence ~ 0 + concentration, data = .))
coef_df_0 <- data_mean_lm_0 %>% 
  do(tidy(.$mod)) %>% 
  bind_rows()


# calculate meal volume and add it to the dataframe: 
data_samples$consumed_volume_0 <- NA
for (i in 1:nrow(data_samples)) {
  if(data_samples[i, "run"] == 1) {
    Slope = coef_df_0$estimate[1]
  } else {Slope = coef_df_0$estimate[2]}
  data_samples[i, "consumed_volume_0"] <- data_samples[i, "fluorescence"] / Slope
}




#### 2.3 Food consumption graph ####
# nice graph based on 0 intercept model to export at 800x600
# Define the labels and their positions
labels_df <- data.frame(treatment = c("high", "mid", "control", "low"),
                        label = c("a", "ab", "ab", "b"),
                        x = c(1, 2, 3, 4),
                        y = c(0.8, 0.8, 0.8, 0.8))
# Plot with labels (exported as )
ggplot(data_samples, aes(x = treatment, y = consumed_volume_0)) +
  geom_boxplot(fill = alpha("grey", 0.5), color = "black", notch = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, color = viridis(1)[1]) +
  labs(x = "Treatments", y = "Honey water consumption [μL]", title = "") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  geom_text(data = labels_df, aes(label = label, x = x, y = y, fontface = "bold"),
            hjust = -0.2, vjust = 0.5, size = 5)

group_by(data_samples, treatment) %>%
  summarise(
    count = n(),
    mean = mean(consumed_volume_0, na.rm = TRUE),
    sd = sd(consumed_volume_0, na.rm = TRUE),
    min = min(consumed_volume_0, na.rm = TRUE), 
    max = max(consumed_volume_0, na.rm = TRUE)
  )


#### 2.4 Fluoprescin Stats ####


# # Load the MASS package
# library(MASS)

# data is clearly not normally distributed. Log transformations did not work to get to a distribution suitable for lmer. 
# As the data is negative exponential distribution we use a glmer negative binomial model with thetha = 1 which coresponds to a negative exp distribution


mod <- glmer.nb(consumed_volume_0 ~ treatment + (1|colony) + (1|petridish) + (1|run), data = data_samples)
summary(mod)
Anova(mod)
cld(emmeans(mod, pairwise ~ "treatment", adjust = "tukey"), Letters = letters)

# 500 ppm ants consumed significantly less compared to controls but no other differences. 




# Food uptake quantification run 2 only. (means that this script should be updated/cleansed up at some point befor submission to peer reviewed journal)
data_run2 <- subset(data_samples, run == 2)

head(data_run2)
boxplot(data_run2$consumed_volume_0 ~ data_run2$treatment)

table(data_run2$treatment)

unique(levels(data_run2$treatment))

for (level in levels(data_run2$treatment)) {
  print(unique(data_run2[data_run2$treatment == level, "petridish"]))
}

# Filter dataset to exclude doubled petri dishes
exclude_petridish_values <- c(2.4, 2.9, 3.4, 3.9)
data_run2_filtered <- subset(data_run2, !(petridish %in% exclude_petridish_values))
table(data_run2_filtered$treatment)
boxplot(data_run2_filtered$consumed_volume_0 ~ data_run2_filtered$treatment)
#
group_by(data_run2_filtered, treatment) %>%
  summarise(
    count = n(),
    mean = mean(consumed_volume_0, na.rm = TRUE),
    sd = sd(consumed_volume_0, na.rm = TRUE),
    min = min(consumed_volume_0, na.rm = TRUE), 
    max = max(consumed_volume_0, na.rm = TRUE)
  )


mod <- glmer.nb(consumed_volume_0 ~ treatment + (1|petridish), data = data_run2_filtered)
summary(mod)
Anova(mod)
cld(emmeans(mod, pairwise ~ "treatment", adjust = "tukey"), Letters = letters)
residuals <- residuals(mod)
plot(residuals)
qqnorm(residuals)
qqline(residuals, col = "red")

labels_df <- data.frame(treatment = c("high", "mid", "control", "low"),
                        label = c("a", "ab", "ab", "b"),
                        x = c(1, 2, 3, 4),
                        y = c(0.8, 0.8, 0.8, 0.8))
# Plot with labels (exported as )
ggplot(data_run2_filtered, aes(x = treatment, y = consumed_volume_0)) +
  geom_boxplot(fill = alpha("grey", 0.5), color = "black", notch = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, color = viridis(1)[1]) +
  labs(x = "Treatments", y = "Honey water consumption [μL]", title = "") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  geom_text(data = labels_df, aes(label = label, x = x, y = y, fontface = "bold"),
            hjust = -0.2, vjust = 0.5, size = 5)


citation("survival")
citation("coxme")
citation("multcomp")









#### 3.  Neonic only survival analysis ####

exp1_data <- read.csv('flupy_survival_test.csv') #read table
#exp1_data$concentration   #### see what happens if we run it as a numerical variable. 
exp1_data$concentration <- factor(exp1_data$concentration) #transform to factor

#### 3.1 Stats neonic only ####
null_model <- coxme ( Surv (time = survival, event = censor) ~ 1                 + ( 1 | petri_dish) , data = exp1_data)
full_model <- coxme ( Surv (time = survival, event = censor) ~ 1 + concentration + ( 1 | petri_dish) , data = exp1_data)
anova(null_model   ,  full_model )
summary(full_model)

summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH"))
letters <- cld(summary(glht( full_model, linfct = mcp (concentration="Tukey")), test=adjusted("BH")))




#### 3.1 Plot neonic only ####

surv_plot <- survfit(Surv (time = survival, event = censor) ~ 1 + concentration , data=exp1_data) #basic plot
#legend <- c('0     - a', '0.5  - ab','1     - ab','5     - a','10   - ab','50   - ab','100  - b','500  - c')
legend <- c('0', '0.5','1','5','10','50','100','500')

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
# remove the grid lines
surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank()) +  
  theme(# Adjust font
    axis.title = element_text(size = 20),          # Font size of axis titles
    axis.text = element_text(size = 14),           # Font size of axis labels
    legend.title = element_text(size = 16),        # Font size of legend title
    legend.text = element_text(size = 14),         # Font size of legend labels
    axis.ticks = element_line(linewidth = 1),           # Size of axis tick marks
    strip.text = element_text(size = 14)           # Font size of facet labels (if applicable)
  )
surv_plot # export eg. as 800*600 














#### 4.  Flupy Fungus combined exposure effects on individual survival ####

getwd()

#read data and adjust data
{
antdata2 <- read.csv('flupy_fungus_individual_survival_data.csv')
antdata2$fungus <- as.factor(antdata2$fungus)
antdata2$petri_dish <- as.factor(antdata2$petri_dish)
antdata2$block <- as.factor(antdata2$block)
# concentration a proper number or factor
antdata2$concentration <- as.numeric(as.character(antdata2$concentration))
antdata2$concentration_factor <- as.factor(antdata2$concentration)
# creating a colony column, random factor must be controlled
antdata2 <- within(antdata2,colony <- substr(petri_dish,2,2))
antdata2$colony <- as.factor(antdata2$colony)

# create a new column "survival_14" and "censor_14 in antdata2 as we want to look at survival effects over 14 days
antdata2$survival_14 <- ifelse(antdata2$survival >= 14, 14, antdata2$survival)
antdata2$censor_14 <- ifelse(antdata2$survival <= 13, 1, 0)
}

antdata2$fungus <- factor(antdata2$fungus , levels=c("S","M"))

antdata2$colony

#### 4.3 FluFu stats ####
# Analyse dataset only for 14 days. 

#baseline model which is then updated to become the full interaction model
{
null_model          <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1                                 + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
flupy_model_f       <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration_factor          + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
fungus_model_f      <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration_factor + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
interaction_model_f <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration_factor * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
anova(null_model, flupy_model_f, fungus_model_f, interaction_model_f)
}

anova(fungus_model_f, interaction_model_f)
Anova(interaction_model_f)


# pairwise comparison
em_means <- emmeans(interaction_model_f, ~ concentration_factor * fungus)
letters <- c("a", "b", "c", "d", "e", "f")
# Compute the CLD (Comparisons of Least Squares Means) using the letters
cld_results <- cld(em_means, Letters = letters)
letters <- cld_results$.group

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

fungus_model_n      <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration + fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
interaction_model_n <- coxme ( Surv (time = survival_14, event = censor_14) ~ 1 + concentration * fungus + (1 | petri_dish) + (1 | colony) + (1 | block), data = antdata2)
anova(fungus_model_n, interaction_model_n)
summary(interaction_model_n)
Anova(interaction_model_n, type="III")
Anova(interaction_model_n)

contrast_matrix <- rbind("slope_fungusS"=c(1,0,0),"slope_fungusM"=c(1,0,1))
summary(glht(interaction_model_n,linfct=contrast_matrix),test=adjusted("BH"))






#### nice interaction plot ####

surv_mod <- survfit(Surv (time = survival_14, event = censor_14) ~ 1 + concentration + fungus, data=antdata2)
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
                             fungus_label = ifelse(fungus == "M", "fungus", "sham"),
                             concentration_label = ifelse(concentration == 0, "0 ppm",
                                                          ifelse(concentration == 5, "5 ppm", "50 ppm")))
letters <- gsub(" ", "", letters)
# Create the survival plot with the updated factor variables
surv_mod <- survfit(Surv(time = survival_14, event = censor_14) ~ 1 + concentration_label + fungus_label, data = antdata_updated)
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
                        legend = c(0.15, 0.25),
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
  geom_text(data = data.frame(x = rep(14.5, length(letters)), y = c(0.77, 0.73, 0.64, 0.51, 0.43, 0.32), label = letters),
            aes(x = x, y = y, label = label),
            color = "black",
            size = 5,
            hjust = 0.5,
            vjust = 0,
            fontface = "bold",
            inherit.aes = FALSE)
surv_plot














#### graph for the poster ####
# antdata2$fungus 
# antdata2$fungus <- factor(antdata2$fungus, levels = c("M", "S"))
# levels(antdata2$fungus) <- c("fungus", "sham")
# 
# antdata_updated <- transform(antdata2,
#                              fungus_label = ifelse(fungus == "M", "fungus", "sham"),
#                              concentration_label = ifelse(concentration == 0, "0 ppm",
#                                                           ifelse(concentration == 5, "5 ppm", "50 ppm")))
# 
# levels(antdata2$fungus)
# 
# color_palette <- c("black", "black", "darkgrey", "darkgrey", "lightgrey", "lightgrey")
# surv_mod <- survfit(Surv (time = survival_14, event = censor_14) ~ 1 + concentration + fungus, data=antdata2)
# # Reverse the linetype mapping for the "fungus" variable
# # antdata2$fungus <- factor(antdata2$fungus, levels = c("Uninfected", "Fungus"))
# surv_plot <- ggsurvplot(surv_mod, data = antdata2,
#                         linetype = "fungus",
#                         size = 1.5,
#                         xlab = 'Time (days)', ylab = 'Proportion Surviving',
#                         palette = color_palette, 
#                         ggtheme = theme_bw(),
#                         xlim = c(0, 14), break.time.by = 2,
#                         censor = FALSE,
#                         conf.int = TRUE,
#                         conf.int.alpha = 0.1, 
#                         legend = "none")
# surv_plot$plot <- surv_plot$plot + theme(panel.grid = element_blank())+  
#   theme(# Adjust font
#     axis.title = element_text(size = 24),          # Font size of axis titles
#     axis.text = element_text(size = 20),           # Font size of axis labels
#     legend.title = element_text(size = 16),        # Font size of legend title
#     legend.text = element_text(size = 14),         # Font size of legend labels
#     axis.ticks = element_line(linewidth = 1),           # Size of axis tick marks
#     strip.text = element_text(size = 14)           # Font size of facet labels (if applicable)
#   ) 
# surv_plot




#### plot #####

surviplot <- survfit(Surv (time = survival_14, event = censor_14) ~ 1 + concentration + fungus, data=antdata2)
aggregate(  censor ~ fungus + concentration, FUN=mean,data=antdata2)
# Choose a viridis palette by name

chosen_colors <- viridis(3, option = "inferno", begin = 0.2, end = 0.95)
color_vector <- rep(chosen_colors, each = 2)
# Create a plot to display the color squares
plot(1:3, rep(1, 3), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0.5, 3.5))
for (i in 1:3) {
  rect(i - 0.4, 0.8, i + 0.4, 1.2, col = chosen_colors[i], border = "white")
}
legend("topright", legend = c("Color 1", "Color 2", "Color 3"), fill = chosen_colors)
par(mar = c(5, 4, 4, 2) + 0.1)


#quicker way of displaying colors: 
library(scales)
show_col(chosen_colors, labels = TRUE, borders = NULL)


{
  surv_plot <- ggsurvplot(surviplot, data = antdata2,
                          pval = FALSE, pval.coord = c(6, 1),
                          linetype = "fungus",
                          lwd = 2,
                          xlab = 'Time (days)', ylab = 'Proportion Surviving',
                          palette = rep(color_vector, 2),
                          ggtheme = theme_bw(),
                          xlim = c(0, 15), break.time.by = 2,
                          panel.labs = list(fungus = c("SHAM", "FUNGUS")), 
                          short.panel.labs = TRUE,
                          panel.labs.font = list(face = "bold.italic", size = 11),
                          censor = FALSE,
                          conf.int = TRUE,
                          conf.int.alpha = 0.2
  )
  surv_plot$plot <- surv_plot$plot + 
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    scale_y_continuous(limits = c(0.2, 1), breaks = seq(0.2, 1, by = 0.2)) +
    theme(
      axis.title = element_text(size = 30),          # Font size of axis titles
      axis.text = element_text(size = 26),           # Font size of axis labels
      axis.ticks = element_line(linewidth = 2),      # Size of axis tick marks
      panel.border = element_rect(linewidth = 2),    # Border around the plot
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  print(surv_plot)
  # Convert ggsurvplot object to a grob
  surv_grob <- ggplotGrob(surv_plot$plot)
  # Save the grob with a transparent background
  ggsave("/Users/gismo/Desktop/survival_plot.png", surv_grob, bg = "transparent", dpi = 320, width = 10, height = 9)
}





