# rm(list = ls())
rm(list = setdiff(ls(), "first_time_use_working_directory"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ###  VITAL - FC1 Analyses   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####
#' This script was written by Nathalie Stroeymeyt and Daniel Schläppi
#' It contains the analysis of the binary food choice experiment conducted prior to the vital tracking experiment
#' Of 16 colonies 2 feeding events were recorded (one week apart) at each the ants were presented with a virus or a control food source
#' Position of virus food sources has been chosen pseudo randomly to avoid learning effects from 1th to 2nd feeding
#' The first data table contains the information on the number of ants present at each of the food sources at any given time 
#' The second table contains annotations of ant feeding durations (Covers only the first feeding event and only until 1h until after the second food source has been discovered (quite some variability)
#' For the second feeding event only the first 5 ants feeding on each food source were recorded.
#' 
#' Structure
#' 1. Prerequisites
#' 2. Part I - Analysis of the initial data set containing the number of ants present on the two food sources
#' 2.1 Data preparation
#' 2.2 Initial data exploration
#' 2.3 Modelling - Overall mean nr of ants per food source & per feeding session
#' 2.4 Modelling - Nr. of ants over time
#' 2.5 Plot nr. of ants over time
#' 2.6 First discovery
#' 2.7 Plot Nr of ants over time with shifted t0
#' 2.8 Main Analysis and randomization of first discovery
#' 3. Part II - Analysis of the second data set on feeding duratoins
#' 3.1 Data preparation
#' 3.2 Mean Duration first 5 feeding events per food source (Analysis + Plot)
#' 3.3 Analysis of completely annotated first feeding session
#' 3.3.1 Average feeding duration per food source
#' 3.3.2 Colony level analysis - Number of feeding events and overall feeding duration
#' 3.3.3 Exploitation of Food sources over time
#' 3.3.4 Plotting predicted exploitation of food and exploitation rate
#' 3.3.5 First feeding
#' 3.3.6 Randomization test



#### 1. Prerequisites ####

# choose directory - set it to the folder containing the data files (in my case vital_fc1)
set_directory <- function() {
  if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") {
    selected_dir <- tcltk::tk_choose.dir(default = "~/", caption = "Select Working Directory")
    if (is.null(selected_dir) || selected_dir == "") {
      cat("No directory selected. Exiting.\n")
      return()
    }
    setwd(selected_dir)
    first_time_use_working_directory <<- getwd()
    cat(crayon::blue(getwd()))
  } else {
    setwd(first_time_use_working_directory)
    cat(crayon::blue(getwd()))
  }
}
set_directory()



# load data
dat_numbers <- read.csv("vital_fc1_data_numbers.csv", header = TRUE)
dat_duration <- read.csv("vital_fc1_data_feedingdurations.csv", header = TRUE)
# remove 6 NA lines (artifact from annotations)
dat_duration <- dat_duration[!is.na(dat_duration$feeding_duration), ]


# libraries
# install.packages("pacman")
pacman::p_load(lubridate, plotrix, scales, car, lme4, Hmisc, 
               dplyr, tidyverse, blmeco, lmtest, lsmeans, lubridate,
               emmeans, multcompView, multcomp, viridis, crayon, 
               e1071, glmmTMB, DHARMa, merTools, tidyr, pheatmap, grid,
               progress, ggplot2)

# functions
sem <- function(x) {sd(x,na.rm=T)/sqrt(length(na.omit(x)))} # standard error of means
source("func_test_norm.R")  # adds test_norm() to the environment (it takes your model as an argument)

# parameters
RUN_ANALYSIS_AND_SIMULATIONS_NR <- TRUE
RUN_ANALYSIS_AND_SIMULATIONS_DURATION <- TRUE



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2. Part I - Analysis of the initial data set: Number of ant present on each of the food sources at a given time since discovery ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### 2.1 data preparation ####

# create time point vector (13 time points from t00, t05...t60)
Time_points_to_include <- 1:13 # full hour 
# Time_points_to_include <- 1:5  # include only first 20 minutes - the visible trend (more ants on virus food more pronounced in the first 20 minutes so it might make sense to look at a subset of the data)
total_iterations <- 500 # number of iterations in randomization loops (use 1000 or 10000 for publication)

# in the loaded data block corresponds to the eight days of video recordings -> needs to be renamed to video_session
# a new variable "block" gets created grouping the 4x4 colonies that were recorded as a group (twice)
# 1-4 instead of 1-8
dat_numbers <- rename(dat_numbers, video_session = block)
block_corrected <- NULL
for (i in 1:nrow(dat_numbers)) {
  block_old           <- dat_numbers[i, "video_session"]
  block_new       <- if(block_old<5) {block_old} else {block_old-4}
  block_corrected <- rbind(block_corrected, data.frame(block_new))
}
dat_numbers$block <- block_corrected$block_new

# create empty variables
dynamic_dat_same_t0_Nb      <- NULL
dynamic_dat_shifted_t0_Nb   <- NULL
dynamic_dat_same_t0_Diff      <- NULL
dynamic_dat_shifted_t0_Diff   <- NULL
# Short bit of info: T0 either defined as the time the first of the two food sources was discovered (same_t0) and or defined separately as the moment of discovery of each of the two food sources (different_t0)
# Nb = the number of ants present on each food, Diff = difference in the number of ants present on the two food sources


# fill up the created data frames
for (i in 1:nrow(dat_numbers)){ # i <- 1
  colony_id                          <- dat_numbers[i,"colony_id"]
  feeding_session                    <- dat_numbers[i,"feeding_session"]
  block                              <- dat_numbers[i,"block"]
  virus_position                     <- dat_numbers[i,"position_virus_corrected"]
  healthy_position                   <- c("left","right")[which(c("left","right")!=virus_position)]
  discovery_time_virus               <- period_to_seconds(hms(dat_numbers[i,paste("time_first_ant_",virus_position,sep="")]))
  discovery_time_healthy             <- period_to_seconds(hms(dat_numbers[i,paste("time_first_ant_",healthy_position,sep="")]))
  
  Nb_ant_virus_same_t0               <- as.numeric(dat_numbers[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(virus_position,1,1)  ,sep="_")])[Time_points_to_include]
  Nb_ant_healthy_same_t0             <- as.numeric(dat_numbers[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(healthy_position,1,1),sep="_")])[Time_points_to_include]
  
  if (discovery_time_virus<discovery_time_healthy){
    first_source_discovered          <- "virus"
    Nb_ant_virus_shifted_t0          <- Nb_ant_virus_same_t0
    Nb_ant_healthy_shifted_t0         <- as.numeric(dat_numbers[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
  }else{
    first_source_discovered          <- "healthy"
    Nb_ant_virus_shifted_t0          <- as.numeric(dat_numbers[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
    Nb_ant_healthy_shifted_t0        <- Nb_ant_healthy_same_t0
  }
  
  time_since_discovery               <- period_to_seconds(hms(dat_numbers[i,c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60")]))[Time_points_to_include]  - period_to_seconds(hms(dat_numbers[i,"time_first_ant"]))
  
  dynamic_dat_same_t0_Nb                <- rbind(dynamic_dat_same_t0_Nb, data.frame (colony_id,
                                                                                     feeding_session,
                                                                                     block,
                                                                                     virus_position,
                                                                                     discovery_time_virus,
                                                                                     discovery_time_healthy,
                                                                                     first_source_discovered,
                                                                                     time=rep(time_since_discovery,2),
                                                                                     status = rep(c("virus","healthy"),each=length(time_since_discovery)),
                                                                                     Nb_ants=c(Nb_ant_virus_same_t0,Nb_ant_healthy_same_t0),
                                                                                     stringsAsFactors = F))
  
  dynamic_dat_shifted_t0_Nb              <- rbind(dynamic_dat_shifted_t0_Nb, data.frame (colony_id,
                                                                                         feeding_session,
                                                                                         block,
                                                                                         virus_position,
                                                                                         discovery_time_virus,
                                                                                         discovery_time_healthy,
                                                                                         first_source_discovered,
                                                                                         time=rep(time_since_discovery,2),
                                                                                         status = rep(c("virus","healthy"),each=length(time_since_discovery)),
                                                                                         Nb_ants=c(Nb_ant_virus_shifted_t0,Nb_ant_healthy_shifted_t0),
                                                                                         stringsAsFactors = F))
  
  dynamic_dat_same_t0_Diff                <- rbind(dynamic_dat_same_t0_Diff, data.frame (colony_id,
                                                                                         feeding_session,
                                                                                         block,
                                                                                         virus_position,
                                                                                         discovery_time_virus,
                                                                                         discovery_time_healthy,
                                                                                         first_source_discovered,
                                                                                         time=time_since_discovery,
                                                                                         status = "virus_minus_healthy",
                                                                                         Delta_Nb_ants=Nb_ant_virus_same_t0-Nb_ant_healthy_same_t0,
                                                                                         stringsAsFactors = F))
  
  dynamic_dat_shifted_t0_Diff              <- rbind(dynamic_dat_shifted_t0_Diff, data.frame (colony_id,
                                                                                             feeding_session,
                                                                                             block,
                                                                                             virus_position,
                                                                                             discovery_time_virus,
                                                                                             discovery_time_healthy,
                                                                                             first_source_discovered,
                                                                                             time=time_since_discovery,
                                                                                             status = "virus_minus_healthy",
                                                                                             Delta_Nb_ants=Nb_ant_virus_shifted_t0-Nb_ant_healthy_shifted_t0,
                                                                                             stringsAsFactors = F))
  
}


#### 2.2 Initial data exploration ####
data <- dynamic_dat_same_t0_Nb
data$time_min <- data$time/60
head(data)
df1 <- group_by(data, colony_id, feeding_session, status) %>% 
  summarise(
    count = n(), 
    mean = mean(Nb_ants), 
    first_source_discovered = first(first_source_discovered)
  )
boxplot(df1$mean ~ df1$status) #plotting the two feeding sessions pooled (two means per colony)
par(mar=c(5.1, 4.5, 4.1, 1.8),
    cex.lab=1.3, 
    cex.axis=1.3, 
    cex.main= 1.3)

boxplot(df1$mean ~ df1$status + df1$feeding_session, main = "mean number of ants over time", ylab = "#ants", xlab = "food source", 
        pars  =  list(xaxt = "n"), ylim = c(0,5.6), cex.lab = 1.3, cex.axis = 1.3, cex.main =1.3)
text(c(1:4), -0.3, labels = c("control", "virus", "control","virus"), pos = 1, xpd = TRUE, cex = 1.3)
text(c(1.5,3.5), 5.4, labels = c("feeding session 1", "feeding session 2"), cex = 1.3)
text(c(1:4), 5, labels = c("a", "c", "b", "c"), font = 2, cex = 1.3) # stats from model below
abline(v=2.5, lty=4)

#### 2.3 Modelling - Overall mean nr of ants per food source & per feeding session  ####

mod <- lmer(mean ~ status*feeding_session + (1|colony_id) + (1|first_source_discovered) , data=df1 )
Anova(mod)
shapiro.test(residuals(mod)) # with 13 time points  not significant so fine (if using fewer data points datatransformation or another model is required)

# test model assumptions
summary(mod)
compareqqnorm(mod)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # just non significant --> assumption of normally distributed residuals is ok
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))
#pairwise differences
marginal = lsmeans(mod, ~ status*feeding_session, data = df1)
CLD = cld(marginal,
          alpha=0.5,
          Letters=letters,
          adjust="BH")
CLD

# # When using fewer time points (e.g. only the first 20 min the above model is should be replaced with glmer becasue of the distribution of residuals)
# mod <- glmer(mean ~ status*feeding_session + (1|colony_id) + (1|first_source_discovered), data=df1 , family = "poisson") 
# Anova(mod)
# #pairwise differences
# marginal = lsmeans(mod, ~ status*feeding_session, data = df1)
# CLD = cld(marginal,
#           alpha=0.5,
#           Letters=letters,
#           adjust="tukey")
# CLD

#### 2.4 Modelling - Nr of ants over time ####
# model glmer with time all time points included instead of means 
data$time <- as.factor(data$time)
data$feeding_session <- as.factor(data$feeding_session)
model <- glmer(Nb_ants ~ status + time + (1|colony_id) + (1|first_source_discovered) + (1|feeding_session), data=data, family = "poisson")
Anova(model)

#### 2.5 Plot nr. of ants over time ####
# mean number of ants over time

data_plot <- data %>% group_by(status, time_min) %>% 
  summarise(
    mean = mean(Nb_ants),
    sem = sem(Nb_ants)
  ) %>% as.data.frame(data_plot)

# smooth or polygon graph (smooth uses the smoothing function in ggplot and the polygon real data)
which_graph <- "polygon" # smooth or polygon
if (which_graph == "smooth") {
ggplot(data = data_plot, aes(x = time_min, y = mean, color = status)) +
  geom_point() +
  geom_smooth() +
  theme_light() +
  scale_color_viridis(end = 0.8,
                      name  ="food source",
                      breaks=c("healthy", "virus"),
                      labels=c("control", "virus"),
                      discrete = TRUE, option = "D") +
  ggtitle("Mean number of ants over time") +
  xlab("time [min]") +
  ylab("mean number of ants") +
  theme(text = element_text(size = 20))
} else if (which_graph == "polygon") {
  ggplot(data = data_plot, aes(x = time_min, y = mean, color = status)) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem), alpha = 0.2) +
    theme_bw() +
    scale_color_viridis(end = 0.8,
                        name  ="food source",
                        breaks=c("healthy", "virus"),
                        labels=c("control", "virus"),
                        discrete = TRUE, option = "D") +
    ggtitle("") +
    xlab("time [min]") +
    ylab("mean number of ants") +
    theme(text = element_text(size = 20))
} else {print("revisit the code ;-)")}


#### 2.6 First discovery ####
# get a new data frame with first discovery times 
data_plot2 <- data %>% group_by(colony_id, feeding_session, status) %>% 
  summarise(
    mean = mean(Nb_ants),
    first_discovered = first(first_source_discovered),
    discovery_time_control = first(discovery_time_healthy), 
    discovery_time_virus = first(discovery_time_virus)
  ) %>% as.data.frame(data_plot2)

data_plot2_wide <-  reshape(data = data_plot2, 
                            idvar=c("colony_id", "feeding_session", "first_discovered"),
                            v.names = "mean",
                            timevar="status",
                            direction = "wide")

# plot mean number of ants and show which of the food sources was discovered first
ggplot(data_plot2_wide) + 
  geom_segment(aes(x = 1, xend = 2, y = mean.healthy, yend = mean.virus, color = first_discovered)) +
  theme_classic() + 
  geom_point(aes(x = 1,  y = mean.healthy)) +
  geom_point(aes(x = 2, y = mean.virus)) +
  scale_x_continuous(name="food source", breaks = c(1,2), labels = c("control", "virus"), limits=c(0.8, 2.2)) +
  scale_color_viridis(end = 0.8,
                      name  ="first discovered",
                      labels=c("control", "virus"),
                      discrete = TRUE, option = "D") +
  ggtitle("Mean number of ants per food pair") +
  ylab("mean number of ants") +
  theme(text = element_text(size = 20))
#maybe something there? 

### histogram of first discovery
n_virus_discovered_first <- sum(data_plot2_wide$first_discovered == "virus")
n_max <- nrow(data_plot2_wide)
t <- table(data_plot2_wide$first_discovered)
barplot(t, main = "", xlab = "food source", ylab = "number of first discoveries", ylim = c(0,26), xaxt="n")
axis(1, at=c(0.7, 1.95), labels=c("control", "virus"))
segments(x0 = 0.7, y0 =24, x1 = 1.9)
text(x = 1.3, y = 25, label = "p = 0.07 (glmer bino)", font = 2, cex = 1.1) #stats see just below

#Is that significantly different from 50/50?
# prop.test(n_virus_discovered_first, n_max, alternative = "two.sided", p = 0.5)
# binom.test(n_virus_discovered_first, n_max, p = 0.5) # test if this is 50/50 selection or not
# these test not suitable because the two feeding sessions are not independant (repeated measurements of the same colonies
# needs to be confirmed with alternative test that takes this into account. 
# ! With the corrected data the virus food might be discovered first more often ! --> needs to be accounted for when analysing the exploitation of the two food sources
data_plot2_wide <- data_plot2_wide %>% 
  mutate(first_discovered_binomial = ifelse(first_discovered == "virus", 1, 0))
mod <-glmer(first_discovered_binomial ~ 1 + (1|feeding_session) + (1|colony_id), family = binomial, data = data_plot2_wide)
summary(mod)

# there might be a trend but it is just not significant... 
ggplot(data_plot2_wide, aes(x = factor(feeding_session), fill = factor(first_discovered_binomial))) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("0" = "#CCEBC5", "1" = "#FBB4AE"), 
                    labels = c("0" = "Control", "1" = "Virus")) +
  labs(x = "Feeding Session", y = "Count of first discovery", fill = "First Discovered") +
  scale_y_continuous(limits = c(0, 18), breaks = c(0, 4, 8, 12, 16)) +
  theme_classic() +
  annotate("text", x = 1.9 , y = 18, label = "p = 0.07 (glmer bino)", 
           hjust = 1.1, vjust = 2, size = 4)

### Do colonies always have the same preference? Decision fidelity? 
## Plot heat map, requires data transformation
heatmap_data <- data_plot2_wide %>%
  dplyr::select(colony_id, feeding_session, first_discovered_binomial) %>%
  pivot_wider(names_from = feeding_session, values_from = first_discovered_binomial) %>%
  tibble::column_to_rownames(var = "colony_id") %>%
  as.matrix()
pheatmap(heatmap_data, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("#CCEBC5", "#FBB4AE"),
         legend = FALSE,
         labels_col = c("Week 1", "Week 2"),
         show_rownames = TRUE,
         show_colnames = TRUE)
# No trend visible, would probably require more data

#### 2.6.1 Time of first discovery #### 
# relabel status
data_plot2_time <- data_plot2 %>%
  rename(food_source = status) %>%
  mutate(food_source = recode(food_source, 
                              "healthy" = "control", 
                              "virus" = "virus"))
# Reshape the 'discovery_time_control' and 'discovery_time_virus' columns into one column discovery time
data_plot2_time <- data_plot2_time %>%
  mutate(discovery_time = ifelse(food_source == "control", 
                                 discovery_time_control, 
                                 discovery_time_virus)) %>%
  dplyr::select(-discovery_time_control, -discovery_time_virus)

# Stats
# time until discovery of each food source
mod <- lmer(sqrt(discovery_time) ~ food_source + (1|colony_id), + (1|feeding_session), data = data_plot2_time)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow=c(1,1))
test_norm(mod)
hist(resid(mod), breaks = 10)
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

# Plot 
ggplot(data_plot2_time, aes(x = food_source, y = discovery_time, fill = food_source)) +
  geom_boxplot(outlier.shape = NA) + #remove outliers to avoid double plotting with the plotting of individual data points.
  geom_jitter(width = 0.3, color = "black", alpha = 0.2) +
  labs(title = "Discovery Time",
       x = "Food Source",
       y = "Time (s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +

  geom_segment(aes(x = 1.1, xend = 1.9, y = 0.9 * max(discovery_time, na.rm = TRUE),
                   yend = 0.9 * max(discovery_time, na.rm = TRUE)),
               color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 0.95 * max(data_plot2_time$discovery_time, na.rm = TRUE),
           label = "p = 0.037 (lmer)", color = "black", size = 3)








#### 2.7 Plot Nr of ants over time with shifted t0 ####
# instead of counting time for both food sources starting when the first one is discovered we now look at the data with discovery time for each food sources separately
# thereby minimizing effect of first discovery (but not removing it!!!) and showing the initial recruitment to the food sources.
data_shifted <- dynamic_dat_shifted_t0_Nb
data_shifted$time_min <- data_shifted$time/60
data_plot3 <- data_shifted %>% group_by(status, time_min) %>% 
  summarise(
    mean = mean(Nb_ants), 
    sem = sem(Nb_ants)
  ) %>% as.data.frame()
ggplot(data = data_plot3, aes(x = time_min, y = mean, color = status)) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem), alpha = 0.2) +
  theme_light() +
  scale_color_viridis(end = 0.8,
                      name  ="food source",
                      breaks=c("healthy", "virus"),
                      labels=c("control", "virus"),
                      discrete = TRUE, option = "D") +
  ggtitle("                       initial recruitment") +
  xlab("time [min]") +
  ylab("mean number of ants")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2.8 Main Analysis initial Data ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

if (RUN_ANALYSIS_AND_SIMULATIONS_NR) { # RUN_ANALYSIS_AND_SIMULATIONS_NR <- TRUE
 for (time_origin in c("same","shifted")){ # time_origin <- "same"
  ### get dynamic data
  dynamic_dat_Nb      <- get(paste("dynamic_dat",time_origin,"t0","Nb",sep="_")) 
  dynamic_dat_Diff    <- get(paste("dynamic_dat",time_origin,"t0","Diff",sep="_")) 
  
  ### calculate means and SE create dynamic_summary_dat_Nb objects
  dynamic_summary_dat_Nb    <- data.frame(as.matrix(aggregate(Nb_ants   ~feeding_session+status+time,function(x)cbind(mean(x),std.error(x)),data=dynamic_dat_Nb)),stringsAsFactors = F)
  names(dynamic_summary_dat_Nb)[which(grepl(".1",names(dynamic_summary_dat_Nb)))] <- gsub(".1","_Mean",names(dynamic_summary_dat_Nb)[which(grepl(".1",names(dynamic_summary_dat_Nb)))])
  names(dynamic_summary_dat_Nb)[which(grepl(".2",names(dynamic_summary_dat_Nb)))] <- gsub(".2","_SE",names(dynamic_summary_dat_Nb)[which(grepl(".2",names(dynamic_summary_dat_Nb)))])
  dynamic_summary_dat_Nb$Nb_ants_Mean <- as.numeric( dynamic_summary_dat_Nb$Nb_ants_Mean)
  dynamic_summary_dat_Nb$Nb_ants_SE   <- as.numeric( dynamic_summary_dat_Nb$Nb_ants_SE)
  dynamic_summary_dat_Nb <- dynamic_summary_dat_Nb[order(dynamic_summary_dat_Nb$time),]
  
  ### calculate means and SE create dynamic_summary_dat_Diff objects
  dynamic_summary_dat_Diff    <- data.frame(as.matrix(aggregate(Delta_Nb_ants   ~feeding_session+status+time,function(x)cbind(mean(x),std.error(x)),data=dynamic_dat_Diff)),stringsAsFactors = F)
  names(dynamic_summary_dat_Diff)[which(grepl(".1",names(dynamic_summary_dat_Diff)))] <- gsub(".1","_Mean",names(dynamic_summary_dat_Diff)[which(grepl(".1",names(dynamic_summary_dat_Diff)))])
  names(dynamic_summary_dat_Diff)[which(grepl(".2",names(dynamic_summary_dat_Diff)))] <- gsub(".2","_SE",names(dynamic_summary_dat_Diff)[which(grepl(".2",names(dynamic_summary_dat_Diff)))])
  dynamic_summary_dat_Diff$Delta_Nb_ants_Mean <- as.numeric( dynamic_summary_dat_Diff$Delta_Nb_ants_Mean)
  dynamic_summary_dat_Diff$Delta_Nb_ants_SE   <- as.numeric( dynamic_summary_dat_Diff$Delta_Nb_ants_SE)
  dynamic_summary_dat_Diff <- dynamic_summary_dat_Diff[order(dynamic_summary_dat_Diff$time),]
  
  ### calculate one mean number of ants for each source, each colony and each feeding session
  summary_dat_Nb      <- aggregate(Nb_ants   ~colony_id+feeding_session+block+virus_position+discovery_time_virus+discovery_time_healthy+first_source_discovered+status,FUN=mean,data=dynamic_dat_Nb)
  
  ### stats
  print( paste("Statistics -",time_origin,"T0"))
  model <- lmer(Nb_ants ~ status*feeding_session + (1|colony_id) + (1|first_source_discovered), data=summary_dat_Nb )
  print(Anova(model))
  print(shapiro.test(residuals(model)))
  
  ### plot  TO IMPROVE: SHOULD BE MEANS AND STANDARD ERRORS RATHR THAN BOXPLOT
  boxplot(Nb_ants~status+feeding_session,data=summary_dat_Nb)
  
  ### NOTE: As feeding session is highly non-significant, we don't need to distinguish between the two in the analysis of randomisation tests
  
  ### Plot Dynamic data For Each Session
  for (session in c(1,2)){ # session <- 1
    ymin <- min( dynamic_summary_dat_Nb$Nb_ants_Mean-dynamic_summary_dat_Nb$Nb_ants_SE )
    ymax <- max( dynamic_summary_dat_Nb$Nb_ants_Mean+dynamic_summary_dat_Nb$Nb_ants_SE )
    
    plot(Nb_ants_Mean ~ time, pch=16, col="red",
         data=dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),],
         ylim=c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin)),bty="l", 
         main =paste("Session",session,"-", time_origin ,"t0"))
    polygon( x= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"time"]
                   ,
                   rev(dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"time"]) ) , 
             y= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_Mean"] - dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_SE"]  
                   ,
                   rev( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_Mean"] + dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_SE"])),
             col=alpha("red",0.5),border=NA
    )
    points(Nb_ants_Mean ~ time, pch=16, col="blue",
           data=dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),])
    polygon(  x= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"time"]
                    ,
                    rev(dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"time"]) ) , 
              y= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_Mean"] - dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_SE"]  
                    ,
                    rev( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_Mean"] + dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_SE"])),
              col=alpha("blue",0.5),border=NA
    )
    legend ("topright",legend=c("virus","healthy"),col=c("red","blue"),pch=16,bty="n")
    
  }
  
  ### Randomisation test for each session
  ### Aim: test whether the results are an artifact from virus food tending to be discovered earlier than healthy food
  ### H0: there is no difference in ant behaviour towards either virus or healthy food, the only relevant point for recruitment is which food was discovered first
  ### In the observed dataset the virus was discovered first a certain number of time
  ### The randomisation will consist in randomly drawing the same number of colonies from the total and saying that in these colonies the first source had virus (irrespective of truth), and in all others the first source was healthy irrespective of truth)
  ### If observed data is an artifact, we should observe a similar difference in number of ants between the two sources in randomised and real data
  
  ### Pre-calculations: nb of colonies that found the virus first in each session, and observed difference in number of ants across all colonies and all times
  nb_discoveries <- aggregate(colony_id~first_source_discovered+feeding_session ,FUN=length, data= summary_dat_Nb[which(summary_dat_Nb$status=="healthy"),])
  observed_Diff    <- aggregate(Delta_Nb_ants   ~ status,FUN=mean,data=dynamic_dat_Diff)
  
  random_data_dynamic_Nb <- NULL
  random_data_Diff       <- NULL 
  
  pb <- progress_bar$new(
    format = "Progress Randomisation: :current/:total [:bar] :percent ETA: :eta",
    total = total_iterations,
    clear = FALSE,
    width = 60
  )
  
  for (i in 1 : total_iterations){ ### randomisation loop # i <- 1
    rand_Nb   <- NULL
    rand_Diff <- NULL
    ### randomise separately for each session
    for (session in c(1,2)){ # session <- 1
      nb_virus       <- nb_discoveries[which(nb_discoveries$first_source_discovered=="virus"&nb_discoveries$feeding_session==session),"colony_id"]
      subset_Nb   <- dynamic_dat_Nb[which(dynamic_dat_Nb$feeding_session==session),]
      subset_Diff <- dynamic_dat_Diff[which(dynamic_dat_Diff$feeding_session==session),]
      colony_list <- unique(subset_Nb$colony_id)
      
      hypothetic_virus_first    <- sample(colony_list,size=nb_virus,replace = F)
      
      ### now create first_source_discovered_RAND column according to hypothetical status
      rand_subset_Nb                                                                                           <- subset_Nb             ; rand_subset_Diff                                                                                             <- subset_Diff  
      rand_subset_Nb$first_source_discovered_RAND                                                              <- "healthy"             ; rand_subset_Diff$first_source_discovered_RAND                                                                <- "healthy" ;
      rand_subset_Nb[which(rand_subset_Nb$colony_id%in%hypothetic_virus_first),"first_source_discovered_RAND"] <- "virus"               ; rand_subset_Diff[which(rand_subset_Diff$colony_id%in%hypothetic_virus_first),"first_source_discovered_RAND"] <- "virus"
      
      ### for Nb: copy status column into a new column status_RAND, and switch the values for those colonies who have been assigned a different first source discovered than in the observed data
      rand_subset_Nb$status_RAND  <- rand_subset_Nb$status                                                                           <- rand_subset_Nb$status 
      rand_subset_Nb[which(rand_subset_Nb$first_source_discovered_RAND!=rand_subset_Nb$first_source_discovered &  rand_subset_Nb$status=="virus"),"status_RAND"] <- "healthy"
      rand_subset_Nb[which(rand_subset_Nb$first_source_discovered_RAND!=rand_subset_Nb$first_source_discovered &  rand_subset_Nb$status=="healthy"),"status_RAND"] <- "virus"
      
      ### for Diff: multiply the rows of colonies who have been assigned a different first source discovered than in the observed data by minus 1
      rand_subset_Diff[which(rand_subset_Diff$first_source_discovered_RAND!=rand_subset_Diff$first_source_discovered),"Delta_Nb_ants"] <- (-1)*rand_subset_Diff[which(rand_subset_Diff$first_source_discovered_RAND!=rand_subset_Diff$first_source_discovered),"Delta_Nb_ants"] 
      
      ### concatenate and store
      rand_Nb   <- rbind(rand_Nb  ,data.frame(RAND=i, feeding_session=session, rand_subset_Nb))
      rand_Diff <- rbind(rand_Diff,data.frame(RAND=i, feeding_session=session, rand_subset_Diff))
    }
    
    ### then we use aggregate to get the mean number of ants at each time point in the hypothetical scenario, across all colonies and both feeding events
    mean_rand_dat_Nb <- aggregate(Nb_ants ~ status_RAND + time, FUN=mean, data=rand_Nb)
    
    ### and we use aggregate to calculate the mean delta nb ant across all times,  all colonies and both feeding events in the hypothetical scenario
    mean_rand_dat_Diff <- aggregate(Delta_Nb_ants ~  1, FUN=mean, data=rand_Diff)
    
    ### and we concatenate
    random_data_dynamic_Nb <- rbind(random_data_dynamic_Nb,data.frame(RAND=i, feeding_session=session, mean_rand_dat_Nb))
    random_data_Diff       <- rbind(random_data_Diff      ,data.frame(RAND=i, feeding_session=session, mean_rand_dat_Diff))
    pb$tick()
    }
  
  ### Plot expected vs observed, Delta ants
  xmin <- min (c(random_data_Diff$Delta_Nb_ants),observed_Diff[,"Delta_Nb_ants"])
  xmax <- max (c(random_data_Diff$Delta_Nb_ants),observed_Diff[,"Delta_Nb_ants"])
  hist(random_data_Diff$Delta_Nb_ants,col=alpha("grey",0.5),xlim=c(xmin,xmax),main=paste("Observed vs Expected,",capitalize(time_origin),"T0"), xlab=expression(paste(Delta, " number of ants")))
  arrows(x0=observed_Diff[,"Delta_Nb_ants"],
         y0=0,
         y1=-10,
         code=1,
         col="black",
         lwd=4,
         length=0.1)
  
  ### P value
  if (observed_Diff[,"Delta_Nb_ants"]>median(random_data_Diff$Delta_Nb_ants)){
    pval <- 2*length(which(random_data_Diff$Delta_Nb_ants>=observed_Diff[,"Delta_Nb_ants"]))/length(random_data_Diff$Delta_Nb_ants)
    mtext(paste("p=",pval),3, cex = 1.3)
  }else{
    pval <- 2*length(which(random_data_Diff$Delta_Nb_ants<=observed_Diff[,"Delta_Nb_ants"]))/length(random_data_Diff$Delta_Nb_ants)
    mtext(paste("p=",pval),3, cex = 1.3)
  }
 }
}

### please note due to the randomization the provided p-value is always a bit different. 
# There is a trend. For shifted t0 it sometimes is significant especially if only looking at the first 20 minutes instead of a full hour.
# for the manuscript lets use 10000 iterations: same_t0 -- p = 0.1044, shifted_t0 -- p = 0.0878









#### 3. PART II - Analysis of annotated Feeding Events i.e. Feeding Duration ####


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3. PART II - Analysis of feeding duration     #### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# For the second feeding event (week 2) onyl the first five feeding events were annotated - data needs to be subsetted accordingly for different parts of the analysis

#### 3.1 Data preparation ####

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


#### 3.2 Mean Duration first 5 feeding events per food source ####

# mean duration of feeding  events based on the two food sources just for the first 5 events per food source
# subset based on event_id (for each event a number saying this was feeding nr X on this food source for left or right)
dat_duration <- dat_duration %>% # extract the number from event_id for new variable
  mutate(event_number = as.numeric(gsub("^[^0-9]*(\\d+).*", "\\1", event_id))) %>% 
  filter(!is.na(event_number))  
dat_duration_sub <- subset(dat_duration, event_number <= 5 )

dat_duration_sub %>%
  group_by(food_source) %>%
  summarize(
    mean_feeding_duration = mean(feeding_duration_seconds, na.rm = TRUE),
    sd_feeding_duration = sd(feeding_duration_seconds, na.rm = TRUE)
  )

# Worth investigating, but feeding duration might not necessarily reflect feeding volume very well... 
# Some ants are full after 30 sec others take 5 min to get full, might depend more on other variables e.g. like consistency of food (which was mash thus not a homogeneous liquid)

# lmm <- lmer(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block) + (1|feeding_session), data = dat_duration_sub)
# summary(lmm)
# Anova(lmm)
#' Not the right model given the distribution of the residuals--> consider changing to glmer with poisson or antother test.
#' We have repeated measures of the same colonies, data not normally distributed --> non-parametric test for repeated measures. 
#' The Friedmann test would be suitable (compares multiple paired groups) alternative to repeated measures Anova.
#' Arrange data and run test: 
mean_feeding_duration <- dat_duration_sub %>%
  group_by(colony_id, feeding_session, food_source) %>%
  summarise(mean_duration = mean(feeding_duration_seconds, na.rm = TRUE)) %>%
  ungroup()
wide_data <- mean_feeding_duration %>%
  pivot_wider(names_from = c(feeding_session, food_source), values_from = mean_duration)
friedman_data <- as.matrix(wide_data[,-1]) # turn into matrix without id column
friedman.test(friedman_data)


### Plot 
ggplot(dat_duration_sub, aes(x = food_source, y = feeding_duration_seconds, fill = food_source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, color = "black", alpha = 0.1) +
  labs(
    title = "Feeding Duration by Food Source (first 5 feeding events only)",
    x = "Food Source",
    y = "Feeding Duration (seconds)"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(dat_duration_sub$feeding_duration_seconds, na.rm = TRUE), 
                 yend = 1.05 * max(dat_duration_sub$feeding_duration_seconds, na.rm = TRUE)),
             color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 1.1 * max(dat_duration_sub$feeding_duration_seconds, na.rm = TRUE), 
           label = "p = 0.54 (Friedmann test)", color = "black", size = 3)

# no difference regarding mean duration but from watching the videos again I am confident that the initial analysis does not do justice to the data (is only scratching the surface)
# Thus, following the above, for the first feeding event week 1 each and every feeding event that happened in the first 60 minutes of discovery of a food got annotated









#### 3.3 Analysis of completely annotated first feeding session ####
dat_duration_first <- subset(dat_duration, feeding_session == 1)

# consider capping any feeding event longer than 5 minutes to limit their impact...
# feeding annotations done visually with start and end being define as the moment mouth parts connected or disconnected with the food source)
# some ants that were feeding for extended periods appeared to become inactive for a while before disconnecting (mouth parts connected but actual food intake might have been stopped - actual food uptake not measurable with this simple set up)

cap_threshold <- 5 * 60  # capping threshold in seconds
dat_duration_first$feeding_duration_seconds_capped <- ifelse( 
  dat_duration_first$feeding_duration_seconds > cap_threshold,
  cap_threshold,
  dat_duration_first$feeding_duration_seconds)
# for now, all following analyses were run with with normal (uncapped) duration, the capping can probably be ignored as its impact is not too big... but it might help to make some graphs nicer... 

#### 3.3.1 Average feeding duration ####
mod <- glmer(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block), data = dat_duration_first, family = "poisson")
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
                  data = dat_duration_first, family = nbinom2)
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

### Plot
ggplot(dat_duration_first, aes(x = food_source, y = feeding_duration_seconds_capped, fill = food_source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, color = "black", alpha = 0.1) +
  labs(title = "Feeding Duration by Food Source",
       x = "Food Source",
       y = "Feeding Duration (seconds)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05 * max(feeding_duration_seconds_capped, na.rm = TRUE), 
                   yend = 1.05 * max(feeding_duration_seconds_capped, na.rm = TRUE)),
               color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 1.1 * max(dat_duration_first$feeding_duration_seconds_capped, na.rm = TRUE), 
           label = "*** (glmer nbinom2)", color = "black", size = 3)

dat_duration_first %>%
  group_by(food_source) %>%
  summarize(
    mean_fd = mean(feeding_duration_seconds_capped, na.rm = TRUE),
    median_fd = median(feeding_duration_seconds_capped, na.rm = TRUE),
    sd_fd = sd(feeding_duration_seconds_capped, na.rm = TRUE)
  ) %>% as.data.frame()




#### 3.3.1.5  Feeding duration in dependence of feeding time ####
# are feeding events different depending on when they happen? Checking for potential effects of food freshness or water content and such on feeding duration

# transform time to seconds
dat_duration_first$feeding_start_seconds <- period_to_seconds(hms(dat_duration_first$feeding_start))
dat_duration_first$feeding_end_seconds <- dat_duration_first$feeding_start_seconds + dat_duration_first$feeding_duration_seconds
dat_duration_first$feeding_end_seconds_capped <- dat_duration_first$feeding_start_seconds + dat_duration_first$feeding_duration_seconds_capped

mod <- lmer(log10(feeding_duration_seconds) ~ food_source*feeding_start_seconds + (1|colony_id), data = dat_duration_first)
Anova(mod)
summary(mod)
test_norm(mod) 
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

mod <- lmer(feeding_start_seconds ~ food_source + (1|colony_id), data = dat_duration_first)
Anova(mod)
summary(mod)
test_norm(mod) 

### plot of feeding durations over time
plot(feeding_duration_seconds_capped ~ feeding_start_seconds, data = dat_duration_first,
     main = "Feeding Duration vs Feeding Start Time",
     xlab = "Feeding Start Time (seconds)",
     ylab = "Feeding Duration (seconds)",
     pch = 16,
     col = rgb(t(col2rgb(viridis::viridis(2)[as.factor(dat_duration_first$food_source)]))/255, 
               alpha = 0.5),
     cex = 0.8)  # Optional: increase size of points for better visibility
abline(lm(feeding_duration_seconds_capped ~ feeding_start_seconds, data = dat_duration_first), col = "red", lwd = 2)
abline(lm(feeding_duration_seconds_capped ~ feeding_start_seconds, 
          data = dat_duration_first[dat_duration_first$food_source == "control", ]), 
       col = "blue", lwd = 2)
abline(lm(feeding_duration_seconds_capped ~ feeding_start_seconds, 
          data = dat_duration_first[dat_duration_first$food_source == "virus", ]), 
       col = "orange", lwd = 2)
abline(v=fixef(mod)["(Intercept)"], lty=4, col = "blue", lwd = 3); abline(v=fixef(mod)["(Intercept)"] + fixef(mod)["food_sourcevirus"], lty=4, col = "orange", lwd = 3)
legend("topright", legend = levels(as.factor(dat_duration_first$food_source)), pch = 16,
       col = rgb(t(col2rgb(viridis::viridis(2)))/255, alpha = 0.5), 
       title = "Food Source")

plot(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration_first,
     main = "Feeding Duration vs Feeding Start Time",
     xlab = "Feeding Start Time (seconds)",
     ylab = "Feeding Duration (seconds)",
     pch = 16,
     col = rgb(t(col2rgb(viridis::viridis(2)[as.factor(dat_duration_first$food_source)]))/255, 
               alpha = 0.5),
     cex = 0.8)  # Optional: increase size of points for better visibility
abline(lm(feeding_duration_seconds ~ feeding_start_seconds, data = dat_duration_first), col = "red", lwd = 2)
abline(lm(feeding_duration_seconds ~ feeding_start_seconds, 
          data = dat_duration_first[dat_duration_first$food_source == "control", ]), 
       col = "blue", lwd = 2)
abline(lm(feeding_duration_seconds ~ feeding_start_seconds, 
          data = dat_duration_first[dat_duration_first$food_source == "virus", ]), 
       col = "orange", lwd = 2)
abline(v=fixef(mod)["(Intercept)"], lty=4, col = "blue", lwd = 3); abline(v=fixef(mod)["(Intercept)"] + fixef(mod)["food_sourcevirus"], lty=4, col = "orange", lwd = 3)
legend("topright", legend = levels(as.factor(dat_duration_first$food_source)), pch = 16,
       col = rgb(t(col2rgb(viridis::viridis(2)))/255, alpha = 0.5), 
       title = "Food Source")





#### 3.3.2 Colony level analysis - Number of feeding events and overall feeding duration ####

dat_summary <- dat_duration_first %>%
  group_by(food_source, colony_id) %>%
  summarise(
    num_feeding_events = n(),
    total_feeding_duration = sum(feeding_duration_seconds, na.rm = TRUE)
  ) %>% as.data.frame()

# number of feeding events per colony per food source 
mod <- lmer(num_feeding_events ~ food_source + (1|colony_id) , data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow = c(1,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #not significant --> good 

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

#### 3.3.2.1 Feeding rate (feeding events per min) ####
### Repeat of this colony level analysis once more but this time capping to 1h of annotations otherwise there is a bias towards the earlier discovered food
#' Note: Annotations were not done for the complete 2h of recordings per colony.
#' Instead, annotations were done until 1h after the second food was discovered.
#' So if for example one food was discovered min 10 and one at min 30 we would get 60 min of annotations until min 90 for the food discovered last
#' However, for the other food that was discovered earlier we would have 90 minutes of annotations.
#' This will thus create a bias if one food source was consistently discoverd earlier/later which now we know was the case.
#' so by looking at only 60 minutes with shifted_discovery we take away some of the bias but there is still the natural bias towards the first discovered food
#' such as pheromone trails and effects of previous feedings (reduced hunger etc) on subsequent feedings

#' --> look at a number of feedings as well as the rate of feeding events during 1h post discovery of each food so rate per min in 60 minutes after a food has been discoverd

dat_duration_first$feeding_annotation_end_seconds <- period_to_seconds(hms(dat_duration_first$feeding_annotation_end))

summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id, food_source) %>%
  summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
subsetted_data <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id", "food_source"))
# filter to only include feeding events starting within the 60 min (3600s) window since discovery (first actual feeding) of the food source of interest (different discovery times)
subsetted_data <- subsetted_data %>%
  mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
  mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>% 
  filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)

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
mod <- lmer(sqrt(num_feeding_events_per_min) ~ food_source + (1|colony_id), data = dat_summary)
summary(mod)
print(Anova(mod, type=3))
test_norm(mod)
compareqqnorm(mod); par(mfrow=c(1,1))
hist(resid(mod))
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

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








#### 3.3.3 Exploitation of Food sources over time #### 
# plot a time curve showing the exploitation of the two food sources
# plot a the rate of exploitation over time 
# create a full grid with seconds, food_source and colonies and add in the the original data:

# transform time to seconds
dat_duration_first$feeding_start_seconds <- period_to_seconds(hms(dat_duration_first$feeding_start))
dat_duration_first$feeding_end_seconds <- dat_duration_first$feeding_start_seconds + dat_duration_first$feeding_duration_seconds_capped



dat_duration_first_modified <- dat_duration_first %>%
  group_by(colony_id) %>%
  mutate(
    first_feeding_virus     = min(feeding_start_seconds[food_source == "virus"], na.rm = TRUE),
    first_feeding_control   = min(feeding_start_seconds[food_source == "control"], na.rm = TRUE),
    first_source_fed_on = ifelse(
      first_feeding_virus < first_feeding_control, "virus",
      ifelse(first_feeding_control < first_feeding_virus, "control", NA_character_))) %>% 
  ungroup() %>% as.data.frame()
# get colony metadata
colony_metadata <- dat_duration_first_modified %>%
  dplyr::select(colony_id, treatment, position_virus, block, 
                first_feeding_virus, first_feeding_control, first_source_fed_on) %>%
  distinct(colony_id, .keep_all = TRUE)

# run the below for shifted or non shifted time points....
for (start_time in c("start_experiment", "discovery_first_food", "shifted_discovery")) { # keep shifted for last...# start_time <- "shifted_discovery"
  if (start_time == "start_experiment") {
    subsetted_data <- dat_duration_first
  }
  if (start_time == "shifted_discovery") { 
    summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id, food_source) %>%
      summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
    subsetted_data <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id", "food_source"))
    # filter to only include feeding events starting within the 60 min (3600s) window since discovery (first actual feeding) of the food source of interest (different discovery times)
    subsetted_data <- subsetted_data %>%
      mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
      mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>% 
      filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)
  }
  if (start_time == "discovery_first_food") {
    summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id) %>%
      summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
    subsetted_data <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id"))
    # filter to only include first hour since discovery of first food. 
    subsetted_data <- subsetted_data %>%
      mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
      mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>%
      filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)
    }
  
  grid <- expand_grid(
    colony_id = unique(subsetted_data$colony_id), 
    food_source = unique(subsetted_data$food_source),
    time = seq(0, max(subsetted_data$feeding_end_seconds, na.rm = TRUE))
  ) %>% as.data.frame()
  
  cumulative_explotation_over_time <- NULL
  
  # loop over colony_id and food_source
  for (colony in unique(grid$colony_id)) { 
    for (food in unique(grid$food_source)) {
      # filter grid and feeding data for the current colony and food source
      current_grid <- grid %>% filter(colony_id == colony, food_source == food)
      current_feedings <- subsetted_data %>% filter(colony_id == colony, food_source == food)
      
      # compute across all rows of current_grid using vectorized operations
      if (nrow(current_feedings) > 0) {
        current_feedings <- current_feedings %>% arrange(feeding_start_seconds) # Sort by start time
        
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
      cumulative_explotation_over_time <- bind_rows(cumulative_explotation_over_time, current_grid) # append to cumulative_exploitation_over_time
    }
  }
  
  # aggregate data and calculate means of duration for plotting
  mean_exploitation_time <- cumulative_explotation_over_time %>%
    group_by(time, food_source) %>%
    summarize(mean_time = mean(cumulated_exploitation_time, na.rm = TRUE),
              sd = sd(cumulated_exploitation_time, na.rm = TRUE), 
              sem = sem(cumulated_exploitation_time))
  title_lab_1 <- paste0("Cumulative Food Exploitation (T0 = ", start_time, ")")
  print(
    ggplot(mean_exploitation_time, aes(x = time, y = mean_time, color = food_source)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean_time - sem, ymax = mean_time + sem, fill = food_source), alpha = 0.2) +
    labs(
      title = title_lab_1,
      x = "Time",
      y = "Mean Cumulated Exploitation Time (seconds)") +
    scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
    theme_bw() +
    theme(legend.title = element_blank())+
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  )
  
  # aggregate data and calculate means of exploitation rate for plotting
  mean_exploitation_rate <- cumulative_explotation_over_time %>%
    group_by(time, food_source) %>%
    summarize(mean_exploitation = mean(exploitation_rate, na.rm = TRUE),
              sd = sd(exploitation_rate, na.rm = TRUE))
  title_lab_2 <- paste0("Exploitation Rate (T0 = ", start_time, ")")
  print(
    ggplot(mean_exploitation_rate, aes(x = time, y = mean_exploitation, color = food_source)) +
      geom_line(size = 1) +
      labs(
        title = title_lab_2,
        x = "Time",
        y = "Mean Cumulated Exploitation Time (seconds)") +
      scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
      theme_bw() +
      theme(legend.title = element_blank())+
      theme(legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  )
  # calculate means and SE
  dynamic_summary_dat_dur <- cumulative_explotation_over_time %>%
    group_by(food_source, time) %>%
    summarize(exploitation_rate_Mean   = mean(exploitation_rate, na.rm = TRUE),
              exploitation_rate_SE = sd(exploitation_rate, na.rm = TRUE) / sqrt(sum(!is.na(exploitation_rate))),
              .groups = 'drop') %>% 
    arrange(time) %>%
    mutate(across(contains("exploitation_rate_Mean"), as.numeric)) %>%
    mutate(across(contains("exploitation_rate_SE"), as.numeric)) %>% as.data.frame()
  # calculate one mean exploitation rate for each source, each colony and each feeding session
  summary_dat_duration <- aggregate(exploitation_rate   ~colony_id+food_source,FUN=mean,data=cumulative_explotation_over_time)
  summary_dat_duration <- summary_dat_duration %>% left_join(colony_metadata, by = "colony_id")
  
  #stats
  print( paste("Statistics - T0 = ",start_time))
  model <- lmer(exploitation_rate ~ food_source + (1|colony_id) + (1|first_source_fed_on), data=summary_dat_duration )
  print(Anova(model))
  pvalue <- round(Anova(model)$`Pr(>Chisq)`,2)
  print(shapiro.test(residuals(model))) # for shifted t0 not significant... 
  plot_data <- summary_dat_duration %>%
    group_by(food_source) %>%
    summarize(
      mean_exploitation_rate = mean(exploitation_rate, na.rm = TRUE),
      se_exploitation_rate = sd(exploitation_rate, na.rm = TRUE) / sqrt(n()), 
      .groups = 'drop') %>% as.data.frame()
  
  # stats on the mean exploitation rate
  title_lab_3 <- paste0("Mean Exploitation Rate (T0 = ", start_time, ")")
  p_val <- paste0("p = ", pvalue, "(lmer)")
  p <- ggplot(plot_data, aes(x = food_source, y = mean_exploitation_rate, fill = food_source)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_exploitation_rate - se_exploitation_rate, 
                      ymax = mean_exploitation_rate + se_exploitation_rate),
                  width = 0.2) +
    labs(
      x = "Food Source",
      y = "Exploitation Rate (Mean ± SE)",
      title = title_lab_3) +
    scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
  geom_segment(aes(x = 1.1, xend = 1.9 , y = 1.05 * max(mean_exploitation_rate, na.rm = TRUE),
                   yend = 1.05 * max(mean_exploitation_rate, na.rm = TRUE)),
               color = "black", size = 0.25) +
    annotate("text", x = 1.5, y = 1.1 * max(plot_data$mean_exploitation_rate, na.rm = TRUE),
             label = p_val, color = "black", size = 3)
  print(p)
}



























#### 3.3.4 Plotting predicted exploitation of food and exploitation rate ####
# use shifted time to model and predict exploitation. 
# Then add predicted values with confidence intervals to plots 
# cumulative exploitation for shifted starting times per food source is already loaded if the code above has been run, if not you need to get it.
# mod <- lmer(cumulated_exploitation_time ~ food_source*time + (1|colony_id), data = cumulative_explotation_over_time)
# cumulative_explotation_over_time %>% as.data.frame()
# cumulative_explotation_over_time$cumulated_exploitation_time[5]
# 
# mod <- lmer(cumulated_exploitation_time ~ food_source*time + (1|colony_id), data = cumulative_explotation_over_time)
# summary(mod)
# Anova(mod, type=3)
# test_norm(mod)
# 
# data_predic <- cumulative_explotation_over_time
# time_range <- max(data_predic$time)-min(data_predic$time)
# 
# new_dat <- expand.grid(
#   time = seq(min(data_predic$time), max(data_predic$time), length.out = 100),
#   food_source = unique(data_predic$food_source)
# )
# 
# new_dat$cumulated_exploitation_time <- predict(mod, new_dat, re.form = NA)  
# mm <- model.matrix(terms(mod),new_dat)
# pvar1 <- diag(mm %*% tcrossprod(vcov(mod),mm))
# cmult <- 1.96 # use 1.96 for a 95% confidence interval
# new_dat <- data.frame(
#   new_dat,
#   ci_lo = new_dat$cumulated_exploitation_time - cmult * sqrt(pvar1),
#   ci_hi = new_dat$cumulated_exploitation_time + cmult * sqrt(pvar1),
#   se_lo = new_dat$cumulated_exploitation_time - sqrt(pvar1),
#   se_hi = new_dat$cumulated_exploitation_time + sqrt(pvar1)
# )
# new_dat_control     <- new_dat[which(new_dat$food_source == "control"),]
# new_dat_virus       <- new_dat[which(new_dat$food_source == "virus"  ),]
# data_predic$food_source <- factor(data_predic$food_source, levels = c("control", "virus"))
# 
# ### plot predicted values ###
# 
# col_control <- "#CCEBC5"
# col_virus <- "#FBB4AE" 
# 
# p_cum_exploitation <- ggplot(data = new_dat, aes(x = time, y = cumulated_exploitation_time)) +
#   geom_line(data = new_dat_control, aes(y = cumulated_exploitation_time), size = 1.2, color = col_control) +
#   geom_line(data = new_dat_virus, aes(y = cumulated_exploitation_time), size = 1.2, color = col_virus) +
#   geom_ribbon(data = new_dat_control, aes(ymin = ci_lo, ymax = ci_hi), fill = alpha(col_control, 0.5), color = NA) + # ribbons for standard error ranges
#   geom_ribbon(data = new_dat_virus, aes(ymin = ci_lo, ymax = ci_hi), fill = alpha(col_virus, 0.5), color = NA) +
#   # Add points for observed data
#   geom_point(alpha = 0.5, size = 0.1, position = position_jitterdodge(jitter.width = 4.2, dodge.width = 9),
#             aes(group = food_source, color = food_source)) +
#   xlab("Time (sec)") +
#   ylab("Cumulated Exploitation Time (sec)") +
#   scale_color_manual(values = c("control" = col_control, "virus" = col_virus)) +
#   theme_bw() +
#   theme(legend.title = element_blank())+
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# print(p_cum_exploitation)
# add points and means as line to plot!

### the combined plot with predicted values and real data does not work yet, but it can wait.

# mean_exploitation_time <- cumulative_explotation_over_time %>%
#   group_by(time, food_source) %>%
#   summarize(mean_time = mean(cumulated_exploitation_time, na.rm = TRUE),
#             sd = sd(cumulated_exploitation_time, na.rm = TRUE), 
#             sem = sem(cumulated_exploitation_time)) %>% as.data.frame()
# 
# combined_plot <- ggplot(mean_exploitation_time, aes(x = time, y = mean_time, color = food_source)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = mean_time - sem, ymax = mean_time + sem, fill = food_source), alpha = 0.2) +
#   geom_line(data = new_dat_control, aes(x = time, y = cumulated_exploitation_time), size = 1.2, color = col_control) +
#   geom_ribbon(data = new_dat_control, aes(x = time, ymin = ci_lo, ymax = ci_hi), fill = alpha(col_control, 0.2), color = NA) +
#   
#   # Plotting actual data for the virus group
#   geom_line(data = new_dat_virus, aes(x = time, y = cumulated_exploitation_time), size = 1.2, color = col_virus) +
#   geom_ribbon(data = new_dat_virus, aes(x = time, ymin = ci_lo, ymax = ci_hi), fill = alpha(col_virus, 0.2), color = NA) +
# 
#   labs(
#     title = title_lab_1,
#     x = "Time",
#     y = "Mean Cumulated Exploitation Time (seconds)"
#   ) +
#   scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # Print the combined plot
# print(combined_plot)


# !not yet sure if that is right

#### ALternative approach? Check in different environment than on the old laptop one... Update that environment to the one I am using here 
#### ALso do the predictions for the rate and plot both cases over the actual data points!!! 

#  †he below might be an al†ernative way of doing the calculations but it crashes R?! ignore for now...
# # Obtain estimated marginal means
# emm <- emmeans(mod, ~ food_source * time, pbkrtest.limit = 121504)
# # Extract predictions and confidence intervals
# emm_df <- as.data.frame(emm)
# ggplot(cumulative_explotation_over_time, aes(x = time, y = cumulated_exploitation_time, color = food_source)) +
#   geom_point() +  # Plot actual data
#   geom_line(data = emm_df, aes(x = time, y = emmean, color = food_source), size = 1) +  # Predicted lines
#   geom_ribbon(data = emm_df, aes(x = time, ymin = lower.CL, ymax = upper.CL, fill = food_source), alpha = 0.2) +  # Confidence intervals
#   labs(title = "Exploitation Over Time by Food Source",
#        x = "Time",
#        y = "Cumulated Exploitation Time") +
#   theme_minimal() +
#   scale_color_manual(values = c("red", "blue")) +  # Adjust colors if needed
#   scale_fill_manual(values = c("red", "blue"))    # Adjust colors if needed



#### 3.3.5 First feeding ####

first_feeding <- dat_duration_first %>% group_by(colony_id, food_source) %>% 
  summarise(time_first_feeding = as.numeric(min(feeding_start_seconds, na.rm = TRUE))) %>% as.data.frame()

# mean time until first feeding per food source
mod <- lmer(log10(time_first_feeding) ~ food_source + (1|colony_id), data = first_feeding)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow=c(1,1))
test_norm(mod)
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

# on average the virus food is feasted upon earlier than the control food source

food_first_eaten_from <- first_feeding %>% group_by(colony_id) %>% 
  arrange(colony_id, time_first_feeding) %>%
  slice(1) %>% dplyr::select(colony_id, food_first_eaten_from = food_source, time_earliest_feeding = time_first_feeding) %>% as.data.frame()
  
t2 <- table(food_first_eaten_from$food_first_eaten_from)
barplot(t2, main = "", xlab = "food source", ylab = "number of first feedings", ylim = c(0,14), xaxt="n")
axis(1, at=c(0.7, 1.95), labels=c("control", "virus"))
segments(x0 = 0.7, y0 =12.5, x1 = 1.9)
text(x = 1.3, y = 13, label = "p = 0.07 (prop.test)", font = 2, cex = 1.1)

n_virus_eaten_first <- sum(food_first_eaten_from$food_first_eaten_from == "virus")
n_max <- nrow(food_first_eaten_from)
binom.test(n_virus_eaten_first, n_max, p = 0.5)
prop.test(n_virus_eaten_first, n_max, alternative = "two.sided", p = 0.5)

# there might be a trend but it is just not significant 
#--> randomization required to disentangle effect of first discovery from a potential preference







#### 3.3.6 Randomization II #### 
#' use randomization to see if difference between the two food sources is caused by first discovery.

# create same variables as used for first randomization written by Nathalie

dat_duration_first <- dat_duration_first %>%
  group_by(colony_id) %>%
  mutate(
    first_feeding_virus     = min(feeding_start_seconds[food_source == "virus"], na.rm = TRUE),
    first_feeding_control   = min(feeding_start_seconds[food_source == "control"], na.rm = TRUE),
    first_source_fed_on = ifelse(
      first_feeding_virus < first_feeding_control, "virus",
      ifelse(first_feeding_control < first_feeding_virus, "control", NA_character_))) %>% 
  ungroup() %>% as.data.frame()

colony_metadata <- dat_duration_first %>%
  dplyr::select(colony_id, treatment, position_virus, block, 
                first_feeding_virus, first_feeding_control, first_source_fed_on) %>%
  distinct(colony_id, .keep_all = TRUE)


if (RUN_ANALYSIS_AND_SIMULATIONS_DURATION) { # RUN_ANALYSIS_AND_SIMULATIONS_DURATION <- TRUE
  ### Get data calculate and plot cumulative exploitation and exploitation rate
  for (start_time in c("shifted_discovery")) { # start_time <- "shifted_discovery" # Loop over different definitions of †0 an other option would be # "start_experiment",  "discovery_first_food", # start_time = "discovery_first_food"
    time_origin <- start_time
    # get subsetted data
    if (start_time == "start_experiment") {
      subsetted_data <- dat_duration_first # complete data set (note: data has only been annotated until sixty minutes after the second food source has been discoverd or until 2h since start of the experiment have elapsed)
    }
    if (start_time == "shifted_discovery") { 
      summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id, food_source) %>%
        summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
      subsetted_data <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id", "food_source"))
      subsetted_data <- subsetted_data %>% # filter to only include feeding events starting within the 60 min (3600s) window since discovery (first actual feeding) of the food source of interest (different discovery times)
        mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
        mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>% 
        filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)
    }
    if (start_time == "discovery_first_food") {
      summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id) %>%
        summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
      subsetted_data <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id"))
      subsetted_data <- subsetted_data %>% # filter to only include first hour since discovery of first food. 
        mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
        mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>%
        filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)
    }
    
    # expand a full grid with second by second rows
    grid <- expand_grid(
      colony_id = unique(subsetted_data$colony_id), 
      food_source = unique(subsetted_data$food_source),
      time = seq(0, max(subsetted_data$feeding_end_seconds, na.rm = TRUE))
    ) %>% as.data.frame()
    
    
    # calculate cumulative exploitation at each second and exploitation rate. 
    cumulative_explotation_over_time <- NULL
    
    for (colony in unique(grid$colony_id)) { # loop over colony_id and food_source
      for (food in unique(grid$food_source)) {
        current_grid <- grid %>% filter(colony_id == colony, food_source == food) # filter grid and feeding data for the current colony and food source
        current_feedings <- subsetted_data %>% filter(colony_id == colony, food_source == food)
        
        # compute across all rows of current_grid using vectorized operations
        if (nrow(current_feedings) > 0) {
          current_feedings <- current_feedings %>% arrange(feeding_start_seconds) # Sort by start time
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
        
        # append to cumulative_exploitation_over_time
        cumulative_explotation_over_time <- bind_rows(cumulative_explotation_over_time, current_grid) 
      }
    }
    
    
    # aggregate data and calculate means of duration for plotting
    mean_exploitation_time <- cumulative_explotation_over_time %>%
      group_by(time, food_source) %>%
      summarize(mean_time = mean(cumulated_exploitation_time, na.rm = TRUE),
                sd = sd(cumulated_exploitation_time, na.rm = TRUE), 
                sem = sem(cumulated_exploitation_time))
    title_lab_1 <- paste0("Cumulative Food Exploitation (T0 = ", start_time, ")")
    print(
      ggplot(mean_exploitation_time, aes(x = time, y = mean_time, color = food_source)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = mean_time - sem, ymax = mean_time + sem, fill = food_source), alpha = 0.2) +
        labs(
          title = title_lab_1,
          x = "Time",
          y = "Mean Cumulated Exploitation Time (seconds)") +
        scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
        theme_bw() +
        theme(legend.title = element_blank())+
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    )
    
    
    # aggregate data and calculate means of exploitation rate for plotting
    mean_exploitation_rate <- cumulative_explotation_over_time %>%
      group_by(time, food_source) %>%
      summarize(mean_exploitation = mean(exploitation_rate, na.rm = TRUE),
                sd = sd(exploitation_rate, na.rm = TRUE))
    title_lab_2 <- paste0("Exploitation Rate (T0 = ", start_time, ")")
    print(
      ggplot(mean_exploitation_rate, aes(x = time, y = mean_exploitation, color = food_source)) +
        geom_line(size = 1) +
        labs(
          title = title_lab_2,
          x = "Time",
          y = "Mean Cumulated Exploitation Time (seconds)") +
        scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
        theme_bw() +
        theme(legend.title = element_blank())+
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    )
    
    
    ### Stats 
    
    # calculate means and SE
    dynamic_summary_dat_dur <- cumulative_explotation_over_time %>%
      group_by(food_source, time) %>%
      summarize(exploitation_rate_Mean   = mean(exploitation_rate, na.rm = TRUE),
                exploitation_rate_SE = sd(exploitation_rate, na.rm = TRUE) / sqrt(sum(!is.na(exploitation_rate))),
                .groups = 'drop') %>% 
      arrange(time) %>%
      mutate(across(contains("exploitation_rate_Mean"), as.numeric)) %>%
      mutate(across(contains("exploitation_rate_SE"), as.numeric)) %>% as.data.frame()
    
    # calculate one mean exploitation rate for each source, each colony and each feeding session
    summary_dat_duration <- aggregate(exploitation_rate   ~colony_id+food_source,FUN=mean,data=cumulative_explotation_over_time)
    summary_dat_duration <- summary_dat_duration %>% left_join(colony_metadata, by = "colony_id")
    
    # model
    print( paste("Statistics - T0 = ",start_time))
    model <- lmer(exploitation_rate ~ food_source + (1|colony_id) + (1|first_source_fed_on), data=summary_dat_duration )
    print(Anova(model))
    print(shapiro.test(residuals(model))) # (consider looking for alternative models because sometimes close. but this is not the test statistics to report so ignore for now
    
    # plot overall mean
    plot_data <- summary_dat_duration %>%
      group_by(food_source) %>%
      summarize(
        mean_exploitation_rate = mean(exploitation_rate, na.rm = TRUE),
        se_exploitation_rate = sd(exploitation_rate, na.rm = TRUE) / sqrt(n()), 
        .groups = 'drop') %>% as.data.frame()
    
    # stats on the mean exploitation rate
    
    p <- ggplot(plot_data, aes(x = food_source, y = mean_exploitation_rate, fill = food_source)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = mean_exploitation_rate - se_exploitation_rate, 
                        ymax = mean_exploitation_rate + se_exploitation_rate),
                    width = 0.2) +
      labs(
        x = "Food Source",
        y = "Exploitation Rate (Mean ± SE)",
        title = "Mean Exploitation Rate by Food Source"
      ) +
      scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
      theme_bw() +
      theme(legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(p)
    
    ### Plot exploitation over time data
    ymin <- min(dynamic_summary_dat_dur$exploitation_rate_Mean - dynamic_summary_dat_dur$exploitation_rate_SE, na.rm = TRUE)
    ymax <- max(dynamic_summary_dat_dur$exploitation_rate_Mean + dynamic_summary_dat_dur$exploitation_rate_SE, na.rm = TRUE)
    plot(exploitation_rate_Mean ~ time, type = "l", col = "red", lwd = 1,  # plot virus line
         xlab = "time [s]",
         data = dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), ],
         ylim = c(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)),
         bty = "l", 
         main = paste("Exploitation Rate -", time_origin, "t0"))
    polygon( # polygon virus
      x = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "time"],
            rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "time"])),
      y = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_Mean"] - dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_SE"],
            rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_Mean"] + dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_SE"])),
      col = alpha("#FBB4AE", 0.5), border = NA)
    lines( # line control food 
      exploitation_rate_Mean ~ time, col = "green", lwd = 1,  
      data = dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), ])
    polygon( #  polygon - control food
      x = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "time"],
            rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "time"])),
      y = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_Mean"] - dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_SE"],
            rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_Mean"] + dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_SE"])),
      col = alpha("#CCEBC5", 0.5), border = NA)
    legend("topright", legend = c("virus", "control"), col = c("red", "green"), lwd = 2, bty = "n")
    
    ### Randomisation test 
    # Aim: test whether the results are an artifact from virus food tending to be discovered earlier than healthy food
    # H0: there is no difference in ant behaviour towards either virus or healthy food, the only relevant point explaining exploitation rate is which food was discovered first
    
    # In the observed dataset the virus was discovered first a certain number of time
    # The randomisation will consist in randomly drawing the same number of colonies from the total and saying that in these colonies the first source had virus (irrespective of truth), and in all others the first source was healthy irrespective of truth)
    # If observed data is an artifact, we should observe similar exploitation rates in the number of ants between the two sources in randomised and real data
    
    ### Pre-calculations: nb of colonies that found the virus first in each session, and observed difference in number of ants across all colonies and all times
    
    
    nb_discoveries <- aggregate(colony_id~first_source_fed_on, FUN=length, data=colony_metadata)
    # mean difference in exploitation rate
    delta_exploitation_rate <- plot_data$mean_exploitation_rate[plot_data$food_source == "virus"] - plot_data$mean_exploitation_rate[plot_data$food_source == "control"]
    observed_Diff <- data.frame(
      status = "virus_minus_control",
      delta_exploitation_rate = delta_exploitation_rate)
    
    summary_dat_duration <- summary_dat_duration %>% left_join(colony_metadata, by = "colony_id")
    
    #initiate randomisation 
    random_data_dynamic_dur <- NULL
    mean_random_data_dynamic_dur <- NULL
    
    pb <- progress_bar$new(
      format = "Progress: :current/:total [:bar] :percent ETA: :eta",
      total = total_iterations,
      clear = FALSE,
      width = 60
    )
    
    for (i in 1 : total_iterations){ ### randomisation loop # i <- 1
      rand_dur   <- NULL
      nb_virus <- nb_discoveries[which(nb_discoveries$first_source_fed_on == "virus"), "colony_id"]
      subset_dur <- cumulative_explotation_over_time %>% left_join(colony_metadata, by = "colony_id") %>% as.data.frame() # which dataset is this  is this
      colony_list <- unique(subset_dur$colony_id)
      hypothetic_virus_first    <- sample(colony_list,size=nb_virus,replace = F)
      
      # now create first_source_discovered_RAND column according to hypothetical status 
      rand_subset_dur                                                                                            <- subset_dur
      rand_subset_dur$first_source_discovered_RAND                                                               <- "control"  # note: it is refering to first food source fed on and not first discoverd (slight difference in the actual data) 
      rand_subset_dur[which(rand_subset_dur$colony_id%in%hypothetic_virus_first),"first_source_discovered_RAND"] <- "virus" 
      
      # copy food_source column into a new column food_source_RAND, and switch the values for those colonies who have been assigned a different first source discovered than in the observed data
      rand_subset_dur$food_source_RAND  <- rand_subset_dur$food_source    
      rand_subset_dur[which(rand_subset_dur$first_source_discovered_RAND!=rand_subset_dur$first_source_fed_on &  rand_subset_dur$food_source=="virus")  , "food_source_RAND"] <- "control"
      rand_subset_dur[which(rand_subset_dur$first_source_discovered_RAND!=rand_subset_dur$first_source_fed_on &  rand_subset_dur$food_source=="healthy"), "food_source_RAND"]    <- "virus"
      
      # concatenate and store
      rand_dur <- rbind(rand_dur  , data.frame(RAND=i, rand_subset_dur))
      
      # get the mean number of ants at each time point in the hypothetical scenario, across all colonies
      mean_rand_dat_dur <- aggregate(exploitation_rate ~ food_source_RAND + time, FUN=mean, data=rand_dur)
      overall_mean_rand_dat_dur <- aggregate(exploitation_rate ~ food_source_RAND, FUN=mean, data=rand_dur)
      
      # concatenate
      random_data_dynamic_dur           <- rbind(random_data_dynamic_dur , data.frame(RAND=i, mean_rand_dat_dur))
      mean_random_data_dynamic_dur      <- rbind(mean_random_data_dynamic_dur , data.frame(RAND=i, overall_mean_rand_dat_dur))
      pb$tick()
    }
    
    mean_random_data_dynamic_dur_DIFF <- mean_random_data_dynamic_dur %>%
      pivot_wider(names_from = food_source_RAND, values_from = exploitation_rate) %>%
      mutate(delta_exploitation_rate = virus - control) %>% as.data.frame()
    
    
    ### Plot expected vs observed Exploitation rate
    xmin <- min (c(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate),observed_Diff[,"delta_exploitation_rate"])
    xmax <- max (c(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate),observed_Diff[,"delta_exploitation_rate"])
    hist(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate,col=alpha("grey",0.5),xlim=c(xmin,xmax),main=paste("Observed vs Expected,",capitalize(time_origin),"T0"), xlab=expression(paste(Delta, " exploitation rate")))
    hist(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate,col=alpha("grey",0.5),xlim=c(xmin,xmax),main=paste("Observed vs Expected,",capitalize(time_origin),"T0"), xlab=expression(paste(Delta, " exploitation rate")))
    arrows(x0=observed_Diff[,"delta_exploitation_rate"],
           y0=0,
           y1=-10,
           code=1,
           col="black",
           lwd=4,
           length=0.1)
    
    ### P value
    if (observed_Diff[,"delta_exploitation_rate"]>median(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate)){
      pval <- 2*length(which(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate>=observed_Diff[,"delta_exploitation_rate"]))/length(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate)
      mtext(paste("p=",pval),3, cex = 1.3)
    }else{
      pval <- 2*length(which(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate<=observed_Diff[,"delta_exploitation_rate"]))/length(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate)
      mtext(paste("p=",pval),3, cex = 1.3)
    }
  }
}

#' Randomisation test: By randomising first discovery we create a normal distribution of random differences between the exploitation of the two food sources 
#' showing the range where we could expect our data to land if first discovery is the driving factor for the difference in exploitation rates
#' If our data is an outlier of that data i.e. p value is significant we can assume that there is an additional factor for example actual preference
#' which drives the observed data away from the randomised data

#### 3.3.7 Probability of Feeding ####
#' Idea: Combine the nb data set and the exploitation data set for the available time points calculate a probability of feeding:
#' I.e. nr of ants feeding over nr of ants present on food source at the time. 
#' Make a polygon plot and do stats to compare if probability of feeding differs for virus and non virus food?
#' Unfortunately this is not true probability of feeding but rather the proportion of ants present at the food which currently feeding at certain points in time.

# create new time variable (seconds since start of the recording) that can be matched between the two data sets.
number_frame <- dynamic_dat_shifted_t0_Nb %>% 
  filter(feeding_session == 1 ) %>% 
  mutate(time_since_start_rec = ifelse(status == "virus", 
                                       time + discovery_time_virus, 
                                       time + discovery_time_healthy),
         status = ifelse(status == "healthy", "control", status)
  ) %>%
  rename(Nb_ants_present = Nb_ants, 
         food_source = status) %>% 
  dplyr::select(colony_id, food_source, Nb_ants_present, time_since_start_rec)


summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id, food_source) %>%
    summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
feeder_frame <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id", "food_source"))

Nb_ants_feeding <- numeric(nrow(number_frame)) #vector to fill 

for (i in 1:nrow(number_frame)) { # i = 7 # for each time point calculate how many ants are feeding at that moment
  colony_id <- number_frame$colony_id[i]
  food_source <- number_frame$food_source[i]
  time_since_start_rec <- number_frame$time_since_start_rec[i]
  Nb_ants_feeding[i] <- sum(feeder_frame$colony_id == colony_id & 
                              feeder_frame$food_source == food_source &
                              feeder_frame$feeding_start_seconds <= time_since_start_rec &
                              feeder_frame$feeding_end_seconds >= time_since_start_rec)
}

number_frame$Nb_ants_feeding <- Nb_ants_feeding

# number_frame <- number_frame %>%
#   rowwise() %>% 
#   mutate(Nb_ants_feeding = sum(feeder_frame$colony_id == colony_id & 
#                                  feeder_frame$food_source == food_source &
#                                  feeder_frame$feeding_start_seconds <= time_since_start_rec &
#                                  feeder_frame$feeding_end_seconds >= time_since_start_rec)) %>%
#   ungroup() %>% as.data.frame()

number_frame <- number_frame %>%
  mutate(feeding_probability = ifelse(Nb_ants_present > 0, Nb_ants_feeding / Nb_ants_present, 0))
number_frame_filtered <- number_frame %>%
  filter(Nb_ants_present !=   0)


# statistical testing
#' binomial models can also handle proportional data (values between 0 and 1), where the dependent variable is a proportion (e.g., a probability).
#' This is done using the logit link for a binomial distribution.

# Create columns for successes and failures from feeding_probability
number_frame_filtered$successes <- round(number_frame_filtered$feeding_probability * number_frame_filtered$Nb_ants_present) # number of ants present and feeding 
number_frame_filtered$failures <- number_frame_filtered$Nb_ants_present - number_frame_filtered$successes # number of ants present but not feeding
# Fit a GLMM using binomial family (with logit link)
mod_glmm <- glmer(cbind(successes, failures) ~ food_source + (1 | colony_id),
                  data = number_frame_filtered, 
                  family = binomial(link = "logit"))
summary(mod_glmm)


ggplot(number_frame_filtered, aes(x = food_source, y = feeding_probability, fill = food_source)) +
  geom_violin(trim = TRUE) + #or just geom_boxplot()
  geom_jitter(width = 0.3, color = "black", alpha = 0.2) +  # Add jittered points
  labs(title = "Feeding Probability of ants present on Food",
       x = "Food Source",
       y = "Feeding Probability") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +  # Custom colors
  geom_segment(aes(x = 1.1, xend = 1.9, y = 0.9 * max(feeding_probability, na.rm = TRUE), 
                   yend = 0.9 * max(feeding_probability, na.rm = TRUE)),
               color = "black", size = 0.25) +  # Line between violins
  annotate("text", x = 1.5, y = 0.95 * max(number_frame_filtered$feeding_probability, na.rm = TRUE), 
           label = "p = 0.98 (lmer)", color = "black", size = 3)  # p-value annotation






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
cat(yellow("Move on to the next script: bead_data_analysis.R"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

