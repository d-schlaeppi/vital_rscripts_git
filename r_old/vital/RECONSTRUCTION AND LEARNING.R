### ### ### ### ### ### ###
### ANALYSIS FOR IUSSI  ###
### ### ### ### ### ### ###

#### prerequisites ####



#load packages
install.packages("pacman")
pacman::p_load(lubridate, plotrix, scales, car, lme4, Hmisc, dplyr, tidyverse, blmeco)

suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(lsmeans))
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(multcompView))
suppressPackageStartupMessages(library(multcomp))
suppressPackageStartupMessages(library(viridis))

#get raw data
setwd("/Users/gismo/Desktop/R/vital") # Daniel

sem <- function(x){sd(x,na.rm=T)/sqrt(length(na.omit(x)))} #defining the function for the standard error of means



#### preparing datasets ####

rm(list=ls())
dat <- read.table("Food_choice_1_data_2021_DS.txt", header = TRUE,stringsAsFactors = F)



# in the loaded dat block corresponds to the eight days of video recordings -> needs to be renamed to video_session
# further, a corrected variable "block" needs to be created" to group the 4x4 colonies that were recorded together
# correct the variable block that goes from 1-8 to 1-4 as it should be
dat <- rename(dat, video_session = block)
block_corrected <- NULL
for (i in 1:nrow(dat)) {
  block_old           <- dat[i, "video_session"]
  block_new       <- if(block_old<5) {block_old} else {block_old-4}
  block_corrected <- rbind(block_corrected, data.frame(block_new))
}
dat$block <- block_corrected$block_new


# create empty objects (to become our data frames which we will use)
dynamic_dat_same_t0_Nb      <- NULL
dynamic_dat_shifted_t0_Nb   <- NULL
dynamic_dat_same_t0_Diff      <- NULL
dynamic_dat_shifted_t0_Diff   <- NULL

# create time point vector (up to 13 time points from t00, t05...t60)
Time_points_to_include <- 1:13 #used for full analysis
Time_points_to_include <- 1:7
Time_points_to_include <- 1:6 
Time_points_to_include <- 1:5 # used for the initial recruitment
#


#### my practice stuff - ignore ####
# fill up objects
# reminder for myself: dat[row,column]
for (i in 1:nrow(dat)) {
  colony_id                             <- dat[i, "colony_id"]
  dynamic_dat_same_t0_Nb                <- rbind(dynamic_dat_same_t0_Nb,      data.frame (colony_id))
  dynamic_dat_shifted_t0_Nb             <- rbind(dynamic_dat_shifted_t0_Nb,   data.frame (colony_id))
  dynamic_dat_same_t0_Diff              <- rbind(dynamic_dat_same_t0_Diff,    data.frame (colony_id))
  dynamic_dat_shifted_t0_Diff           <- rbind(dynamic_dat_shifted_t0_Diff, data.frame (colony_id))
}
dynamic_dat_same_t0_Nb
dynamic_dat_shifted_t0_Nb
dynamic_dat_same_t0_Diff
dynamic_dat_shifted_t0_Diff

dynamic_dat_same_t0_Nb      <- NULL
dynamic_dat_shifted_t0_Nb   <- NULL
dynamic_dat_same_t0_Diff      <- NULL
dynamic_dat_shifted_t0_Diff   <- NULL


#### filling the new data sets ####
for (i in 1:nrow(dat)){
  colony_id                          <- dat[i,"colony_id"]
  feeding_session                    <- dat[i,"feeding_session"]
  block                              <- dat[i,"block"]
  virus_position                     <- dat[i,"pos_corrected"]
  healthy_position                   <- c("left","right")[which(c("left","right")!=virus_position)]
  discovery_time_virus               <- period_to_seconds(hms(dat[i,paste("time_first_ant_",virus_position,sep="")]))
  discovery_time_healthy             <- period_to_seconds(hms(dat[i,paste("time_first_ant_",healthy_position,sep="")]))
  
  Nb_ant_virus_same_t0               <- as.numeric(dat[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(virus_position,1,1)  ,sep="_")])[Time_points_to_include]
  Nb_ant_healthy_same_t0             <- as.numeric(dat[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(healthy_position,1,1),sep="_")])[Time_points_to_include]
  
  if (discovery_time_virus<discovery_time_healthy){
    first_source_discovered          <- "virus"
    Nb_ant_virus_shifted_t0          <- Nb_ant_virus_same_t0
    Nb_ant_healthy_shifted_t0         <- as.numeric(dat[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
  }else{
    first_source_discovered          <- "healthy"
    Nb_ant_virus_shifted_t0          <- as.numeric(dat[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
    Nb_ant_healthy_shifted_t0        <- Nb_ant_healthy_same_t0
  }
  time_since_discovery               <- period_to_seconds(hms(dat[i,c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60")]))[Time_points_to_include]                           - period_to_seconds(hms(dat[i,"time_first_ant"]))
  
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

#we created the dataframes -> lets explore, learn, and do some prelim analysis
data <- dynamic_dat_same_t0_Nb
data$time_min <- data$time/60
#

#### comparison means at both food sources (feeding events separately) ####

head(data)
df1 <- group_by(data, colony_id, feeding_session, status) %>% 
  summarise(
    count = n(), 
    mean = mean(Nb_ants), 
    first_source_discovered = first(first_source_discovered)
  )


#plotting the two feeding sessions pooled (two means per colony)
boxplot(df1$mean ~ df1$status)
#boxplot feeding session separately
boxplot(df1$mean ~ df1$status + df1$feeding_session)

#### boxplots comparison 4 groups for presentation (based on 13 time points) ####

par(mar=c(5.1, 4.5, 4.1, 1.8),
    cex.lab=1.3, 
    cex.axis=1.3, 
    cex.main= 1.3)

boxplot(df1$mean ~ df1$status + df1$feeding_session, main = "mean number of ants over time", ylab = "#ants", xlab = "food source", 
        pars  =  list(xaxt = "n"), ylim = c(0,5.6), cex.lab = 1.3, cex.axis = 1.3, cex.main =1.3)
text(c(1:4), -0.3, labels = c("control", "virus", "control","virus"), pos = 1, xpd = TRUE, cex = 1.3)
text(c(1.5,3.5), 5.4, labels = c("feeding session 1", "feeding session 2"), cex = 1.3)
text(c(1:4), 5, labels = c("a", "bc", "ab", "c"), font = 2, cex = 1.3)
abline(v=2.5, lty=4)

#### model - overall mean nr of ants per food source & per feeding session ####

#model lmer --> used for IUSSI presentation even though the model assumptions might be violated with the smaller data set when not all the time points are included
mod <- lmer(mean ~ status*feeding_session + (1|colony_id) + (1|first_source_discovered), data=df1 )
Anova(mod)
shapiro.test(residuals(mod)) # with 13 time points just not significant 
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
          adjust="tukey")
CLD

# ev as a glmer if we compare smaller data sets?
mod <- glmer(mean ~ status*feeding_session + (1|colony_id) + (1|first_source_discovered), data=df1 , family = "poisson")
Anova(mod)
summary(mod)
#pairwise differences
marginal = lsmeans(mod, ~ status*feeding_session, data = df1)
CLD = cld(marginal,
          alpha=0.5,
          Letters=letters,
          adjust="tukey")
CLD





#### model glmer with time all time points included instead of means ####
head(dynamic_dat_same_t0_Nb)
head(data)
data$time <- as.factor(data$time)
data$feeding_session <- as.factor(data$feeding_session)

model <- glmer(Nb_ants ~ status + time + (1|colony_id) + (1|first_source_discovered) + (1|feeding_session), data=data, family = "poisson")
Anova(model)
estimates <- fixef(model)
estimates






#### plot time series with standard error of means or with geom smooth for a nice outline ####

data_plot <- data %>% group_by(status, time_min) %>% 
  summarise(
    mean = mean(Nb_ants),
    sem = sem(Nb_ants)
  )
data_plot <- as.data.frame(data_plot)

p <- ggplot(data = data_plot, aes(x = time_min, y = mean, color = status)) +
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
print(p)

#test with theme classic

### same graph but with different colors and with polygons instead of geom smooth



p <- ggplot(data = data_plot, aes(x = time_min, y = mean, color = status)) +
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
print(p)





#### graph on first discovery bias ####
# How the different feeding events looked like
data_plot2 <- data %>% group_by(colony_id, feeding_session, status) %>% 
  summarise(
    mean = mean(Nb_ants),
    first_discovered = first(first_source_discovered),
    discovery_time_control = first(discovery_time_healthy), 
    discovery_time_virus = first(discovery_time_virus)
  )
data_plot2 <- as.data.frame(data_plot2)
data_plot2_wide <-  reshape(data = data_plot2, 
                            idvar=c("colony_id", "feeding_session", "first_discovered"),
                            v.names = "mean",
                            timevar="status",
                            direction = "wide"
                            )
head(data_plot2)


p <- ggplot(data_plot2_wide) + 
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
print(p)



data_plot4 <- data %>% group_by(colony_id, feeding_session) %>% 
  summarise(
    mean = mean(Nb_ants),
    first_discovered = first(first_source_discovered),
    discovery_time_control = first(discovery_time_healthy), 
    discovery_time_virus = first(discovery_time_virus)
  )
data_plot4 <- as.data.frame(data_plot4)
head(data_plot4)
str(data_plot4)

boxplot(discovery_time_control, discovery_time_virus, data = data_plot4)
boxplot(data_plot4$discovery_time_control, data_plot4$discovery_time_virus)
# --> probably no difference at all which means we stick to the story with first discoveries for now. 


#### simple histogramm of first discovery ####
head(data_plot2_wide)
t <- table(data_plot2_wide$first_discovered)
barplot(t, main = "", xlab = "food source", ylab = "number of first discoveries", ylim = c(0,22), xaxt="n")
axis(1, at=c(0.7, 1.95), labels=c("control", "virus"))
segments(x0 = 0.7, y0 =21, x1 = 1.9)
text(x = 1.3, y = 21.5, label = "ns", font = 2, cex = 1.3)

pars  =  list(xaxt = "n")

#Is that significantly different from 50/50?
prop.test(12, 32, alternative = "two.sided")
binom.test(11, 32) 



#### time series plot for the first 20 min with shifted t0 to make ####
data_shifted <- dynamic_dat_shifted_t0_Nb
data_shifted$time_min <- data_shifted$time/60


data_plot3 <- data_shifted %>% group_by(status, time_min) %>% 
  summarise(
    mean = mean(Nb_ants), 
    sem = sem(Nb_ants)
  )
data_plot3 <- as.data.frame(data_plot3)
p <- ggplot(data = data_plot3, aes(x = time_min, y = mean, color = status)) +
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
  ylab("mean number of ants") +
  theme(text = element_text(size = 20))
print(p)

#not used:
#  geom_line() + 
#position dodge for point ant line
#geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0,
#              position=position_dodge(1)) +










#### randomisation test and last histogramm ####
#was done in nathalies scirpt