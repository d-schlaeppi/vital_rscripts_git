### ### ### ### ### ### ### ### ### ###
### ### ###  FC1 Analyses   ### ### ###
### ### ### ### ### ### ### ### ### ###

rm(list=ls())

#### prerequisites ####
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lme4)) #lmer
suppressPackageStartupMessages(library(car)) #contains levene test
suppressPackageStartupMessages(library(lmtest)) #bptest(model) -  heteroscedasticity
suppressPackageStartupMessages(library(blmeco)) #contains compareqqnorm  (multiple qq boxplots)
suppressPackageStartupMessages(library(ggplot2))

setwd("/Users/gismo/Desktop/R/vital") # homeoffice mac
dat <- read.table("food_choice_1_data_2021_DS.txt", header = TRUE)
dat$line_nr <- c(1:32) # add a row number
sem <- function(x){sd(x,na.rm=T)/sqrt(length(na.omit(x)))} #defining the function for the standard error of means
names(dat)
head(dat)

# variables and explanation
# experimental_id - colony identifier (not used)
# colony_id - second colony identifier (used)
# position_virus - where the virus was positioned seen from the nest entrance
# treatment - based on virus position in the two feeding events: right right, left left, right left, left right
# treatment_nr - just a numeric version of the treatment (needs to be translated to factor if used)
# video_position - Which of the four video positions (4 recording stations in a single incubator) was used for the colony recording
# feeding session - First or second feeding event for each colony
# block - Always for colonies were recorded at the same time as a block. Each block corresponds to one (of eight) recording session with four colonies (please note: slightly special definitin as each colony is present in two blocks (e.g.1&5) 
#         which is different from the four blocks as shown in the lab meeting where the two feeding events of each colony are in the same block (1-4)
# camera - which of the four cameras was used to record the colony (almost corresponds to video-position)
# position - in the first two recording sessions only two cameras were used and the l-left & r-right stand for the positoin of the colony in the videofootage
# time_first_ant - time since start of the recording until the first ant arrives at one of the two food sources (defined as t00)
# side_of_first_ant - which food source (left or right) was discovered first
# second_side - which food source (left or right) was discovered second
# t00; t05; t10....t60 - the times at which the number of ants was recorded in five min intervals (starting from t00)
# t00_l; t00_r; t05_l; t05_r; ... t60_r - number of ants recorded on the left and right food source at t00-t60
# time_first_ant_left - time the first ant is on the left food source
# time_first_ant_right - time the first ant is on the right food source
# second - which food source (left or right) was discovered second = second_side
# trl00; trl15; trl20...trl60 - Time shifted 5 min recording intervalls strarting from the time shifted t00 which is when the second food source was discovered
# rl00, rl15, rl20...rl60 - number of ants recorded on the food source discovered second at the time shifted t00-t60


#### To Dos ####

# Analyse first discovery as binomial data by translating control/virus into 1 or 0 and then fit a glmer wiht binomial data structure with colony as random factor
# Fit a logistic curve over each feeding event... and then compare curves
# Individual plot for each feeding event 
# make the time series graph with first and second discovered!
# randomization tests? explore what this is
# Ev. explore cumulative sum instead of mean? 
# Analyse the two feeding events separately



#### First Data Exploration ####
head(dat)

# create data frame (tibble) containing means 
df <-  group_by(dat, line_nr)  %>%
  summarise(
    count = n(),
    mean_second = sum(rl00, rl05, rl10, rl15, rl20, rl25, rl30, rl35, rl40, rl45, rl50, rl55, rl60)/13,
    mean_first = if(second == "l") {
      sum(t00_r, t05_r, t10_r, t15_r, t20_r, t25_r, t30_r, t35_r, t40_r, t45_r, t50_r, t55_r, t60_r)/13
    } else {
      sum(t00_l, t05_l, t10_l, t15_l, t20_l, t25_l, t30_l, t35_l, t40_l, t45_l, t50_l, t55_l, t60_l)/13},
    first_side = if(second == "l") {"right"} else {"left"},
    right_mean = if(first_side == "right") {mean_first} else {mean_second},
    left_mean = if(first_side == "left") {mean_first} else {mean_second},
    virus_position = position_virus, #change back to position corrected if required! Check out the video file!!! Block 1 cam 1 position left (1)
    virus_mean = if(virus_position == "right") {right_mean} else {left_mean},
    control_mean = if(virus_position == "right") {left_mean} else {right_mean},
    diff = control_mean - virus_mean, 
    first_discovered = if(virus_position == "right" & first_side == "left" | virus_position == "left" & first_side == "right") {"control_first"} else {"virus_first"}
  )

df %>% print(n = Inf)
df %>% head()


# How the different feeding events looked like
ggplot(df) + 
  geom_segment(aes(x = 1, xend = 2, y = control_mean, yend = virus_mean, colour = first_discovered)) + 
  theme_classic() + 
  geom_point(aes(x = 1,  y = control_mean)) +
  geom_point(aes(x = 2, y = virus_mean)) +
  scale_x_discrete(
    breaks = c("1", "2"),
    labels = c("control", "virus"),
    limits = c(1, 2)
  ) + 
  labs(y = "average number of ants per foodsource", x = "food source") + 
  scale_color_manual(values=c("darkblue", "green")) +
  ggtitle("Mean (t00-t60) number of ants for each pair of food sources")

#### time series plots ####

# normal 
# create data frame that contains the number of ants present at the control and at the virus food source (based on left or right)
data_plot <-  group_by(dat, line_nr)  %>%
  summarise(
    c_00 =  if(pos_corrected == "right") {t00_l} else {t00_r}, 
    c_05 =  if(pos_corrected == "right") {t05_l} else {t05_r},
    c_10 =  if(pos_corrected == "right") {t10_l} else {t10_r},
    c_15 =  if(pos_corrected == "right") {t15_l} else {t15_r},
    c_20 =  if(pos_corrected == "right") {t20_l} else {t20_r},
    c_25 =  if(pos_corrected == "right") {t25_l} else {t25_r},
    c_30 =  if(pos_corrected == "right") {t30_l} else {t30_r},
    c_35 =  if(pos_corrected == "right") {t35_l} else {t35_r},
    c_40 =  if(pos_corrected == "right") {t40_l} else {t40_r},
    c_45 =  if(pos_corrected == "right") {t45_l} else {t45_r},
    c_50 =  if(pos_corrected == "right") {t50_l} else {t50_r},
    c_55 =  if(pos_corrected == "right") {t55_l} else {t55_r},
    c_60 =  if(pos_corrected == "right") {t60_l} else {t60_r},
    v_00 =  if(pos_corrected == "left") {t00_l} else {t00_r}, 
    v_05 =  if(pos_corrected == "left") {t05_l} else {t05_r},
    v_10 =  if(pos_corrected == "left") {t10_l} else {t10_r},
    v_15 =  if(pos_corrected == "left") {t15_l} else {t15_r},
    v_20 =  if(pos_corrected == "left") {t20_l} else {t20_r},
    v_25 =  if(pos_corrected == "left") {t25_l} else {t25_r},
    v_30 =  if(pos_corrected == "left") {t30_l} else {t30_r},
    v_35 =  if(pos_corrected == "left") {t35_l} else {t35_r},
    v_40 =  if(pos_corrected == "left") {t40_l} else {t40_r},
    v_45 =  if(pos_corrected == "left") {t45_l} else {t45_r},
    v_50 =  if(pos_corrected == "left") {t50_l} else {t50_r},
    v_55 =  if(pos_corrected == "left") {t55_l} else {t55_r},
    v_60 =  if(pos_corrected == "left") {t60_l} else {t60_r},
  )
number_of_ants <- as.vector(as.matrix(data_plot[,c(2:27)]))
plot_frame <- data.frame(number_of_ants)
times <- c(0,5,10,15,20,25,30,35,40,45,50,55,60)
plot_frame$food_bowl <- c(rep(1:32, times = 13), rep(33:64, times = 13))
plot_frame$times <- rep(times, each = 32, times = 2)
plot_frame$treatment <- rep(c("control","virus"), each = 16*26)
plot_frame$id <- rep(dat$colony_id, times = 26)
plot_frame$feeding_session <- rep(c("feeding session 1","feeding session 2"), each = 16, times = 26)



#plot time series with standard error of means
x <- plot_frame %>% group_by(treatment, times) %>%  
  summarise(
    mean = mean(number_of_ants), 
    sd = sd(number_of_ants), 
    sem = sem(number_of_ants)
  )
x <- as.data.frame(x)

p <- ggplot(data = x, aes(x = times, y = mean, color = treatment)) + 
  geom_point(position=position_dodge(1)) + 
  geom_smooth() +
  geom_line(position=position_dodge(1)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(1)) +
  ggtitle("Mean number of ants over time") +
  xlab("time (t00-t60)") +
  ylab("# ants (mean ± sem)")
print(p)

#plot time series with standard deviation

p <- ggplot(data = x, aes(x = times, y = mean, color = treatment)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(1))
print(p)




#plot time series for each colony
par(mfrow=c(4,8))
for (i in c(1:32)) {
  plot(plot_frame$number_of_ants[plot_frame$food_bowl==i] ~ plot_frame$times[plot_frame$food_bowl==i], 
       ylim = c(0, max(plot_frame$number_of_ants[plot_frame$food_bowl==i], plot_frame$number_of_ants[plot_frame$food_bowl==i+32])),
       col = "darkblue", xlab = "time", ylab = "# ants per food source", 
       main = paste(plot_frame$id[i], c("--") ,plot_frame$feeding_session[i])) 
  lines(plot_frame$number_of_ants[plot_frame$food_bowl==i] ~ plot_frame$times[plot_frame$food_bowl==i], col = "darkblue")
  points(plot_frame$number_of_ants[plot_frame$food_bowl==i+32] ~ plot_frame$times[plot_frame$food_bowl==i+32], col = "green") 
  lines(plot_frame$number_of_ants[plot_frame$food_bowl==i+32] ~ plot_frame$times[plot_frame$food_bowl==i+32], col = "green")
  }
par(mfrow=c(1,1))





#### time shifted time series plots ####
# with the second values for starting from the discovery of the second food source
# to do adjust graph of time shifted data so that the scale on the y axis fits the one of the non-timeshifted series

data_plot2 <-  group_by(dat, line_nr)  %>%
  summarise(
    c_00 = if(pos_corrected == "right") {if(second == "r") {t00_l} else {rl00} } else {if(second == "r") {rl00} else {t00_r}},
    c_05 = if(pos_corrected == "right") {if(second == "r") {t05_l} else {rl05} } else {if(second == "r") {rl05} else {t05_r}},
    c_10 = if(pos_corrected == "right") {if(second == "r") {t10_l} else {rl10} } else {if(second == "r") {rl10} else {t10_r}},
    c_15 = if(pos_corrected == "right") {if(second == "r") {t15_l} else {rl15} } else {if(second == "r") {rl15} else {t15_r}},
    c_20 = if(pos_corrected == "right") {if(second == "r") {t20_l} else {rl20} } else {if(second == "r") {rl20} else {t20_r}},
    c_25 = if(pos_corrected == "right") {if(second == "r") {t25_l} else {rl25} } else {if(second == "r") {rl25} else {t25_r}},
    c_30 = if(pos_corrected == "right") {if(second == "r") {t30_l} else {rl30} } else {if(second == "r") {rl30} else {t30_r}},
    c_35 = if(pos_corrected == "right") {if(second == "r") {t35_l} else {rl35} } else {if(second == "r") {rl35} else {t35_r}},
    c_40 = if(pos_corrected == "right") {if(second == "r") {t40_l} else {rl40} } else {if(second == "r") {rl40} else {t40_r}},
    c_45 = if(pos_corrected == "right") {if(second == "r") {t45_l} else {rl45} } else {if(second == "r") {rl45} else {t45_r}},
    c_50 = if(pos_corrected == "right") {if(second == "r") {t50_l} else {rl50} } else {if(second == "r") {rl50} else {t50_r}},
    c_55 = if(pos_corrected == "right") {if(second == "r") {t55_l} else {rl55} } else {if(second == "r") {rl55} else {t55_r}},
    c_60 = if(pos_corrected == "right") {if(second == "r") {t60_l} else {rl60} } else {if(second == "r") {rl60} else {t60_r}},
    v_00 = if(pos_corrected == "right") {if(second == "r") {rl00} else {t00_r} } else {if(second == "r") {t00_l} else {rl00}},
    v_05 = if(pos_corrected == "right") {if(second == "r") {rl05} else {t05_r} } else {if(second == "r") {t05_l} else {rl05}},
    v_10 = if(pos_corrected == "right") {if(second == "r") {rl10} else {t10_r} } else {if(second == "r") {t10_l} else {rl10}},
    v_15 = if(pos_corrected == "right") {if(second == "r") {rl15} else {t15_r} } else {if(second == "r") {t15_l} else {rl15}},
    v_20 = if(pos_corrected == "right") {if(second == "r") {rl20} else {t20_r} } else {if(second == "r") {t20_l} else {rl20}},
    v_25 = if(pos_corrected == "right") {if(second == "r") {rl25} else {t25_r} } else {if(second == "r") {t25_l} else {rl25}},
    v_30 = if(pos_corrected == "right") {if(second == "r") {rl30} else {t30_r} } else {if(second == "r") {t30_l} else {rl30}},
    v_35 = if(pos_corrected == "right") {if(second == "r") {rl35} else {t35_r} } else {if(second == "r") {t35_l} else {rl35}},
    v_40 = if(pos_corrected == "right") {if(second == "r") {rl40} else {t40_r} } else {if(second == "r") {t40_l} else {rl40}},
    v_45 = if(pos_corrected == "right") {if(second == "r") {rl45} else {t45_r} } else {if(second == "r") {t45_l} else {rl45}},
    v_50 = if(pos_corrected == "right") {if(second == "r") {rl50} else {t50_r} } else {if(second == "r") {t50_l} else {rl50}},
    v_55 = if(pos_corrected == "right") {if(second == "r") {rl55} else {t55_r} } else {if(second == "r") {t55_l} else {rl55}},
    v_60 = if(pos_corrected == "right") {if(second == "r") {rl60} else {t60_r} } else {if(second == "r") {t60_l} else {rl60}},
  )


number_of_ants2 <- as.vector(as.matrix(data_plot2[,c(2:27)]))
plot_frame2 <- data.frame(number_of_ants2)
plot_frame2$food_bowl <- c(rep(1:32, times = 13), rep(33:64, times = 13))
plot_frame2$times <- rep(times, each = 32, times = 2)
plot_frame2$treatment <- rep(c("control","virus"), each = 16*26)

x2 <- plot_frame2 %>% group_by(treatment, times) %>%  
  summarise(
    mean = mean(number_of_ants2), 
    sd = sd(number_of_ants2), 
    sem = sem(number_of_ants2)
  )
x2 <- as.data.frame(x2)

p <- ggplot(data = x2, aes(x = times, y = mean, color = treatment)) + 
  geom_point(position=position_dodge(1)) + 
  geom_line(position=position_dodge(1)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(1))+ 
  ggtitle("Mean number of ants over time (timeshifted") +
  xlab("time (t00-t60)") +
  ylab("# ants (mean ± sem)")
print(p)


#### Further preliminary data explorations ####
# general look at left and right (independant of food source)
boxplot(df$left_mean, df$right_mean, main = "left vs right mean", xlab = "foodsource", ylab = "# ants (mean)", names = c("left", "right"))

# potential difference between virus and control? 

max(df$virus_mean)

boxplot(df$control_mean, df$virus_mean, main = "control vs virus (mean)", xlab = "foodsource", ylab = "# ants (mean)", names = c("control", "virus")) 
t.test(df$virus_mean, df$control_mean, paired = TRUE) # simplest test indicates that the two means are different - however not independant samples so another test will be required
boxplot(df$diff, main = "mean difference between control and virus", xlab = "Δ mean", ylab = "# ants (mean)", names = c())
abline(0, 0, lty = 2, lwd = 2)


par(mfrow=c(1,2))
boxplot(df$control_mean, df$virus_mean, main = "mean number of ants per food", xlab = "food source", ylab = "#ants", pars  =  list(xaxt = "n"))
axis(1, at=c(1,2), labels = FALSE)
text(c(1,2), -0.2, labels = c("control", "virus"), srt = 45, pos = 1, xpd = TRUE)
boxplot(df$diff, main = "Mean Δants of food pairs(control-virus)", xlab = "mean difference", ylab = "Δ #ants")
abline(0, 0, lty = 2, lwd = 2)
par(mfrow=c(1,1))


# Exploring the bias towards the first side discovered
boxplot(df$diff ~ df$first_side)

boxplot((df$virus_mean+df$control_mean) ~ df$first_discovered)

boxplot(df$virus_mean ~ df$first_discovered) # If the virus was discovered first the overall response of the ants seems to be stronger than when the control food source is discovered first?
boxplot(df$control_mean ~ df$first_discovered)

boxplot(df$diff ~ df$first_discovered) #if virus first discovered a  bias towards virus, if control first then a slight trend towards control but less strong?
abline(0,0)

# test for effect of first discovered!
t <- table(df$first_discovered)
barplot(t, main = "Frequency first discovered")
#Is that significantly different from 50/50?
prop.test(11, 32, alternative = "two.sided")
binom.test(11, 32) 
# to do: find correct way to do this with glmer (binomial response 1 or 0)



#### Model - Mean nr of Ants over time (60 min) on virus and control food ####

#create suitable data frame containting the right resonse variable (means) and the required random factors etz
treatment_model <- c(rep("control", 32), rep("virus", 32))
mean_nr_ants <- c(df$control_mean, df$virus_mean)
colony_id <- rep(dat$colony_id, 2)
position_virus <- rep(dat$position_virus, 2)
treatment_expdesign <- rep(dat$treatment, 2)
block <- as.factor(rep(dat$block, 2))
feeding_session <- as.factor(rep(dat$feeding, 2))
camera_position <- as.factor(rep(dat$camera, 2))
first_discovered <- rep(df$first_discovered, 2)
unit <- rep(1:32, 2) #do I need this somehow? Or is it accounted for if I say that colony_ID is a random factor

df_model <- data.frame(treatment_model, mean_nr_ants, colony_id, position_virus, treatment_expdesign, block, feeding_session, camera_position, first_discovered, unit)
str(df_model)
head(df_model)

boxplot(df_model$mean_nr_ants ~ df_model$treatment_model*df_model$first_discovered, main = "mean #ants ~ foodsource and first discovery", ylab = "#ants", xlab = "food source and discovery", pars  =  list(xaxt = "n"))
text(c(1:4), -0.2, labels = c("control, control first","control, virus first", "virus, control first","virus, virus first"), srt = 20, pos = 1, xpd = TRUE)
# Bias towards first discovered food source stronger in virus?

# simple linear model (mod1) - quick first look before including potential random factors
boxplot(df_model$mean_nr_ants ~ df_model$treatment_model)
mod1 <- lm(mean_nr_ants ~ treatment_model, data = df_model)
summary(mod1)

# of course model 1 is not ok because we need to take into account pseudo replication and check the random factors
#model: mean nr of ants control side vs virus side, random effects - colony id, treatment, position, block 1-4, run 1/2, position of virus

#potential random effects
boxplot(df_model$mean_nr_ants ~ df_model$treatment_expdesign) 
summary(lm(df_model$mean_nr_ants ~ df_model$treatment_expdesign)) 
boxplot(df_model$mean_nr_ants ~ df_model$position_virus) 
summary(lm(df_model$mean_nr_ants ~ df_model$position_virus)) 
boxplot(df_model$mean_nr_ants ~ df_model$block)
abline(mean(df_model$mean_nr_ants), 0)
summary(lm(df_model$mean_nr_ants ~ df_model$block)) 
boxplot(df_model$mean_nr_ants ~ df_model$feeding_session)
abline(mean(df_model$mean_nr_ants), 0)
summary(lm(df_model$mean_nr_ants ~ df_model$feeding_session)) 
boxplot(df_model$mean_nr_ants ~ df_model$camera_position) 
summary(lm(df_model$mean_nr_ants ~ df_model$camera_position))
boxplot(df_model$mean_nr_ants ~ df_model$first_discovered)
summary(lm(df_model$mean_nr_ants ~ df_model$first_discovered))

#model containing all random factors
mod2 <- lmer(mean_nr_ants ~ treatment_model + (1|first_discovered) + (1|colony_id) + (1|position_virus) + (1|treatment_expdesign) + (1|block) + (1|feeding_session) + (1|camera_position) + (1|unit), data = df_model)
summary(mod2)
Anova(mod2)
estimates <- fixef(mod2)
estimates
#looks nice. However - warning of singular fit because random effects structure is too complex to be suported by the data -> remove random factors that are not requiredd
# of the random factors feeding session explains most of the variance.

# check individual parts of the model to decide which random factors to discard
# baseline <- lmer(mean_nr_ants ~ 1 + (1|colony_id), data = df_model)
# m1 <- update(baseline, .~. + treatment_model)
# m2 <- update(m1 , .~. + (1|treatment_expdesign))
# m3 <- update(m1 , .~. + (1|block))
# m4 <- update(m1 , .~. + (1|feeding_session))
# m5 <- update(m1 , .~. + (1|camera_position))
# m6 <- update(m1 , .~. + (1|position_virus))
# anova(baseline, m1, m2, m3, m4, m5, m6)

baseline <- lmer(mean_nr_ants ~ 1 + (1|colony_id), data = df_model)
m1 <- update(baseline, .~. + treatment_model)
m2 <- update(m1 , .~. + feeding_session)
m3 <- update(m1 , .~. + block)
m4 <- update(m1 , .~. + treatment_expdesign)
m5 <- update(m1 , .~. + camera_position)
m6 <- update(m1 , .~. + position_virus)
m7 <- update(m1 , .~. + first_discovered)
m8 <- update(m1 , .~. + unit)
anova(baseline, m1, m2, m3, m4, m5, m6, m7, m8)
# decide to keep colony id and feeding session for the final model 



#### final model ####
mod3 <- lmer(mean_nr_ants ~ treatment_model + (1|colony_id) + (1|first_discovered), data = df_model) # + (1|unit) would be great to include unit as well, but then the model will be overfitted?
summary(mod3)
Anova(mod3)
estimates <- fixef(mod3)
estimates

# test model assumtions
compareqqnorm(mod3)
par(mfrow=c(1,1))
leveneTest(residuals(mod3) ~ df_model$treatment_model) # non-significant --> homogenity of variance ok
boxplot(residuals(mod3) ~ df_model$treatment_model)
aov_residuals <- residuals(object = mod3)
shapiro.test(x = aov_residuals) # non significant --> ok ##plot(mod, 2)#normality
par(mfrow=c(2,2))
plot(mod3)
scatter.smooth(fitted(mod3),resid(mod3)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod3), main="normal QQ-plot, residuals") 
qqline(resid(mod3))  # qq of residuals
scatter.smooth(fitted(mod3), sqrt(abs(resid(mod3))))  # homogeneity of variance
par(mfrow=c(1,1))
# qq of random effects
qq_re <- ranef(mod3)$colony_id[,1] 
qqnorm(qq_re)
qqline(qq_re)
boxplot(resid(mod3)~df_model$colony_id)
qq_re <- ranef(mod3)$first_discovered[,1]
qqnorm(qq_re)
qqline(qq_re)
boxplot(resid(mod3)~df_model$first_discovered)




mod4 <- lmer(mean_nr_ants ~ treatment_model*feeding_session + (1|colony_id) + (1|first_discovered), data = df_model)
Anova(mod4)

compareqqnorm(mod4)
par(mfrow=c(1,1))
leveneTest(residuals(mod3) ~ df_model$treatment_model) # non-significant --> homogenity of variance ok
shapiro.test(x = aov_residuals) # non significant --> ok ##plot(mod, 2)#normality
par(mfrow=c(2,2))
scatter.smooth(fitted(mod3),resid(mod3)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod3), main="normal QQ-plot, residuals") 
qqline(resid(mod3))  # qq of residuals
scatter.smooth(fitted(mod3), sqrt(abs(resid(mod3))))  # homogeneity of variance
par(mfrow=c(1,1))
