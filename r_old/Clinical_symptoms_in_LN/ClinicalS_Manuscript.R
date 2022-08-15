### ### ### ### ### ### ### ### ### ### ### ### ### 
### Clinical symthoms of ABPV in lasius niger   ###
### ### ### ### ### ### ### ### ### ### ### ### ### 

# set working directory and load data
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/Clinical_symptoms_in_LN")     
col <- read.table("table1.txt", header = TRUE)  
in_col <- read.table("table2.txt", header = TRUE)
mov <- read.table("speed.txt", header = TRUE)
ANV <- read.table("ANVlong_full.txt", header = TRUE)
virus_mov <- read.table("20190226_ANV.txt", header = TRUE)
mortality <- read.table("mortality.txt", header = TRUE)


# load packages
install.packages(c("lmerTest", "chron", "lme4", "arm", "car", "blmeco"))
library(lmerTest)
library(chron) #contains the funciton times
library(lme4)  #contains lmer for linear mixed effect modelling
library(arm)   #contains the function sim
library(car)   #contains the function Anova
library(blmeco) #contains compareqqnorm  (multiple qq boxplots)


#### worker mortality ####

head(mortality)
attach(mortality)

plot(worker_mortality ~ treatment)
mean(worker_mortality)
sd(worker_mortality)
shapiro.test(worker_mortality) #significant -> not normally distributed --> non parametric test --> Mann-Whitney U-test
wilcox.test(worker_mortality ~ treatment, exact=FALSE)

detach(mortality)

#### Virus titers ####

### DVW 
queens <- subset(col, sample.1 == "queen")
positive_queens <- subset(queens, PN_assignment_DWVb == 1)
shapiro.test(log10(positive_queens$vcps_DWVb))
mean(log10(positive_queens$vcps_DWVb))
sd(log10(positive_queens$vcps_DWVb))
quantile(log10(positive_queens$vcps_DWVb))

workers <- subset(col, sample.1 == "pooled_workers")
positive_workers <- subset(workers, PN_assignment_DWVb == 1)
shapiro.test(log10(positive_workers$vcps_DWVb))
mean(log10(positive_workers$vcps_DWVb))
sd(log10(positive_workers$vcps_DWVb))

# indiviudal workers
individual_workers <- subset(in_col, sample == "worker")
positive_worker <- subset(individual_workers, PN_assignment_DWVb == 1)
shapiro.test(log10(positive_worker$vcps_DWVb))
mean(log10(positive_worker$vcps_DWVb))
sd(log10(positive_worker$vcps_DWVb))

# eggs


# feeding pupae
fp <- subset(in_col, sample == "feeding_pupae")
shapiro.test(log10(fp$vcps_DWVb))
mean(log10(fp$vcps_DWVb))
sd(log10(fp$vcps_DWVb))

### ABPV
queens <- subset(queens, treatment == "virus")
positive_queens <- subset(queens, PN_assignment_ABPV == 1)
shapiro.test(log10(positive_queens$vcps_ABPV))
mean(log10(positive_queens$vcps_ABPV))
sd(log10(positive_queens$vcps_ABPV))
quantile(log10(positive_queens$vcps_ABPV))

control_queen <- log10(3.49E+06)

positive_workers <- subset(workers, PN_assignment_ABPV == 1)
shapiro.test(log10(positive_workers$vcps_ABPV))
mean(log10(positive_workers$vcps_ABPV))
sd(log10(positive_workers$vcps_ABPV))

# indiviudal workers
shapiro.test(log10(individual_workers$vcps_ABPV))
mean(log10(individual_workers$vcps_ABPV))
sd(log10(individual_workers$vcps_ABPV))

# eggs

# feeding pupae 
shapiro.test(log10(fp$vcps_ABPV))
mean(log10(fp$vcps_ABPV))
sd(log10(fp$vcps_ABPV))



#### Movement Analyses #### 


#### prerequisites #### 
dat <- mov

dat$treat1 <- as.factor(dat$treat1)
dat$treat1 <- factor(dat$treat1, levels = c("control", "low", "high"))
dat$treat2 <- as.factor(dat$treat2)

#response variables: initial speed (ini_speed), time_active, av_speed, overall movment --> sum mov
#random factor: colonies and day (and eventually time)

#transform time into Minutes only
dat$time2 <- times(dat$time)
dat$time_min <- 60 * hours(dat$time2) + minutes(dat$time2)
dat$time_min

# create a subset (half the dataset to exclude all samples from the neonic treatment)
dat <- subset(dat, treat1 == "control")





#### overall movement ####
par(mfrow=c(1,1))
par(oma=c(1,1,0,0))
boxplot(sum_move ~ date, data=dat)  #day might have an effect and will be included as random 
plot(dat$sum_move ~ dat$time_min)
abline(lm(dat$sum_move ~ dat$time_min))              # time of day not relevant 
hist(dat$sum_move)
boxplot(dat$sum_move ~ dat$treat2) # plot data with boxplots to have a first idea of the data

mod <- lmer(sum_move ~  treat2 + (1|colony) + (1|date), data=dat, REML=FALSE)
summary(mod)
Anova(mod)

#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$colony[,1]) # qq of random effects
qqline(ranef(mod)$colony[,1])
qqnorm(ranef(mod)$date[,1]) # The varince of date is so small that in the qq-plot it collapses to zero --> it is a bug and it can be assumed to be close to zero
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~dat$date)

# simulate
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
head(newdat)
Xmat <- model.matrix(~treat2, data=newdat)    
fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 
newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat

# Graph: 
col <- c("grey78", "grey78")
par(mfrow=c(1,1))
par(oma=c(0,0,0,0))
par(mar =c(5.1, 5.1, 4.1, 2.1))
boxplot(dat$sum_move ~ dat$treat2, border = "white", 
        main = "", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
axis(1, at = c(1,2), labels=c("no virus", "virus"), cex.axis = 1.5)
title(ylab="overall movment [~cm]", line=3.3, cex.lab=1.5)
stripchart(sum_move ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19), col= col,
           add = TRUE) 
#plot posteriors
points(c(1,2), newdat$fit, pch = 19, cex=2)
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5)

#Calculate the effectsize respectively how much the virus affects movement
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))
percentile <- ecdf(full)
percentile(0)











#### time active active_time --> time incative #### 
par(mfrow=c(1,1))
dat$incative <- 120-dat$active_time
dat$inactive_time <- round(dat$inactive_time, digits = 0)

hist(dat$inactive)
boxplot(dat$inactive ~ dat$treat2)

boxplot(dat$inactive_time ~ dat$date)             
plot(dat$inactive_time ~ dat$time_min)
abline(lm(dat$inactive_time ~ dat$time_min))         


mod <- glmer(inactive_time ~ treat2 + (1|colony) + (1|date), data = dat, family = "poisson")
summary(mod)
Anova(mod)
#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$colony[,1]) # qq of random effects
qqline(ranef(mod)$colony[,1])
qqnorm(ranef(mod)$date[,1]) # The varince of date is so small that in the qq-plot it collapses to zero --> it is a bug and it can be assumed to be close to zero
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~dat$date)
compareqqnorm(mod)


#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)
apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat1<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
head(newdat1)
Xmat <- model.matrix(~treat2, data=newdat1)    
fitmat <- matrix(ncol=nsim, nrow=nrow(newdat1))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 
newdat1$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat1$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat1$fit <- Xmat%*%fixef(mod)
newdat1

#Graph: 
par(mfrow=c(1,1))
par(oma=c(1,1,0,0))
boxplot(dat$inactive_time ~ dat$treat2, border = "white", 
        main = "time inactive", ylab = "time inactive [s]", xlab = "treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
axis(1, at = c(1,2), labels=c("no virus", "virus"), cex.axis = 1.5)
title(ylab="time inactive [s]", line=3, cex.lab=1.5)
stripchart(inactive_time ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE) 
#plot posteriors
points(c(1,2), newdat1$fit, pch = 19, cex=1)
segments(c(1,2), newdat1$lower, c(1,2), newdat1$upper, lwd = 5)

#Calculate the effectsize respectively how much the virus affects movement
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))







#### average speed av_speed #### 
par(mfrow=c(1,1))
hist(dat$av_speed, breaks = 7)
boxplot(dat$av_speed ~ dat$treat2)

boxplot(dat$av_speed ~ dat$date)              #day might have an effect and will be included as a random factor
plot(dat$av_speed ~ dat$time_min)      
abline(lm(dat$av_speed ~ dat$time_min))     #time of day  not relevant 

mod <- lmer(av_speed ~  treat2 + (1|colony) + (1|date), data=dat, REML=FALSE)
summary(mod)
Anova(mod)

#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$colony[,1]) # qq of random effects
qqline(ranef(mod)$colony[,1])
qqnorm(ranef(mod)$date[,1]) 
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~dat$date)

nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat2<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
head(newdat2)
Xmat <- model.matrix(~treat2, data=newdat2)      ####second part of model without of random factors
fitmat <- matrix(ncol=nsim, nrow=nrow(newdat2))

for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 
newdat2$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat2$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat2$fit <- Xmat%*%fixef(mod)
newdat

#Graph
par(mfrow=c(1,1))
par(oma=c(1,1,0,0))
boxplot(dat$av_speed ~ dat$treat2, border = "white", 
        main = "", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n")
axis(1, at = c(1,2), labels=c("no virus", "virus"))
title(ylab="average speed [~cm/s]", line=3)
stripchart(av_speed ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19), col= col,
           add = TRUE, cex= 0.8) 
#plot posteriors
points(c(1,2), newdat2$fit, pch = 19, cex=2)
segments(c(1,2), newdat2$lower, c(1,2), newdat2$upper, lwd = 5)


#Calculate the effectsize respectively how much the virus affects movement
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))









#### initial speed ####
dat$ini_speed <- dat$ini_speed/10    #adapt initial speed so it almost reflects cm/second
hist(dat$ini_speed, breaks = 7)
boxplot(dat$ini_speed ~ dat$treat2) # plot data with boxplots to have a first idea of the data
boxplot(ini_speed ~ date, data=dat)      #day might have an effect and will be included as a random factor 
plot(dat$ini_speed ~ dat$time_min)
abline(lm(dat$ini_speed ~ dat$time_min)) # time of day does not really seem to be relevant and will not be included in the model


mod <- lmer(ini_speed ~  treat2 + (1|colony) + (1|date), data=dat, REML=FALSE)
summary(mod)
Anova(mod)


#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$colony[,1]) # qq of random effects
qqline(ranef(mod)$colony[,1])
qqnorm(ranef(mod)$date[,1]) # The varince of date is so small that in the qq-plot it collapses to zero --> it is a bug and it can be assumed to be close to zero
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~dat$date)

#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))


#fitted values with 95% credible intervals
newdat3<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
head(newdat3)

Xmat <- model.matrix(~treat2, data=newdat3)      #second part of model without of random factors
fitmat <- matrix(ncol=nsim, nrow=nrow(newdat3))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat3$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat3$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat3$fit <- Xmat%*%fixef(mod)
newdat3

#Graph
par(mfrow=c(1,1))


boxplot(dat$ini_speed ~ dat$treat2, border = "white", 
        main = "", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n")
axis(1, at = c(1,2), labels=c("no virus", "virus"))
title(ylab="initial speed [~cm/s]", line=2.5)
stripchart(ini_speed ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19), col= col,
           add = TRUE, cex=0.8) 
#plot posteriors
points(c(1,2), newdat3$fit, pch = 19, cex=2)
segments(c(1,2), newdat3$lower, c(1,2), newdat3$upper, lwd = 5)

#Calculate the effectsize respectively how much the virus affects movement
str(fitmat)
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))






#### Colonysize increase ####

ANV <- subset(ANV, treatment_1 =="control") #remove neonic treatment
levels(ANV$treatment_2)[levels(ANV$treatment_2)=="control"] <- "no virus"
ANV <- subset(ANV, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")
ANV_reduced <- subset(ANV, id!="C8V")

# initial colony size 
date1 <- subset(ANV, date == "16.07.2017")
date2 <- subset(ANV, date == "01.10.2017")
date1_reduced <- subset(ANV_reduced, date == "16.07.2017")
date2_reduced <- subset(ANV_reduced, date == "01.10.2017")

boxplot(date1$adults ~ date1$treatment_2)
shapiro.test(date1$adults) # non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date1$adults, date1$treatment_2) # nonsignificant --> homoheneity of variance can be assumed
t.test(date1$adults ~ date1$treatment_2) 
mean(date1$adults[date1$treatment_2=="virus"])
sd(date1$adults[date1$treatment_2=="virus"])
mean(date1$adults[date1$treatment_2=="no virus"])
sd(date1$adults[date1$treatment_2=="no virus"])

boxplot(date1$pupae ~ date1$treatment_2)
shapiro.test(date1$pupae)# non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date1$pupae, date1$treatment_2) # nonsignificant --> homoheneity of variance can be assumed
t.test(date1$pupae ~ date1$treatment_2) 

mean(date1$pupae[date1$treatment_2=="virus"])
sd(date1$pupae[date1$treatment_2=="virus"])
mean(date1$pupae[date1$treatment_2=="no virus"])
sd(date1$pupae[date1$treatment_2=="no virus"])


# colony size at end 
boxplot(date2$adults ~ date2$treatment_2)
shapiro.test(date2$adults) # non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date2$adults, date2$treatment_2) # nonsignificant --> homoheneity of variance can be assumed
t.test(date2$adults ~ date2$treatment_2) 

boxplot(date2$pupae ~ date2$treatment_2)
shapiro.test(date2$pupae) # non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date2$pupae, date2$treatment_2) # nonsignificant --> homogeneity of variance can be assumed
t.test(date2$pupae ~ date2$treatment_2) 

mean(date2$pupae[date2$treatment_2=="virus"])
sd(date2$pupae[date2$treatment_2=="virus"])
mean(date2$pupae[date2$treatment_2=="no virus"])
sd(date2$pupae[date2$treatment_2=="no virus"])

# with the outlier removed
boxplot(date2_reduced$adults~date2_reduced$treatment_2)
shapiro.test(date2_reduced$adults) # non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date2_reduced$adults, date2_reduced$treatment_2) # nonsignificant --> homoheneity of variance can be assumed
t.test(date2_reduced$adults ~ date2_reduced$treatment_2) 

boxplot(date2_reduced$pupae~date2_reduced$treatment_2)
shapiro.test(date2_reduced$pupae) # non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date2_reduced$pupae, date2_reduced$treatment_2) # nonsignificant --> homoheneity of variance can be assumed
t.test(date2_reduced$pupae ~ date2_reduced$treatment_2) 




# difference in the colonysize 
date1$diff_adults <- date2$adults-date1$adults
date1$diff_adults
boxplot(date1$diff_adults ~ date1$treatment_2, xlab= "treatments", ylab = "nr. of adult ants", main = "")
shapiro.test(date1$diff_adults) # non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date1$diff_adults, date2$treatment_2) # nonsignificant --> homoheneity of variance can be assumed
t.test(date1$diff_adults ~ date2$treatment_2) 
mean(date1$diff_adults[date1$treatment_2=="virus"])
sd(date1$diff_adults[date1$treatment_2=="virus"])
mean(date1$diff_adults[date1$treatment_2=="no virus"])
sd(date1$diff_adults[date1$treatment_2=="no virus"])

mod <- lmer(diff_adults ~ treatment_2 + (1|adults), data = date1, REML = FALSE) # singular fit warning... 
summary(mod)
Anova(mod)


date1$diff_pupae <- date2$pupae - date1$pupae
date1$diff_pupae
boxplot(date1$diff_pupae ~ date1$treatment_2, xlab= "treatments", ylab = "nr. of adult ants", main = "")
shapiro.test(date1$diff_pupae) # non significant --> Nullhypothesis that Data is normally distributed can be accepted
leveneTest(date1$diff_pupae, date2$treatment_2) # nonsignificant --> homoheneity of variance can be assumed
t.test(date1$diff_pupae ~ date2$treatment_2) 


# queen weight / weight adults: 
ANV_final <- subset(ANV, day == 432)
head(ANV_final)
controls <- subset(ANV_final, treatment_2 == "no virus")
virus <- subset(ANV_final, treatment_2 == "virus")

boxplot(ANV_final$weight_queen ~ ANV_final$treatment_2)
shapiro.test(ANV_final$weight_queen) # non significant --> nullhypothesis of normally distributed data can be assumed --> t-test.
t.test(ANV_final$weight_queen ~ ANV_final$treatment_2)

mean(controls$weight_queen)
sd(controls$weight_queen)
mean(virus$weight_queen)
sd(virus$weight_queen)

boxplot(ANV_final$weight_20adults ~ ANV_final$treatment_2)
shapiro.test(ANV_final$weight_20adults) # non significant --> nullhypothesis of normally distributed data can be assumed --> t-test.
t.test(ANV_final$weight_20adults ~ANV_final$treatment_2)

mean(controls$weight_20adults)
sd(controls$weight_20adults)
mean(virus$weight_20adults)
sd(virus$weight_20adults)


#### Graphs for manuscript ####

# 1x4 graph crediabl intervals:  

par(mfrow=c(1,4))
par(oma=c(0,0,0,0))
par(mar =c(5.1, 5.1, 3.1, 2.1))
labels <- c("No virus", "Virus")

boxplot(dat$sum_move ~ dat$treat2, border = "white", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Distance [~cm]", line=3.5, cex.lab=1.5)
stripchart(sum_move ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19),  col= col, cex = 1,
           add = TRUE) 
points(c(1,2), newdat$fit, pch = 19, cex=1.5)
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 3.5)
text(c(1:2), 410, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('A', side=3, line=0.5, at=-0, cex = 1.5, font = 2)

boxplot(dat$inactive_time ~ dat$treat2, border = "white", 
        ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Time inactive [s]", line=3, cex.lab=1.5)
stripchart(inactive_time ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19), col= col, cex = 1,
           add = TRUE) 
points(c(1,2), newdat1$fit, pch = 19, cex=1.5)
segments(c(1,2), newdat1$lower, c(1,2), newdat1$upper, lwd = 3.5)
mtext('B', side=3, line=0.5, at=-0, cex = 1.5, font = 2)
text(c(1:2), 87.5, labels = c("a", "a"), font = 2, cex = 1.5)

boxplot(dat$av_speed ~ dat$treat2, border = "white", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Average speed [~cm/s]", line=3.5, cex.lab=1.5)
stripchart(av_speed ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19), col= col,
           add = TRUE, cex = 1) 
points(c(1,2), newdat2$fit, pch = 19, cex=1.5)
segments(c(1,2), newdat2$lower, c(1,2), newdat2$upper, lwd = 3.5)
mtext('C', side=3, line=0.5, at=-0, cex = 1.5, font = 2)
text(c(1:2), 3.55, labels = c("a", "b"), font = 2, cex = 1.5)

boxplot(dat$ini_speed ~ dat$treat2, border = "white", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Initial speed [~cm/s]", line=3, cex.lab=1.5)
stripchart(ini_speed ~ treat2, data = dat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19), col= col,
           add = TRUE, cex = 1) 
points(c(1,2), newdat3$fit, pch = 19, cex=1.5)
segments(c(1,2), newdat3$lower, c(1,2), newdat3$upper, lwd = 3.5)
mtext('D', side=3, line=0.5, at=-0, cex = 1.5, font = 2)
text(c(1:2), 5.93, labels = c("a", "b"), font = 2, cex = 1.5)

# 2x2 boxplot graph with Overall movement, time inactive, overall movement, initial speed.

par(mfrow=c(2,2))

boxplot(dat$sum_move ~ dat$treat2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Distance [~cm]", line=3.5, cex.lab=1.5)
text(c(1:2), 410, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('A', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)

boxplot(dat$inactive_time ~ dat$treat2, border = "black", 
        ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Time inactive [s]", line=3, cex.lab=1.5)
mtext('B', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)
text(c(1:2), 87.5, labels = c("a", "a"), font = 2, cex = 1.5)

boxplot(dat$av_speed ~ dat$treat2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Average speed [~cm/s]", line=3.5, cex.lab=1.5)
mtext('C', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)
text(c(1:2), 3.55, labels = c("a", "b"), font = 2, cex = 1.5)

boxplot(dat$ini_speed ~ dat$treat2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="Initial speed [~cm/s]", line=3, cex.lab=1.5)
mtext('D', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)
text(c(1:2), 5.93, labels = c("a", "b"), font = 2, cex = 1.5)



# 2x2 Graph with initial colony size and difference in the number of adults+pupae

par(mfrow=c(2,2))

boxplot(date1$adults ~ date1$treatment_2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="# Adults", line=3.5, cex.lab=1.5)
text(c(1:2), 53.5, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('A', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)

boxplot(date1$pupae ~ date1$treatment_2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="# Pupae", line=3.5, cex.lab=1.5)
text(c(1:2), 10, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('B', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)

boxplot(date1$diff_adults ~ date1$treatment_2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab=expression(paste(Delta, " Adults")), line=3.5, cex.lab=1.5)
text(c(1:2), 126, labels = c("a", "b"), font = 2, cex = 1.5)
mtext('C', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)

boxplot(date1$diff_pupae ~ date1$treatment_2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", ylim= c(-11,17),
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab=expression(paste(Delta, " Pupae")), line=3.5, cex.lab=1.5)
text(c(1:2), 16, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('D', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)

# 1x2 Graph colonysize, just the adults: 
par(mfrow=c(1,2))
boxplot(date1$adults ~ date1$treatment_2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab="# Adults", line=3.5, cex.lab=1.5)
text(c(1:2), 53.5, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('A', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)
boxplot(date1$diff_adults ~ date1$treatment_2, border = "black", 
        main = "", ylab = "", xlab = "Treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, las=1)
axis(1, at = c(1,2), labels=labels, cex.axis = 1.5)
title(ylab=expression(paste(Delta, " # Adults")), line=3.5, cex.lab=1.5)
text(c(1:2), 126, labels = c("a", "b"), font = 2, cex = 1.5)
mtext('B', side=3, line=0.5, at=0.3, cex = 1.5, font = 2)



#### Virus titers in the movement experiment  ####
head(virus_mov)

#subsets: Only no-neonic treatments and only those for which virus analyses have been done. 
data <- subset(virus_mov, treatment_1 == "control")
data <- subset(data, virus_analyses == "y")
data

#split dataset for queens and workers
workers <- subset(data, sample =="worker")
queens  <- subset(data, sample =="queen")

boxplot(workers$virus_titer~workers$treatment_2)
boxplot(queens$virus_titer~queens$treatment_2)

#significant difference between treatments with regards of virus titers in queens and workers  

shapiro.test(workers$virus_titer) #significant -> not normally distributed -> Mann whitney U test
wilcox.test(workers$virus_titer ~ workers$treatment_2)

shapiro.test(queens$virus_titer) #significant -> not normally distributed -> Mann whitney U test
wilcox.test(queens$virus_titer ~ queens$treatment_2)

#(unfortunetely one queen so much positive that it overrules the treatment effect?)
# insted compare the percentages of samples positive with a fisher test (exact and ok for small samplesizes of not normally distributed data)



workers$virus_qualitative[workers$treatment_2 == "control"]
workers$virus_qualitative[workers$treatment_2 == "virus"]
a <- c(8,3)
b <- c(0,4)
chi <- data.frame(a,b)
fisher.test(chi)

queens$virus_qualitative[queens$treatment_2 == "control"]
queens$virus_qualitative[queens$treatment_2 == "virus"]

queens$id

a <- c(5,1)
b <- c(3,7)
chi <- data.frame(a,b)
fisher.test(chi)


# correlation of initial_speed and virus load within virus treatment
treatmentgroup <- subset(workers, treatment_2 == "virus")
treatmentgroup$id

plot(treatmentgroup$virus_titer)
head(mov)

a <- c(mean(mov$ini_speed[mov$id == "C1V"]),
  mean(mov$ini_speed[mov$id == "C4V"]),
  mean(mov$ini_speed[mov$id == "C5V"]),
  mean(mov$ini_speed[mov$id == "C11V"]),
  mean(mov$ini_speed[mov$id == "C13V"]),
  mean(mov$ini_speed[mov$id == "C17V"]),
  mean(mov$ini_speed[mov$id == "C18V"]))

cor <- data.frame(a,treatmentgroup$virus_titer)
cor.test(cor$a, cor$treatmentgroup.virus_titer, method = "spearman") ### no correlation and a much to small samplesize









