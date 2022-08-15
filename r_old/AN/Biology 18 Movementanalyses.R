


### analysing ant movement depending on virus/neonic for biology 18 ###

setwd("H:/DS-Bees/R/ANV")
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV")
install.packages("chron")
install.packages("lmerTest")
install.packages("arm")
install.packages("car")
library(lmerTest)
library(chron)
library(lme4)
library(arm)
library(car)

#### Virus in the samples ####
vir <- read.table("Lniger_vcps.txt", header = TRUE)
head(vir)

ABPV <- subset(vir, virus == "ABPV") 
DWV <- subset(vir, virus == "DWV")

###Values
#DWV
colonies <- subset(DWV, type=="colony" & assignment==1  )
summary(colonies)

queens <- subset(DWV, type=="queen" & assignment==1  )
summary(queens)

fp <- subset(DWV, type=="feeding_pupae" & assignment==1)
summary(fp)

sa <- subset(DWV, type=="single_ant" & assignment==1)
summary(sa)

#ABPV
colonies <- subset(ABPV, type=="colony" & assignment==1 )
summary(colonies)
colonies <- subset(ABPV, type=="queen" & assignment==1 & sample== "C7")
summary(colonies)


queens <- subset(ABPV, type=="queen" & assignment==1 & sample!= "Q7")
summary(queens)
queens <- subset(ABPV, type=="queen" & assignment==1 & sample== "Q7")
summary(queens)

fp <- subset(APBV, type=="feeding_pupae" & assignment==1)
summary(fp)

sa <- subset(ABPV, type=="single_ant" & assignment==1)
summary(sa)








#### Effects of the viruses ####

dat <- read.table("speed.txt", header = TRUE)
head(dat)
str(dat)

dat$treat1 <- as.factor(dat$treat1)
dat$treat1 <- factor(dat$treat1, levels = c("control", "low", "high"))
dat$treat2 <- as.factor(dat$treat2)


#response variables: initial speed (inispeed), time_active, avspeed movment nrgrooming
#random factor: colonies and day
#transform time into Minutes only
dat$time2 <- times(dat$time)
dat$time_min <- 60 * hours(dat$time2) + minutes(dat$time2)
dat$time_min




##### And now only for the controls that did not have neonics for biology 18  #####

#### initial speed ####

cdat <- subset(dat, treat1 == "control") 
cdat$ini_speed <- cdat$ini_speed/10


boxplot(ini_speed ~ date, data=dat)  #day might have an effect and will be included as random 
boxplot(ini_speed ~ date, data=cdat)


plot(cdat$ini_speed ~ cdat$time_min)
mod_time <- lm(cdat$ini_speed ~ cdat$time_min)
abline(mod_time)              # time of day does not really seem to be relevant 

hist(cdat$ini_speed, breaks = 7)
par(mfrow=c(1,1))
boxplot(cdat$ini_speed ~ cdat$treat2) # plot data with boxplots to have a first idea of the data
mod <- lmer(ini_speed ~  treat2 + time_min + (1|colony) + (1|date), data=cdat, REML=FALSE)
mod
summary(mod)
str(mod)
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
boxplot(resid(mod)~cdat$date)

#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)

apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
newdat$time_min=mean(cdat$time_min)
head(newdat)

Xmat <- model.matrix(~treat2 + time_min, data=newdat)      ####second part of model without of random factors

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat

#create a super nice Graph: 
par(mar=c(5,6,5,2))
par(mfrow=c(1,1))
par(oma=c(1,1,0,0))
par(cex=1)

col <- c("grey78", "grey78")
boxplot(cdat$ini_speed ~ cdat$treat2, border = "white", 
        main = "initial speed [10s]", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
axis(1, at = c(1,2), labels=c("no virus", "virus"), cex.axis = 1.5)
title(ylab="initial speed [~cm/s]", line=2.5, cex.lab=1.5)

stripchart(ini_speed ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE) 
#plot posteriors
points(c(1,2), newdat$fit, pch = 19, cex=2, col = c("darkolivegreen3", "red"))
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5, col = c("darkolivegreen3", "red"))
         
#Calculate the effectsize respectively how much the virus affects movement
fitmat[,1]


str(fitmat)
novirus <- apply(fitmat[1,],2, mean)
  virus <- apply(fitmat[1,],2, mean)
  
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))







#### overall movement ####
boxplot(sum_move ~ date, data=dat)  #day might have an effect and will be included as random 
boxplot(sum_move ~ date, data=cdat)


plot(cdat$sum_move ~ cdat$time_min)
mod_time <- lm(cdat$sum_move ~ cdat$time_min)
abline(mod_time)              # time of day does not really seem to be relevant 

hist(cdat$sum_move, breaks = 7)

par(mfrow=c(1,1))
boxplot(cdat$sum_move ~ cdat$treat2) # plot data with boxplots to have a first idea of the data

mod <- lmer(sum_move ~  treat2 + time_min + (1|colony) + (1|date), data=cdat, REML=FALSE)
mod
summary(mod)
str(mod)
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
boxplot(resid(mod)~cdat$date)

#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)

apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
newdat$time_min=mean(cdat$time_min)
head(newdat)

Xmat <- model.matrix(~treat2 + time_min, data=newdat)      ####second part of model without of random factors

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat

#create a super nice Graph: 

par(mfrow=c(1,1))
par(oma=c(1,1,0,0))

boxplot(cdat$sum_move ~ cdat$treat2, border = "white", 
        main = "overall movement [120s]", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
axis(1, at = c(1,2), labels=c("no virus", "virus"), cex.axis = 1.5)
title(ylab="overall movment [~cm]", line=3.3, cex.lab=1.5)

stripchart(sum_move ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE) 


#plot posteriors
points(c(1,2), newdat$fit, pch = 19, cex=2, col = c("darkolivegreen3", "red"))
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5, col = c("darkolivegreen3", "red"))

#Calculate the effectsize respectively how much the virus affects movement
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))









#### average speed av_speed #### 

plot(cdat$av_speed ~ cdat$time_min)
mod_time <- lm(cdat$av_speed ~ cdat$time_min)
abline(mod_time)              # time of day does not really seem to be relevant 

hist(cdat$av_speed, breaks = 7)

par(mfrow=c(1,1))
boxplot(cdat$av_speed ~ cdat$treat2) # plot data with boxplots to have a first idea of the data

mod <- lmer(av_speed ~  treat2 + time_min + (1|colony) + (1|date), data=cdat, REML=FALSE)
mod
summary(mod)
str(mod)
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
boxplot(resid(mod)~cdat$date)

#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)

apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
newdat$time_min=mean(cdat$time_min)
head(newdat)

Xmat <- model.matrix(~treat2 + time_min, data=newdat)      ####second part of model without of random factors

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat

#create a super nice Graph: 

par(mfrow=c(1,1))
par(oma=c(1,1,0,0))

boxplot(cdat$av_speed ~ cdat$treat2, border = "white", 
        main = "average speed", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
axis(1, at = c(1,2), labels=c("no virus", "virus"), cex.axis=1.5)
title(ylab="average speed [~cm/s]", line=3, cex.lab=1.5)


stripchart(av_speed ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE) 


#plot posteriors
points(c(1,2), newdat$fit, pch = 19, cex=2, col = c("darkolivegreen3", "red"))
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5, col = c("darkolivegreen3", "red"))

#Calculate the effectsize respectively how much the virus affects movement
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))









#### time active active_time #### 
head(cdat)

plot(cdat$active_time ~ cdat$time_min)
mod_time <- lm(cdat$active_time ~ cdat$time_min)
abline(mod_time)              # time of day does not really seem to be relevant 

hist(cdat$active_time, breaks = 7)

par(mfrow=c(1,1))
boxplot(cdat$active_time ~ cdat$treat2) # plot data with boxplots to have a first idea of the data

mod <- lmer(active_time ~  treat2 + time_min + (1|colony) + (1|date), data=cdat, REML=FALSE)
mod
summary(mod)
str(mod)
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
boxplot(resid(mod)~cdat$date)

#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)

apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
newdat$time_min=mean(cdat$time_min)
head(newdat)

Xmat <- model.matrix(~treat2 + time_min, data=newdat)      ####second part of model without of random factors

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat

#create a super nice Graph: 

par(mfrow=c(1,1))
par(oma=c(1,1,0,0))

boxplot(cdat$active_time ~ cdat$treat2, border = "white", 
        main = "time active", ylab = "time active [s]", xlab = "treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
axis(1, at = c(1,2), labels=c("no virus", "virus"), cex.axis = 1.5)
title(ylab="time active [s]", line=3, cex.lab=1.5)

stripchart(active_time ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE) 


#plot posteriors
points(c(1,2), newdat$fit, pch = 19, cex=2, col = c("darkolivegreen3", "red"))
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5, col = c("darkolivegreen3", "red"))

#Calculate the effectsize respectively how much the virus affects movement
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))







