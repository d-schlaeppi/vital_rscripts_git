setwd("H:/DS-Bees/R/ANV")
library(lme4)

dat <- read.table("speed.txt", header = TRUE)
head(dat)
str(dat)

summary(dat$treat1)

dat$treat1 <- as.factor(dat$treat1)
dat$treat2 <- as.factor(dat$treat2)

dat.c <- subset(dat, treat1=="control")





#transform time into Minutes only
install.packages("chron")
library(chron) 
dat.c$time2 <- times(dat.c$time)
dat.c$time_min <- 60 * hours(dat.c$time2) + minutes(dat.c$time2)
dat.c$time_min


hist(dat.c$ini_speed)
boxplot(dat.c$ini_speed ~ dat.c$treat2) 

mod <- lmer(ini_speed ~  treat2 + time_min + (1|colony) + (1|date), data=dat, REML=FALSE)
mod
summary(mod)
anova(mod)


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

??sim

### baeysian

nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)

apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat<-expand.grid(
  treat1=factor(c('control','low','high'), levels=levels(dat$treat1)),
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
newdat$time_min=mean(dat$time_min)
head(newdat)

Xmat <- model.matrix(~treat1 + treat2 + treat1:treat2 + time_min, data=newdat)      ####second part of model without of random factors

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat



