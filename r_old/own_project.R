#########################################################################################
### Baeysian Stats course ###############################################################
### own project #########################################################################
#########################################################################################


getwd()
setwd("/Users/gismo/Desktop/stat_course/R")
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV")
library(lme4)

dat <- read.table("speed.txt", header = TRUE)
head(dat)
str(dat)

dat$treat1 <- as.factor(dat$treat1)
dat$treat1 <- factor(dat$treat1, levels = c("control", "low", "high"))
dat$treat2 <- as.factor(dat$treat2)


#response variables: initial speed (inispeed), time_active, avspeed movment nrgrooming
#random factor: colonies and day

#transform time into Minutes only
install.packages("chron")
library(chron) 
dat$time2 <- times(dat$time)
dat$time_min <- 60 * hours(dat$time2) + minutes(dat$time2)
dat$time_min

plot(dat$inispeed ~ dat$date) #day might have an effect and will be included as random foactor
plot(dat$inispeed ~ dat$time_min)
mod_time <- lm(dat$inispeed ~ dat$time_min)
abline(mod_time)              # time of day does not really seem to be relevant 

#### For initial speed ####
hist(dat$ini_speed, breaks = 7)

par(mfrow=c(1,1))
boxplot(dat$ini_speed ~ dat$treat1 + dat$treat2) # plot data with boxplots to have a first idea of the data

mod <- lmer(ini_speed ~ treat1 + treat2 + treat1:treat2 + time_min + (1|colony) + (1|date), data=dat, REML=FALSE)
mod
summary(mod)
str(mod)

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

#create a super nice Graph: 

par(mfrow=c(1,1))
par(oma=c(1,0,0,0))

col <- c("darkolivegreen1", "darkorange3", "darkred")
boxplot(dat$ini_speed ~ dat$treat1 + dat$treat2, border = "white", 
        main = "Initial speed", ylab = "initial movment speed [~cm/s]", xlab = "Treatment", 
        family = c("mono"), xaxt="n")
axis(1, at = c(1.2,2,2.8,4.2,5,5.8), labels=c("","no virus","","","virus",""))

stripchart(inispeed ~ treat1 + treat2, data = dat, 
           vertical = TRUE, at = c(1.2,2,2.8,4.2,5,5.8) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col = alpha(col, 0.8),
           add = TRUE) 

#plot posteriors
points(c(1.2,2,2.8,4.2,5,5.8), newdat$fit, pch = 19, cex=1.5)
segments(c(1.2,2,2.8,4.2,5,5.8), newdat$lower, c(1.2,2,2.8,4.2,5,5.8), newdat$upper, lwd = 3)

#add legends 
legend(0.3, 13.4, c("control", "low", "high"), col = col, title = "neonicotinoid treatment", pch = c(19) , bty = "o", bg="transparent")
legend(0.7, 13.2, c("","",""), col = col, title = "", pch = c(17) , bty = "n", bg="transparent")
#add graph description

mtext("Real data of the initial speed plotted with the predicted mean and its 95% credible interval for each treatment", side = 1, outer = TRUE)

str(fitmat)
head(fitmat)
fitmat[,1]

#Calculate the effectsize of the non-virus / virus treatments (credible interfavals for differences in the mean)




novirus <- apply(fitmat[c(1,2,3),] , 2, mean)


virus <- apply(fitmat[c(4,5,6),], 2, mean)

full <- novirus-virus
quantile(full, prob = c(0.025, 0.5, 0.975))


#Calculate differences in means between the neonic treatments

control <- apply(fitmat[c(1,4),],2,mean)
low <- apply(fitmat[c(2,5),],2,mean)
high <- apply(fitmat[c(3,6),],2,mean)

cl <- control-low
  quantile(cl, prob = c(0.025, 0.5, 0.975))
ch <- control-high
  quantile(ch, prob = c(0.025, 0.5, 0.975))
lh <- low-high
  quantile(lh, prob = c(0.025, 0.5, 0.975))

#we have now some 95 confidence intervalls of the possible effect size of the neonicotinoid treatment that range from -8/8 to -11/4
#based on these assumtions there seems to be no change speed of biological relevance due to the neonicotinoid treatmetns

  
  
#### For overall movement ####
  
head(dat)
hist(dat$movment)
  
par(mfrow=c(1,1))
boxplot(dat$movment ~ dat$treat1 + dat$treat2) # plot data with boxplots to have a first idea of the data
  
mod <- glmer(movment ~ treat1 + treat2 + treat1:treat2 + (1|colony), data=dat, family = poisson, REML=FALSE)  #add  + (1|date) + (1|time) as soon as the variables are converted to the right form
mod
summary(mod)
str(mod)
  
#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$colony[,1]) # qq of random effects
qqline(ranef(mod)$colony[,1])
  
#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)
  
apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))
  
#fitted values with 95% credible intervals
newdat<-expand.grid(
treat1=factor(c('control','low','high'), levels=levels(dat$treat1)),
treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
  
Xmat <- model.matrix(~treat1 + treat2 + treat1:treat2, data=newdat)      ####second part of model without of random factors
  
fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 
  
newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat
  
#create a super nice Graph: 
par(mfrow=c(1,1))
boxplot(dat$movment ~ dat$treat1 + dat$treat2, border = "white")
stripchart(movment ~ treat1 + treat2, data = dat, 
          vertical = TRUE, method = "jitter", 
          pch = c(19,19,19,17,17,17), col = alpha(c("orange", "red", "darkred"), 0.5),
          add = TRUE) 
  
points(1:6, newdat$fit, pch = 19, cex=1.5)
segments(1:6, newdat$lower, 1:6, newdat$upper, lwd = 3)
  
str(fitmat)
head(fitmat)
fitmat[,1]











#### For average speed ####






#### Effect of the Treatments on colony Development ####
head(dat)

datc <- subset(dat, final == "year2" )


datc <- read.table("Year2.txt", header = TRUE)
head(datc)
str(dat)

summary(datc$treatment_1)
summary(datc$treatment_2)


datc$treatment_1 <- as.factor(datc$treatment_1)
datc$treatment_1 <- factor(datc$treatment_1, levels = c("control", "low", "high"))
datc$treatment_2 <- as.factor(datc$treatment_2)
datc

datc$treatment_1


datc$workforce <- datc$pupae+datc$adults
datc$workforce
mean(datc$workforce)

boxplot(workforce ~ treatment_1 + treatment_2, data=datc)
boxplot(workforce ~ treatment_2 + treatment_1, data=datc)

#exclude colonies with queenfailure:

datc <- subset(datc, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")
summary(datc$treatment_1)

# The chosen model might be  a two way ANOVA --> FIT  (with interaction?) 

# modc <- glm(workforce ~ treatment_2 + treatment_1 +treatment_2:treatment_1, data=datc, family=poisson)  # compare to discussion in file: Eurobee.R

modc <- lm(workforce ~ treatment_2 + treatment_1 +treatment_2:treatment_1, data=datc)
modc
summary(modc)
summary(modc)$sigma
Anova(modc)



hist(datc$workforce, breaks = 12)
shapiro.test(datc$workforce)

#check the model assumptions
par(mfrow=c(2,2))
plot(modc)      
acf(resid(modc))
datc$treat21 <- factor(paste(datc$treat2, datc$treat1))
datc$treat21
plot(resid(modc)~datc$treat21)

#bayesian inference using sim
# add predicted values (group means)

newdatc<-expand.grid(
  treatment_2=factor(c('control','virus'),levels=levels(datc$treatment_2)),
  treatment_1=factor(c('control','low','high'), levels=levels(datc$treatment_1)))
newdatc$fit <- predict(modc, newdata=newdatc)
nsim <- 2000
bsimc <- sim(modc, n.sim=nsim)
str(bsimc)
fitmat <- matrix(ncol=nsim, nrow=nrow(newdatc))
Xmat <- model.matrix(formula(modc)[c(1,3)], data=newdatc)
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsimc@coef[i,]
newdatc$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdatc$upper <- apply(fitmat, 1, quantile, prob=0.975)

newdatc


#create a super nice Graph: 
par(mfrow=c(1,1))
col2 <- c("darkolivegreen1","darkolivegreen1", "darkorange3","darkorange3", "darkred", "darkred")
boxplot((datc$workforce) ~ datc$treatment_2 + datc$treatment_1, border = "white", 
        main = "Effects on colonysize", ylab = "# of adults and pupae", xlab = "treatment (c=control, v=virus)", 
        family = c("mono"), xaxt="n")

boxplot(log(datc$workforce) ~ datc$treatment_2 + datc$treatment_1, border = "white", 
        main = "Effect on colonysize", ylab = "# of adults and pupae", xlab = "treatment (c=control, v=virus)", 
        family = c("mono"), xaxt="n")

x <- c("c","v") 
x <- rep(x, times=3)
axis(1, at = c(1.2, 1.8, 3.2, 3.8, 5.2,5.8), labels=x)

stripchart(log(datc$workforce) ~ datc$treatment_2 + datc$treatment_1, data = datc, 
           vertical = TRUE, at = c(1.2, 1.8, 3.2, 3.8, 5.2,5.8) ,method = "jitter", 
           pch = c(19,17,19,17,19,17), col = col2,
           add = TRUE) 

stripchart((datc$workforce) ~ datc$treatment_2 + datc$treatment_1, data = datc, 
           vertical = TRUE, at = c(1.2, 1.8, 3.2, 3.8, 5.2,5.8) ,method = "jitter", 
           pch = c(19,17,19,17,19,17), col = col2,
           add = TRUE) 


#plot posteriors and add legends
points(c(1.2, 1.8, 3.2, 3.8, 5.2,5.8)+0.15, newdatc$fit, pch = 19, cex=1.5)
segments(c(1.2, 1.8, 3.2, 3.8, 5.2,5.8)+0.15, newdatc$lower, c(1.2, 1.8, 3.2, 3.8, 5.2,5.8)+0.15, newdatc$upper, lwd = 3)

col <- c("darkolivegreen1", "darkorange3", "darkred")
legend(5, 350, c("control", "low", "high"), col = col, title = "neonicotinoid treatment", pch = c(19) , bty = "n", bg="transparent")
legend(5.25, 349, c("","",""), col = col, title = "", pch = c(17) , bty = "n", bg="transparent")


legend(4.5, 6, c("control", "low", "high"), col = col, title = "neonicotinoid treatment", pch = c(19) , bty = "n", bg="transparent")
legend(4.75, 6, c("","",""), col = col, title = "", pch = c(17) , bty = "n", bg="transparent")

#Calculate differences in means between the non-virus / virus treatments
novirus <- c(fitmat[1,],fitmat[3,],fitmat[5,])
virus <- c(fitmat[2,],fitmat[4,],fitmat[6,])
quantile((novirus-virus), prob = c(0.025, 0.5, 0.975))

# The effect size of virus on the colony development is with a 95% a decrease between 4 and 28 Individuals 
#2.5%       50%     97.5% 
#-18.84831  13.14856  54.88614

 x <- apply(fitmat[c(1,3,5),],2, mean)
 y <- apply(fitmat[c(2,4,6),],2, mean)
 quantile((x-y), prob = c(0.025, 0.5, 0.975))

#Calculate differences between the neonicotinoid treatments: 
control <-  c(fitmat[1,], fitmat[2,])
low <-  c(fitmat[3,], fitmat[4,])
high <-  c(fitmat[5,], fitmat[6,])

x <- apply(fitmat[c(1,2),],2, mean)
y <- apply(fitmat[c(3,4),],2, mean)
z <- apply(fitmat[c(5,6),],2, mean)

cl <- x-y
quantile(cl, prob = c(0.025, 0.5, 0.975))
ch <- x-z
quantile(ch, prob = c(0.025, 0.5, 0.975))
lh <- y-z
quantile(lh, prob = c(0.025, 0.5, 0.975))




