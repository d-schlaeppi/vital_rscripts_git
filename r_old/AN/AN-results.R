### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### AN - Ants & Neonics ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### prerequisites ####
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV")
setwd("/Users/gismo/Desktop/R/AN")

#read tables
# 
#ANV <- read.table("20190226_ANV.txt", header = TRUE) #mostly included in dat_all
#dat <- read.table("ANVlong_full.txt", header = TRUE) # mostly included in dat_all
# survival <- read.table("ANV_queen-survival.txt", header = TRUE) # now included in all dat. 
# DLS <- read.table("lifestage_duration.txt", header = TRUE) # now included in dat_all
# neonic <- read.table("ANV_Neonic_in_workers_and_queens.txt", header = TRUE) # now included in all_dat

dat_all <- read.table <- read.table("AN_data.txt", header = TRUE)
dat_all$treatment <- factor(dat_all$treatment, levels(dat_all$treatment)[c(1,3,2)])

#load libraries --> these libraries need to be cited
library(survival) # contains kalan meier plot function
library(car) #contains levene's test
library(arm) #contains the sim function for bayesian inference
library(dunn.test) #contains the dunn test
library(nlme) #lme
library(emmeans) #lsmeans to compare pairwise differences after the modeling using bonferroni corrections
library(blmeco) #contains compareqqnorm  (multiple qq boxplots)
library(lme4) #lmer
library(lsmeans) #cld
library(multcomp) #cld 2
library(multcompView) #for the CLD functions

#### Cumulative survival of queens in % ####

#mean survival time accross the treatments
tapply(dat_all$day_of_death[dat_all$y2_survival ==1], dat_all$treatment[dat_all$y2_survival==1], mean)

plot(survfit(Surv(dat_all$day_of_death, dat_all$y2_survival)~1)) #plots first basic plot with 95 confidence intervalls for all the data
fit=survfit(Surv(dat_all$day_of_death, dat_all$y2_survival)~dat_all$treatment)
plot(fit) 
summary(fit)

#plot cumulative queen survival
library(survminer)
ggsurvplot(fit, data = dat_all, pval = TRUE)
plot(fit, col=c(1:3), xlab="days since queens were placed in chambers", ylab="cumulative queen survival (%) ", lty = c(1,2,3), lwd = 2, yaxt="n", ylim=c(0.6, 1), cex.axis = 1.5, cex.lab = 1.5)
axis(2, at = c(0, 0.20, 0.40, 0.60, 0.80, 1), labels=c("0","20","40","60","80","100") , cex.axis = 1.5)
legend(12, 0.71, legend = c("0.0 ?g thiamethoxam/L", "4.5 ?g thiamethoxam/L", "30 ?g thiamethoxam/L"),  col=(1:3), lty= c(1,2,3) , lwd=2, cex = 1.5, bty = "n")

#calculate potential statistical difference between the groups
survdiff(Surv(day_of_death, y2_survival) ~ treatment, data = dat_all)
mod <- survdiff(Surv(day_of_death, y2_survival) ~ treatment, data=dat_all)
summary(mod)
mod



#### Number of Queens producing adults - Is there a sig. difference in the number of queens producing adults between the treatment? ####

dat_all$adult_production <- as.factor(dat_all$adult_production)
plot(dat_all$adult_production ~ dat_all$treatment)
chisq.test(dat_all$adult_production, dat_all$treatment)



#### Duration of different Lifestages ####

#define standard error function #3 not needed here... 
sem <- function(x){sd(x,na.rm=T)/sqrt(length(na.omit(x)))}

#eggs
plot(dat_all$duration_eggs ~ dat_all$treatment, main = "duration egg development", xlab= "treatments", ylab = "duration [days]", cex.axis = 1.5, cex.lab= 1.5, cex.main= 1.5)
require(dplyr)
group_by(dat_all, treatment) %>%
  summarise(
    count = n(),
    N = sum(!is.na(duration_pupae)),
    mean = mean(duration_eggs, na.rm = TRUE),
    sd = sd(duration_eggs, na.rm = TRUE)
  )
tapply(dat_all$duration_eggs, dat_all$treatment, quantile, na.rm = TRUE)
#statistical testing: 
shapiro.test(dat_all$duration_eggs) #shapiro ist signifikant --> Nullhypothese (Daten sind Normalverteilt muss verworfen werden)
leveneTest(dat_all$duration_eggs, dat_all$treatment) # not significant --> Nullhypthose f?r equal variance of the data can be assumed
kruskal.test(duration_eggs ~ treatment, data = dat_all, na.action=na.omit)


#larvae
plot(dat_all$duration_larva ~ dat_all$treatment, main = "duration larva development", xlab= "treatments", ylab = "duration [days]", cex.axis = 1.5, cex.lab= 1.5, cex.main= 1.5)
group_by(dat_all, treatment) %>%
  summarise(
    count = n(),
    N = sum(!is.na(duration_pupae)),
    mean = mean(duration_larva, na.rm = TRUE),
    sd = sd(duration_larva, na.rm = TRUE)
  )
tapply(dat_all$duration_larva, dat_all$treatment, quantile, na.rm = TRUE)
#statistical testing: 
shapiro.test(dat_all$duration_larva) #shapiro ist signifikant --> Nullhypothese (Daten sind Normalverteilt muss verworfen werden)
leveneTest(dat_all$duration_larva, dat_all$treatment) # not significant --> Nullhypthose f?r equal variance of the data can be assumed
kruskal.test(duration_larva ~ treatment, data = dat_all, na.action=na.omit)

#pupae
plot(dat_all$duration_pupae ~ dat_all$treatment, main = "duration pupae development", xlab= "treatments", ylab = "duration [days]", cex.axis = 1.5, cex.lab= 1.5, cex.main= 1.5)
group_by(dat_all, treatment) %>%
  summarise(
    count = n(),
    N = sum(!is.na(duration_pupae)),
    mean = mean(duration_pupae, na.rm = TRUE),
    sd = sd(duration_pupae, na.rm = TRUE)
  )
tapply(dat_all$duration_pupae, dat_all$treatment, quantile, na.rm = TRUE)
#statistical testing: 
shapiro.test(dat_all$duration_pupae) #shapiro ist nicht signifikant --> Nullhypothese (Daten sind Normalverteilt)  kann angenommen werden 
leveneTest(dat_all$duration_pupae, dat_all$treatment) # not significant --> Nullhypthose of equal variance of the data can be assumed
mod <- aov(duration_pupae ~ treatment, data = dat_all, na.action=na.omit)
summary(mod)

#complete ontogenesis
plot(dat_all$duration_adults ~ dat_all$treatment)
tapply(dat_all$duration_adults, dat_all$treatment, mean, na.rm = TRUE)
tapply(dat_all$duration_adults, dat_all$treatment, sd, na.rm = TRUE)
tapply(dat_all$duration_adults, dat_all$treatment, quantile, na.rm=TRUE)
#statistical testing: 
shapiro.test(dat_all$duration_adults)
leveneTest(dat_all$duration_adults, dat_all$treatment)
kruskal.test(duration_adults ~ treatment, data = dat_all, na.action=na.omit)

#preoviposition
plot(dat_all$duration_preoviposition ~ dat_all$treatment)
tapply(dat_all$duration_preoviposition, dat_all$treatment, mean, na.rm = TRUE)
tapply(dat_all$duration_preoviposition, dat_all$treatment, sd, na.rm = TRUE)
tapply(dat_all$duration_preoviposition, dat_all$treatment, quantile, na.rm = TRUE)

shapiro.test(dat_all$duration_preoviposition)
leveneTest(dat_all$duration_preoviposition, dat_all$treatment)
kruskal.test(duration_preoviposition ~ treatment, data = dat_all, na.action = na.omit)


#### Final Colony size after season 1 and season 2 ####
#### START - Final colony size year 1 + 2 (final counts) #####

# First graph of colonysize for y1 and y2 to get an idea on the data
boxplot(dat_all$y2_adults ~ dat_all$treatment)
boxplot(dat_all$y1_adults ~ dat_all$treatment)

### Y1 adults ###
mod <- lm(y1_adults ~ treatment, data = dat_all)
par(mfrow=c(2,2))
summary(mod)
anova(mod) # no significant differences
fit <- aov(y1_adults ~ treatment, data=dat_all) #same same but differet
summary(fit)
compareqqnorm(mod)
plot(mod)
acf(resid(mod))
#normality
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # ok
#heteroskedastity 
bptest(mod) #non significant --> ok
#homogenity of variance 
leveneTest(mod) #non-significant = ok 


# Descriptive statistics Y1 - Mean and SD values
tapply(dat_all$y1_adults, dat_all$treatment, mean, na.rm=T)
tapply(dat_all$y1_adults, dat_all$treatment, sd  , na.rm=T)
shapiro.test(dat_all$y1_adults)
tapply(dat_all$y1_adults, dat_all$treatment, quantile, na.rm=T)


#Y2 adults
mod <- lm(y2_adults ~ treatment, data=dat_all)
par(mfrow=c(2,2))
plot(mod)
acf(resid(mod))
summary(mod)
anova(mod) 
a1 <- aov(mod)
TukeyHSD(a1)

#Mean and SD values
tapply(dat_all$y2_adults, dat_all$treatment, mean, na.rm=T)
tapply(dat_all$y2_adults, dat_all$treatment, sd  , na.rm=T)
shapiro.test(dat_all$y2_adults) #signifikant --> reject null hypothesis of normal distribution --> quantiles and not sd
tapply(dat_all$y2_adults, dat_all$treatment, quantile, na.rm=T)


#### Graph Colonysize before the overwinterings #### 

#nicer boxplot (1200 x 880) / or as pdf (9x7 inches)
par(mfrow=c(1,2))
label <- c("Control", "Low", "High")
boxplot(dat_all$y1_adults ~ dat_all$treatment, main = "Week 13", xlab = "Treatments", ylab = "Number of workers", ylim=c(0, 26), las=1, xaxt = "n")
axis(1, at=1:3, labels = label)
text(c(1:3), 25, labels = c("a", "a", "a"), font = 2)
mtext('A', side=3, line=1.5, at=0, cex = 1.5, font = 2)
boxplot(dat_all$y2_adults ~ dat_all$treatment, main = "Week 64", xlab = "Treatments", ylab = "Number of workers", ylim=c(120, 420), las=1, xaxt = "n" )
axis(1, at=1:3, labels = label)
text(c(1:3), 410, labels = c("a", "b", "b"), font = 2)
mtext('B', side=3, line=1.5, at=0, cex = 1.5, font = 2)

#### Year 2 Bayesian Inference to see the effect size #### 
# is not in the manuscript and not updated with dat_all, so it might not work!
mod <- lm(y2_adults ~ treatment, data = dat_all)
mod
summary(mod)
anova(mod)
a1 <- aov(mod) #post hoc test to see pairwise differences
TukeyHSD(a1)
#check the model assumptions
par(mfrow=c(2,2))
plot(mod)
acf(resid(mod))
final$treat21 <- factor(paste(final$treatment_2, final$treatment_1))
final$treat21
plot(resid(mod)~final$treatment_1)

#bayesian inference using sim
modc <- mod
newdatc<-expand.grid(
  treatment=factor(c('control','low','high'), levels=levels(final$treatment)))
newdatc$fit <- predict(modc, newdata=newdatc)
nsim <- 2000
bsimc <- sim(modc, n.sims=2000)
str(bsimc)
fitmat <- matrix(ncol=nsim, nrow=nrow(newdatc))
Xmat <- model.matrix(formula(modc)[c(1,3)], data=newdatc)

for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsimc@coef[i,]
newdatc$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdatc$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdatc

#create a nice Graph (alternative to the boxplots)
par(mfrow=c(1,1))
col2 <- c("darkolivegreen1","darkorange3", "darkred")
boxplot((final$workforce) ~ final$treatment_1, border = "white", 
        main = "Effects on colonysize", ylab = "# of adults and pupae", xlab = "treatments", 
        family = c("mono"), xaxt="n")
x <- c("control","low","high") 
axis(1, at = c(1,2,3), labels=x)
stripchart(workforce ~ treatment_1, data = final, 
           vertical = TRUE, at = c(1,2,3) , method = "jitter", 
           pch = c(19,17,19,17,19,17), col = col2,
           add = TRUE) 
#plot posteriors and add legends
points(c(1,2,3)+0.15, newdatc$fit, pch = 19, cex=1.5)
segments(c(1,2,3)+0.15, sort(newdatc$lower, decreasing = TRUE), c(1,2,3)+0.15, sort(newdatc$upper, decreasing = TRUE), lwd = 3)
 
#Calculate differences in means between the 3 treatments
control_low <- fitmat[1,]-fitmat[2,]
quantile(control_low, prob = c(0.025, 0.5, 0.975))

control_high <- fitmat[1,]-fitmat[3,]
quantile(control_high, prob = c(0.025, 0.5, 0.975))

low_high <- fitmat[2,]-fitmat[3,]
quantile(low_high, prob = c(0.025, 0.5, 0.975))



#### Effects on colony developement y1) with bayesian inference (not used in paper) ####
head(final1)
boxplot(final1$adults ~ final1$treatment_1)
#create a simple linar model
modc <- lm(adults ~ treatment_1, data=final1)
par(mfrow=c(2,2))
plot(modc)
acf(resid(modc))
final1$treat21 <- factor(paste(final1$treatment_2, final1$treatment_1))
final1$treat21
plot(resid(modc)~final1$treatment_1)

newdatc<-expand.grid(
  treatment_1=factor(c('control','low','high'), levels=levels(final1$treatment_1)))
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
col2 <- c("darkolivegreen1","darkorange3", "darkred")
boxplot((final1$adults) ~ final1$treatment_1, border = "white", 
        main = "Effects on colonysize year 1", ylab = "# of adults and pupae", xlab = "treatments", 
        family = c("mono"), xaxt="n")
x <- c("control","low","high") 
axis(1, at = c(1,2,3), labels=x)
stripchart(adults ~ treatment_1, data = final1, 
           vertical = TRUE, at = c(1,2,3) , method = "jitter", 
           pch = c(19,17,19,17,19,17), col = col2,
           add = TRUE) 

#plot posteriors and add legends
points(c(1,2,3)+0.15, newdatc$fit, pch = 19, cex=1.5)
segments(c(1,2,3)+0.15, newdatc$lower, c(1,2,3)+0.15, newdatc$upper, lwd = 3)

#Calculate differences in means between the 3 treatments
control_low <- fitmat[1,]-fitmat[2,]
quantile(control_low, prob = c(0.025, 0.5, 0.975))
control_high <- fitmat[1,]-fitmat[3,]
quantile(control_high, prob = c(0.025, 0.5, 0.975))
low_high <- fitmat[2,]-fitmat[3,]
quantile(low_high, prob = c(0.025, 0.5, 0.975))


#### Colonydevelopement - Differences in available brood ####
# no longer in the mansuscript. Was used to see the course of the development over the year. 
# Plot the development graphs for year 1
#Eggs, Larva, Pupae, Adults Year I
sem<-function(x){sd(x,na.rm=T)/sqrt(length(na.omit(x)))}

pointcolor<-c("gray20","grey45","grey70")
linecolor<-c("grey25","grey60","grey80")
label <- c("Control", "Low", "High")
layout(matrix(1:4, 2, 2, byrow=T))
dev.off()

LnData <- subset(dat, year == 1)
LnData <- subset(LnData, treatment_2 == "control")
#get rid of day 43 which does not seem to have been counted
LnData <- subset(LnData, day != 43)
table(LnData$treatment_1)
summary(LnData$treatment_2)
LnData$treatment_1 <- relevel(LnData$treatment_1, ref ="low")
LnData$treatment_1 <- relevel(LnData$treatment_1, ref ="control")

op <- par(mfrow = c(2,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
par(mfrow = c(2,2))
par(mar=c(4, 4.1, 3.1, 1.1))

#Eggs #
semE_pergroup<-tapply(LnData$eggs,list(LnData$day,LnData$treatment_1),sem)
semE_vector<-c(as.numeric(semE_pergroup[,1]),as.numeric(semE_pergroup[,2]),as.numeric(semE_pergroup[,3]))
meanE_pergroup<-tapply(LnData$eggs,list(LnData$day,LnData$treatment_1),mean,na.rm=T)
day <- c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,18,19,20,22,23,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,41,44,45,46,47,49,50,53,54,56,58,60,62,64,66,68,70,72,74,76)
plot(day, meanE_pergroup[,1], pch=16, las = 1, xlab = "Time [d]", ylab = "Number of eggs", col = pointcolor[1],
     main = "Eggs", ylim=c(0,90), xlim = c(0, 76), cex.lab = 1.2, cex.axis = 1.2)
points(day, meanE_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanE_pergroup[,3],pch=16,col=pointcolor[3])
legend(60, 90, label, col = pointcolor, pch = 16, title = "Treatments", cex = 1.2)
meanE_pergroup
max(meanE_pergroup[,3])
for(i in 1:3){
  arrows(x0=day,x1=day,
         y0=meanE_pergroup[,i]-semE_pergroup[,i],y1=meanE_pergroup[,i]+semE_pergroup[,i],
         angle=90,length=0.025,code=3,col=linecolor[i])}
for(i in 1:3){lines(day, meanE_pergroup[,i], col = pointcolor[i])}
mtext('A', side=3, line=1, at=-10, cex = 1.5, font = 2)

#Larvae #
semL_pergroup<-tapply(LnData$larva,list(LnData$day,LnData$treatment_1),sem)
semL_vector<-c(as.numeric(semL_pergroup[,1]),as.numeric(semL_pergroup[,2]),as.numeric(semL_pergroup[,3]))
meanL_pergroup<-tapply(LnData$larva,list(LnData$day,LnData$treatment_1),mean,na.rm=T)
plot(day, meanL_pergroup[,1],pch=16, las = 1, xlab = "Time [d]", ylab = "Number of larvae", col = pointcolor[1],
     main = "Larvae", ylim=c(0,50), xlim = c(15,76), cex.lab = 1.2, cex.axis = 1.2)
points(day, meanL_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanL_pergroup[,3],pch=16,col=pointcolor[3])
legend(15, 50, label, col = pointcolor, pch = 16, title = "Treatments", cex = 1.2)
for(i in 1:3){arrows(x0=day,x1=day,
                     y0=meanL_pergroup[,i]-semL_pergroup[,i],y1=meanL_pergroup[,i]+semL_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}
for(i in 1:3){lines(day, meanL_pergroup[,i], col = pointcolor[i])}
mtext('B', side=3, line=1, at=6.5, cex = 1.5, font = 2)

#Pupae #
semP_pergroup<-tapply(LnData$pupae,list(LnData$day,LnData$treatment_1),sem)
semP_vector<-c(as.numeric(semP_pergroup[,1]),as.numeric(semP_pergroup[,2]),as.numeric(semP_pergroup[,3]))
meanP_pergroup<-tapply(LnData$pupae,list(LnData$day,LnData$treatment_1),mean,na.rm=T)
plot(day, meanP_pergroup[,1],pch=16, las = 1, xlab = "Time [d]", ylab = "Number of pupae", col = pointcolor[1],
     main = "Pupae", ylim=c(0,20), xlim = c(24,76), cex.lab = 1.2, cex.axis = 1.2)
points(day, meanP_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanP_pergroup[,3],pch=16,col=pointcolor[3])
legend(65, 20, label, col = pointcolor, pch = 16, title = "Treatments", cex = 1.2)
for(i in 1:3){arrows(x0=day,x1=day,
                     y0=meanP_pergroup[,i]-semP_pergroup[,i],y1=meanP_pergroup[,i]+semP_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}
for(i in 1:3){lines(day, meanP_pergroup[,i], col = pointcolor[i])}
mtext('C', side=3, line=1, at=17, cex = 1.5, font = 2)

#Adults #
semA_pergroup<-tapply(LnData$adults,list(LnData$day,LnData$treatment_1),sem)
semA_vector<-c(as.numeric(semA_pergroup[,1]),as.numeric(semA_pergroup[,2]),as.numeric(semA_pergroup[,3]))
meanA_pergroup<-tapply(LnData$adults,list(LnData$day,LnData$treatment_1),mean,na.rm=T)
plot(day, meanA_pergroup[,1],pch=16, las = 1, xlab = "Time [d]", ylab = "Number of workers", col = pointcolor[1],
     main = "Workers", ylim=c(0,18), xlim = c(40,76), cex.lab = 1.2, cex.axis = 1.2)
points(day, meanA_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanA_pergroup[,3],pch=16,col=pointcolor[3])
legend(40, 18, label, col = pointcolor, pch = 16, title = "Treatments", cex = 1.2)
for(i in 1:3){arrows(x0=day,x1=day,
                     y0=meanA_pergroup[,i]-semA_pergroup[,i],y1=meanA_pergroup[,i]+semA_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}
for(i in 1:3){lines(day, meanA_pergroup[,i], col = pointcolor[i])}
mtext('D', side=3, line=1, at=35, cex = 1.5, font = 2)

par(mar=c(5.1, 4.1, 4.1, 2.1))




#### peak day Analyses  ####
# also no longer in the manuscript because the during the season data was not shown. 
#peak day --> peakday for eggs is at day 20
sum(LnData$eggs[which(LnData$day == 18)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 19)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 20)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 22)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 23)], na.rm = TRUE)

# --> the peak day for eggs is measurementday nr. 20 
egg <- subset(LnData, day == 20)
head(egg)
shapiro.test(egg$eggs) # shapiro ist nicht signifikant --> Nullhypothese, dass daten Normal verteilt sind kann angenommen werden --> ANova oder LM
mod <- aov(eggs ~ treatment_1, data = egg, na.action=na.omit)
summary(mod)
plot(mod)
mod <- lm(eggs ~ treatment_1, data = egg)
summary(mod)

tapply(egg$eggs, egg$treatment_1, mean, na.rm=T)
tapply(egg$eggs, egg$treatment_1, sd, na.rm=T)
tapply(egg$eggs, egg$treatment_1, quantile, na.rm=T)
# there is no significant difference for eggs at their peak day in year 1

###diference at the end of Season one (Day 76)
egg <- subset(LnData, day == 76)
head(egg)
boxplot(egg$eggs ~ egg$treatment_1)
shapiro.test(egg$eggs) # shapiro ist  signifikant --> Nullhypothese, dass daten Normal verteilt muss verworfeb werden --> Kruskal wallis oder glm
leveneTest(egg$eggs, egg$treatment_1) # Nicht signifikant --> Nullhypothese von gleicher varianz kann angenommen werden
kruskal.test(egg$eggs ~ egg$treatment_1)

tapply(egg$eggs, egg$treatment_1, mean, na.rm=T)
tapply(egg$eggs, egg$treatment_1, sd, na.rm=T)
tapply(egg$eggs, egg$treatment_1, quantile, na.rm=T)



### Larva

sum(LnData$larva[which(LnData$day == 25)], na.rm = TRUE)
sum(LnData$larva[which(LnData$day == 26)], na.rm = TRUE)
sum(LnData$larva[which(LnData$day == 27)], na.rm = TRUE)
sum(LnData$larva[which(LnData$day == 28)], na.rm = TRUE)
sum(LnData$larva[which(LnData$day == 29)], na.rm = TRUE)

sum(LnData$larva[which(LnData$day == 70)], na.rm = TRUE)
sum(LnData$larva[which(LnData$day == 72)], na.rm = TRUE)
sum(LnData$larva[which(LnData$day == 74)], na.rm = TRUE)
sum(LnData$larva[which(LnData$day == 76)], na.rm = TRUE)


larva <- subset(LnData, day == 27)
larva2 <- subset(LnData, day == 76)
head(larva)
head(larva2)

shapiro.test(larva$larva) # shapiro ist  signifikant --> Nullhypothese, dass Daten Normal verteilt sind muss verworfen werden --> Kruskal wallis oder glm

shapiro.test(larva2$larva) # shapiro ist  signifikant --> Nullhypothese, dass Daten Normal verteilt sind muss verworfen werden --> Kruskal wallis oder glm
leveneTest(larva2$larva, larva2$treatment_1) # Nicht signifikant --> Nullhypothese von gleicher varianz kann angenommen werden
kruskal.test(larva ~ treatment_1, data = larva2)
plot(larva2$larva ~ larva2$treatment_1)

tapply(larva2$larva, larva2$treatment_1, mean, na.rm=T)
tapply(larva2$larva, larva2$treatment_1, sd, na.rm=T)
tapply(larva2$larva, larva2$treatment_1, quantile, na.rm=T)

#peak day for pupae
sum(LnData$pupae[which(LnData$day == 45)], na.rm = TRUE)
sum(LnData$pupae[which(LnData$day == 46)], na.rm = TRUE)
sum(LnData$pupae[which(LnData$day == 47)], na.rm = TRUE)
sum(LnData$pupae[which(LnData$day == 49)], na.rm = TRUE)
sum(LnData$pupae[which(LnData$day == 50)], na.rm = TRUE)

#peakday is day 47
pupae <- subset(LnData, day == 47)
shapiro.test(pupae$pupae) # shapiro ist  signifikant --> Nullhypothese, dass Daten Normal verteilt sind muss verworfen werden --> Kruskal wallis oder glm
leveneTest(pupae$pupae, pupae$treatment_1) # Nicht signifikant --> Nullhypothese von gleicher varianz kann angenommen werden

kruskal.test(pupae ~ treatment_1, data = pupae)
plot(pupae$pupae ~ pupae$treatment_1)

tapply(pupae$pupae, pupae$treatment_1, mean, na.rm=T)
tapply(pupae$pupae, pupae$treatment_1, sd, na.rm=T)
tapply(pupae$pupae, pupae$treatment_1, quantile, na.rm=T)

#pupae end of season 1
pupae <- subset(LnData, day == 76)
shapiro.test(pupae$pupae) # shapiro ist  signifikant --> Nullhypothese, dass Daten Normal verteilt sind muss verworfen werden --> Kruskal wallis oder glm
leveneTest(pupae$pupae, pupae$treatment_1) # Nicht signifikant --> Nullhypothese von gleicher varianz kann angenommen werden
kruskal.test(pupae ~ treatment_1, data = pupae)
plot(pupae$pupae ~ pupae$treatment_1)
tapply(pupae$pupae, pupae$treatment_1, quantile, na.rm=T)




#### final count larva, pupae eggs ####

#eggs
boxplot(dat_all$y2_eggs ~ dat_all$treatment)
shapiro.test(dat_all$y2_eggs)
leveneTest(dat_all$y2_eggs, dat_all$treatment)
kruskal.test(y2_eggs ~ treatment, data = dat_all)
tapply(dat_all$y2_eggs, dat_all$treatment, mean, na.rm=T)
tapply(dat_all$y2_eggs, dat_all$treatment, sd, na.rm=T)
tapply(final$eggs, final$treatment_1, quantile, na.rm=T)
tapply(dat_all$y2_eggs, dat_all$treatment, quantile, na.rm=TRUE)

#larva
boxplot(dat_all$y2_larva ~ dat_all$treatment)
shapiro.test(dat_all$y2_larva) # significant
leveneTest(dat_all$y2_larva, dat_all$treatment) #  signifikant --> Nullhypothese von gleicher varianz muss verworfen werden
kruskal.test(y2_larva ~ treatment, data = dat_all) #signifikant --> post hoc test f??r Kruskalwallis: Dunn-Bonnferroni
dunn.test(dat_all$y2_larva, dat_all$treament, method = "bonferroni")
x <- dat_all$y2_larva[!is.na(dat_all$y2_larva)]
g <- factor(rep(1:3, c(8,8,8)),
            labels = c("control", "low", "high"))
dunn.test(x, g, method = "bonferroni")
tapply(dat_all$y2_larva, dat_all$treatment, mean, na.rm=TRUE)
tapply(dat_all$y2_larva, dat_all$treatment, sd, na.rm=TRUE)
tapply(dat_all$y2_larva, dat_all$treatment, quantile, na.rm=TRUE)

#pupae
boxplot(dat_all$y2_pupae ~ dat_all$treatment)
shapiro.test(dat_all$y2_pupae) #significant
leveneTest(dat_all$y2_pupae, dat_all$treatment) #  non signifikant --> Nullhypothese von gleicher varianz kann angenommen werden
kruskal.test(y2_pupae ~ treatment, data = dat_all)
tapply(dat_all$y2_pupae, dat_all$treatment, mean, na.rm=TRUE)
tapply(dat_all$y2_pupae, dat_all$treatment, sd, na.rm = TRUE)
tapply(dat_all$y2_pupae, dat_all$treatment, quantile, na.rm=TRUE)


#Pupae at peak day # no longer in the manuscript
Y2 <- subset(ANV, year=="2")
Y2 <- subset(Y2, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")

#376,381,396,404,411,417,426,433,439,446,453
sum(Y2$pupae[which(Y2$days_since_collection == 376)], na.rm = TRUE)
sum(Y2$pupae[which(Y2$days_since_collection == 381)], na.rm = TRUE)
sum(Y2$pupae[which(Y2$days_since_collection == 396)], na.rm = TRUE)
sum(Y2$pupae[which(Y2$days_since_collection == 404)], na.rm = TRUE)
sum(Y2$pupae[which(Y2$days_since_collection == 411)], na.rm = TRUE)
sum(Y2$pupae[which(Y2$days_since_collection == 417)], na.rm = TRUE)
sum(Y2$pupae[which(Y2$days_since_collection == 426)], na.rm = TRUE)
sum(Y2$pupae[which(Y2$days_since_collection == 433)], na.rm = TRUE)

# Peakday for pupae is 396
pupae <- subset(Y2, days_since_collection == 396)
shapiro.test(pupae$pupae) # shapiro ist nicht signifikant --> Nullhypothese, dass Daten Normal verteilt sind kann angenommen werden
leveneTest(pupae$pupae, pupae$treatment_1) # Nicht signifikant --> Nullhypothese von gleicher varianz kann angenommen werden
plot(pupae$pupae ~ pupae$treatment_1)
mod <- aov(pupae ~ treatment_1, data = pupae, na.action=na.omit)
summary(mod)
plot(mod)

pupaec <- subset(pupae, treatment_1 == "control")
pupael <- subset(pupae, treatment_1 == "low")
pupaeh <- subset(pupae, treatment_1 == "high")

mean(pupaec$pupae, na.rm = TRUE)
sd(pupaec$pupae, na.rm = TRUE)
mean(pupael$pupae, na.rm = TRUE)
sd(pupael$pupae, na.rm = TRUE)
mean(pupaeh$pupae, na.rm = TRUE)
sd(pupaeh$pupae, na.rm = TRUE)
# pupae nirgens signigikant






#### weight of workers and queens ####
boxplot(dat_all$weight_queens ~ dat_all$treatment)
boxplot(dat_all$weight_adults ~ dat_all$treatment)

#workers
shapiro.test(dat_all$weight_adults)  # nicht signifikant und daher nehmen kann normalverteilung angenommen werden. 
leveneTest(dat_all$weight_adults, dat_all$treatment) # nicht signifikant --> Nullhypothese von gleicher Varianz kann angenommen werden. 
fligner.test(dat_all$weight_adults, dat_all$treatment) #Test for Homoskedacity, nicht signifikant --> Nullhypothese kann angenommen werden --> ANOVA 
mod <- lm(dat_all$weight_adults ~ dat_all$treatment)
summary(mod)
anova(mod) # signifikante differenz zwischen den drei gruppen
#pairwise post hoc tests
#bonferroni 
pairwise.t.test(dat_all$weight_adults, dat_all$treatment, p.adj = "bonf")
#Tukey HSD
a1 <- aov(mod)
TukeyHSD(a1)

#queens
shapiro.test(dat_all$weight_queen) # nicht signifikant und daher nehmen kann normalverteilung angenommen werden. 
leveneTest(dat_all$weight_queen, dat_all$treatment) # nicht signifikant --> Nullhypothese von gleicher Varianz kann angenommen werden. 
fligner.test(dat_all$weight_queen, dat_all$treatment) #Test for Homoskedacity, nicht signifikant --> Nullhypothese kann angenommen werden --> ANOVA
mod <- lm(dat_all$weight_queen ~ dat_all$treatment)
anova(mod) # keine signifikante differenz zwischen den drei gruppen
summary(mod)

#summary stats
require(dplyr)
group_by(dat_all, treatment) %>%
  summarise(
    count = n(),
    N_workers = sum(!is.na(weight_adults)),
    N_queens = sum(!is.na(weight_queens)),
    mean_workers = mean(weight_adults, na.rm = TRUE),
    sd_workers = sd(weight_adults, na.rm = TRUE),
    mean_queens = mean(weight_queens, na.rm = TRUE),
    sd_queens = sd(weight_queens, na.rm = TRUE)
  )


#### Graph Bodymass of Queens and Workers ####
#more elaborate boxplot #export as a 9x7 inch pdf
par(mar=c(4.1, 4.1, 3.1, 1.1))
par(oma = c(0,0,0,0)+0.5)

par(mfrow=c(1,2))
label <- c("Control", "Low", "High")
boxplot(dat_all$weight_queen ~ dat_all$treatment, main = "Queens", xlab = "Treatments", ylab="Weight [mg]", ylim=c(22, 47), las=1, xaxt = "n")
axis(1, at=1:3, labels = label)
text(c(1,2,3), 45.6, labels = c("a", "a", "a"), font = 2)
mtext('A', side=3, line=1.5, at=0, cex = 1.5, font = 2)
boxplot(dat_all$weight_adults ~dat_all$treatment, main = "Workers", xlab = "Treatments", ylab="Weight [mg]", ylim=c(15, 25), las=1, xaxt = "n")
axis(1, at=1:3, labels = label)
text(c(1,2,3), 24.5, labels = c("a", "ab", "b"), font = 2)
mtext('B', side=3, line=1.5, at=0, cex = 1.5, font = 2)
par(mfrow=c(1,1))




#### Neonicotinoid residues in ants ####

# new
dat_all$tw <- dat_all$tw_thiamethoxam_workers_dryconcentration
dat_all$cw <- dat_all$cw_clothianidin_workers_dryconcentration
dat_all$tq <- dat_all$tq_thiamethoxam_queens_dryconcentration
dat_all$cq <- dat_all$cq_clothianidin_queens_dryconcentration

# 1. Residue concentrations compared across the treatments and castes
#needed: treatment (1), tiam + cloth, for both castes, colony identity
dat_reduced <- subset(dat_all, y2_survival == 0) #take only the colonies that survived until the end of the experiment

table(dat_reduced$treatment) 
#create new dataframe with workers and queens after each other. 
caste_data <- data.frame(sample = c(1:24, 1:24), 
                         identity = rep(dat_reduced$identity, times = 2), 
                         treatment = rep(dat_reduced$treatment, times = 2),
                         caste = c(rep("worker", times = length(dat_reduced$identity)), rep("queen", times = length(dat_reduced$identity))),
                         thiamethoxam = c(dat_reduced$tw, dat_reduced$tq), 
                         clothianidin = c(dat_reduced$cw, dat_reduced$cq)
)

#dry weight makes most sense, dry weight, fresh weight and per abdomen --> look all similar, but per dry weight makes most sense
boxplot(caste_data$thiamethoxam ~ caste_data$caste*caste_data$treatment)
boxplot(caste_data$clothianidin ~ caste_data$caste*caste_data$treatment)

# descriprive stats
require(dplyr)
group_by(caste_data, treatment, caste) %>%
  summarise(
    mean_thiamethoxam = mean(thiamethoxam, na.rm = TRUE),
    sd_thiamethoxam = sd(thiamethoxam, na.rm = TRUE),
    mean_clothianidin = mean(clothianidin, na.rm = TRUE),
    sd_clothianidin = sd(clothianidin, na.rm = TRUE)
  )


#### absolute concentration (for gaetan, not in the manuscript) ####

dat$THX_conc_absolute <- dat$conc_thiamethoxam/dat$nr_abdomens
dat$cloth_conc_absolute <- dat$conc_clothinianidin/dat$nr_abdomens

boxplot(dat$THX_conc_absolute ~ dat$sample*dat$treatment_1)
boxplot(dat$cloth_conc_absolute ~ dat$sample*dat$treatment_1)

boxplot(dat$conc_T_abd ~ dat$sample*dat$treatment_1)
boxplot(dat$conc_C_abd ~ dat$sample*dat$treatment_1) 

tapply(dat$conc_T_abd[dat$sample=="worker"], dat$treatment_1[dat$sample=="worker"], mean)
tapply(dat$conc_T_abd[dat$sample=="worker"], dat$treatment_1[dat$sample=="worker"], sd)
tapply(dat$conc_T_abd[dat$sample=="queen"], dat$treatment_1[dat$sample=="queen"], mean)
tapply(dat$conc_T_abd[dat$sample=="queen"], dat$treatment_1[dat$sample=="queen"], sd)

tapply(dat$conc_C_abd[dat$sample=="worker"], dat$treatment_1[dat$sample=="worker"], mean)
tapply(dat$conc_C_abd[dat$sample=="worker"], dat$treatment_1[dat$sample=="worker"], sd)
tapply(dat$conc_C_abd[dat$sample=="queen"], dat$treatment_1[dat$sample=="queen"], mean)
tapply(dat$conc_C_abd[dat$sample=="queen"], dat$treatment_1[dat$sample=="queen"], sd)

#### Thiamethoxam modeling LMER ####

### drop the controls, as they all were zero
# The controls were all zero, which shows that there was no detactebla contamination prior to the experiment
# The difference between the treatments and the controls is measured with the treatment means of the model (even if controls are not in the model)

caste_reduced <- subset(caste_data, treatment != "control")
caste_reduced$treatment <- droplevels(caste_reduced)$treatment
table(caste_reduced$treatment)

#thiamethoxam
caste_reduced$T_rec <- (1/caste_reduced$thiamethoxam)

# baseline <- lmer(T_rec ~ (1|identity), data = caste_reduced, REML = FALSE)
baseline <- lmer(T_rec ~ (1|identity), data = caste_reduced, REML = FALSE)
sample_M <- update(baseline, .~. + caste)
treatment_M <- update(sample_M, .~. + treatment)
interaction_M <- update(treatment_M, .~. + caste:treatment)
anova(baseline, sample_M, treatment_M, interaction_M)
mod <- interaction_M
summary(mod)
estimates <- 1/fixef(mod)
estimates

leveneTest(residuals(mod) ~ dat_reduced$treatment_1*dat_reduced$sample)
boxplot(residuals(mod) ~ dat_reduced$treatment_1*dat_reduced$sample)

mod <- interaction_M
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
#normality
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> assumption of normally distributed residuals is ok

marginal <- lsmeans(mod,  ~ caste * treatment, adjust="tukey")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# mean and confidence interval limits can be backtransformed, but not standard deviations
# https://www.cienciasinseso.com/en/tag/reciprocal-transformation/



#### Clothianidin modeling LMER ####

# reciprocal transformation
x <- caste_reduced %>% filter(clothianidin != 0) %>% summarize(min = min(clothianidin))
caste_reduced$clothianidin[caste_reduced$clothianidin == 0] <- 0.5*x[1,1]
caste_reduced$C_rec <- (1/caste_reduced$clothianidin)

baseline <- lmer(C_rec ~ 1|identity, data = caste_reduced, REML = FALSE)
sample_M <- update(baseline, .~. + caste)
treatment_M <- update(sample_M, .~. + treatment)
interaction_M <- update(treatment_M, .~. + caste*treatment)
anova(baseline, sample_M, treatment_M, interaction_M)

model <- interaction_M
leveneTest(residuals(model) ~ caste_reduced$treatment*caste_reduced$caste)
boxplot(residuals(model) ~ dat_reduced$treatment_1*dat_reduced$sample)
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance
compareqqnorm(model)

marginal <- lsmeans(model,  ~ sample * treatment_1, adjust="tukey")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD




#### CLothianidin/Thiamethoxam ratio ####

caste_reduced <- subset(caste_data, treatment != "control")
caste_reduced$treatment <- droplevels(caste_reduced)$treatment
table(caste_reduced$treatment)

caste_ratio <- subset(caste_reduced, clothianidin > 0 & thiamethoxam > 0)
caste_ratio$ratio <- caste_ratio$clothianidin/caste_ratio$thiamethoxam

table(caste_ratio$caste, caste_ratio$treatment)
boxplot(caste_ratio$ratio ~ caste_ratio$caste)

caste_ratio$treatment <- droplevels(caste_ratio)$treatment
caste_ratio$identity <- droplevels(caste_ratio)$identity
levels(caste_ratio$identity)
summary(caste_ratio$identity)
x <- data.frame(ratio = caste_ratio$ratio, caste = caste_ratio$caste, identity = caste_ratio$identity)

baseline <- lmer(ratio ~ (1|identity), data = x)
model <- lmer(ratio ~ caste + (1|identity), data = x)
anova(baseline, model)

plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance
compareqqnorm(model)
leveneTest(residuals(model) ~ ratio$sample)
boxplot(residuals(model) ~ ratio$sample)


shapiro.test(x$ratio[x$caste == "queen"])
mean(x$ratio[x$caste == "queen"], na.rm=TRUE)
sd(x$ratio[x$caste == "queen"], na.rm=TRUE)
quantile(x$ratio[x$caste == "queen"])

shapiro.test(x$ratio[x$caste == "worker"])
mean(x$ratio[x$caste == "worker"], na.rm=TRUE)
sd(x$ratio[x$caste == "worker"], na.rm=TRUE)
quantile(x$ratio[x$caste == "worker"])



#### Graph Neonicotinoid residues for paper ####

# that one is not updated. 

par(mar=c(4.5,4,2,2))
par(mfrow=c(2,2))

# Thiamethoxam

boxplot(dat_reduced$conc_T_dry ~ dat_reduced$sample + dat_reduced$treatment_1, las=1,    
        ylab = "ng thiamethoxam / g dryweight", 
        xlab = "Treatments", 
        col  = c("white","grey"), 
        xaxt ="n",
        ylim = c(0.1, 1.6), 
        at = c(0.7, 1.7, 3.3, 4.3), 
        cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2, lwd = 1
)
label <- c("Low", "High")
axis(1, at = c(1.2,3.8), labels=label, cex.axis = 1.2)

legend(0.1, 1.42, legend = c("Queens", "Workers"),  c("white","grey"), bty = "n", cex=1.2)
text(c(0.7, 1.7, 3.3, 4.3), 1.52, labels = c("a", "a", "b", "b"),  cex = 1.2, font = 2)
mtext('A', side=3, line=0.5, at=-0.5, cex = 1.5, font = 2)


# Clothianidin
boxplot(dat_reduced$conc_C_dry ~ dat_reduced$sample + dat_reduced$treatment_1, las=1, 
        ylab = "ng clothianidin / g dryweight", 
        xlab = "Treatments", 
        col = c("white","grey"), 
        xaxt="n", 
        ylim = c(0, 4), 
        at = c(0.7, 1.7, 3.3, 4.3),
        cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2, lwd = 1
)

axis(1, at = c(1.2,3.8), labels=label, cex.axis = 1.2)
legend(0.1, 3.5, legend = c("Queens", "Workers"),  c("white","grey"), bty = "n", cex=1.2)
text(c(0.7, 1.7, 3.3, 4.3), 3.8, labels = c("a", "ab", "b", "c"),  cex = 1.2, font = 2)
mtext('B', side=3, line=0.5, at=-0.5, cex = 1.5, font = 2)

# ratio
boxplot(ratio$ratio ~ ratio$sample, las=1,
        ylab = "c[clothianidin] / c[thiamethoxam]", 
        xlab = "Caste",
        xaxt="n", 
        col = c("white","white"), 
        cex.lab = 1.2, cex.main = 1, cex.axis = 1.2, lwd = 1
)
axis(1, at = c(1,2), labels=c("Queens", "Workers"), cex.axis = 1.2)
text(c(1,2), 8, labels = c("a", "b"),  cex = 1.2, font = 2)
mtext('C', side=3, line=0.5, at=0.2, cex = 1.5, font = 2)



#### citations #### 
citation()
R.Version()
citation("survival") # contains kalan meier plot function
citation("arm") #contains the sim function for bayesian inference
citation("dunn.test") #contains the dunn test
citation("lme4") #lmer



#### weight differences

weight <- read.table("weight.txt", header = TRUE)
head(weight)

AN <- subset(weight, t2 == 1)
head(AN)

AN$fw_abdomen_w <- AN$fresh_weight_workers/AN$nr_Abdomens
AN$dw_addomen_w <- AN$dry_weight_workers/AN$nr_Abdomens

boxplot(AN$fw_abdomen_w~AN$t1)
boxplot(AN$dw_addomen_w~AN$t1)

boxplot(AN$fresh_weight_queens~AN$t1)
boxplot(AN$dry_weight_queens~AN$t1)


summary(aov(AN$fw_abdomen_w ~ AN$t1))
summary(aov(AN$dw_addomen_w ~ AN$t1))







#### Graphs in Manuscript ####
par(mar=c(4.1, 4.1, 3.1, 1.1),
    oma = c(0,0,0,0)+0.5,
    cex.lab=1.5, 
    cex.axis=1.5, 
    cex.main= 1.5,
    mfrow=c(1,2))
label <- c("Control", "Low", "High")

#Bodymass export as 10x8 inch pdfs 
boxplot(dat_all$weight_queens ~ dat_all$treatment, main = "Queens", xlab = "Treatments", ylab="Mass [mg]", ylim=c(22, 47), las=1, xaxt = "n")
axis(1, at=1:3, labels = label)
text(c(1,2,3), 45.6, labels = c("a", "a", "a"), font = 2, cex = 1.5)
mtext('a', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(dat_all$weight_adults ~dat_all$treatment, main = "Workers", xlab = "Treatments", ylab="Mass [mg]", ylim=c(15, 25), las=1, xaxt = "n")
axis(1, at=1:3, labels = label)
text(c(1,2,3), 24.5, labels = c("a", "ab", "b"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)


# Colony size
boxplot(dat_all$y1_adults ~ dat_all$treatment, main = "Week 13", xlab = "Treatments", ylab = "Number of workers", ylim=c(0, 26), las=1, xaxt = "n")
axis(1, at=1:3, labels = label)
text(c(1:3), 25, labels = c("a", "a", "a"), font = 2, cex = 1.5)
mtext('a', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(dat_all$y2_adults ~ dat_all$treatment, main = "Week 64", xlab = "Treatments", ylab = "", ylim=c(120, 420), las=1, xaxt = "n", yaxt = "n")
title(ylab="Number of workers",mgp=c(3,2,0))
axis(2, cex.axis = 1.25, las = 1)
axis(1, at=1:3, labels = label)
text(c(1:3), 410, labels = c("a", "b", "b"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)

# Neonic-Residues (export 8x8 inch pedf)
par(mfrow = c(2,2))
boxplot(caste_reduced$thiamethoxam ~ caste_reduced$caste + caste_reduced$treatment, las=1,    
        ylab = "ng thiamethoxam / g dryweight", 
        xlab = "Treatments", 
        col  = c("white","grey"), 
        xaxt ="n",
        ylim = c(0.1, 1.6), 
        at = c(0.7, 1.7, 3.3, 4.3), 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, lwd = 1)
label <- c("Low", "High")
axis(1, at = c(1.2,3.8), labels=label, cex.axis = 1.5)
text(c(0.7, 1.7, 3.3, 4.3), 1.52, labels = c("a", "a", "b", "b"),  cex = 1.5, font = 2)
mtext('a', side=3, line=0.5, at=-0.5, cex = 2, font = 2)

boxplot(caste_reduced$clothianidin ~ caste_reduced$caste + caste_reduced$treatment, las=1, 
        ylab = "ng clothianidin / g dryweight", 
        xlab = "Treatments", 
        col = c("white","grey"), 
        xaxt="n", 
        ylim = c(0, 4), 
        at = c(0.7, 1.7, 3.3, 4.3),
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, lwd = 1)
axis(1, at = c(1.2,3.8), labels=label, cex.axis = 1.5)
text(c(0.7, 1.7, 3.3, 4.3), 3.8, labels = c("a", "ab", "b", "c"),  cex = 1.5, font = 2)
mtext('b', side=3, line=0.5, at=-0.5, cex = 2, font = 2)

boxplot(caste_ratio$ratio ~ caste_ratio$caste, las=1,
        ylab = "c [clothianidin] / c [thiamethoxam]", 
        xlab = "Caste",
        xaxt="n", 
        col = c("white","white"), 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, lwd = 1)
axis(1, at = c(1,2), labels=c("Queens", "Workers"), cex.axis = 1.5)
text(c(1,2), 8, labels = c("a", "b"),  cex = 1.5, font = 2)
mtext('c', side=3, line=0.5, at=0.2, cex = 2, font = 2)