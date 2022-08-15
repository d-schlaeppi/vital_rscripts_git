
#########################################
###########    ANV for real   ###########
#########################################

#### Introduction ####

setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV")
dat <- read.table("ANVlong_full.txt", header = TRUE)

# load necessary librarys

library(lme4)
library(plyr)
library(lmerTest)
library(chron)
install.packages("arm")
library(arm)
library(car)
library(survival) # contains kalan meier plot function
library(nortest) #lilliefors test for normality
library(reshape2) #useful for data transformation
library(plyr) #to easily rename column names etc.
library(ggplot2) #not used in the end (I think)
library(car)
library(MASS)
library(multcomp)
library(tidyverse)
library(nlme)

if(!require(psych)){install.packages("psych")}
if(!require(nlme)){install.packages("nlme")}
if(!require(car)){install.packages("car")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(rcompanion)){install.packages("rcompanion")}

#Citing R Studio
RStudio.Version()

#möglichkeiten für die Analysen: Bis zum 16.07. mit allen Kolonien! Last time point, time series ohne last time point. 

#### Datenaufbereitung ####


# relevel treatment 1 to control low high 
levels(dat$treatment_1)
dat$treatment_1 <- relevel(dat$treatment_1, ref ="low")
dat$treatment_1 <- relevel(dat$treatment_1, ref ="control")

# Schneller look am 16.07 wo die treatments gesplitted wurden
cutoff <- subset(dat, date=="16.07.2017")
head(cutoff)
boxplot(cutoff$adults ~ cutoff$treatment_1)
#noch kein effekt auf koloniegrösse ersichtlich... 






#### Duration of different Lifestages ####

DLS <- read.table("lifestage_duration.txt", header = TRUE)
DLS <- subset(DLS, treatment_2 == "control")
DLS$treatment_1 <- relevel(DLS$treatment_1, ref ="low")
DLS$treatment_1 <- relevel(DLS$treatment_1, ref ="control")
head(DLS)

#define standard error function
sem <- function(x){sd(x,na.rm=T)/sqrt(length(na.omit(x)))}

#eggs
plot(DLS$duration_eggs ~ DLS$treatment_1, main = "duration egg development", xlab= "treatments", ylab = "duration [days]", cex.axis = 1.5, cex.lab= 1.5, cex.main= 1.5)
tapply(DLS$duration_eggs,DLS$treatment_1,mean,na.rm=T)
tapply(DLS$duration_eggs,DLS$treatment_1,sem)
#statistical testing: 
shapiro.test(DLS$duration_eggs) #shapiro ist signifikant --> Nullhypothese (Daten sind Normalverteilt muss verworfen werden)
levene.test(DLS$duration_eggs, DLS$treatment_1) # not significant --> Nullhypthose für equal variance of the data can be assumed
kruskal.test(duration_eggs ~ treatment_1, data = DLS, na.action=na.omit)  

#larvae
plot(DLS$duration_larva ~ DLS$treatment_1, main = "duration larva development", xlab= "treatments", ylab = "duration [days]", cex.axis = 1.5, cex.lab= 1.5, cex.main= 1.5)
tapply(DLS$duration_larva,DLS$treatment_1,mean,na.rm=T)
tapply(DLS$duration_larva,DLS$treatment_1,sem)
#statistical testing: 
shapiro.test(DLS$duration_larva) #shapiro ist signifikant --> Nullhypothese (Daten sind Normalverteilt muss verworfen werden)
levene.test(DLS$duration_larva, DLS$treatment_1) # not significant --> Nullhypthose für equal variance of the data can be assumed
kruskal.test(duration_larva ~ treatment_1, data = DLS, na.action=na.omit)

#pupae
plot(DLS$duration_pupae ~ DLS$treatment_1, main = "duration pupae development", xlab= "treatments", ylab = "duration [days]", cex.axis = 1.5, cex.lab= 1.5, cex.main= 1.5)
tapply(DLS$duration_pupae,DLS$treatment_1,mean,na.rm=T)
tapply(DLS$duration_pupae,DLS$treatment_1,sem)
#statistical testing: 
shapiro.test(DLS$duration_pupae) #shapiro ist nicht signifikant --> Nullhypothese (Daten sind Normalverteilt)  kann angenommen werden 
levene.test(DLS$duration_pupae, DLS$treatment_1) # not significant --> Nullhypthose für equal variance of the data can be assumed
mod <- aov(duration_pupae ~ treatment_1, data = DLS, na.action=na.omit)
summary(mod)

#complete ontogenesis
DLS$duration_adults <- DLS$Day_Adults-DLS$Day_Eggs 
plot(DLS$duration_adults ~ DLS$treatment_1)
tapply(DLS$duration_adults,DLS$treatment_1,mean,na.rm=T)
tapply(DLS$duration_adults,DLS$treatment_1,sem)
#statistical testing: 
shapiro.test(DLS$duration_adults) #shapiro ist  signifikant --> Nullhypothese (Daten sind Normalverteilt)  muss verworfen werden 
levene.test(DLS$duration_adults, DLS$treatment_1) # not significant --> Nullhypthose für equal variance of the data can be assumed
kruskal.test(duration_adults ~ treatment_1, data = DLS, na.action=na.omit)

#preoviposition
plot(DLS$Day_Eggs ~ DLS$treatment_1)
tapply(DLS$Day_Eggs,DLS$treatment_1,mean,na.rm=T)
tapply(DLS$Day_Eggs,DLS$treatment_1,sem)
shapiro.test(DLS$Day_Eggs) #shapiro ist  signifikant --> Nullhypothese (Daten sind Normalverteilt)  muss verworfen werden 
levene.test(DLS$Day_Eggs, DLS$treatment_1) # not significant --> Nullhypthose für equal variance of the data can be assumed
kruskal.test(Day_Eggs ~ treatment_1, data = DLS, na.action=na.omit)





#### Cumulative survival of queens in % ####
#kaplan meier oder so in %
survival <- read.table("ANV_queen-survival.txt", header = TRUE)
levels(survival$treatment_1)
survival$treatment_1 <- relevel(survival$treatment_1, ref ="low")
survival$treatment_1 <- relevel(survival$treatment_1, ref ="control")
levels(survival$treatment_1)
head(survival)

#Data ANV only
ANV <- subset(survival, treatment_2=="control")
table(ANV$treat)

#mean survival time accross the treatments
tapply(ANV$deathday_2[ANV$sy2==1], ANV$treatment_1[ANV$sy2==1], mean)
plot(survfit(Surv(ANV$deathday_2, ANV$sy2)~1)) #plots first basic plot with 95 confidence intervalls for all the data
fit=survfit(Surv(ANV$deathday_2, ANV$sy2)~ANV$treatment_1)
plot(fit) 

#plot cumulative queen survival
plot(fit, col=c(1:3), xlab="days since queens were placed in chambers", ylab="cumulative queen survival (%) ", lty = c(1,2,3), lwd = 2, yaxt="n", ylim=c(0.6, 1), cex.axis = 1.5, cex.lab = 1.5)
axis(2, at = c(0, 0.20, 0.40, 0.60, 0.80, 1), labels=c("0","20","40","60","80","100") , cex.axis = 1.5)
legend(12, 0.71, legend = c("0.0 µg thiamethoxam/L", "4.5 µg thiamethoxam/L", "30 µg thiamethoxam/L"),  col=(1:3), lty= c(1,2,3) , lwd=2, cex = 1.5, bty = "n")

#calculate potential statistical difference between the groups
mod <- survdiff(Surv(deathday_2, sy2) ~ treatment_1, data=ANV)
summary(mod)
mod



#### Number of Queens producing adults - Is there a sig. difference in the number of queens producing adults between the treatment? ####
adult_production <- read.table("adult_production.txt", header = TRUE)
a_prod <- subset(adult_production, treatment_2 == "control")
head(a_prod)

a_prod$treatment_1 <- relevel(a_prod$treatment_1, ref ="low")
a_prod$treatment_1 <- relevel(a_prod$treatment_1, ref ="control")
plot(a_prod$adult_production~a_prod$treatment_1)

m <- matrix( c(9,10,8,1,0,2), byrow =  T, nrow = 2)
chisq.test(m)

chisq.test(a_prod$adult_production, a_prod$treatment_1)




#### COLONY DEVELOPEMENT ####


#calculate days of analysis
# eggs year 1
# larva year 1
# pupae year 1
# pupae y2  
# see lifestage developement over time?






#### Time series analysis, IM MOMENT niccht die Analyse Variante #### 

ANV <- subset(dat, treatment_2=="control") #nur Kolonien ohne Virustreatment
library(nlme); ?corClasses
head(ANV)


#year1 + 2
model <- glmmPQL(adults ~ treatment_1, random=~1|id, data=ANV, family = poisson, correlation=corCAR1(form = ~measurement_y2|id))  # compare model performance (if they fit to modelassumptions) for both lmer with normal distributed data and for poisson, as we clearly have count data). 
summary(model)
model <- lme(adults ~ treatment_1, random=~1|id, correlation=corCAR1(form = ~measurement_y2|id), data=ANV, method="REML", na.action = na.exclude) #na.action = na.exclude
summary(model)

#only year 2 
model <- glmmPQL(adults ~ treatment_1, random=~1|id, data=year1, family = poisson, correlation=corCAR1(form = ~measurement_y2|id))  # compare model performance (if they fit to modelassumptions) for both lmer with normal distributed data and for poisson, as we clearly have count data). 
summary(model)
model <- lme(adults ~ treatment_1, random=~1|id, correlation=corCAR1(form = ~measurement_y2|id), data=year1, method="REML", na.action = na.exclude) #na.action = na.exclude
summary(model)

#time series analysis do not seem to turn out significant (probably to low sample size with an overlap at the beginning...)

last2months <- subset(ANV, measurement_y2 > 11)
model <- glmmPQL(adults ~ treatment_1, random=~1|id, data=last2months, family = poisson, correlation=corCAR1(form = ~measurement_y2|id))  # compare model performance (if they fit to modelassumptions) for both lmer with normal distributed data and for poisson, as we clearly have count data). 
summary(model)
model <- lme(adults ~ treatment_1, random=~1|id, correlation=corCAR1(form = ~measurement_y2|id), data=last2months, method="REML", na.action = na.exclude) #na.action = na.exclude
summary(model)


last2months <- subset(ANV, measurement_y2 > 12)
model <- glmmPQL(log10(adults) ~ treatment_1, random=~1|id, data=last2months, family = poisson, correlation=corCAR1(form = ~measurement_y2|id))  # compare model performance (if they fit to modelassumptions) for both lmer with normal distributed data and for poisson, as we clearly have count data). 
summary(model)
model <- lme(log10(adults) ~ treatment_1, random=~1|id, correlation=corCAR1(form = ~measurement_y2|id), data=last2months, method="REML", na.action = na.exclude) #na.action = na.exclude
summary(model)

plot(model)

par(mfrow=c(2,2))
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance

shapiro.test(log10(last2months$adults)) #data not normaly distributed 
shapiro.test(log10(ANV$adults))




##### START - Final colony size year 1 + 2 (real counts) #####

# create new subsets of dat --> exclude the virus treated groups and only look at the final counts for y1 and y2
ANV <- subset(dat, treatment_2=="control")
final <- subset(ANV, final=="year2")
final1 <- subset(ANV, final=="year1")
head(final)
summary(final$treatment_1)
summary(final$treatment_2)

#exclude failure queens for the two years 
final1 <- subset(final1, id!= "C9C" & id!= "H4C" & id!="H11V")
final <- subset(final, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")

summary(final$treatment_1)
summary(final$treatment_2)

final$treatment_1 <- as.factor(final$treatment_1)
final$treatment_1 <- factor(final$treatment_1, levels = c("control", "low", "high"))
final$treatment_2 <- as.factor(final$treatment_2)


final$workforce <- final$pupae+final$adults # create a new variable the looks at the total workforce available to the colony in a short time --> will mostlikely be neglected later, but might be interesting
final$workforce
mean(final$workforce)
boxplot(pupae ~ treatment_1, data=final)
boxplot(adults ~ treatment_1, data=final)
boxplot(workforce ~ treatment_1, data=final) # not really different and thus neglected (it is at the end of the season and thus not a lot of pupae are in the colonies anyways)

# First graph of colonysize for y1 and y2
boxplot(final$adults ~ final$treatment_1)
boxplot(final1$adults ~ final1$treatment_1)




#### weight of workers and queens ####
boxplot(final$weight_queen ~ final$treatment_1)
boxplot(final$weight_20adults ~final$treatment_1)

#workers
shapiro.test(final$weight_20adults)  # nicht signifikant und daher nehmen kann normalverteilung angenommen werden. 
leveneTest(final$weight_20adults, final$treatment_1) # nicht signifikant --> Nullhypothese von gleicher Varianz kann angenommen werden. 
fligner.test(final$weight_20adults, final$treatment_1) #Test for Homoskedacity, nicht signifikant --> Nullhypothese kann angenommen werden --> ANOVA 
mod <- lm(final$weight_20adults ~ final$treatment_1)
summary(mod)
anova(mod) # signifikante differenz zwischen den drei gruppen
#pairwise post hoc tests
#bonferroni 
pairwise.t.test(final$weight_20adults, final$treatment_1, p.adj = "bonf")
#Tukey HSD
a1 <- aov(mod)
TukeyHSD(a1)


#queens
shapiro.test(final$weight_queen)  # nicht signifikant und daher nehmen kann normalverteilung angenommen werden. 
leveneTest(final$weight_queen, final$treatment_1) # nicht signifikant --> Nullhypothese von gleicher Varianz kann angenommen werden. 
fligner.test(final$weight_queen, final$treatment_1) #Test for Homoskedacity, nicht signifikant --> Nullhypothese kann angenommen werden --> ANOVA 
mod <- lm(final$weight_queen ~ final$treatment_1)
anova(mod) # keine signifikante differenz zwischen den drei gruppen
summary(mod)


#pairwise post hoc tests
#bonferroni 
pairwise.t.test(final$weight_queen, final$treatment_1, p.adj = "bonf")
#Tukey HSD
a1 <- aov(mod)
TukeyHSD(a1)


control <- subset(final, treatment_1 == "control")
low     <- subset(final, treatment_1 == "low")
high    <- subset(final, treatment_1 == "high")
control1 <- subset(final1, treatment_1 == "control")
low1     <- subset(final1, treatment_1 == "low")
high1    <- subset(final1, treatment_1 == "high")



mean(final$weight_queen)
sd(final$weight_queen, na.rm = TRUE)

mean(control$weight_queen)
sd(control$weight_queen, na.rm = TRUE)
mean(low$weight_queen)
sd(low$weight_queen, na.rm = TRUE)
mean(high$weight_queen)
sd(high$weight_queen, na.rm = TRUE)

mean(control$weight_20adults)
sd(control$weight_20adults, na.rm = TRUE)
mean(low$weight_20adults)
sd(low$weight_20adults, na.rm = TRUE)
mean(high$weight_20adults)
sd(high$weight_20adults, na.rm = TRUE)




#more elaborate boxplot

par(mfrow=c(1,2))
boxplot(final$weight_queen ~final$treatment_1, main = "queen weight", xlab = "treatments", ylab="weight [mg]", ylim=c(22, 47) )
text(c(1,2,3), 45.6, labels = c("a", "a", "a"), font = 2)
mtext('A', side=3, line=1.5, at=0, cex = 1.5, font = 2)
boxplot(final$weight_20adults ~final$treatment_1, main = "weight of 20 workers", xlab = "treatments", ylab="weight [mg]", ylim=c(15, 25) )
text(c(1,2,3), 24.5, labels = c("a", "ab", "b"), font = 2)
mtext('B', side=3, line=1.5, at=0, cex = 1.5, font = 2)




#### Year I + Year 2 Final Count Adults ####

boxplot(final1$adults ~ final1$treatment_1)
boxplot(final$adults ~ final$treatment_1)

#Y1 adults
mod <- lm(adults ~ treatment_1, data=final1)
par(mfrow=c(2,2))
plot(mod)      
acf(resid(mod))

summary(mod)
anova(mod) # keine signifikante differenz zwischen den drei gruppen --> auch kein post hoc test der differenz nötig. 

#Mean and SD values
mean(control1$adults)
sd(control1$adults, na.rm = TRUE)
mean(low1$adults)
sd(low1$adults, na.rm = TRUE)
mean(high1$adults)
sd(high1$adults, na.rm = TRUE)


#Y2 adults
mod <- lm(adults ~ treatment_1, data=final)
par(mfrow=c(2,2))
plot(mod)      
acf(resid(mod))

summary(mod)
anova(mod) 

a1 <- aov(mod)
TukeyHSD(a1)

#Mean and SD values
mean(control$adults)
sd(control$adults, na.rm = TRUE)
mean(low$adults)
sd(low$adults, na.rm = TRUE)
mean(high$adults)
sd(high$adults, na.rm = TRUE)

#nicer boxplot
par(mfrow=c(1,2))
boxplot(final1$adults ~ final1$treatment_1, main = "number of workers year 1", xlab = "treatments", ylab = "nr. of workers", ylim=c(0, 26))
text(c(1,2,3), 25, labels = c("a", "a", "a"), font = 2)
mtext('A', side=3, line=1.5, at=0, cex = 1.5, font = 2)
boxplot(final$adults ~ final$treatment_1, main = "number of workers year 2", xlab = "treatments", ylab = "nr. of workers", ylim=c(120, 420) )
text(c(1,2,3), 410, labels = c("a", "b", "b"), font = 2)
mtext('B', side=3, line=1.5, at=0, cex = 1.5, font = 2)

#### Year 2 Bayesian Inference to see the effect size #### 

mod <- lm(adults ~ treatment_1, data = final)
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
  treatment_1=factor(c('control','low','high'), levels=levels(final$treatment_1)))
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

#create a super nice Graph: 
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
#aufgrund von releveling etz werden die segmente falsch im graphen geprintet --> wir ordenen den vektor neu an. 

#Calculate differences in means between the 3 treatments
control_low <- fitmat[1,]-fitmat[2,]
quantile(control_low, prob = c(0.025, 0.5, 0.975))

control_high <- fitmat[1,]-fitmat[3,]
quantile(control_high, prob = c(0.025, 0.5, 0.975))

low_high <- fitmat[2,]-fitmat[3,]
quantile(low_high, prob = c(0.025, 0.5, 0.975))









#### Effects on colony developement y1) ####

head(final1)
boxplot(final1$adults ~ final1$treatment_1)

summary(final1$treatment_1)
summary(final1$treatment_2)

#exclude failure queens
final1 <- subset(final1, id!= "C9C" &  id!= "H4C" & id!="H11V")

final1$treatment_1 <- as.factor(final1$treatment_1)
final1$treatment_1 <- factor(final1$treatment_1, levels = c("control", "low", "high"))
final1$treatment_2 <- as.factor(final1$treatment_2)

final1$workforce <- final1$pupae+final1$adults
final1$workforce
mean(final1$workforce)
boxplot(workforce ~ treatment_1, data=final1)
boxplot(adults ~ treatment_1, data=final1)

#create a simple linar model

modc <- lm(workforce ~ treatment_1, data=final1)

modc
summary(modc)
Anova(modc)

final1$treatment_1 <- relevel(final1$treatment_1, ref ="low")
modc <- lm(adults ~ treatment_1, data=final1)
summary(modc)


final1$treatment_1 <- relevel(final1$treatment_1, ref ="control")

#hist(datc$workforce, breaks = 12)
#shapiro.test(datc$workforce)

#check the model assumptions
par(mfrow=c(2,2))
plot(modc)      
acf(resid(modc))
final1$treat21 <- factor(paste(final1$treatment_2, final1$treatment_1))
final1$treat21
plot(resid(modc)~final1$treatment_1)

#bayesian inference using sim
# add predicted values (group means)

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
boxplot((final1$workforce) ~ final1$treatment_1, border = "white", 
        main = "Effects on colonysize year 1", ylab = "# of adults and pupae", xlab = "treatments", 
        family = c("mono"), xaxt="n")

x <- c("control","low","high") 
axis(1, at = c(1,2,3), labels=x)

stripchart(workforce ~ treatment_1, data = final1, 
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











#### Plot Growth curves Year 1####

#get rid of day 43 which does not seem to have been counted
ANV <- subset(ANV, day != 43)

control <- subset(ANV, treatment_1 == "high" & year==1)
low<- subset(ANV, treatment_1 == "low" & year==1)
high <- subset(ANV, treatment_1 == "control" & year==1)



#a simple way to plot some of the things with boxplots
#boxplot(high$eggs~high$day, xlab = "Days", ylab = "Eggs", main = "Eggs per Day (Treatment 'High')", las = 1)
#however: confidence intervalls are nicer

#calculate confidence intervals for y1: 

#High

q1HE<-tapply(high$eggs, high$day, quantile, probs = 0.25, na.rm = T)
q2HE<-tapply(high$eggs, high$day, quantile, probs = 0.5, na.rm = T)
q3HE<-tapply(high$eggs, high$day, quantile, probs = 0.75, na.rm = T)

q1HL<-tapply(high$larva, high$day, quantile, probs = 0.25, na.rm = T)
q2HL<-tapply(high$larva, high$day, quantile, probs = 0.5, na.rm = T)
q3HL<-tapply(high$larva, high$day, quantile, probs = 0.75, na.rm = T)

q1HP<-tapply(high$pupae, high$day, quantile, probs = 0.25, na.rm = T)
q2HP<-tapply(high$pupae, high$day, quantile, probs = 0.5, na.rm = T)
q3HP<-tapply(high$pupae, high$day, quantile, probs = 0.75, na.rm = T)

q1HA<-tapply(high$adults, high$day, quantile, probs = 0.25, na.rm = T)
q2HA<-tapply(high$adults, high$day, quantile, probs = 0.5, na.rm = T)
q3HA<-tapply(high$adults, high$day, quantile, probs = 0.75, na.rm = T)


#Low
q1LE<-tapply(low$eggs, low$day, quantile, probs = 0.25, na.rm = T)
q2LE<-tapply(low$eggs, low$day, quantile, probs = 0.5, na.rm = T)
q3LE<-tapply(low$eggs, low$day, quantile, probs = 0.75, na.rm = T)

q1LL<-tapply(low$larva, low$day, quantile, probs = 0.25, na.rm = T)
q2LL<-tapply(low$larva, low$day, quantile, probs = 0.5, na.rm = T)
q3LL<-tapply(low$larva, low$day, quantile, probs = 0.75, na.rm = T)

q1LP<-tapply(low$pupae, low$day, quantile, probs = 0.25, na.rm = T)
q2LP<-tapply(low$pupae, low$day, quantile, probs = 0.5, na.rm = T)
q3LP<-tapply(low$pupae, low$day, quantile, probs = 0.75, na.rm = T)

q1LA<-tapply(low$adults, low$day, quantile, probs = 0.25, na.rm = T)
q2LA<-tapply(low$adults, low$day, quantile, probs = 0.5, na.rm = T)
q3LA<-tapply(low$adults, low$day, quantile, probs = 0.75, na.rm = T)


#Control
q1CE<-tapply(control$eggs, control$day, quantile, probs = 0.25, na.rm = T)
q2CE<-tapply(control$eggs, control$day, quantile, probs = 0.5, na.rm = T)
q3CE<-tapply(control$eggs, control$day, quantile, probs = 0.75, na.rm = T)

q1CL<-tapply(control$larva, control$day, quantile, probs = 0.25, na.rm = T)
q2CL<-tapply(control$larva, control$day, quantile, probs = 0.5, na.rm = T)
q3CL<-tapply(control$larva, control$day, quantile, probs = 0.75, na.rm = T)

q1CP<-tapply(control$pupae, control$day, quantile, probs = 0.25, na.rm = T)
q2CP<-tapply(control$pupae, control$day, quantile, probs = 0.5, na.rm = T)
q3CP<-tapply(control$pupae, control$day, quantile, probs = 0.75, na.rm = T)

q1CA<-tapply(control$adults, control$day, quantile, probs = 0.25, na.rm = T)
q2CA<-tapply(control$adults, control$day, quantile, probs = 0.5, na.rm = T)
q3CA<-tapply(control$adults, control$day, quantile, probs = 0.75, na.rm = T)


#Creating necessary subsets & shortcuts
LnC_ls <- subset(control, control$day >= 18 & control$day <= 25)
LnL_ls <- subset(low, low$day >= 18 & low$day <= 25)
LnH_ls <- subset(high, high$day >= 18 & high$day <= 25)

CEq1 <- tapply(LnC_ls$eggs, LnC_ls$day, quantile, probs = 0.25, na.rm = T)
CEq2 <- tapply(LnC_ls$eggs, LnC_ls$day, quantile, probs = 0.5, na.rm = T)
CEq3 <- tapply(LnC_ls$eggs, LnC_ls$day, quantile, probs = 0.75, na.rm = T)

LEq1 <- tapply(LnL_ls$eggs, LnL_ls$day, quantile, probs = 0.25, na.rm = T)
LEq2 <- tapply(LnL_ls$eggs, LnL_ls$day, quantile, probs = 0.5, na.rm = T)
LEq3 <- tapply(LnL_ls$eggs, LnL_ls$day, quantile, probs = 0.75, na.rm = T)

HEq1 <- tapply(LnH_ls$eggs, LnH_ls$day, quantile, probs = 0.25, na.rm = T)
HEq2 <- tapply(LnH_ls$eggs, LnH_ls$day, quantile, probs = 0.5, na.rm = T)
HEq3 <- tapply(LnH_ls$eggs, LnH_ls$day, quantile, probs = 0.75, na.rm = T)



#### TIME SERIES PLOTt Year 1####
#Eggs, Larva, Pupae, Adults Year I

sem<-function(x){sd(x,na.rm=T)/sqrt(length(na.omit(x)))}
pointcolor<-c("green4","brown3","orange2")
pointcolor<-c("gray20","grey45","grey70")
linecolor<-c("grey25","grey60","grey80")
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

par(mar=c(4.1, 4.1, 3.1, 1.1))


#Eggs #
semE_pergroup<-tapply(LnData$eggs,list(LnData$day,LnData$treatment_1),sem)
semE_vector<-c(as.numeric(semE_pergroup[,1]),as.numeric(semE_pergroup[,2]),as.numeric(semE_pergroup[,3]))
meanE_pergroup<-tapply(LnData$eggs,list(LnData$day,LnData$treatment_1),mean,na.rm=T)

day <- c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,18,19,20,22,23,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,41,44,45,46,47,49,50,53,54,56,58,60,62,64,66,68,70,72,74,76)

plot(day, meanE_pergroup[,1], pch=16, las = 1, xlab = "" , ylab = "# eggs", col = pointcolor[1],
     main = "Number of eggs over time ", ylim=c(0,90), xlim = c(0, 76))
points(day, meanE_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanE_pergroup[,3],pch=16,col=pointcolor[3])
legend(60, 90, c("control", "low", "high"), col = pointcolor, pch = 16, title = "treatments")

meanE_pergroup
max(meanE_pergroup[,3])

for(i in 1:3){
  arrows(x0=day,x1=day,
         y0=meanE_pergroup[,i]-semE_pergroup[,i],y1=meanE_pergroup[,i]+semE_pergroup[,i],
         angle=90,length=0.025,code=3,col=linecolor[i])}

for(i in 1:3){lines(day, meanE_pergroup[,i], col = pointcolor[i])}

mtext('A', side=3, line=1, at=0, cex = 1.5, font = 2)



#Larvae #
semL_pergroup<-tapply(LnData$larva,list(LnData$day,LnData$treatment_1),sem)
semL_vector<-c(as.numeric(semL_pergroup[,1]),as.numeric(semL_pergroup[,2]),as.numeric(semL_pergroup[,3]))
meanL_pergroup<-tapply(LnData$larva,list(LnData$day,LnData$treatment_1),mean,na.rm=T)

plot(day, meanL_pergroup[,1],pch=16, las = 1,xlab = "", ylab = "# larvae", col = pointcolor[1],
     main = "Number of larva over time", ylim=c(0,50), xlim = c(15,76))
points(day, meanL_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanL_pergroup[,3],pch=16,col=pointcolor[3])
legend(15, 50, c("control", "low", "high"), col = pointcolor, pch = 16, title = "treatments")

for(i in 1:3){arrows(x0=day,x1=day,
                     y0=meanL_pergroup[,i]-semL_pergroup[,i],y1=meanL_pergroup[,i]+semL_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}

for(i in 1:3){lines(day, meanL_pergroup[,i], col = pointcolor[i])}

mtext('B', side=3, line=1, at=15, cex = 1.5, font = 2)

#Pupae #
semP_pergroup<-tapply(LnData$pupae,list(LnData$day,LnData$treatment_1),sem)
semP_vector<-c(as.numeric(semP_pergroup[,1]),as.numeric(semP_pergroup[,2]),as.numeric(semP_pergroup[,3]))
meanP_pergroup<-tapply(LnData$pupae,list(LnData$day,LnData$treatment_1),mean,na.rm=T)

plot(day, meanP_pergroup[,1],pch=16, las = 1, xlab = "time [d]", ylab = "# pupae", col = pointcolor[1],
     main = "Number of pupae over time", ylim=c(0,20), xlim = c(24,76))
points(day, meanP_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanP_pergroup[,3],pch=16,col=pointcolor[3])
legend(65, 20, c("control", "low", "high"), col = pointcolor, pch = 16, title = "treatments")

for(i in 1:3){arrows(x0=day,x1=day,
                     y0=meanP_pergroup[,i]-semP_pergroup[,i],y1=meanP_pergroup[,i]+semP_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}

for(i in 1:3){lines(day, meanP_pergroup[,i], col = pointcolor[i])}

mtext('C', side=3, line=1, at=24, cex = 1.5, font = 2)

#Adults #
semA_pergroup<-tapply(LnData$adults,list(LnData$day,LnData$treatment_1),sem)
semA_vector<-c(as.numeric(semA_pergroup[,1]),as.numeric(semA_pergroup[,2]),as.numeric(semA_pergroup[,3]))
meanA_pergroup<-tapply(LnData$adults,list(LnData$day,LnData$treatment_1),mean,na.rm=T)

plot(day, meanA_pergroup[,1],pch=16, las = 1, xlab = "time [d]", ylab = "# workers", col = pointcolor[1],
     main = "Number of workers over time", ylim=c(0,18), xlim = c(40,76))
points(day, meanA_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanA_pergroup[,3],pch=16,col=pointcolor[3])
legend(40, 18, c("control", "high", "low"), col = pointcolor, pch = 16, title = "treatments")

for(i in 1:3){arrows(x0=day,x1=day,
                     y0=meanA_pergroup[,i]-semA_pergroup[,i],y1=meanA_pergroup[,i]+semA_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}

for(i in 1:3){lines(day, meanA_pergroup[,i], col = pointcolor[i])}

mtext('D', side=3, line=1, at=40, cex = 1.5, font = 2)

par(mar=c(5.1, 4.1, 4.1, 2.1))


#### TIME SERIES PLOTt Year 2####
# Pupae, Adults Year 2


#### year 2 ####
head(ANV)
Y2 <- subset(ANV, year=="2")
head(Y2)

control <- subset(ANV, treatment_1 == "high" & year==2)
low<- subset(ANV, treatment_1 == "low" & year==2)
high <- subset(ANV, treatment_1 == "control" & year==2)

#Adults #

# semA_pergroup<-tapply(Y2$adults,list(Y2$days_since_collection,Y2$treatment_1),sem)
# semA_vector<-c(as.numeric(semA_pergroup[,1]),as.numeric(semA_pergroup[,2]),as.numeric(semA_pergroup[,3]))
# meanA_pergroup<-tapply(Y2$adults,list(Y2$days_since_collection,Y2$treatment_1),mean,na.rm=T)
# day <- c(271,285,314,321,328,335,343,352,358,366,376,381,396,404,411,417,426,433,439,446,453,458,467,473,474)

# plot(day, meanA_pergroup[,1],pch=17, las = 1, xlab = "#days since collection", ylab = "#adults", ylim=c(0,350), xlim = c(min(day), max(day)), 
#      col = pointcolor[1], main = "number of adults over time")
# points(day, meanA_pergroup[,2],pch=17, col = pointcolor[2])
# points(day, meanA_pergroup[,3],pch=17,col=pointcolor[3])
# legend(270, 350, c("control", "high", "low"), col = pointcolor, pch = 17, title = "treatments")
# 
# for(i in 1:3){arrows(x0=day, x1=day,
#                      y0=meanA_pergroup[,i]-semA_pergroup[,i],y1=meanA_pergroup[,i]+semA_pergroup[,i],
#                      angle=90,length=0.05,code=3,col=linecolor[i])}
# 
# for(i in 1:3){lines(day, meanA_pergroup[,i], col = pointcolor[i])}
#with excluded final count 
Y2 <- subset(Y2, days_since_collection != 474)
day <- c(271,285,314,321,328,335,343,352,358,366,376,381,396,404,411,417,426,433,439,446,453,458,467,473)

semA_pergroup<-tapply(Y2$adults,list(Y2$days_since_collection,Y2$treatment_1),sem)
semA_vector<-c(as.numeric(semA_pergroup[,1]),as.numeric(semA_pergroup[,2]),as.numeric(semA_pergroup[,3]))
meanA_pergroup<-tapply(Y2$adults,list(Y2$days_since_collection,Y2$treatment_1),mean,na.rm=T)

par(mfrow=c(2,1))
plot(day, meanA_pergroup[,1],pch=16, las = 1, xlab = "time [d]", ylab = "# adults", col = pointcolor[1], ylim=c(0,160), xlim = c(min(day), max(day)),
     main = "Number of adults over time")
points(day, meanA_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanA_pergroup[,3],pch=16,col=pointcolor[3])
legend(270, 150, c("control", "low", "high"), col = pointcolor, pch = 16, title = "treatments")
for(i in 1:3){arrows(x0=day, x1=day,
                     y0=meanA_pergroup[,i]-semA_pergroup[,i],y1=meanA_pergroup[,i]+semA_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}
for(i in 1:3){lines(day, meanA_pergroup[,i], col = pointcolor[i])}

mtext('A', side=3, line=2, at=251, cex = 1.5, font = 2)

#pupae
semA_pergroup<-tapply(Y2$pupae,list(Y2$measurement_y2,Y2$treatment_1),sem)
semA_vector<-c(as.numeric(semA_pergroup[,1]),as.numeric(semA_pergroup[,2]),as.numeric(semA_pergroup[,3]))
meanA_pergroup<-tapply(Y2$pupae,list(Y2$measurement_y2,Y2$treatment_1),mean,na.rm=T)

plot(day, meanA_pergroup[,1],pch=16, las = 1, xlab = "time [d]", ylab = "# pupae", col = pointcolor[1], ylim=c(0,50), xlim = c(min(day), max(day)), 
     main = "Pupae Development Y2")
points(day, meanA_pergroup[,2],pch=16, col = pointcolor[2])
points(day, meanA_pergroup[,3],pch=16,col=pointcolor[3])
legend(270, 45, c("control", "high", "low"), col = pointcolor, pch = 16, title = "treatments")
for(i in 1:3){arrows(x0=day,x1=day,
                     y0=meanA_pergroup[,i]-semA_pergroup[,i],y1=meanA_pergroup[,i]+semA_pergroup[,i],
                     angle=90,length=0.05,code=3,col=linecolor[i])}
for(i in 1:3){lines(day, meanA_pergroup[,i], col = pointcolor[i])}
mtext('B', side=3, line=2, at=251, cex = 1.5, font = 2)

#### final count larva, pupae eggs ####

boxplot(final$eggs ~ final$treatment_1)
boxplot(final$larva ~ final$treatment_1)
boxplot(final$pupae ~ final$treatment_1)



#### peak days  ####


#peak day --> peakday for eggs is at day 20
sum(LnData$eggs[which(LnData$day == 18)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 19)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 20)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 22)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 23)], na.rm = TRUE)

#check between treatments
sum(LnData$eggs[which(LnData$day == 18 & LnData$treatment_1 == "control")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 19 & LnData$treatment_1 == "control")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 20 & LnData$treatment_1 == "control")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 22 & LnData$treatment_1 == "control")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 23 & LnData$treatment_1 == "control")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 18 & LnData$treatment_1 == "low")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 19 & LnData$treatment_1 == "low")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 20 & LnData$treatment_1 == "low")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 22 & LnData$treatment_1 == "low")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 23 & LnData$treatment_1 == "low")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 18 & LnData$treatment_1 == "high")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 19 & LnData$treatment_1 == "high")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 20 & LnData$treatment_1 == "high")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 22 & LnData$treatment_1 == "high")], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 23 & LnData$treatment_1 == "high")], na.rm = TRUE)
#peak day would be different if looked at only for each of the treatments: control 20 low 19 high 22


sum(LnData$eggs[which(LnData$day == 40)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 47)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 50)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 51)], na.rm = TRUE)
sum(LnData$eggs[which(LnData$day == 48)], na.rm = TRUE)
#second peak day over all would be at day 47 



# --> the peak day for eggs is measurementday nr. 20 
egg <- subset(LnData, day == 20)
head(egg)
shapiro.test(egg$eggs) # shapiro ist nicht signifikant --> Nullhypothese, dass daten Normal verteilt sind kann angenommen werden --> ANova oder LM

mod <- aov(eggs ~ treatment_1, data = egg, na.action=na.omit)
summary(mod)
plot(mod)

mod <- lm(eggs ~ treatment_1, data = egg)
summary(mod)

# there is no significant difference for eggs at their peak day in year 1

#peak day --> peakday for larva is at day 27 and at day 76

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
hist(larva$larva)

shapiro.test(larva2$larva) # shapiro ist  signifikant --> Nullhypothese, dass Daten Normal verteilt sind muss verworfen werden --> Kruskal wallis oder glm
hist(larva2$larva)

kruskal.test(larva ~ treatment_1, data = larva)

kruskal.test(larva ~ treatment_1, data = larva2)

mod <- glm(eggs ~ treatment_1, data = larva)
summary(mod)

# there is no significant difference for eggs at their peak days in year 1

















#### Effects of Neonicotinoids on Locomotion behavior ####

speed <- read.table("speed.txt", header = TRUE)
head(speed)
str(dat)

speed$treat1 <- as.factor(speed$treat1)
speed$treat1 <- factor(speed$treat1, levels = c("control", "low", "high"))
speed$treat2 <- as.factor(speed$treat2)

#ANV only 
speed <- subset(speed, treat2 =="control")

#response variables: initial speed (inispeed), time_active, avspeed movment nrgrooming
#random factor: colonies and day

#transform time into Minutes only

speed$time2 <- times(speed$time)
speed$time_min <- 60 * hours(speed$time2) + minutes(speed$time2)
speed$time_min

plot(speed$ini_speed ~ speed$date) #day might have an effect and will be included as random foactor
plot(speed$ini_speed ~ speed$time_min)
mod_time <- lm(speed$ini_speed ~ speed$time_min)
abline(mod_time)              # time of day does not really seem to be relevant 

dat <- speed


#### For initial speed ####
hist(dat$ini_speed, breaks = 7)

head(speed)

par(mfrow=c(1,1))
boxplot(speed$ini_speed ~ speed$treat1) # plot data with boxplots to have a first idea of the data
boxplot(speed$av_speed ~ speed$treat1) 
boxplot(speed$active_time ~ speed$treat1)
boxplot(speed$sum_move ~ speed$treat1) 
boxplot(speed$nr_stops ~ speed$treat1)
boxplot(speed$drops ~ speed$treat1)

#step 1 standardize explanatory variables to mean zero and a Sd of 1 to make it easyer to compare effect sizes use comand sclae(
# dat$variable_new <- scale(dat$variable)



mod <- glm(speed$nr_stops ~ speed$treat1, family=poisson)


mod <- lmer(ini_speed ~ treat1 + time_min + (1|colony) + (1|date), data=dat, REML=FALSE)
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
library(extrafont)
fonts()

par(mfrow=c(1,1))
par(oma=c(1,0,0,0))

col <- c("darkolivegreen1", "darkorange3", "darkred")
boxplot(dat$inispeed ~ dat$treat1 + dat$treat2, border = "white", 
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







##### Duration of lifestages #####

DLS <- read.table("Life Stages.txt", header = TRUE)
DLS$Treatment <- as.factor(DLS$treatment)

#duration
boxplot(DLS$D_EggStage ~ DLS$Treatment)
boxplot(DLS$D_LarvaeStage ~ DLS$Treatment)
boxplot(DLS$D_PupaeStage ~ DLS$Treatment)


#day
boxplot(DLS$Day_Eggs ~ DLS$Treatment)
boxplot(DLS$Day_Larvae ~ DLS$Treatment)
boxplot(DLS$Day_Pupae ~ DLS$Treatment)
boxplot(DLS$`Day_Pupae without cocoon` ~ DLS$Treatment)
boxplot(DLS$Day_Adults ~ DLS$Treatment)


##### Repeat the results from Lars #####

dat <- read.table("danielschlaeppi_AN-tableforLars_2018.txt", header = TRUE)
dat$treatment_lars <- as.factor(dat$treatment_lars)
levels(dat$treatment_lars)
levels(dat$treatment_lars) <- c("control","low","high")
dat$treatment_yearspecific <- as.factor(dat$treatment_yearspecific)
summary(dat$treatment_lars)
summary(dat$treatment_yearspecific)

finals <- subset(dat, measurement_day==76 | measurement_day==101)
final1 <- subset(dat, measurement_day==76)
final2 <- subset(dat, measurement_day==101)

boxplot(finals$eggs ~ finals$treatment_yearspecific, ylab = "eggs", xlab = "treatments over 2 years clhclh")
hist(finals$eggs, xlab = "eggs", ylab= "frequency")
hist(final1$eggs, xlab = "eggs", ylab= "frequency")
hist(final2$eggs, xlab = "eggs", ylab= "frequency")

shapiro.test(finals$eggs) #not normally distributed 
shapiro.test(final1$eggs) #not normally distributed 
shapiro.test(final2$eggs) #not normally distributed



#full model
mod <- glmer(eggs ~ treatment_lars * year + (1|year), data = dat, family = poisson)
summary(mod)

scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$colony[,1]) # qq of random effects
qqline(ranef(mod)$colony[,1])


# separate days
boxplot(finals$eggs ~ finals$treatment_yearspecific)

shapiro.test(final1$eggs) #not normally distributed 
kruskal.test(final1$eggs ~ final1$treatment_lars)

shapiro.test(final2$eggs) #not normally distributed 
kruskal.test(final2$eggs ~ final1$treatment_lars)

mod <- glm(final2$eggs ~ final2$treatment_lars)
summary(mod)

plot(mod)

summary(final2$treatment_lars)


### simple analyses at end of year 2#

#eggs
shapiro.test(final2$eggs) #significant not normally distributed 
hist(final2$eggs)
boxplot(final2$eggs~final2$treatment_lars, main = "eggs", xlab= "treatments", ylab= "#eggs", cex.axis= 1.5, cex.main = 1.5, cex.lab=1.5)
mod <- glm(final2$eggs~final2$treatment_lars, family="poisson")
summary(mod)
plot(mod)

final2$treatment_lars <- relevel(final2$treatment_lars, ref = "low")
mod <- glm(final2$eggs~final2$treatment_lars, family="poisson")
summary(mod)
final2$treatment_lars <- relevel(final2$treatment_lars, ref = "control")

#larva
shapiro.test(final2$larva) #significant not normally distributed 
hist(final2$larva)
boxplot(final2$larva~final2$treatment_lars, main = "larva", xlab= "treatments", ylab= "#larva", cex.axis= 1.5, cex.main = 1.5, cex.lab=1.5)
mod <- glm(final2$larva~final2$treatment_lars, family="poisson")
summary(mod)
plot(mod)

final2$treatment_lars <- relevel(final2$treatment_lars, ref = "low")
mod <- glm(final2$larva~final2$treatment_lars, family="poisson")
summary(mod)
final2$treatment_lars <- relevel(final2$treatment_lars, ref = "control")


#pupae
shapiro.test(final2$pupae) #significant not normally distributed 
hist(final2$pupae)
boxplot(final2$pupae~final2$treatment_lars, main = "pupae", xlab= "treatments", ylab= "#pupae", cex.axis= 1.5, cex.main = 1.5, cex.lab=1.5)
mod <- glm(final2$pupae~final2$treatment_lars, family="poisson")
summary(mod)
plot(mod)

final2$treatment_lars <- relevel(final2$treatment_lars, ref = "low")
mod <- glm(final2$pupae~final2$treatment_lars, family="poisson")
summary(mod)
final2$treatment_lars <- relevel(final2$treatment_lars, ref = "control")


#adults 
levels(finals$treatment_yearspecific) <- c("control y1","low y1","high y1","control y2","low y2","high y2")
boxplot(finals$adults~finals$treatment_yearspecific, main = "adults", xlab= "treatments", ylab= "#adults", cex.axis= 1.5, cex.main = 1.5, cex.lab=1.5)
head(finals)
mod <- glmer(adults ~ treatment_lars * year + (1|colony), data=finals, family = poisson)
summary(mod)

boxplot(final1$adults~final1$treatment_lars, main = "adults year 1", xlab= "treatments", ylab= "#adults", cex.axis= 1.5, cex.main = 1.5, cex.lab=1.5)
boxplot(final2$adults~final2$treatment_lars, main = "adults year 2", xlab= "treatments", ylab= "#adults", cex.axis= 1.5, cex.main = 1.5, cex.lab=1.5)

shapiro.test(final1$adults) #significant not normally distributed 
hist(final1$adults)
mod <- glm(final1$adults~final1$treatment_lars, family="poisson")
summary(mod)
plot(mod)

final1$treatment_lars <- relevel(final1$treatment_lars, ref = "low")
mod <- glm(final1$adults~final1$treatment_lars, family="poisson")
summary(mod)
final1$treatment_lars <- relevel(final1$treatment_lars, ref = "control")
  
shapiro.test(final2$adults) #significant not normally distributed 
hist(final2$adults)
mod <- glm(final2$adults~final2$treatment_lars, family="poisson")
summary(mod)
plot(mod)

final2$treatment_lars <- relevel(final2$treatment_lars, ref = "low")
mod <- glm(final2$adults~final2$treatment_lars, family="poisson")
summary(mod)
final2$treatment_lars <- relevel(final2$treatment_lars, ref = "control")
















###############################################################################
### Neonics residues in Ants ##################################################
###############################################################################

NIC <- read.table("ANV_Neonic_in_workers_and_queens.txt", header = TRUE)
head(NIC)

#as.factor() at the moments the treatments are already in the right state ()
NIC$treatment <- as.factor(NIC$treatment)
NIC$TREATMENT_LARS <- as.factor(NIC$TREATMENT_LARS)


# relevel treatment 1 to control low high 
levels(NIC$treatment_1)
NIC$treatment_1 <- relevel(NIC$treatment_1, ref ="low")
NIC$treatment_1 <- relevel(NIC$treatment_1, ref ="control")
str(NIC)

#remove the one half of the data that was assigned to the virus treatments
AN <- subset(NIC, treatment_2=="control")
table(AN$treatment)
levels(AN$treatment)

#only samples that were tested and from the colonies that lasted to the end / exclude failure colonies
#AN <- subset(AN, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")

AN <- subset(AN, final == "y")
table(AN$treatment_1)

workers <- subset(AN, sample=="worker")
queens <- subset(AN, sample =="queen")

control_w <- subset(workers, treatment_1 == "control")
low_w <- subset(workers, treatment_1 == "low")
high_w <- subset(workers, treatment_1 == "high") 
control_q <- subset(queens, treatment_1 == "control")
low_q <- subset(queens, treatment_1 == "low")
high_q <- subset(queens, treatment_1 == "high")


par(mfrow=c(1,1))

boxplot(workers$conc_T_dry ~ workers$treatment_1)
boxplot(workers$conc_T_fresh ~ workers$treatment_1)
boxplot(workers$conc_T_abd ~ workers$treatment_1)
boxplot(workers$conc_C_dry ~ workers$treatment_1)
boxplot(workers$conc_C_fresh ~ workers$treatment_1)
boxplot(workers$conc_C_abd ~ workers$treatment_1)

boxplot(queens$conc_T_dry ~ queens$treatment_1)
boxplot(queens$conc_T_fresh ~ queens$treatment_1)
boxplot(queens$conc_T_abd ~ queens$treatment_1)
boxplot(queens$conc_C_dry ~ queens$treatment_1)
boxplot(queens$conc_C_fresh ~ queens$treatment_1)
boxplot(queens$conc_C_abd ~ queens$treatment_1)


# basically if per fresh weight, dry weight or per abdomen does not really matter --> per dry weight seems to make the most sense 
#create a nice 4x4 graph
par(mfrow=c(2,2))
boxplot(workers$conc_T_dry ~ workers$treatment_1)
boxplot(workers$conc_C_dry ~ workers$treatment_1)
boxplot(queens$conc_T_dry ~ queens$treatment_1)
boxplot(queens$conc_C_fresh ~ queens$treatment_1)


#remove outlier if necessary
# If it is obvious that the outlier is due to incorrectly entered or measured data, you should drop the outlier. 
# If the outlier does not change the results but does affect assumptions, you may drop the outlier.  But note that in a footnote of your paper.
# More commonly, the outlier affects both results and assumptions.  In this situation, it is not legitimate to simply drop the outlier.  You may run the analysis both with and without it, but you should state in at least a footnote the dropping of any such data points and how the results changed.
# If the outlier creates a significant association, you should drop the outlier and should not report any significance from your analysis.
# luckily the outlier was a mistake from the datasheet :-) 

par(mfrow=c(1,1))

mod <- lm(workers$nr_adults ~ workers$conc_T_fresh)
plot(workers$nr_adults ~ workers$conc_T_fresh)
abline(mod)
mod <- lm(queens$nr_adults ~ queens$conc_T_fresh)
plot(queens$nr_adults ~ queens$conc_T_fresh)
abline(mod)
mod <- lm(workers$nr_adults ~ workers$conc_C_fresh)
plot(workers$nr_adults ~ workers$conc_C_fresh)
abline(mod)
mod <- lm(queens$nr_adults ~ queens$conc_C_fresh)
plot(queens$nr_adults ~ queens$conc_C_fresh)
abline(mod)




#plot(workers$nr_adults ~ queens$nr_adults)
#x <- workers$nr_adults
#y <- queens$nr_adults
#writeClipboard(as.character(x)) # copy paste a vetor into exel
#writeClipboard(as.character(y))
#plot(workers$conc_C_fresh ~ queens$conc_T_fresh)

# correlation test
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r








