
#### Eurobee numbers ####

setwd("H:/DS-Bees/R/ANV")
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV")
dat <- read.table("ANVlong_full.txt", header = TRUE)


# load necessary librarys
library("lme4")
library("plyr")
library(lmerTest)
library(chron)
library(lme4)
library(arm)
library(car)


### fancy function to reset graphical parameters of par to its original values ###

#use the function by typing reset_par()




# Colonylevel effects of Virus? # =Difference in Adults and Pupae between first and last feeding. 16.07.17-01.10.2017
# full dataset for understanding the effect

head(dat)
dif <- subset(dat, date == "16.07.2017" | date == "01.10.2017" | date =="16.10.2017")
final <- subset(dat, final == "year2")

# Exclude the colonies with queen failure (dead or never developped brood)
# Colonies to excluce: C9C, C12C C15V, C16V, L3C, L10V, L14C, H4C, H11V, H18V, H19C

dif <- subset(dif, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")
final <- subset(final, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")


#Correct the number of Adults on date 1 compared to final one based on a simple factor
# eg. estimate date1 135, final count 300 --> 300/135 = 2.222222 

date1 <- subset(dif, date == "16.07.2017")
date2 <- subset(dif, date == "01.10.2017")
date3 <- subset(dif, date == "16.10.2017")

date1$adults_c <- round(date1$adults*(final$adults/date3$adults))
date2$adults_c <- round(date2$adults*(final$adults/date3$adults))
date3$adults_c <- round(date3$adults*(final$adults/date3$adults))

date1$final_adults <- final$adults

#create the difference between time variables (development of adults and pupae)
date1$dif_adults <- date2$adults-date1$adults
date1$dif_pupae  <- date2$pupae-date1$pupae
date1$dif_adults2 <- date3$adults-date1$adults
date1$dif_pupae2 <- date3$pupae-date1$pupae
date1$dif_adults_c <- date3$adults_c-date1$adults_c


#Explore the data
boxplot(date1$dif_adults ~ date1$treat)
boxplot(date1$dif_adults2 ~ date1$treat)
boxplot(date1$dif_adults_c ~ date1$treat)

boxplot(date1$dif_pupae ~ date1$treat)
boxplot(date1$adults ~ date1$treat)
boxplot(date1$pupae ~ date1$treat)
boxplot(date2$adults ~ date1$treat)
boxplot(date2$pupae ~ date1$treat)


### There seems to be a colony effect on the number of adults that emerged during that time and the new developement of pupae --> less brood
# Graph only for Eurobee: Colony level effect only for Neoniccontrols

noneonic <- subset(date1, treatment_1 == "control")
boxplot(noneonic$dif_adults_c ~ noneonic$treat)
boxplot(noneonic$final_adults ~ noneonic$treat)
boxplot(noneonic$dif_pupae ~ noneonic$treat)

#### create a model and a nice graph based on bayesian inference to see wether or not these are significant differences ####

#The variables that might influence the difference in number of ants between the two dates: Virustreatment, colonysize at beginning, number of pupae at begingng
plot(noneonic$dif_adults_c ~ noneonic$pupae)
plot(noneonic$dif_adults_c ~ noneonic$adults_c)
boxplot(noneonic$dif_adults_c ~ noneonic$treat)
boxplot(noneonic$adults_c ~ noneonic$treat)
boxplot(noneonic$adults ~ noneonic$treat)

hist(noneonic$dif_adults_c, breaks = 7)


#see if there is an inital significant difference in adult number at the beginning: 
shapiro.test(noneonic$adults) # data is normaly distributed as shapirotest ist 0.14 and with that the nullhypothesis that the data is Normally distributed is not rejected
t.test(noneonic$adults_c~noneonic$treat)  #or glm as data is count data and thus follows a poisson distribution ??? or can it be treated the same? 

#es ist richtig, dass du Zähldaten hast und keine kontinuierliche Messungen. Aus diesem Grunde wäre die erste Modellwahl ein Poissonmodell. 
#Das Poissonmodell macht aber sehr strenge Annahmen über die Verteilung, inbesondere nimmt es an, dass die Varianz gleich dem Mittelwert ist. 
#In vielen biologischen Daten ist die Varianz aber höher als der Mittelwert. Das ist auch in deinen Daten so. Das sieht man am Plot der Residuen gegen die Leverages. 
#Einge Beobachtungen haben einen starken Einfluss (Cooks Distanz >0.5). Das deutet relativ deutlich darauf hin, dass deine Ameisen stärker streuen als es eine Poissonverteilung annimmt. 
#Die Poissonverteilung nähert sich einer Normalverteilung an, wenn der Mittelwert hoch ist. Das ist schon bei Mittelwerten über 20 der Fall. 
#Deine Mittelwerte sind alle höher und deshalb kannst du gut eine Normalverteilung annehmen. 
#Die Normalverteilung hat den grossen Vorteil, dass die Varianz aus den Daten geschätzt wird und deshalb sind die Vertrauensintervalle zuverlässiger als im Poissonmodell, 
#wo sie eben zu schmal sind, weil eine zu kleine Varianz angenommen wurde.

#Kein Modell ist perfekt, aber in deinem Fall liefert das Poissonmodell eine Resultat, das eine zu hohe Sicherheit vortäuscht. 
#Das Modell mit Normalverteilung liefert zwar krumme Zahlen (reale Zahlen anstatt Zähldaten), wenn du Vorhersagen machen möchtest, aber dafür werden die Vertrauensintervalle
#basierend auf der  aus den Daten geschätzten Varianz berechnet. Die Residuen bieten keinen Grund zur Sorge, deshalb kannst du sehr gut  mit der Normalverteilung fahren. 
#Wenn du ganz sicher sein möchtest, kontrollierst du nochmals den Residuenplot, wo die Wurzel der absoluten Residuen gegen den Anpassungswert geplottet werden. 
#In Zähldaten kann es sein, dass die Varianz sich mit dem Mittelwert verändert. Dafür korrigiert das normale lineare Modell nicht. 

# --> Das bedeutet wir nehmen normal lm mit Normalverteilung sofern die Assumptinsplots gut aussehen. 


mod <- lm(adults_c~treat, data=noneonic)
mod <- lm(adults~treat, data=noneonic)
qqplot(mod)
summary(mod)

#--> initially there is no significant difference in colony size if normal count is used  --> however there is a significant difference if 
# corrected collony size is used as the factor that gets applied is diffrent as it is based on the final counts which differ due to the treatments

boxplot((final$adults/date3$adults) ~ date3$treat) # to see that the factor of replication is indeed not the same over the treatments

boxplot(noneonic$dif_adults_c ~ noneonic$treat)
par(mfrow=c(1,1))



mod <- lm(dif_adults_c ~ treat + adults + pupae, data = noneonic) 

# mod <- glm(dif_adults_c ~ treat + adults + pupae, data = noneonic, family="poisson") 

mod
summary(mod)str(mod)
Anova(mod)


mean(noneonic$dif_adults_c)
sd(noneonic$dif_adults_c)
hist(noneonic$dif_adults_c, breaks = 10)



#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
#qqnorm(ranef(mod)$adults_c[,1]) # qq of random effects
#qqline(ranef(mod)$adults_c[,1])
#qqnorm(ranef(mod)$pupae[,1]) # The varince of date is so small that in the qq-plot it collapses to zero --> it is a bug and it can be assumed to be close to zero
#qqline(ranef(mod)$pupae[,1])


#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)


apply(bsim@coef, 2,quantile, prob=c(0.025,0.5,0.975))

noneonic$treat <- as.factor(noneonic$treat)

#fitted values with 95% credible intervals
newdat<-expand.grid(treat=factor(c('1','2'),levels=levels(noneonic$treat))) 

newdat$adults=mean(noneonic$adults)
newdat$pupae=mean(noneonic$pupae)
head(newdat)



Xmat <- model.matrix(~treat + adults + pupae, data=newdat)      ####second part of model without of random factors

??Xmat

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@coef[i,] 

newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%coef(mod)
newdat

#create a super nice Graph: 
par(mar=c(5,6,5,2))
par(mfrow=c(1,1))
par(oma=c(1,1,0,0))
par(cex=1)

col <- c("grey78", "grey60")
boxplot((noneonic$dif_adults_c) ~ noneonic$treat, border = "white", 
        main = "Difference in Nr. of Adults between t1 and t2", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
axis(1, at = c(1,2), labels=c("no virus", "virus"), cex.axis = 1.5)
title(ylab="Number of Adults", line=2.5, cex.lab=1.5)

stripchart(dif_adults_c ~ treat, data = noneonic, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE) 
#plot posteriors
points(c(1,2), newdat$fit, pch = 19, cex=2)
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5)

#Calculate the effectsize respectively how much the virus affects colonydevelopement
str(fitmat)
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))




##### Presentation Journal Club #######

head(date3)
boxplot(date1$adults ~ date1$treatment_2, main = "adults, 16.07.2017", xlab = "treatments", ylab = "number of adults", cex.main=2, cex.lab=1.5, cex.axis=1.5)
shapiro.test(date1$adults) # significant --> not normally distributed --> man whitney test
wilcox.test(adults ~treatment_2, data=date1)  #not significant

boxplot(date3$adults ~ date3$treatment_2, main = "adults, 16.10.2017", xlab = "treatments", ylab = "number of adults", cex.main=2, cex.lab=1.5, cex.axis=1.5)
shapiro.test(date3$adults) # nonsignificant --> normally distributed --> students t-test
t.test(adults ~ treatment_2, data=date3)  #not significant


boxplot(date1$pupae ~ date1$treatment_2, main = "pupae, 16.07.2017", xlab = "treatments", ylab = "number of pupae", cex.main=2, cex.lab=1.5, cex.axis=1.5)
shapiro.test(date1$pupae) # nonsignificant -->  normally distributed --> students t test
wilcox.test(pupae ~treatment_2, data=date1)  #not significant


boxplot(date3$pupae ~ date3$treatment_2, main = "pupae, 16.10.2017", xlab = "treatments", ylab = "number of pupae", cex.main=2, cex.lab=1.5, cex.axis=1.5)
shapiro.test(date3$pupae) # significant --> not normally distributed --> man whitney wilcox test
wilcox.test(pupae ~ treatment_2, data=date3)


boxplot(date1$dif_adults2 ~ date1$treatment_2, main = "change in number of adults", xlab = "treatments", ylab = "difference in # of adults", cex.main=2, cex.lab=1.5, cex.axis=1.5)
shapiro.test(date1$dif_adults2) # knapp non significant -->  normally distributed --> students t test
t.test(dif_adults2 ~ treatment_2, data=date1) # knapp nicht signifikant

boxplot(date1$dif_pupae2 ~ date1$treatment_2, main = "change in number of pupae", xlab = "treatments", ylab = "difference in # of pupa", cex.main=2, cex.lab=1.5, cex.axis=1.5) 
shapiro.test(date1$dif_pupae2) #significant -->  not normally distributed --> wilcox test
t.test(dif_pupae2 ~ treatment_2, data=date1) #signifikant









#### Colony level effects --> increase in Colony size  ####


dat <- read.table("ANVlong_full.txt", header = TRUE)
date1 <- subset(dat, date == "16.07.2017")
date2 <- subset(dat, date == "01.10.2017")
date1 <- subset(date1, treatment_1=="control") #(without neonic treatments)
date2 <- subset(date2, treatment_1=="control")

date1$diff_adults <- date2$adults-date1$adults
date1$diff_pupae <- date2$pupae-date1$pupae

date1$diff_adults

par(cex = 1.8)
levels(date1$treatment_2)[levels(date1$treatment_2)=="control"] <- "no virus"
boxplot(date1$diff_adults ~ date1$treatment_2, xlab= "treatments", ylab = "nr. of adult ants", main = "")


boxplot(date1$diff_pupae ~ date1$treatment_2)
boxplot(date1$pupae ~ date1$treatment_2)
boxplot(date1$adults ~ date1$treatment_2)
boxplot(date2$adults ~ date2$treatment_2)
boxplot(date2$pupae ~ date2$treatment_2)



# t test for adults
shapiro.test(date1$diff_adults) #non significant --> Ho: Normaldistributed data can be assumed
hist(date1$diff_adults)
t.test(date1$diff_adults ~ date1$treatment_2)

# t test for pupae
shapiro.test(date1$diff_pupae) #non significant --> Ho: Normaldistributed data can be assumed
hist(date1$diff_pupae)
t.test(date1$diff_pupae ~ date1$treatment_2)

# significant difference for the difference in numbers of adults in the two groups 







#### Movement analyses for graphs of eurobee poster ####

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


dat <- read.table("speed.txt", header = TRUE)

dat$treat1 <- as.factor(dat$treat1)
dat$treat1 <- factor(dat$treat1, levels = c("control", "low", "high"))
dat$treat2 <- as.factor(dat$treat2)

#response variables: initial speed (inispeed), time_active, avspeed movment nrgrooming
#random factor: colonies and day

#transform time into Minutes only
dat$time2 <- times(dat$time)
dat$time_min <- 60 * hours(dat$time2) + minutes(dat$time2)
dat$time_min

##### only for the controls that did not have neonics for biology 18  #####
cdat <- subset(dat, treat1 == "control") 

#### initial speed ####
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
Anova(mod) #we used this p valiue


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
newdat1<-expand.grid(
  treat2=factor(c('control','virus'),levels=levels(dat$treat2))) 
newdat1$time_min=mean(cdat$time_min)
head(newdat1)

Xmat <- model.matrix(~treat2 + time_min, data=newdat1)      ####second part of model without of random factors
fitmat <- matrix(ncol=nsim, nrow=nrow(newdat1))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat1$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat1$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat1$fit <- Xmat%*%fixef(mod)
newdat1

#create a super nice Graph: 
par(mar=c(5,6,5,2))
par(mfrow=c(1,2))
par(oma=c(1,1,0,0))
par(cex=1.5)

col <- c("grey78", "grey78")
boxplot(cdat$ini_speed ~ cdat$treat2, border = "white", 
        main = "", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n")
axis(1, at = c(1,2), labels=c("no virus", "virus"))
title(ylab="initial speed [~cm/s]", line=2.5)

stripchart(ini_speed ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE, cex=0.8) 

#plot posteriors
points(c(1,2), newdat1$fit, pch = 19, cex=2, col = c("black", "black"))
segments(c(1,2), newdat1$lower, c(1,2), newdat1$upper, lwd = 5, col = c("black", "black"))

#Calculate the effectsize respectively how much the virus affects movement
str(fitmat)
head(fitmat)

novirus <- apply(fitmat[1,],2, mean)
virus <- apply(fitmat[1,],2, mean)
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
        main = "", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n")
axis(1, at = c(1,2), labels=c("no virus", "virus"))
title(ylab="average speed [~cm/s]", line=3)

stripchart(av_speed ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE, cex= 0.8) 

#plot posteriors
points(c(1,2), newdat$fit, pch = 19, cex=2, col = c("darkolivegreen3", "red"))
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5, col = c("darkolivegreen3", "red"))

#Calculate the effectsize respectively how much the virus affects movement
full <- fitmat[1,] - fitmat[2,]
quantile(full, prob = c(0.025, 0.5, 0.975))





### triplegraph with average speed followed by initial speed

par(mar=c(4,4,1,1))
par(mfrow=c(1,3))
par(oma=c(1,1,0,0))
par(cex=2)

boxplot(cdat$av_speed ~ cdat$treat2, border = "white", 
        main = "", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n")
axis(1, at = c(1,2), labels=c("no virus", "virus"))
title(ylab="average speed [~cm/s]", line=3)
stripchart(av_speed ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE, cex=0.8) 
points(c(1,2), newdat$fit, pch = 19, col = c("black", "black"))
segments(c(1,2), newdat$lower, c(1,2), newdat$upper, lwd = 5, col = c("black", "black"))


boxplot(cdat$ini_speed ~ cdat$treat2, border = "white", 
        main = "", ylab = "", xlab = "treatments", 
        family = c("mono"), xaxt="n")
axis(1, at = c(1,2), labels=c("no virus", "virus"))
title(ylab="initial speed [~cm/s]", line=2.5)
stripchart(ini_speed ~ treat2, data = cdat, 
           vertical = TRUE, at = c(1,2) ,method = "jitter", 
           pch = c(19,19,19,17,17,17), col= col,
           add = TRUE, cex=0.8) #???
points(c(1,2), newdat1$fit, pch = 19, col = c("black", "black"))
segments(c(1,2), newdat1$lower, c(1,2), newdat1$upper, lwd = 5, col = c("black", "black"))


dat <- read.table("ANVlong_full.txt", header = TRUE)
date1 <- subset(dat, date == "16.07.2017")
date2 <- subset(dat, date == "01.10.2017")
date1 <- subset(date1, treatment_1=="control") #(without neonic treatments)
date2 <- subset(date2, treatment_1=="control")
date1$diff_adults <- date2$adults-date1$adults
date1$diff_adults
levels(date1$treatment_2)[levels(date1$treatment_2)=="control"] <- "no virus"
boxplot(date1$diff_adults ~ date1$treatment_2, xlab= "treatments", ylab = "nr. of adult ants", main = "")
