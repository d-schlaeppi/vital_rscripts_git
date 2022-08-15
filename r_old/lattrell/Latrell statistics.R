#set working directory and load table

getwd()
setwd("H:/DS-Bees/R/lattrell")
dat <- read.table("stattable.txt", header = TRUE)
names(dat)

#load packages
library("nlme")           # for linear models
library("lme4")           # linear mixed effect models
library("ggplot2")        # for nice graphs
library(plotrix)


# sample = number of pcr probe, vcps = viral copy per sample, name = lattrells sample names, week = time of sampling, id = colony or pupae number


#get ab idea of the data
str(dat)



attach(dat)
table(treat,week)
hist(vcps)
hist(log(vcps))
hist(log10(vcps))
detach(dat)

qqplott <- qplot(sample = log10(dat$vcps), stat="qq")
qqplott

dat.ana <- subset(dat, treat!=4)
c <- subset(dat, treat==3) 
t12 <- subset(dat, treat==1 | treat==2) 
t4 <- subset(dat, treat==4)

qqplott <- qplot(sample = log10(c$vcps), stat="qq")
qqplott

qqplott <- qplot(sample = log10(t12$vcps), stat="qq")
qqplott


qqplott <- qplot(sample = log10(t4$vcps), stat="qq")
qqplott

normal.lm = lm(log10(vcps) ~ treat, data = dat.ana)
normal.stdres = rstandard(normal.lm)
qqnorm (normal.stdres, ylab="Standardized Residuals",  xlab="Normal Scores", main="Normaldistribution?")
qqline(normal.stdres)

shapiro.test(log10(dat.ana$vcps))
hist(log10(dat.ana$vcps), breaks = 15, freq = FALSE)
curve(dnorm(x, mean=mean(log10(dat.ana$vcps)), sd=sd(log10(dat.ana$vcps)), add=TRUE, col=blue, lwd=2)) 

# vcps seems to be normally distributed with logtransformation

w1<- subset(dat, week==1)
hist(log(w1$vcps))
hist(log10(w1$vcps))

#subset nur mit treatments um viral copies zu sehen (quantile)
t12 <- subset(dat, treat==1 | treat==2) 
quantile(t12$vcps)
quantile(c$vcps)
quantile(t4$vcps)

# model: vcps depends on the treatment, there are three treatment = 3 (1 control, 2 real treatments (treatment 4 excluted as it is only a positive control)
# week as random factor? in the full model yes, but not if weeks are looked at separately (thus no need to look at them as a time series with related data)

#subeset whithout treatment 4

dat.t <- subset(dat, treat !=4)
hist(log(dat.t$vcps))
hist(log10(dat.t$vcps))

dat.t$treat <- as.factor(dat.t$treat)
mfull <- lmer(log10(vcps) ~ factor(treat) + (1|week), data=dat.t)
summary(mfull)


attach(dat.t)
dat.t$treat3 <- relevel(dat.t$treat, ref="3")
mfull3 <- lmer(log10(vcps)~factor(treat3) + (1|week), data=dat.t)
summary(mfull3)

dat.t$treat1 <- relevel(dat.t$treat, ref="1")
mfull1 <- glmer(log10(vcps)~factor(treat1) + (1|week), data=dat.t)
summary(mfull1)
detach(dat.t)

### full model did not show any significant differences -> thus we analyze the weeks separately


hist(log10(dat$vcps))
two #DIFFERENCES FOR EACH WEEK: 
w1<- subset(dat.t, week==1)
w2<- subset(dat.t, week==2)
w3<- subset(dat.t, week==3)
w4<- subset(dat.t, week==4)


# week 1
m1 <- lm(log10(vcps) ~ factor(treat), dat=w1)
summary(m1)
predict.lm(m1, se.fit=T, interval="confidence")

glmm1 <- glmer(log10(vcps) ~ factor(treat)+ (1|id), data=w1) 


w1$treat1 <- relevel(factor(w1$treat), ref="3")
m1 <- lm(log10(vcps)~factor(treat1), dat=w1)
summary(m1)

# 1a, 2a, 3b

#week 2 
m2 <- lm(log(vcps)~factor(treat), dat=w2)
summary(m2)
predict.lm(m2, se.fit=T, interval="confidence")

w2$treat1 <- relevel(w2$treat, ref="1")
m2 <- lm(log(vcps)~factor(treat1), dat=w2)
summary(m2)
# 1a, 2a, 3b

#week 3
m3 <- lm(log(vcps)~factor(treat), dat=w3)
summary(m3)
predict.lm(m3, se.fit=T, interval="confidence")

w3$treat1 <- relevel(w3$treat, ref="1")
m3 <- lm(log(vcps)~factor(treat1), dat=w3)
summary(m3)

#week 4 
m4 <- lm(log(vcps)~factor(treat), dat=w4)
summary(m4)
predict.lm(m4, se.fit=T, interval="confidence")



###
# Create a nice graph with the predicted values for the weekly models

dat.plot <- read.table("plottable.txt", header = TRUE) 
dat.plot$treat <- as.factor(dat.plot$treat) 

#### in ggplot  we stopped in the middle and did it without ggplot
#finalplot <- ggplot(dat.plot, aes( week, logvcps, color=x, group = treat))
#finalplot + geom_point() + geom_errorbar(aes(ymin=lwr,ymax=upr, position=logvcps), width=0.1) + geom_line(aes(color=x))
#x <- c("black", "blue", "yellow","black", "blue", "yellow","black", "blue", "yellow","black", "blue", "yellow")




t1 <- subset(dat.plot, dat.plot$treat == 1)
t1$treat

t2 <- subset(dat.plot, dat.plot$treat == 2)
t2$treat

t3 <- subset(dat.plot, dat.plot$treat == 3)
t3$treat

par(mar=c(5,5,2,2))

plotCI(dat.plot$week, dat.plot$logvcps, ui = dat.plot$upr, li = dat.plot$lwr, xlab = "week", ylab="LN [predicted genomic DWV copies / ant] ", 
       cex.lab = 2, cex.axis= 2, cex.main=4, pch=19, col=c("red", "blue", "black"), lwd = 2, xaxt="n")
axis(side=1, at=seq(0,4,by=1), cex.axis=2)
lines(t1$week, t1$logvcps, add=TRUE, col="red", lwd = 3 )
lines(t2$week, t2$logvcps, add=TRUE, col="blue", lwd = 3)
lines(t3$week, t3$logvcps, add=TRUE, col="black", lwd = 3)
legend(3.1,23.5, c("treatment 1", "treatment 2", "control"), box.lwd=0, bty='n', bg="transparent", pch=19, col=c("red", "blue", "black"), cex=2)
text(1,14, c("***"), cex=2)
text(2,14, c("***"), cex=2)
text(3,14, c("***"), cex=2)
text(4,14, c("***"), cex=2)

####
# New graph with patricks table: 

dat.plot <- read.table("plottablog10.txt", header = TRUE) 
dat.plot$treat <- as.factor(dat.plot$treat)


t1 <- subset(dat.plot, dat.plot$treat == 1)
t1$treat
t2 <- subset(dat.plot, dat.plot$treat == 2)
t2$treat
t3 <- subset(dat.plot, dat.plot$treat == 3)
t3$treat

plotCI(dat.plot$week, dat.plot$logvcps, ui = dat.plot$upr, li = dat.plot$lwr, xlab = "week", ylab="log [predicted genomic DWV copies / ant] ", 
       cex.lab = 2, cex.axis= 2, cex.main=4, pch=19, col=c("red", "blue", "black"), lwd = 2, xaxt="n")
axis(side=1, at=seq(9,12,by=1), cex.axis=2)
lines(t1$week, t1$logvcps, add=TRUE, col="red", lwd = 3 )
lines(t2$week, t2$logvcps, add=TRUE, col="blue", lwd = 3)
lines(t3$week, t3$logvcps, add=TRUE, col="black", lwd = 3)
legend(11,10.5, c("treatment 1", "treatment 2", "control"), box.lwd=0, bty='n', bg="transparent", pch=19, col=c("red", "blue", "black"), cex=2)
text(9,6, c("***"), cex=2)
text(10,6, c("***"), cex=2)
text(11,6, c("***"), cex=2)
text(12,6, c("***"), cex=2)

############################################################
### Check if threre are changes of viral titers over time ##
############################################################


timemodel <- lmer(log10(vcps) ~ factor(week) + (1|treat), data=dat.t)
summary(timemodel)
t12 <- subset(dat, treat==1 | treat==2) 
timemodelt12 <- lmer(log10(vcps) ~ factor(week) + (1|treat), data=t12)
summary(timemodelt12)

t1 <- subset(dat, treat==1) 
t2 <- subset(dat, treat==2) 
c <- subset(dat, treat==3) 
t12


mt12 <- lm(log10(vcps) ~ factor(week), dat=t12)
summary(mt12)

t12$week2 <- relevel(factor(t12$week), ref="2")
mt12.2 <- lm(log10(vcps) ~ factor(week2), dat=t12)
summary(mt12.2)

mfull3 <- lmer(log10(vcps)~factor(treat3) + (1|week), data=dat.t)

mt1 <- lm(log10(vcps) ~ factor(week), dat=t1)
summary(mt1)

mt2 <- lm(log10(vcps) ~ factor(week), dat=t2)
summary(mt2)


mc <- lm(log10(vcps) ~ factor(week), dat=c)
summary(mc)



