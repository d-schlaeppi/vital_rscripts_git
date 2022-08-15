#####################
## Lattrell Stats ###
#####################

getwd()
setwd("H:/DS-Bees/R/lattrell")

dat <- read.table("stattable.txt", header = TRUE)
names(dat)

#load packages
library("nlme")           # for linear models
library("lme4")           # linear mixed effect models
library("ggplot2")        # for nice graphs
library("plotrix")

# sample = number of pcr probe, vcps = viral copy per sample, name = lattrells sample names, week = time of sampling, id = colony or pupae number

#create a new logtransformed vcps variable
dat$log_vcps <- log10(dat$vcps) 
hist(dat$log_vcps)
dat$treat <- as.factor(dat$treat)  #treatment as factor


#create subsets
c   <- subset(dat, treat==3) 
t1  <- subset(dat, treat==1)
t2  <- subset(dat, treat==2) 
t12 <- subset(dat, treat==1 | treat==2) 
t4  <- subset(dat, treat==4)
a   <- subset(dat, treat!=4)

#viral copies in the treatments
quantile(c$vcps)
quantile(t1$vcps)
quantile(t2$vcps)
quantile(t4$vcps)


###MODELING###



###Create a full model to predict viral copies per sample --> how to read the output of this model? 

mfull <- lmer(log_vcps ~ treat + week + treat*week + (1|id), data=a)
summary(mfull)

qqnorm(resid(mfull))
qqline(resid(mfull))



### individual analyses as we are not shure how to read and interpret the full model
#DIFFERENCES FOR EACH WEEK: 

w1<- subset(a, week==1)
w2<- subset(a, week==2)
w3<- subset(a, week==3)
w4<- subset(a, week==4)





### Differences of the tree treatments with regard of viral copies during each weak ###

### week 1 ###
m1 <- lm(log_vcps ~ factor(treat), dat=w1)
summary(m1)

qqnorm(resid(m1))   #visual check for normality of the model residues --> despite two ouliers the residues fit a normal distribution quite nicely!
qqline(resid(m1))

w1$treat1 <- relevel(factor(w1$treat), ref="3")
m1.2 <- lm(log_vcps ~ factor(treat1), dat=w1)
summary(m1.2)

### week 2 ###

m2 <- lm(log_vcps ~ factor(treat), dat=w2)
summary(m2)
qqnorm(resid(m2))
qqline(resid(m2)) #one outlier rest is fine --> assumtion of normal distribution ok. 

w2$treat1 <- relevel(factor(w2$treat), ref="3")
m2.1 <- lm(log10(vcps)~factor(treat1), dat=w2)
summary(m2.1)

### week 3 ###

m3 <- lm(log_vcps ~ factor(treat), dat=w3)
summary(m3)
qqnorm(resid(m3))
qqline(resid(m3)) 

w3$treat1 <- relevel(factor(w3$treat), ref="3")
m3.1 <- lm(log10(vcps)~factor(treat1), dat=w3)
summary(m3.1)

### week 4 ###

m4 <- lm(log_vcps ~ factor(treat), dat=w4)
summary(m4)
qqnorm(resid(m4))
qqline(resid(m4)) 

w4$treat1 <- relevel(factor(w4$treat), ref="3")
m4.1 <- lm(log10(vcps)~factor(treat1), dat=w4)
summary(m4.1)


#predicted model values for a graph: 
predict.lm(m1, se.fit=T, interval="confidence")
predict.lm(m2, se.fit=T, interval="confidence")
predict.lm(m3, se.fit=T, interval="confidence")
predict.lm(m4, se.fit=T, interval="confidence")

#values are saved in an separate excel to make a graph of the predicted values: table

dat.plot <- read.table("plottable.txt", header = TRUE) 
dat.plot$treat <- as.factor(dat.plot$treat) 

###
# Create graph with the predicted values for the weekly models
###

dat.plot <- read.table("plottablog10.txt", header = TRUE) 
dat.plot$treat <- as.factor(dat.plot$treat)

pt1 <- subset(dat.plot, dat.plot$treat == 1)
pt2 <- subset(dat.plot, dat.plot$treat == 2)
pt3 <- subset(dat.plot, dat.plot$treat == 3)

plotCI(dat.plot$week, dat.plot$logvcps, ui = dat.plot$upr, li = dat.plot$lwr, xlab = "week", ylab="log [predicted genomic DWV copies / ant] ", 
       cex.lab = 2, cex.axis= 2, cex.main=4, pch=19, col=c("red", "blue", "black"), lwd = 2, xaxt="n")
axis(side=1, at=seq(9,12,by=1), cex.axis=2)
lines(pt1$week, pt1$logvcps, add=TRUE, col="red", lwd = 3 )
lines(pt2$week, pt2$logvcps, add=TRUE, col="blue", lwd = 3)
lines(pt3$week, pt3$logvcps, add=TRUE, col="black", lwd = 3)
legend(11,10.5, c("treatment 1", "treatment 2", "control"), box.lwd=0, bty='n', bg="transparent", pch=19, col=c("red", "blue", "black"), cex=2)
text(9,6, c("***"), cex=2)
text(10,6, c("***"), cex=2)
text(11,6, c("***"), cex=2)
text(12,6, c("***"), cex=2)

###
### Is there a change of viral titers over time in the two treatments?
###

#treatment 1: 
m.time1 <- lm(log_vcps ~ factor(week), dat=t1)
summary(m.time1)
qqnorm(resid(m.time1))
qqline(resid(m.time1)) 

t1$week1 <- relevel(factor(t1$week), ref="2")
m.time1.1 <- lm(log10(vcps)~factor(week1), dat=t1)
summary(m.time1.1)

t1$week2 <- relevel(factor(t1$week), ref="3")
m.time1.2 <- lm(log10(vcps)~factor(week2), dat=t1)
summary(m.time1.2)

#treatment 1: 
m.time2 <- lm(log_vcps ~ factor(week), dat=t2)
summary(m.time2)
qqnorm(resid(m.time2))
qqline(resid(m.time2)) 

t2$week1 <- relevel(factor(t2$week), ref="2")
m.time2.1 <- lm(log10(vcps)~factor(week1), dat=t2)
summary(m.time2.1)

t2$week2 <- relevel(factor(t2$week), ref="3")
m.time2.2 <- lm(log10(vcps)~factor(week2), dat=t2)
summary(m.time2.2)





