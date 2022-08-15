#####################
## Lattrell Stats ###
#####################
setwd("C:/Users/Patrick/Desktop/Uni/Bachelorarbeit")
dat <- read.table("Stattable_final - Kopie.txt", header = TRUE, dec=",")
Stattable_final_Kopie <- read.table("Stattable_final - Kopie.txt", header = TRUE, dec=",")
str(Stattable_final_Kopie)

t1_2_C<-subset(Stattable_final_Kopie, treat !=4)
show(t1_2_C)
t1_2t<-subset(t1_2_C, treat !=3)
t1_C<-subset(t1_2_C, treat !=2)
tC<-subset(t1_C, treat !=1)






#### NEW SD calculations DS ####

setwd("H:/DS-Bees/R/lattrell")                                    # set your own WD and read in the table with the added variable log_vcps which is the already
dat_new <- read.table("Stattable_final.txt", header = TRUE)       # log10 transformed number of viral copies  per sample

names(dat_new)
str(dat_new)
dat_new$treat <- as.factor(dat_new$treat)

###treatment 1

#all weeks
t1 <- subset(dat_new, treat == "1")
min(t1$vcps)
max(t1$vcps)
mean(t1$log_vcps)
sd(t1$log_vcps)        

#week 1
t1w1 <- subset(t1, week == "1")
min(t1w1$vcps)
max(t1w1$vcps)
mean(t1w1$log_vcps)
sd(t1w1$log_vcps)





#### Your old file ####







#
logmeanC= log10(mean(tC$copies.ant))
logmeanC

# mean + sd of t1
t1<-subset(t1_2t, treat !=2)
logmean1=log10(mean(t1$copies.ant))
logmean1
sd1= log10(sd(t1$copies.ant))
sd1
#
t1_123<-subset(t1, week !=4)
t1_12<-subset(t1_123, week !=3)
# mean + sd of week 1 t1
t1_1<-subset(t1_12, week !=2)
logmean1_1=log10(mean(t1_1$copies.ant))
logmean1_1
sd1_1= sd(t1_1$copies.ant)
log10(sd1_1)
# mean + sd of week 2 t1
t1_2<-subset(t1_12, week !=1)
logmean1_2=log10(mean(t1_2$copies.ant))
logmean1_2
sd1_2= log10(sd(t1_2$copies.ant))
sd1_2
#
t1_234<-subset(t1, week !=1)
t1_34<-subset(t1_234, week !=2)
#mean + sd of week 3 t1
t1_3<-subset(t1_34, week !=4)
logmean1_3=log10(mean(t1_3$copies.ant))
logmean1_3
sd1_3= log10(sd(t1_3$copies.ant))
sd1_3
#mean + sd of week 4 t1
t1_4<-subset(t1_34,week !=3)
logmean1_4=log10(mean(t1_4$copies.ant))
logmean1_4
sd1_4= log10(sd(t1_4$copies.ant))
sd1_4
#mean + sd of t2
t2<-subset(t1_2t, treat !=1)
logmean2=log10(mean(t2$copies.ant))
logmean2
sd2= log10(sd(t2$copies.ant))
sd2
#
t2_123<-subset(t2, week !=4)
t2_12<-subset(t2_123, week !=3)
# mean + sd of week 1 t2
t2_1<-subset(t2_12, week !=2)
logmean2_1=log10(mean(t2_1$copies.ant))
logmean2_1
sd2_1= log10(sd(t2_2$copies.ant))
sd2_1
# mean + sd of week 2 t2
t2_2<-subset(t2_12, week !=1)
logmean2_2=log10(mean(t2_2$copies.ant))
logmean2_2
sd2_2= log10(sd(t2_2$copies.ant))
sd2_2
#
t2_234<-subset(t2, week !=1)
t2_34<-subset(t2_234, week !=2)
# mean + sd of week 3 t2
t2_3<-subset(t2_34, week !=4)
logmean2_3=log10(mean(t2_3$copies.ant))
logmean2_3
sd2_3= log10(sd(t2_3$copies.ant))
sd2_3
# mean + sd of week 4 t2
t2_4<-subset(t2_34,week !=3)
logmean2_4=log10(mean(t2_4$`copies/ant`))
logmean2_4
sd2_4= log10(sd(t2_4$`copies/ant`))
sd2_4

# mean + sd of t4 (bee puppae feed to the ants)
t4_234<-subset(Stattable_final_Kopie, treat !=1)
t4_34<-subset(t4_234, treat !=2)
t4<-subset(t4_34, treat !=3)
logmean4=log10(mean(t4$`copies/ant`))
logmean4
sd4=sd(t4$`copies/ant`)
logsd4=log10(sd4)
logsd4
vcps4=t4$`copies/ant`
log10(sd(vcps4))

setwd("C:/Users/Patrick/Desktop/Uni/Bachelorarbeit")

dat <- read.table("Stattable_final - Kopie.txt", header = TRUE, dec = ",")
names(dat)

show(dat)
#load packages
library("nlme")           # for linear models
library("lme4")           # linear mixed effect models
library("ggplot2")        # for nice graphs
library("plotrix")

# sample = number of pcr probe, vcps = viral copy per sample, name = lattrells sample names, week = time of sampling, id = colony or pupae number

#create a new logtransformed vcps variable
dat$log_copies.ant <- log10(dat$copies.ant) 
hist(dat$log_copies.ant)
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

dat.plot <- read.table("plotTab_1.txt", header = TRUE) 
dat.plot$treat <- as.factor(dat.plot$treat)

dat.plot$week <- c(1,1,1,2,2,2,3,3,3,4,4,4)
dat.plot$treat <- c(1,2,3,1,2,3,1,2,3,1,2,3)
###
# Create graph with the predicted values for the weekly models
###

dat.plot <- read.table("plottablog10.txt", header = TRUE) 
dat.plot$treat <- as.factor(dat.plot$treat)

names(dat.plot)

pt1 <- subset(dat.plot, dat.plot$treat == 1)
pt2 <- subset(dat.plot, dat.plot$treat == 2)
pt3 <- subset(dat.plot, dat.plot$treat == 3)

attach(dat.plot)

plotCI(week, fit, ui = Upr, li = lwr)

plotCI(dat.plot$week, dat.plot$fit, ui = dat.plot$Upr, li = dat.plot$lwr, xlab = "week", ylab="log [predicted genomic DWV copies / ant] ", 
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

#treatment 2: 
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





