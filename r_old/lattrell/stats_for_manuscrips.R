#####################
## Lattrell Stats ###
#####################


#### Final stats ####

#load packages
library(lme4) #lmer
library(car) #contains levene's test
library(blmeco) #contains compareqqnorm  (multiple qq boxplots)
library(emmeans) #lsmeans to compare pairwise differences after the modeling using bonferroni corrections
library(arm) #contains the sim function for bayesian inference

setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/lattrell")
setwd("/Users/gismo/Desktop/R/lattrell") # homeoffice mac
dat <- read.table("Stattable_final.txt", header = TRUE)   

names(dat)

# sample = numbering during, vcps = viral copy per sample, name = lattrells sample names, week = time of sampling (1-4 = week 13-16), id = colony or pupae number, treat = treatment
head(dat)
str(dat)
dat$treat <- as.factor(dat$treat)
levels(dat$treat) <- c("Treatment 1", "Treatment 2", "Control", "Feeding Pupae")
dat$treat
table(dat$treat)

dat$log_vcps <- dat$log_vcps_new


#Descriptive statisti(positive control)
#titre detection limit

tapply(dat$log_vcps, dat$treat, shapiro.test) #all non signigicant --> Nullhypothesis of normaldistribution can be assumed 
tapply(dat$log_vcps, dat$treat, mean, na.rm=T)
tapply(dat$log_vcps, dat$treat, sd, na.rm=T)
tapply(dat$log_vcps, dat$treat, min, na.rm=T)
tapply(dat$log_vcps, dat$treat, max, na.rm=T)

#create subsets: 
c <- subset(dat, treat=="Control")
t1 <- subset(dat, treat=="Treatment 1") 
t2 <- subset(dat, treat=="Treatment 2")
t12 <- subset(dat, treat=="Treatment 1" | treat=="Treatment 2")

boxplot(c$log_vcps_new ~ c$week)
abline(h, 0)
boxplot(t1$log_vcps_new ~ t1$week)
boxplot(t2$log_vcps_new ~ t2$week)


#more descriptive stats 
tapply(c$log_vcps, c$week, mean, na.rm=T)
tapply(c$log_vcps, c$week, sd, na.rm=T)
tapply(t1$log_vcps, t1$week, mean, na.rm=T)
tapply(t1$log_vcps, t1$week, sd, na.rm=T)
tapply(t2$log_vcps, t2$week, mean, na.rm=T)
tapply(t2$log_vcps, t2$week, sd, na.rm=T)


#exculde feeding pupae from further analyses, because they are just a positive control 
dat <- subset(dat, treat != "Feeding Pupae")
dat$treat <- droplevels(dat)$treat
dat$id <- droplevels(dat)$id
dat$id = factor(dat$id,levels(dat$id)[c(1, 8:15, 2:7)])
levels(dat$id)


#create a linear mixed effect model 
baseline_M <- lmer(log_vcps ~ (1|id), data = dat, REML = FALSE)
baseline_M <- lmer(log_vcps_new ~ (1|id), data = dat, REML = FALSE)

week_M <- update(baseline_M, .~. + week)
treat_M <- update(week_M, .~. + treat)
interaction_M <- update(treat_M, .~. + week:treat)
anova(baseline_M, week_M, treat_M, interaction_M) 


#check model assumptions
mod <- interaction_M
plot(mod)
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])

#model output
marginal <- lsmeans(mod,  ~ week + treat + week:treat , adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD
lsmeans(mod, pairwise~treat*week, adjust="bonferroni")
summary(mod)



# predict values form model
mod <- lmer(log_vcps ~ week + treat + week:treat + (1|id), data=dat, REML=FALSE)
mod <- lmer(log_vcps_new ~ week + treat + week:treat + (1|id), data=dat, REML=FALSE)

#bayesian inference using sim
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
str(bsim)
apply(bsim@fixef, 2,quantile, prob=c(0.025,0.5,0.975))

#fitted values with 95% credible intervals
newdat<-expand.grid(
  week=seq(1:4),
  treat=factor(c('Treatment 1','Treatment 2','Control'), levels=levels(dat$treat))
  )
head(newdat)
Xmat <- model.matrix(~ week + treat + week:treat, data=newdat)      ####second part of model without of random factors

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))

for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,] 

newdat$lower <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upper <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- Xmat%*%fixef(mod)
newdat


#create a nice graph  
par(mfrow=c(1,1))
col <- c("gray70", "gray20", "black")
pt1 <- subset(newdat, newdat$treat == "Treatment 1")
pt2 <- subset(newdat, newdat$treat == "Treatment 2")
pt3 <- subset(newdat, newdat$treat == "Control")

library("gplots")
plotCI(newdat$week, newdat$fit, ui = newdat$upper, li = newdat$lower, xlab = "Time [Week]", ylab = "Predicted genomic copies of DWV [log]", 
       las = 1, col = c("white"), xaxt = "n")
legend(x=3.2, y=10.2, box.lty=0, legend = c("Treatment 1", "Treatment 2", "Control"), pch=c(1,2,19), col="black",  bg="transparent" )
lim <- par("usr")

h <- log10(4036.127537)
abline(h,0, lty=2)

#grey box
#rect(0, lim[3], lim[4], h, border = "grey90", col = "grey90")
#box()
axis(1, at = c(1,2,3,4), labels=c("13","14","15","16"))

#plot posteriors
points(c(1:4), pt1$fit, pch = 1)
segments(c(1:4), pt1$lower, c(1:4), pt1$upper)
points(c(1:4), pt2$fit, pch = 2)
segments(c(1:4), pt2$lower, c(1:4), pt2$upper)
points(c(1:4), pt3$fit, pch = 19)
segments(c(1:4), pt3$lower, c(1:4), pt3$upper)

lines(pt1$week, pt1$fit, add=TRUE, col="black")
lines(pt2$week, pt2$fit, add=TRUE, col="black")
lines(pt3$week, pt3$fit, add=TRUE, col="black")


# Alternative graph with real data instead of predicted value

head(dat)
pointcolor<-c("gray20","grey45","grey70")
linecolor<-c("grey25","grey60","grey80")
label <- c("Treatment 1", "Treatment 2","Control")


SD_per_treat <- tapply(dat$log_vcps_new, list(dat$week,dat$treat), sd)
SD_vector <- c(as.numeric(SD_per_treat[,1]), as.numeric(SD_per_treat[,2]), as.numeric(SD_per_treat[,3])) 
mean_per_treat <- tapply(dat$log_vcps_new, list(dat$week,dat$treat), mean, na.rm=T)
week <- c(13:16)

plot(week, mean_per_treat[,1], las = 1, xlab = "Time [week]", ylab = "Genomic copies of DWV [log]", pch = 0, xaxt = "n",
     ylim=c(3,10), cex.lab = 1.5, cex.axis = 1.5,  bty="n", cex = 2)
box(lwd = 2)
axis(1, c(13:16), cex.axis = 1.5)
points(week, mean_per_treat[,2], pch=1, cex = 2)
points(week, mean_per_treat[,3], pch=2, cex = 2)
legend(15.3, 10.3, legend = c("Treatment 1", "Treatment 2", "Control"), cex = 1.5, bg = "transparent", pch = c(0:2), box.lty=0)
for(i in 1:3){
  arrows(x0=week,x1=week,
         y0=mean_per_treat[,i]-SD_per_treat[,i],y1=mean_per_treat[,i]+SD_per_treat[,i],
         angle=90,length=0.05,code=3, lwd = 2)}
for(i in 1:2){lines(week, mean_per_treat[,i], lwd = 3, lty = i) }
lines(week, mean_per_treat[,3], lwd = 3, lty = 4)
h <- log10(4036.127537)
abline(h,0, lty=3, lwd = 2)


#citations

citation()
RStudio.Version()
citation("lme4")
citation("survival")
citation("dunn.test")
citation("arm")



















