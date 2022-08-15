########################################################################
###How to Plot Survival Curves (Kaplan-Meier survival curve analysis)###
########################################################################

#factor is treatment which has 6 levels
#levels: Not Injected, PBS, solution1, solution2, solution3, solution4
#survivorship object is Surv(time to death in days, censoring indicator(=1 if event [death] occured,) # 0 = still alive)



dat <- read.table(file.choose(), header = T)
attach(dat)

dat.run1 <- subset(dat, run==1 & deathday!=0)
dat.run2 <- subset(dat, run==2 & deathday!=0)
  
detach(dat)

#get a basic feel for the data throug summary stats and histogram
#balanced between subject desing

table(dat$treatment)

#mean survival time accross the treatments
tapply(deathday[status==1], treatment[status==1], mean)


#load package 
library(survival) # contains kalan meier plot function

plot(survfit(Surv(dat$deathday, dat$status)~1)) #plots first basic plot with 95 confidence intervalls for all the data

fit=survfit(Surv(dat$deathday, dat$status)~dat$treatment)
plot(fit) 

# add labels and colours... looks already qite nice
plot(fit, col=c(1:6), xlab="days after injection", ylab="survival probability", lwd = 2)
legend(0.4, 0.4, col=(1:6), lwd=2)


#get the confidence intervalls out of the model (6 blocks, one for each treatment)
summary(fit)

#median time of death 
fit

#are there differences between the survival curves respectively effects of the treatments (log-rank test most popular  rho=0, generalize wilcoxon test rho=1)
#tests: Ho - no differences in treatments, H1 - difference between at least one group
#if p-value is bigger than 0.05 we do not reject the null hypothesis, "if p is low Null must GO" -> H1(Arbeitshypothese) is true

survdiff(Surv(dat$death,dat$status)~dat$treatment, rho=0)
survdiff(Surv(dat$death,dat$status)~dat$treatment, rho=1)

#To compare the means use an ANOVA, but survival function gives much more information! 



### first run
  
plot(survfit(Surv(dat.run1$deathday, dat.run1$status)~1)) 

fit=survfit(Surv(dat.run1$deathday, dat.run1$status)~dat.run1$treatment)
plot(fit) 

plot(fit, col=c(1:6), xlab="days after injection", ylab="survival probability", lwd = 2, main="Run 1")
legend(0.4, 0.4, c("NI", "PBS", "S1", "S2", "S3", "S4"),  col=(1:6), lwd=2)


### second run: 

plot(survfit(Surv(dat.run2$deathday, dat.run2$status)~1)) 

fit2=survfit(Surv(dat.run2$deathday, dat.run2$status)~dat.run2$treatment)
plot(fit2) 

plot(fit2, col=c(1:8), xlab="days after injection", ylab="survival probability", lwd = 2, main="Run 2")
legend(0.4, 0.4, c("NI 2", "PBS 2", "BS1.2", "S1.2"),  col=(1:6), lwd=2)
summary(fit)





