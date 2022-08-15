### ANT COUNT CORRECTION  ####

### TEST ###

setwd("H:/DS-Bees/R/ANV")
dat <- read.table("test.txt", header = TRUE)

dat
dat$t <- c((rep("C", 10)), rep("L", 10), rep("H", 11))
dat$t <- as.factor(dat$t)
dat$t <- factor(dat$t, levels = c("C", "L", "H"))


dat$diff <- dat$real - dat$estimate
dat$diff


### easy way -> defining a proportion based on real and count data and apply the factor ### 

dat$factor <- dat$real/dat$estimate

dat$corrected21 <- dat$count1 * dat$factor
dat$corrected22 <- dat$count2 * dat$factor
dat$corrected23 <- dat$count3 * dat$factor
dat$corrected24 <- dat$estimate * dat$factor


allcounts <- c(dat$count1,dat$count2,dat$count3,dat$estimate)
allcounts
allcorect <- c(dat$corrected21, dat$corrected22, dat$corrected23, dat$corrected24)

allcorect

plot(allcounts, ylim=c(0, 420), col="blue", pch=19)
points(allcorect, col= "red", pch=1)
points(c(94:124), y= dat$real,  , col= "red", pch=19)



### try to do it with a model ###
# problem with lm: the intercept through zero is not a legitable option... thus a quadratic model would be 

fit <- lm(dat$diff~dat$real)
summary(fit)
fit$coefficients

dat$corrected11 <- -13.0456161 + dat$count1*0.4503032 
dat$corrected12 <- -13.0456161 + dat$count2*0.4503032 
dat$corrected13 <- -13.0456161 + dat$count3*0.4503032 
dat$corrected14 <- -13.0456161 + dat$real*0.4503032 

dat$corrected11
plot(dat$corrected11, ylim=c(-13,2))
points(dat$count1)


plot(dat$corrected12, ylim=c(-13,29))
points(dat$count2)

max(dat$count2)

dat$corrected13
dat$corrected14
  
  
plot(dat$diff~dat$real)
abline(fit)

dat$corrected1 <- predict(fit,dat$count1)


# try with a quadratic model: 

dat$real2 <- dat$real^2
qfit <- lm(dat$diff ~ dat$real+dat$real2)
predicted <- predict(qfit,list(dat$real,
                               dat$real2) )

lines(dat$real,predicted)







