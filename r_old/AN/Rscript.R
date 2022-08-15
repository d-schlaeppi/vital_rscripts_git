#####################
##       ANV      ###
#####################

getwd()
setwd("H:/DS-Bees/R/ANV")

dat <- read.table("first_look.txt", header = TRUE)
names(dat)

dat$workforce <- dat$pupae+dat$adults


plot(dat$neonic, dat$adults)
plot(dat$neonic, dat$pupae)
plot(dat$neonic, dat$workforce)

m1 <- lm(workforce ~ factor(neonic), dat=dat)
summary(m1)
qqnorm(resid(m1))  
qqline(resid(m1))


plot(dat$virus, dat$adults)
plot(dat$virus, dat$pupae)
plot(dat$virus, dat$workforce)

plot(dat$treat, dat$adults)
plot(dat$treat, dat$pupae)
plot(dat$treat, dat$workforce)


