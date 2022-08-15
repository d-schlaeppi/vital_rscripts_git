#####################
##       ANV      ###
#####################

getwd()
setwd("H:/DS-Bees/R/ANV")

dat <- read.table("ANV_Y2_long.txt", header = TRUE)
dat <- read.table("first_look2.txt", header = TRUE)
names(dat)

dat$treat1 <- as.factor(dat$treat1)
dat$treat2 <- as.factor(dat$treat2)
dat$id <- as.factor(dat$id)

dat$workforce <- dat$pupae+dat$adults

lm <- dat

lm <- subset(dat, measurement == 21) #last measurement
complete.cases(lm)  #see which rows have NAs in it or are complete 
lm$eggs <- NULL  #remove an entire column
lm$larva <- NULL
lm <- na.omit(lm)
lm$treat1

control <- subset(lm, treat1=="control")
low <- subset(lm, treat1=="low")
high <- subset(lm, treat1=="high")


plot(lm$treat1, lm$adults)
plot(lm$treat1, lm$pupae)
plot(lm$treat1, lm$workforce)

m1 <- lm(adults ~ factor(treat1), dat=lm)
summary(m1)
qqnorm(resid(m1))  
qqline(resid(m1))

lm$treat1 <- relevel(lm$treat1, ref="low")


plot(lm$treat2, lm$adults)
plot(lm$treat2, lm$pupae)
plot(lm$treat2, lm$workforce)

plot(control$treat2, control$adults)
plot(low$treat2, low$adults)
plot(high$treat2, high$adults)

plot(control$treat2, control$pupae)
plot(low$treat2, low$pupae)
plot(high$treat2, high$pupae)





