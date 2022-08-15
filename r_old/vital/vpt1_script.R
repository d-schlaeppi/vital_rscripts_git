#### vpt1 ####

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lme4)) #lmer
suppressPackageStartupMessages(library(car)) #contains anova
suppressPackageStartupMessages(library(blmeco)) #contains compareqqnorm  (multiple qq boxplots)

#### prerequisites ####
setwd("/Users/gismo/Desktop/R/vital") # homeoffice mac
dat <- read.table("vpt1.txt", header = TRUE)
dat$colony <- as.factor(dat$colony)
head(dat)


# sum number of deaht per treatment
require(dplyr)
group_by(dat, colony) %>%
  summarise(
    sum_dead = sum(death), 
    sum_adults_all  = sum(death) + workers[day_since_start == 123]
  )

# controls less productive???

### graph - incidents over time 
mean_incidents_pg <- tapply(dat$incidents,list(dat$day_since_start, dat$treatment), mean, na.rm=TRUE)
sem_pg<-tapply(dat$incidents,list(dat$day_since_start,dat$treatment),sem)
sd_pg<-tapply(dat$incidents,list(dat$day_since_start,dat$treatment),sd)
sem_vector <- c(as.numeric(sem_pg[,1]),as.numeric(sem_pg[,2]))
sd_vector <- c(as.numeric(sd_pg[,1]),as.numeric(sd_pg[,2]))

day <- dat$day_since_start[dat$treatment == "control"]
pointcolor <- c("grey45","grey80")
linecolor  <- c("grey80","grey45")
label <- c("Virus", "Control")
sem<-function(x){sd(x,na.rm=T)/sqrt(length(na.omit(x)))}

plot(day, mean_incidents_pg[,2], pch=16, las = 1, xlab = "Time [d]", ylab = "Number of incidents", col = pointcolor[1],
     main = "Incidents", cex.lab = 1.2, cex.axis = 1.2, ylim=c(0,28))
points(day, mean_incidents_pg[,1], pch=16, col = pointcolor[2])
for(i in 1:2){lines(day, mean_incidents_pg[,i], col = linecolor[i])}
legend(2, 27, label, col = pointcolor, pch = 16, title = "Treatments", cex = 1.2)
for(i in 1:2){
  arrows(x0=day,x1=day,
         y0=mean_incidents_pg[,i]-sem_pg[,i],y1=mean_incidents_pg[,i]+sem_pg[,i],
         angle=90,length=0.025,code=3,col=linecolor[i])}

### graph - deaths over time + clinical symtoms over time
mean_symptoms_pg <- tapply(dat$sum_symptoms, list(dat$day_since_start, dat$treatment), mean, na.rm=TRUE)
sem_pg<-tapply(dat$sum_symptoms,list(dat$day_since_start, dat$treatment), sem)
sem_vector <- c(as.numeric(sem_pg[,1]),as.numeric(sem_pg[,2]))
plot(day, mean_symptoms_pg[,2], pch=16, las = 1, xlab = "Time [d]", ylab = "Number of workers showing CS", col = pointcolor[1],
     main = "Clinical Symptoms", cex.lab = 1.2, cex.axis = 1.2, ylim=c(0,22), xlim = c(60,124))
points(day, mean_symptoms_pg[,1], pch=16, col = pointcolor[2])
for(i in 1:2){lines(day, mean_symptoms_pg[,i], col = linecolor[i])}
legend(60, 21, label, col = pointcolor, pch = 16, title = "Treatments", cex = 1.2)
for(i in 1:2){
  arrows(x0=day,x1=day,
         y0=mean_symptoms_pg[,i]-sem_pg[,i],y1=mean_symptoms_pg[,i]+sem_pg[,i],
         angle=90,length=0.025,code=3,col=linecolor[i])}

mean_clinical_pg <- tapply(dat$death, list(dat$day_since_start, dat$treatment), mean, na.rm=TRUE)
sem_pg<-tapply(dat$death,list(dat$day_since_start, dat$treatment), sem)
sem_vector <- c(as.numeric(sem_pg[,1]),as.numeric(sem_pg[,2]))
plot(day, mean_death_pg[,2], pch=16, las = 1, xlab = "Time [d]", ylab = "Number of dead workers", col = pointcolor[1],
     main = "Dead workers", cex.lab = 1.2, cex.axis = 1.2, ylim=c(0,22), xlim = c(60,124))
points(day, mean_death_pg[,1], pch=16, col = pointcolor[2])
for(i in 1:2){lines(day, mean_death_pg[,i], col = linecolor[i])}
legend(60, 21, label, col = pointcolor, pch = 16, title = "Treatments", cex = 1.2)
for(i in 1:2){
  arrows(x0=day,x1=day,
         y0=mean_death_pg[,i]-sem_pg[,i],y1=mean_death_pg[,i]+sem_pg[,i],
         angle=90,length=0.025,code=3,col=linecolor[i])}

### histogramm of different clinical symptoms
vec <- c(sum(dat$slow), sum(dat$unresposive), sum(dat$posture_abnormal), sum(dat$moving_zigzag_limping), sum(dat$trembly_hinking),
         sum(dat$troubles), sum(dat$impaired_locomotion_pronounced), sum(dat$partial_paralysis), sum(dat$paralysis), sum(dat$newly_emerged_struggling))
names(vec) <- names(dat)[34:43]

par(mar = c(15.1, 4.1, 4.1, 2.1))
barplot(vec, ylab = "#", main = "clinical symptoms", las= 2)

?barplot

par(mar = c(5.1, 4.1, 4.1, 2.1))


### Overall difference between treatments 
boxplot(dat$incidents ~ dat$treatment, main = "incidents/measurement - treatment")
boxplot(dat$incidents ~ dat$colony, main = "incidents/measurement - colony")

# subset only starting from 63 when the first samples were collected and incidents started to rise. 
dat_short <- subset(dat, dat$day_since_start >= 63)
head(dat_short)

boxplot(dat_short$incidents ~ dat_short$treatment, main = "incidents/measurement - treatment")
boxplot(dat_short$incidents ~ dat_short$colony, main = "incidents/measurement - colony")

boxplot(dat_short$incidents ~ dat_short$treatment, xlab = "treatments", ylab = "nr of incidents per recording", main = "incidents")
boxplot(dat_short$death ~ dat_short$treatment, xlab = "treatments", ylab = "nr of dead ants per recording", main = "mortality")
boxplot(dat_short$sum_symptoms ~ dat_short$treatment, xlab = "treatments", ylab = "nr of ants displaying CS per recording", main = "clinical symptoms")

mod <- glmer(incidents ~ treatment + (1|colony) + (1|day_since_start), data = dat_short, family = "poisson")
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

# test model assumtions
compareqqnorm(mod)
leveneTest(residuals(mod) ~ dat_short$treatment) # non-significant --> homogenity of variance ok
boxplot(residuals(mod) ~ dat_short$treatment)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # normality of residues is not fullfilled... but we use poisson distribution anyways?
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

mod <- glmer(death ~ treatment + (1|colony) + (1|day_since_start), data = dat_short, family = "poisson")
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

mod <- glmer(sum_symptoms ~ treatment + (1|colony) + (1|day_since_start), data = dat_short, family = "poisson")
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

### colors ###
#subset without controls and only from the day the white ants were painted
colors <- subset(dat_short, day_since_start >= 84)
colors <- subset(colors, treatment == "virus")

# mortality per colony 
death_yellow <- tapply(colors$death_yellow, list(colors$day_since_start), sum)
death_blue <- tapply(colors$death_blue, list(colors$day_since_start), sum)
death_white <- tapply(colors$death_white, list(colors$day_since_start), sum) 
death_black <- tapply(colors$death_black, list(colors$day_since_start), sum)
df1 <- data.frame(matrix(unlist(death_yellow), nrow=length(death_yellow), byrow=TRUE),stringsAsFactors=FALSE)
df2 <- data.frame(matrix(unlist(death_blue), nrow=length(death_blue), byrow=TRUE),stringsAsFactors=FALSE)
df3 <- data.frame(matrix(unlist(death_white), nrow=length(death_white), byrow=TRUE),stringsAsFactors=FALSE)
df4 <- data.frame(matrix(unlist(death_black), nrow=length(death_black), byrow=TRUE),stringsAsFactors=FALSE)
clinic_yellow <- tapply(colors$clinic_yellow, list(colors$day_since_start), sum)
clinic_blue <- tapply(colors$clinic_blue, list(colors$day_since_start), sum)
clinic_white <- tapply(colors$clinic_white, list(colors$day_since_start), sum) 
clinic_black <- tapply(colors$clinic_black, list(colors$day_since_start), sum)
df5 <- data.frame(matrix(unlist(clinic_yellow), nrow=length(clinic_yellow), byrow=TRUE),stringsAsFactors=FALSE)
df6 <- data.frame(matrix(unlist(clinic_blue), nrow=length(clinic_blue), byrow=TRUE),stringsAsFactors=FALSE)
df7 <- data.frame(matrix(unlist(clinic_white), nrow=length(clinic_white), byrow=TRUE),stringsAsFactors=FALSE)
df8 <- data.frame(matrix(unlist(clinic_black), nrow=length(clinic_black), byrow=TRUE),stringsAsFactors=FALSE)


df <- data.frame(matrix(nrow = 4*nrow(df1), ncol = 4))
df$X1 <- rep(unique(colors$day_since_start), 4)
df$X2 <- factor(rep(c("yellow", "blue", "white", "black"), each = 19), levels = c("yellow", "blue", "white", "black"))
df$X3 <- c(df1[,1], df2[,1], df3[,1], df4[,1])
df$X4 <- c(df5[,1], df6[,1], df7[,1], df8[,1])
names(df) <- c("day", "color", "dead", "clinic")


boxplot(df$dead ~ df$color, main = "mortality", ylab = "number of dead ants per recording", xlab = "treatments")
boxplot(df$clinic ~ df$color, main = "clinical symptoms", ylab = "number of dead ants per recording", xlab = "treatments")

mod <- glm(dead ~ color, data = df, family = "poisson")
Anova(mod)
summary(mod)
lsmeans(mod, pairwise ~ color, adjust = "tukey")

mod <- glm(clinic ~ color, data = df, family = "poisson")
Anova(mod)
summary(mod)
lsmeans(mod, pairwise ~ color, adjust = "tukey")

data <- subset(dat_short, day == 123)

boxplot()
