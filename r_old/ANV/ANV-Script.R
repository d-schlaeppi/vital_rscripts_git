### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ANV - Ants, Neonics and Viruses - Interaction of stressors  ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### load libraries ####
suppressPackageStartupMessages(library(survival))# contains kalan meier plot function
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(chron))
suppressPackageStartupMessages(library(car)) #contains levene's test + Anova
#library(arm) #contains the sim function for bayesian inference
#library(dunn.test) #contains the dunn test
suppressPackageStartupMessages(library(nlme)) #lme
suppressPackageStartupMessages(library(emmeans)) #lsmeans to compare pairwise differences after the modeling using bonferroni corrections
suppressPackageStartupMessages(library(blmeco)) #contains compareqqnorm  (multiple qq boxplots)
suppressPackageStartupMessages(library(lme4)) #lmer
suppressPackageStartupMessages(library(lsmeans)) #cld
suppressPackageStartupMessages(library(multcomp)) #cld 2
suppressPackageStartupMessages(library(glmmTMB))
suppressPackageStartupMessages(library(multcompView)) #for the CLD functions
suppressPackageStartupMessages(library(lmtest)) #bptest(model) -  heteroscedasticity
suppressPackageStartupMessages(library(survminer)) #contains ggsurvplot()

#### prerequisites ####
# setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV") #office computer
setwd("/Users/gismo/Desktop/R/ANV") # homeoffice mac
dat <- read.table("ANV_all_Data.txt", header = TRUE)
mov <- read.table("speed.txt", header = TRUE)

# preparation of Dataframe --> dat 
head(dat)
str(dat)
dat$treatment_1 <- factor(dat$treatment_1, levels(dat$treatment_1)[c(1,3,2)])
dat$treatment <- as.factor(dat$treatment)
dat$treatment <- revalue(dat$treatment, c("1"="control-control", "2"="control-virus","3" = "low-control","4" = "low-virus","5" = "high-control","6" = "high-virus" ))
dat$adult_production <- as.factor(dat$adult_production)
dat$adult_production <- revalue(dat$adult_production, c("1" = "yes", "0" = "no"))



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Cumulative queen survival in percent ####

#mean survival time accross the treatments
tapply(dat$day_of_death[dat$y2_survival_status==1], dat$treatment[dat$y2_survival_status==1], mean)
plot(dat$day_of_death ~ dat$treatment)
plot(survfit(Surv(dat$day_of_death, dat$y2_survival_status)~1)) #plots first basic plot with 95 confidence intervalls for all the data
fit=survfit(Surv(dat$day_of_death, dat$y2_survival_status)~dat$treatment) 
plot(fit) 
summary(fit)
#plot cumulative queen survival
ggsurvplot(fit, data = dat, pval = TRUE) #modify graoh to put it into supplementary material
#calculate potential statistical difference between the groups
mod <- survdiff(Surv(day_of_death, y2_survival_status) ~ treatment, data=dat)
mod





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Number of Queens producing adults - Is there a sig. difference in the number of queens producing adults between the treatment? ####

plot(dat$adult_production ~ dat$treatment)
chisq.test(dat$adult_production, dat$treatment)

#generalised logistic mixed model!!!!!!!!!!! (ev ordinal family, Bernulli ja / nein 1/0
#logistic regression

mod <- glm(dat$adult_production ~ dat$treatment_1*dat$treatment_2, family = binomial)
summary(mod)
Anova(mod)

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1:treatment_2, data = dat, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# test model assumtions
leveneTest(mod) #non-significant = ok #homogenity of variance (plot(mod,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #non significant --> assumption of normally distributed residuals is ok ##plot(mod, 2)#normality
bptest(mod) #non significant --> assumtions of heteroskedasticity not violated ok 
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
par(mfrow=c(1,1))
summary(mod)

chisq.test(dat$adult_production, dat$treatment_1)
#m <- matrix(c(9,9,10,10,8,8,1,1,0,0,2,2), byrow = T, nrow = 2) #is the same as transforming it to a matrix 
#chisq.test(m)
# no significant effect ??? Is there the potential for another analysis because there might be a chance that the neonic treatment (1) has an effect on the production of adults. 







### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Bodymass of Queens and Workers ####


### workers ###
boxplot(dat$weight_adults ~ dat$treatment)

#stat summary:
require(dplyr)
group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(weight_adults)), 
    mean = mean(weight_adults, na.rm = TRUE),
    sd = sd(weight_adults, na.rm = TRUE),
    shapiro = shapiro.test(dat$weight_adults)$p,
    levene = leveneTest(dat$weight_adults, dat$treatment)$`Pr(>F)`[1],
    q_0 = quantile(weight_adults, probs = (0), na.rm = TRUE),
    q_25 = quantile(weight_adults, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(weight_adults, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(weight_adults, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(weight_adults, probs = (1), na.rm = TRUE),
  )

#model
mod <- lm(weight_adults ~ treatment_1*treatment_2, data = dat)
summary(mod)
Anova(mod)
#pairwise differences
lsmeans(mod, pairwise ~ treatment_1*treatment_2, data = dat, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1*treatment_2, data = dat)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# test model assumtions
leveneTest(mod) #non-significant = ok #homogenity of variance (plot(mod,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #non significant --> assumption of normally distributed residuals is ok ##plot(mod, 2)#normality
bptest(mod) #non significant --> assumtions of heteroskedasticity not violated ok 
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
par(mfrow=c(1,1))



### queens
boxplot(dat$weight_queens ~ dat$treatment)


#stat summary:
require(dplyr)
group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(weight_queens)), 
    mean = mean(weight_queens, na.rm = TRUE),
    sd = sd(weight_queens, na.rm = TRUE),
    shapiro = shapiro.test(dat$weight_queens)$p,
    levene = leveneTest(dat$weight_queens, dat$treatment)$`Pr(>F)`[1],
    q_0 = quantile(weight_queens, probs = (0), na.rm = TRUE),
    q_25 = quantile(weight_queens, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(weight_queens, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(weight_queens, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(weight_queens, probs = (1), na.rm = TRUE),
  )


#model
mod <- lm(weight_queens ~ treatment_1*treatment_2, data = dat)
summary(mod)
Anova(mod)

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1:treatment_2, data = dat, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# test model assumtions
leveneTest(mod) #non-significant = ok #homogenity of variance (plot(mod,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #non significant --> assumption of normally distributed residuals is ok ##plot(mod, 2)#normality
bptest(mod) #non significant --> assumtions of heteroskedasticity not violated ok 
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
par(mfrow=c(1,1))





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Colonysize after year I ####


#Y1 only for treatment 1 because 2 not yet established, y2 for all 6 treatments
# no need for random factor beceause there are no repeated measures and conditions the same for all apart from the treatment

boxplot(dat$y1_adults ~ dat$treatment_1)
hist(dat$y1_adults)
shapiro.test(dat$y1_adults) #significant -> null hypothesis that data is normally distributed must be rejected
leveneTest(dat$y1_adults, dat$treatment_1) #non significant --> nullhypothesis of equal variances is not rejected
kruskal.test(dat$y1_adults ~ dat$treatment_1, na.action = na.omit)

boxplot(dat$y1_eggs ~ dat$treatment_1) #the one big outlier is from colony c8v which had a delayed development
shapiro.test(dat$y1_eggs) #significant -> null hypothesis that data is normally distributed must be rejected
leveneTest(dat$y1_eggs, dat$treatment_1) #non significant --> nullhypothesis of equal variances is not rejected
kruskal.test(dat$y1_eggs ~ dat$treatment_1, na.action = na.omit)

boxplot(dat$y1_larva ~ dat$treatment_1)
shapiro.test(dat$y1_larva) #significant -> null hypothesis that data is normally distributed must be rejected
leveneTest(dat$y1_larva, dat$treatment_1) #non significant --> nullhypothesis of equal variances is not rejected
kruskal.test(dat$y1_larva ~ dat$treatment_1, na.action = na.omit)

boxplot(dat$y2_pupae ~ dat$treatment_1)
shapiro.test(dat$y2_pupae) #significant -> null hypothesis that data is normally distributed must be rejected
leveneTest(dat$y2_pupae, dat$treatment_1) #non significant --> nullhypothesis of equal variances is not rejected
kruskal.test(dat$y2_pupae ~ dat$treatment_1, na.action = na.omit)



#no differences between the three treatments y1


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Colony after year II - adults and brood ####


### adults ###
boxplot(dat$y_adults ~ dat$treatment)
mod <- lm(y_adults ~ treatment_1*treatment_2, data = dat)
Anova(mod)
summary(mod)

#check model assumptions
plot(mod)
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
par(mfrow=c(1,1))
#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")


### eggs ###
boxplot(dat$y2_eggs ~ dat$treatment)
mod <- lm(y2_eggs ~ treatment_1*treatment_2, data = dat)
Anova(mod)
summary(mod)

#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
par(mfrow=c(1,1))


### larva ###
boxplot(dat$y2_larva ~ dat$treatment)
mod <- lm(y2_larva ~ treatment_1*treatment_2, data = dat)
Anova(mod)
summary(mod)

#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
par(mfrow=c(1,1))

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")


### pupae ###
boxplot(dat$y2_pupae ~ dat$treatment)
dat$y2_pupae_trans <- log10(dat$y2_pupae+0.1)
boxplot(dat$y2_pupae_trans ~ dat$treatment)

pupae_anova <- aov(y2_pupae_trans ~ treatment_1*treatment_2, data = dat)
summary(pupae_anova)
#multiple comparisons
TukeyHSD(pupae_anova)

#test anova validity: 
#homogenity of variance 
plot(pupae_anova, 1)
leveneTest(y2_pupae ~ treatment_1*treatment_2, data = dat) # non significant and thus ok. 
#normality
plot(pupae_anova, 2)
aov_residuals <- residuals(object = pupae_anova)
shapiro.test(x = aov_residuals) #non significant --> assumption of normally distributed residuals is ok for log transformed data 

#calculate summary statistics
require(dplyr)
group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    mean = mean(y2_pupae, na.rm = TRUE),
    sd = sd(y2_pupae, na.rm = TRUE), 
    non_NA_count = sum(!is.na(y2_pupae))
  )






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### behavioural assay ####

head(mov) 
#resposne variables - active_time, initial_speed, average_speed, overall_movement
plot(mov$drops~mov$treatment)
plot(mov$number_of_stops~mov$treatment)
plot(mov$active_time~mov$treatment)
plot(mov$overall_movement~mov$treatment)
plot(mov$average_speed~mov$treatment)
plot(mov$initial_speed~mov$treatment)

str(mov)
mov$treatment_1 <- factor(mov$treatment_1, levels(mov$treatment_1)[c(1,3,2)])
mov$treatment <- as.factor(mov$treatment)
mov$treatment <- revalue(mov$treatment, c("1"="control-control", "2"="control-virus","3" = "low-control","4" = "low-virus","5" = "high-control","6" = "high-virus" ))
mov$time <- times(mov$time)
mov$time_minutes <- 60*hours(mov$time) + minutes(mov$time)

### overall movement ###
par(mfrow = c(1,1))
boxplot(mov$overall_movement ~ mov$date) #day might have and effect and will be included in the model as random factor
plot(mov$overall_movement ~ mov$time_minutes)
abline(lm(mov$overall_movement ~ mov$time_minutes))              # time of day not relevant
hist(mov$overall_movement)
boxplot(mov$overall_movement ~ mov$treatment)

baseline <- lmer(overall_movement ~ (1|date), data=mov, REML=FALSE)
colony_M <- lmer(overall_movement ~ (1|id) + (1|date), data=mov, REML=FALSE)
treatment_1_M <- lmer(overall_movement ~ treatment_1 + (1|id) + (1|date), data=mov, REML=FALSE)
treatment_2_M <- lmer(overall_movement ~ treatment_1 + treatment_2 + (1|id) + (1|date), data=mov, REML=FALSE)
interaction_M <- lmer(overall_movement ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|date), data=mov, REML=FALSE)
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)


mod <- interaction_M
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

#check model assumptions
leveneTest(residuals(mod) ~ mov$treatment_1*mov$treatment_2)
boxplot(residuals(mod) ~ mov$treatment_1*mov$treatment_2)
plot(mod)

par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$date[,1])
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~mov$date)
compareqqnorm(mod)
par(mfrow=c(1,1))

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "bonferroni")

boxplot(mov$overall_movement ~ mov$treatment_2*mov$treatment_1)

### active_time - time_inactive ###
hist(mov$time_inactive)
mov$time_inactive <- 121-mov$active_time
hist(mov$time_inactive)
mov$time_inactive <- log10(mov$time_inactive)
boxplot(mov$time_inactive ~ mov$treatment)
hist(mov$time_inactive)

boxplot(mov$time_inactive ~ mov$date) #day might have and effect and will be included in the model as random factor
plot(mov$time_inactive ~ mov$time_minutes)
abline(lm(mov$time_inactive ~ mov$time_minutes)) # time of day seems not relevant and will be neglected in the model

baseline <- lmer(time_inactive ~ (1|date), data=mov, REML = FALSE)
colony_M <- lmer(time_inactive ~ (1|id) + (1|date), data=mov, REML = FALSE)
treatment_1_M <- lmer(time_inactive ~ treatment_1 + (1|id) + (1|date), data=mov, REML = FALSE)
treatment_2_M <- lmer(time_inactive ~ treatment_1 + treatment_2 + (1|id) + (1|date), data=mov, REML = FALSE)
interaction_M <- lmer(time_inactive ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|date), data=mov, REML = FALSE)
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)

mod <- interaction_M
mod

Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

#check model assumptions
leveneTest(residuals(mod) ~ mov$treatment_1*mov$treatment_2) #not significant --> homogeneity of variance can be assumed
boxplot(residuals(mod) ~ mov$treatment_1*mov$treatment_2)
plot(mod)

par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$date[,1])
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~mov$date)
compareqqnorm(mod)
par(mfrow=c(1,1))

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "bonferroni")



### average speed ###
boxplot(mov$average_speed ~ mov$date) #day might have and effect and will be included in the model as random factor
plot(mov$average_speed ~ mov$time_minutes)
abline(lm(mov$average_speed ~ mov$time_minutes))              # time of day not relevant
hist(mov$average_speed)
boxplot(mov$average_speed ~ mov$treatment)

baseline <- lmer(average_speed ~ (1|date), data=mov, REML=FALSE)
colony_M <- lmer(average_speed ~ (1|id) + (1|date), data=mov, REML=FALSE)
treatment_1_M <- lmer(average_speed ~ treatment_1 + (1|id) + (1|date), data=mov, REML=FALSE)
treatment_2_M <- lmer(average_speed ~ treatment_1 + treatment_2 + (1|id) + (1|date), data=mov, REML=FALSE)
interaction_M <- lmer(average_speed ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|date), data=mov, REML=FALSE)
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)

mod <- interaction_M
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

#check model assumptions
leveneTest(residuals(mod) ~ mov$treatment_1*mov$treatment_2)
boxplot(residuals(mod) ~ mov$treatment_1*mov$treatment_2)
plot(mod)

par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$date[,1])
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~mov$date)
compareqqnorm(mod)
par(mfrow=c(1,1))

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "bonferroni")



### initial speed ###
boxplot(mov$initial_speed ~ mov$date) #day might have and effect and will be included in the model as random factor
plot(mov$initial_speed ~ mov$time_minutes)
abline(lm(mov$initial_speed ~ mov$time_minutes))              # time of day not relevant
hist(mov$initial_speed)
boxplot(mov$initial_speed ~ mov$treatment)

baseline <- lmer(initial_speed ~ (1|date), data=mov, REML=FALSE)
colony_M <- lmer(initial_speed ~ (1|id) + (1|date), data=mov, REML=FALSE)
treatment_1_M <- lmer(initial_speed ~ treatment_1 + (1|id) + (1|date), data=mov, REML=FALSE)
treatment_2_M <- lmer(initial_speed ~ treatment_1 + treatment_2 + (1|id) + (1|date), data=mov, REML=FALSE)
interaction_M <- lmer(initial_speed ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|date), data=mov, REML=FALSE)
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)

mod <- interaction_M
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

#check model assumptions
leveneTest(residuals(mod) ~ mov$treatment_1*mov$treatment_2)
boxplot(residuals(mod) ~ mov$treatment_1*mov$treatment_2)
plot(mod)

par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$date[,1])
qqline(ranef(mod)$date[,1])
boxplot(resid(mod)~mov$date)
compareqqnorm(mod)
par(mfrow=c(1,1))
#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")







### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Neonicotinoid residues ####

dat$tw <- dat$tw_thiamethoxam_workers_dryconcentration
dat$cw <- dat$cw_clothianidin_workers_dryconcentration
dat$tq <- dat$tq_thiamethoxam_queens_dryconcentration
dat$cq <- dat$cq_clothianidin_queens_dryconcentration

# 1. Residue concentrations compared across the treatments and castes
#needed: treatment, treatment_1, treatment_2, tiamet + clotianidin for both castes, colony identity


#create virus titer per mg bodyweight to compare them between castes and species. 

caste_data <- data.frame(sample =c(1:60, 1:60), 
                         identity = rep(dat$identity, times = 2), 
                         treatment = rep(dat$treatment, times = 2), 
                         treatment_1 = rep(dat$treatment_1, times = 2), 
                         treatment_2 = rep(dat$treatment_2, times = 2), 
                         caste = c(rep("w", times = length(dat$identity)), rep("q", times = length(dat$identity))),
                         thiamethoxam = c(dat$tw, dat$tq), 
                         clothianidin = c(dat$cw, dat$cq), 
                         virus_titre = c(dat$virus_copies_per_worker, dat$virus_copies_per_queen)
)

par(mfrow = c(2,2))
boxplot(dat$tq ~ dat$treatment)
boxplot(dat$tw ~ dat$treatment)
boxplot(dat$cq ~ dat$treatment)
boxplot(dat$cw ~ dat$treatment)

#thiamethoxam
#reciprocal transformation
min(caste_data$thiamethoxam[caste_data$thiamethoxam != 0], na.rm = TRUE)
caste_data$thiamethoxam[caste_data$thiamethoxam == 0] <- 0.5*min(caste_data$thiamethoxam[caste_data$thiamethoxam != 0], na.rm = TRUE)
caste_data$thiamethoxam_rec <- (1/caste_data$thiamethoxam) #reciprocal transformation
caste_data$thiamethoxam_rec


#baseline <- lmer(thiamethoxam ~ (1|identity), data = caste_data, REML = FALSE) #model assumptions not fullfilled
baseline <- lmer(thiamethoxam_rec ~ (1|id), data = caste_data)
caste_M <- update(baseline, .~. + caste)
treatment_1_M <- update(caste_M, .~. + treatment_1)
treatment_2_M <- update(treatment_1_M, .~. + treatment_2)
interaction_M <- update(treatment_2_M, .~. + treatment_1:treatment_2 )
anova(baseline, caste_M, treatment_1_M, treatment_2_M, interaction_M)

mod <- interaction_M
summary(mod)
estimates <- 1/fixef(mod)
estimates

mod <- interaction_M
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(model))))  # homogeneity of variance
compareqqnorm(mod)






### 1.4 Thiamethoxam residues in workers TW ###

boxplot(dat$tw_thiamethoxam_workers_dryconcentration ~ dat$treatment)
#stat summary:
require(dplyr)
group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    mean = mean(tw_thiamethoxam_workers_dryconcentration, na.rm = TRUE),
    sd = sd(tw_thiamethoxam_workers_dryconcentration, na.rm = TRUE), 
    non_NA_count = sum(!is.na(tw_thiamethoxam_workers_dryconcentration))
  )

dat$tw_log <- log10(dat$tw+0.001)

group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    mean = mean(tw_log, na.rm = TRUE)
  )

### drop the controls, as they all were zero -> shows that there was no detecteble contamination prior to the experiment
# The difference between the treatments and the controls is measured with the treatment means of the model (even if controls are not in the model)
dat_reduced <- subset(dat, treatment_1 != "control")
dat_reduced$treatment_1 <- droplevels(dat_reduced)$treatment_1
dat_reduced$treatment <- droplevels(dat_reduced)$treatment
plot(dat_reduced$tw ~ dat_reduced$treatment)

#model with original data and exponential distribution or because model assumptions violated when analyzed normally or logtransformed and one outlier removed)
#alternative: non-parametric test for interaction? not a good option
#mod <- glm(tw+0.01 ~ treatment_1*treatment_2, data = dat_reduced, family = Gamma), summary(mod, dispersion=1) #to have an actual exponential distribution and not just one of many gamma distributions
mod <- lm(tw_log ~ treatment_1*treatment_2, data = dat_reduced)
mod <- lm(tw_log ~ treatment_1*treatment_2, data = dat)
Anova(mod)
summary(mod)

#check model assumptions
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
plot(mod)
par(mfrow=c(1,1))

#heteroskedastity 
library(lmtest)
bptest(mod) #non significant --> ook

#normality
plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #non significant --> assumption of normally distributed residuals is ok

#homogenity of variance 
plot(mod, 1)
leveneTest(mod) #significant = not ok? 

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey") #or TukeyHSD(aov(mod))
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD









### clothianidin residues in Workers CW ###

boxplot(dat$cw ~ dat$treatment)
#stat summary:
require(dplyr)
group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    mean = mean(cw, na.rm = TRUE),
    sd = sd(cw, na.rm = TRUE), 
    non_NA_count = sum(!is.na(cw))
  )

dat$cw_log <- log10(dat$cw+0.01)
plot(dat$cw_log ~ dat$treatment)

# drop the controls, as they all were zero
dat_reduced <- subset(dat, treatment_1 != "control")
dat_reduced$treatment_1 <- droplevels(dat_reduced)$treatment_1
dat_reduced$treatment <- droplevels(dat_reduced)$treatment
plot(dat_reduced$cw_log ~ dat_reduced$treatment)

mod <- lm(cw_log ~ treatment_1*treatment_2, data = dat_reduced)
Anova(mod)
summary(mod)

#check model assumptions
compareqqnorm(mod)
par(mfrow=c(2,2))
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
plot(mod)
par(mfrow=c(1,1))

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey") #or TukeyHSD(aov(mod))
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat_reduced)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

#modelassumtions
#homogenity of variance 
plot(mod, 1)
leveneTest(mod) #non-significant = ok 
#normality
plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # significant --> assumption of normally distributed residuals is violated


### thiamethoxam residues in Queens TQ ###
boxplot(dat$tq ~ dat$treatment)
#stat summary:
require(dplyr)
group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    mean = mean(tq, na.rm = TRUE),
    sd = sd(tq, na.rm = TRUE), 
    non_NA_count = sum(!is.na(tq)),
    non_zero_count = sum(!is.na(cq))-sum(cq==0, na.rm = TRUE)
  )



### clothianidin residues in Queens CQ ###
boxplot(dat$cq ~ dat$treatment)
#stat summary:
require(dplyr)
group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    mean = mean(cq, na.rm = TRUE),
    sd = sd(cq, na.rm = TRUE), 
    non_NA_count = sum(!is.na(cq)),
    non_zero_count = sum(!is.na(cq))-sum(cq==0, na.rm = TRUE)
  )






#comparison of the castes


#ratio comparison 








### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Virus titers ####

### Virus detection treshold
vdt_old <- 1325	#Virus detection treshold
vdt <- 5543
vdt_log <- log10(vdt)

#variable set up
head(dat)
dat$q_vc <- dat$virus_copies_per_queen
dat$q_vc_log <- log10(dat$q_vc)
#dat$q_vc_log_2 <- log10(dat$q_vc + 0.5*min(dat$q_vc, na.rm = TRUE))   # no longer necessary
dat$w_vc <- dat$virus_copies_per_worker
dat$w_vc_log <- log10(dat$w_vc) #+ 0.5*min(dat$w_vc[dat$w_vc != 0], na.rm = TRUE)) 

### look at data ###

#stat summary:
require(dplyr)
worker_data <- group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(w_vc)),
    positive_w = sum(!is.na(w_vc)) - sum(pos_neg_workers=="neg", na.rm = TRUE),
    negative_w = sum(pos_neg_workers=="neg", na.rm = TRUE),
    mean_titer_w = mean(w_vc, na.rm = TRUE),
    logmean_w = log10(mean_titer_w),
    sd_virustiter_w = sd(w_vc, na.rm = TRUE)
  )
worker_data

queen_data <- group_by(dat, treatment_1, treatment_2) %>%
      summarise(
        count = n(),
    N_Q = sum(!is.na(q_vc)),
    positive_q = sum(!is.na(q_vc)) - sum(pos_neg_queens=="neg", na.rm = TRUE),
    negative_q = sum(pos_neg_queens=="neg", na.rm = TRUE),
    mean_titer_q = mean(q_vc, na.rm = TRUE),
    logmean_q = log10(mean_titer_q),
    sd_titer_q = sd(q_vc, na.rm = TRUE),
  )
queen_data 

quantiles_workers <- group_by(dat, treatment_1, treatment_2) %>%
  summarise(
    q_0 = quantile(w_vc_log, probs = (0), na.rm = TRUE),
    q_25 = quantile(w_vc_log, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(w_vc_log, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(w_vc_log, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(w_vc_log, probs = (1), na.rm = TRUE),
  )
quantiles_workers

quantiles_queens <- group_by(dat, treatment_1, treatment_2) %>% summarise(
    q_0 = quantile(q_vc_log, probs = (0), na.rm = TRUE),
    q_25 = quantile(q_vc_log, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(q_vc_log, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(q_vc_log, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(q_vc_log, probs = (1), na.rm = TRUE),
  )
quantiles_queens

plot(dat$q_vc_log, dat$w_vc_log)
boxplot(dat$q_vc_log ~ dat$treatment)
abline(vdt_log, 0)
boxplot(dat$w_vc_log ~ dat$treatment)
abline(log10(vdt_old), 0)


#virus titer depending on treatment (for each caste)

# virus titer quees
boxplot(dat$q_vc_log ~ dat$treatment)
boxplot(dat$q_vc_log ~ dat$treatment_2*dat$treatment_1)
boxplot(dat$q_vc_log ~ dat$treatment_1*dat$treatment_2)

mod <- lm(q_vc_log ~ treatment_1*treatment_2, data = dat)
summary(mod)
Anova(mod)

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey") #or TukeyHSD(aov(mod))
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

#check model assumptions
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
compareqqnorm(mod)
par(mfrow=c(1,1))
#normality
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> assumption of normally distributed residuals is ok
#heteroskedastity 
bptest(mod) #non significant --> ok
#homogenity of variance 
leveneTest(mod) #non-significant = ok 




# virus titer workers 
boxplot(dat$w_vc_log ~ dat$treatment_1*dat$treatment_2)
boxplot(dat$w_vc_log ~ dat$treatment_2*dat$treatment_1)

library(ARTool) #The Aligned Rank Transform for nonparametric factorial analyses using only ANOVA procedures
dat_noNA <- subset(dat, dat$w_vc_log != "NA")

m = art(w_vc_log ~ treatment_1 * treatment_2, data=dat_noNA) # uses a linear mixed model
anova(m)

#Pairwise comparison within one treatment
library(emmeans)
emmeans(artlm(m, "treatment_1"), pairwise ~ treatment_1)
emmeans(artlm(m, "treatment_2"), pairwise ~ treatment_2)
library(phia)
testInteractions(artlm(m, "treatment_1:treatment_2"), pairwise=c("treatment_1","treatment_2"), adjustment="holm") # i do not fully get how to interprete this


# An good option for cross-factor pairwise comparisons is to use either nonparametric Mann-Whitney U tests or Wilcoxon signed-rank tests on the original data.
# Use the former when subjects were in only one of the conditions being compared (i.e., between-subjects), and use the latter when subjects were in both of the conditions being compared (i.e., within-subjects).
# For example, if the interaction between X1 and X2 is statistically significant from the art and anova calls above, and if X1 now, for simplicity, has levels a and b, and X2 has levels d and e, then we can do:

# cross-factor pairwise comparisons using Mann-Whitney U tests. Note each comparison assumes subjects were in only one of the compared conditions.
uncorrected_pairwise_p_values <- c(
  cc_vs_cv = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "control" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  cc_vs_lc = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "low" & dat$treatment_2 == "control",]$w_vc_log)$p.value,
  cc_vs_lv = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "low" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  cc_vs_hc = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "control",]$w_vc_log)$p.value,
  cc_vs_hv = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  cv_vs_lc = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "virus",]$w_vc_log,   dat[dat$treatment_1 == "low" & dat$treatment_2 == "control",]$w_vc_log)$p.value,
  cv_vs_lv = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "virus",]$w_vc_log,   dat[dat$treatment_1 == "low" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  cv_vs_hc = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "virus",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "control",]$w_vc_log)$p.value,
  cv_vs_hv = wilcox.test(dat[dat$treatment_1 == "control" & dat$treatment_2 == "virus",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  lc_vs_lv = wilcox.test(dat[dat$treatment_1 == "low" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "low" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  lc_vs_hc = wilcox.test(dat[dat$treatment_1 == "low" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "control",]$w_vc_log)$p.value,
  lc_vs_hv = wilcox.test(dat[dat$treatment_1 == "low" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  lv_vs_hc = wilcox.test(dat[dat$treatment_1 == "low" & dat$treatment_2 == "virus",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "control",]$w_vc_log)$p.value,
  lv_vs_hv = wilcox.test(dat[dat$treatment_1 == "low" & dat$treatment_2 == "virus",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "virus",]$w_vc_log)$p.value,
  hc_vs_hv = wilcox.test(dat[dat$treatment_1 == "high" & dat$treatment_2 == "control",]$w_vc_log,   dat[dat$treatment_1 == "high" & dat$treatment_2 == "virus",]$w_vc_log)$p.value
  )
uncorrected_pairwise_p_values
# correct for multiple comparisons using Holm's sequential Bonferroni procedure (Holm 1979)
p.adjust(uncorrected_pairwise_p_values, method="holm")



#### citations ####
#citing R Studio
#RStudio.Version()
#citation()
#R.Version()
#citation("survival") # contains kalan meier plot function
#citation("arm") #contains the sim function for bayesian inference
#citation("dunn.test") #contains the dunn test
#citation("lme4") #lmer
