#####################################################################
### iUSSI analyses ##################################################
#####################################################################

library(dplyr)
library (multcomp)
library(lme4)
library(nlme)
library(plyr)
library(lmerTest)
library(chron)
library(arm)
library(car)
library(survival) # contains kalan meier plot function
library(nortest) #lilliefors test for normality
library(reshape2) #useful for data transformation
library(plyr) #to easily rename column names etc.
library(ggplot2) #not used in the end (I think)
library(MASS)
library(multcomp)
library(tidyverse)
if(!require(psych)){install.packages("psych")}
if(!require(nlme)){install.packages("nlme")}
if(!require(car)){install.packages("car")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(all.effect)){install.packages("all.effect")}
if(!require(emmeans)){install.packages("emmeans")}
library(emmeans)
library(sjPlot)
library(sjmisc)


install.packages("lme4")


setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/iUSSI")
dat <- read.table("ANV_Neonic_in_workers_and_queens.txt", header = TRUE)
head(dat)
dat$treatment <- as.factor(dat$treatment)
dat$TREATMENT_LARS <- as.factor(dat$TREATMENT_LARS)

# relevel treatment 1 to control low high 
levels(dat$treatment_1)
dat$treatment_1 <- factor(dat$treatment_1, levels(dat$treatment_1)[c(1,3,2)])


#all data used -> virus treatments stay included
#only samples that were tested and from the colonies that survived until the end of the experiment
dat <- subset(dat, final == "y")
table(dat$treatment_1)
dat <- subset(dat, sample != "blank")  #exclude blanks and drop the sample level so it does not appear anymore
dat$sample <- droplevels(dat)$sample
table(dat$sample)



#dry weight makes most sense! 
boxplot(dat$conc_T_dry ~ dat$sample*dat$treatment_1)
boxplot(dat$conc_C_dry ~ dat$sample*dat$treatment_1)

workers <- subset(dat, sample == "worker")
queens <- subset(dat, sample == "queen")


#test with linear mixed effect models 

attach(dat)

#thiamethoxam  - First model with no interaction... looks fine and qq of residuals is ok. 
model <- lmer(conc_T_dry ~ sample + treatment_1 + (1|id), data = dat)
summary(model)
anova(model)

plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance
marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD

#Clothianidin
model <- lmer(conc_C_dry ~ sample + treatment_1 + (1|id), data = dat)
summary(model)
anova(model)

plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance

marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="tukey")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD
detach(dat)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Ratio from thiamethoxam to clothianidin ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Because of the many zeros it makes no sense to calculate the variable relation... But we can just calculate the means for the treatments and calculate the relation afterwards
# control is ignored

#create subsets
WC <- subset(workers, treatment_1 == "control")
QC <- subset(queens, treatment_1 == "control")
WL <- subset(workers, treatment_1 == "low")
QL <- subset(queens, treatment_1 == "low")
WH <- subset(workers, treatment_1 == "high")
QH <- subset(queens, treatment_1 == "high")
WHL<- subset(workers, treatment_1 == "high" | treatment_1 == "low")
QHL<- subset(queens, treatment_1 == "high" | treatment_1 == "low")


# low t/c workers and queens
TCRL_workers <- mean(WL$conc_T_dry)/mean(WL$conc_C_dry)
TCRL_workers
TCRL_queens  <- mean(QL$conc_T_dry)/mean(QL$conc_C_dry)
TCRL_queens

# high t/c workers and queens
TCRH_workers <- mean(WH$conc_T_dry)/mean(WH$conc_C_dry)
TCRH_workers
TCRH_queens  <- mean(QH$conc_T_dry)/mean(QH$conc_C_dry)
TCRH_queens

# low and high workers and queens
TCRLH_workers <- mean(WHL$conc_T_dry)/mean(WHL$conc_C_dry)
TCRLH_workers
TCRLH_queens <- mean(QHL$conc_T_dry)/mean(QHL$conc_C_dry)
TCRLH_queens 

### Ratio only for high treatment (no zeros) so we can calculate it for all samples and do stats for it ###

WH$conc_C_dry
QH$conc_C_dry

QH$conc_C_dry[QH$conc_C_dry == 0] <- NA
WH$ratio <- WH$conc_T_dry/WH$conc_C_dry
QH$ratio <- QH$conc_T_dry/QH$conc_C_dry
boxplot()
boxplot(QH$ratio)

#neu für den ganzen Datensatz und dann alle Zeors mit NA replacen für statistik comparison with glmer
dat$c <- dat$conc_C_dry
dat$c[dat$c == 0] <- NA
dat$c
dat$t <- dat$conc_T_dry
dat$t[dat$t == 0] <- NA
dat$t
dat$ratio_TC <- dat$t/dat$c
dat$ratio_CT <- dat$c/dat$t #the nicer one

boxplot(dat$ratio_CT ~ dat$sample)


#for the model take only the complete cases
dat_model <- subset(dat, ratio_CT != "NA")
head(dat_model)
levels(dat_model$id)
dat_model$id <- droplevels(dat_model)$id
levels(dat_model$id)

dat_model_reduced <- subset(dat_model, id == "L19C" | id == "L20V" |id == "H1V" |id == "H2C" |id == "H3V" |id == "H5V" | id == "H7C" |id == "H9V" | id == "H10C" |id == "H13V" |id == "H14C" |id == "H15V" | id == "H16C" | id == "H20C")
boxplot(dat_model_reduced$ratio_CT ~ dat_model_reduced$sample)
model <- glmer(ratio_CT ~ sample + 1|id,  data = dat_model_reduced, family = "Gamma") #ev we should also include the treatment
summary(model)
anova(model)
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance

#model super complicated and with only a subset of values and so... 
#simple comparison of two means t test if normally distributed
shapiro.test(dat$ratio_CT)
#signifikant und offensichtlich nicht normal verteilt --> wilcoxtest
wilcox.test(dat$ratio_CT ~ dat$sample)

mean(dat$ratio_CT[dat$sample == "queen"], na.rm=TRUE)
sd(dat$ratio_CT[dat$sample == "queen"], na.rm=TRUE)

mean(dat$ratio_CT[dat$sample == "worker"], na.rm=TRUE)
sd(dat$ratio_CT[dat$sample == "worker"], na.rm=TRUE)


#### ratio Graph for Journalclub talk ####
#Thiamethoxam
par(mar=c(5,5,5,1))
boxplot(dat$ratio_CT ~ dat$sample,
        main = "clothianidin-thiamethoxam ratio", 
        ylab = "clothianidin [ng/g] / thiamethoxam [ng/g]", 
        xlab = "caste", 
        col = c("white","white"), 
        cex.lab = 1.5, cex.main = 3, cex.axis = 1.5, lwd = 1
)

text(c(1,2), 8, labels = c("a", "b"),  cex = 1.5, font = 2)


#idea: Queens have relation to thiamethoxam concentration lower clothianidin conc. than workers due to better detoxabilities



### drop the controls, as they all were zero anyways despite one... to make the graphs look better

dat_reduced <- subset(dat, treatment_1 != "control")
table(dat_reduced$treatment_1)
dat_reduced$treatment_1 <- droplevels(dat_reduced)$treatment_1
table(dat_reduced$treatment_1)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### new statistics with interaction effects #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

attach(dat_reduced)

#thiamethoxam
#lme
baseline <- lme(conc_T_dry ~ 1, random = ~1|id, data = dat_reduced, method = "ML")
sample_M <- update(baseline, .~. + sample)
treatment_M <- update(sample_M, .~. + treatment_1)
interaction_M <- update(treatment_M, .~. + sample*treatment_1 )
anova(baseline, sample_M, treatment_M, interaction_M)

model <- interaction_M
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance

marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD

#lmer
baseline <- lmer(conc_T_dry ~ 1|id, data = dat_reduced)
sample_M <- update(baseline, .~. + sample)
treatment_M <- update(sample_M, .~. + treatment_1)
interaction_M <- update(treatment_M, .~. + sample*treatment_1 )
anova(baseline, sample_M, treatment_M, interaction_M)
model <- interaction_M
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance
marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD


#glmer for thiamethoxam
baseline <- glmer(conc_T_dry ~ 1|id, data = dat_reduced, family = Gamma(link=log))  #gamma funktioniert nicht?
sample_M <- update(baseline, .~. + sample)
treatment_M <- update(sample_M, .~. + treatment_1)
interaction_M <- update(treatment_M, .~. + sample*treatment_1 )
anova(baseline, sample_M, treatment_M, interaction_M)
model <- interaction_M
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance
marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD


# glmer for clothianidin
dat_reduced$conc_C_dry_plus1 <- dat_reduced$conc_C_dry + 1
dat_reduced$conc_C_dry_plus1
dat_reduced$conc_C_dry_log <- log10(dat_reduced$conc_C_dry_plus1)
dat_reduced$conc_C_dry_log
baseline <- lmer(conc_C_dry_log ~ 1|id, data = dat_reduced)  
sample_M <- update(baseline, .~. + sample)
treatment_M <- update(sample_M, .~. + treatment_1)
interaction_M <- update(treatment_M, .~. + sample*treatment_1 )
anova(baseline, sample_M, treatment_M, interaction_M)
model <- interaction_M
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance
marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD

baseline <- glmer(conc_C_dry_plus1 ~ 1|id, data = dat_reduced, family = Gamma(link=log))  #gamma funktioniert nicht?
sample_M <- update(baseline, .~. + sample)
treatment_M <- update(sample_M, .~. + treatment_1)
interaction_M <- update(treatment_M, .~. + sample*treatment_1 )
anova(baseline, sample_M, treatment_M, interaction_M)
model <- interaction_M
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance
marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Final analyzes on the poster ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# glmer would be ok as we could use a gamma distribution with the log link for a nice model fit and good residual distributions 
# however the same can be achieved by using the simpler lmer models with normaly distributed residuals after log-transformation (with plus 1 so that log10 von 1 = 0)


# lmer for Thiamethoxam

dat_reduced$conc_T_dry_plus1 <- dat_reduced$conc_T_dry + 1
dat_reduced$conc_T_dry_plus1
dat_reduced$conc_T_dry_log <- log10(dat_reduced$conc_T_dry_plus1)
dat_reduced$conc_T_dry_log

baseline <- lmer(conc_T_dry_log ~ 1|id, data = dat_reduced)  
sample_M <- update(baseline, .~. + sample)
treatment_M <- update(sample_M, .~. + treatment_1)
interaction_M <- update(treatment_M, .~. + sample*treatment_1 )
anova(baseline, sample_M, treatment_M, interaction_M)

model <- interaction_M
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance

marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD

# lmer for Clothianidin

dat_reduced$conc_C_dry_plus1 <- dat_reduced$conc_C_dry + 1
dat_reduced$conc_C_dry_plus1
dat_reduced$conc_C_dry_log <- log10(dat_reduced$conc_C_dry_plus1)
dat_reduced$conc_C_dry_log

baseline <- lmer(conc_C_dry_log ~ 1|id, data = dat_reduced)  
sample_M <- update(baseline, .~. + sample)
treatment_M <- update(sample_M, .~. + treatment_1)
interaction_M <- update(treatment_M, .~. + sample*treatment_1 )
anova(baseline, sample_M, treatment_M, interaction_M)

model <- interaction_M
plot(model)
scatter.smooth(fitted(model),resid(model)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(model), main="normal QQ-plot, residuals") 
qqline(resid(model))  # qq of residuals
scatter.smooth(fitted(model), sqrt(abs(resid(model))))  # homogeneity of variance

marginal <- lsmeans(model,  ~ sample + treatment_1, adjust="bonferroni")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="bonferroni")
CLD

detach(dat_reduced)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
####  graphs for poster ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# drop the controls, as they all were zero anyways despite one... to make the graphs look better

dat_reduced <- subset(dat, treatment_1 != "control")
table(dat_reduced$treatment_1)
dat_reduced$treatment_1 <- droplevels(dat_reduced)$treatment_1
table(dat_reduced$treatment_1)

#### Graph for poster ####
par(mfrow=c(1,1))
par(mar=c(5,5,5,1))

#Thiamethoxam
boxplot(dat_reduced$conc_T_dry ~ dat_reduced$sample + dat_reduced$treatment_1, 
        main = "thiamethoxam", 
        ylab = "ng thiamethoxam / g dryweight", 
        xlab = "neonicotinoid treatment", 
        col = c("white","grey"), 
        xaxt="n", 
        at = c(0.7, 1.7, 3.3, 4.3), 
        cex.lab = 1.5, cex.main = 3, cex.axis = 1.5, lwd = 1
)

x <- c("low","high") 
axis(1, at = c(1.2,3.8), labels=x, cex.axis = 1.5)

legend(0.1, 2.5, legend = c("queens", "workers"),  c("white","grey"), bty = "n", cex=1.5)
text(c(0.6, 1.6, 3.2, 4.2), 2, labels = c("a", "a", "a", "b"),  cex = 1.5, font = 2)

#Clothianidin
boxplot(dat_reduced$conc_C_dry ~ dat_reduced$sample + dat_reduced$treatment_1, 
        main = "clothianidin", 
        ylab = "ng clothianidin / g dryweight", 
        xlab = "neonicotinoid treatment", 
        col = c("white","grey"), 
        xaxt="n", 
        at = c(0.7, 1.7, 3.3, 4.3),
        cex.lab = 1.5, cex.main = 3, cex.axis = 1.5, lwd = 1
)

x <- c("low","high") 
axis(1, at = c(1.2,3.8), labels=x, cex.axis = 1.5)
legend(0.1, 3.6, legend = c("queens", "workers"),  c("white","grey"), bty = "n", cex=1.5)
text(c(0.6, 1.6, 3.2, 4.2), 2.9, labels = c("a", "b", "ab", "c"),  cex = 1.5, font = 2)







