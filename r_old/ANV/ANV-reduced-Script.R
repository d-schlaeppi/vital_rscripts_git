### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ANV - Ants, Neonics and Viruses - Interaction of stressors  ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### load libraries ####
suppressPackageStartupMessages(library(survival))# contains kalan meier plot function
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(chron))
suppressPackageStartupMessages(library(car))
#library(car) #contains levene's test
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
suppressPackageStartupMessages(library(ARTool)) #The Aligned Rank Transform for nonparametric factorial analyses using only ANOVA procedures
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(effsize)


#### prerequisites ####
# setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV") #office computer
setwd("/Users/gismo/Desktop/R/ANV") # homeoffice mac
dat <- read.table("ANV_all_Data.txt", header = TRUE)
mov <- read.table("speed.txt", header = TRUE)



# preparation of Dataframe --> dat 
head(dat)
str(dat)
dat_r <- subset(dat, treatment_1 != "low") #remove low because it does not really add valuable information to the study
dat_r$treatment_1 <- droplevels(dat_r)$treatment_1
dat_r$treatment <- droplevels(dat_r)$treatment
dat_r$identity <- droplevels(dat_r)$identity
dat_r$treatment <- as.factor(dat_r$treatment)
dat_r$treatment <- revalue(dat_r$treatment, c("1"="Control", "2"="Virus","5" = "Neonic","6" = "Interaction"))
dat_r$tr <- dat_r$treatment
dat_r$tr <- relevel(dat_r$tr, "Neonic") 
dat_r$tr <- relevel(dat_r$tr, "Control")
table(dat_r$tr)

#dat_r$treatment <- revalue(dat_r$treatment, c("control-control"="Control", "control-virus"="Virus","high-control" = "Neonicotinoid","high-virus" = "Mixed" ))
dat_r$treatment_1 <- revalue(dat_r$treatment_1, c("control" = "Control", "high" = "Neonicotinoid"))
dat_r$treatment_2 <- revalue(dat_r$treatment_2, c("control" = "Control", "virus" = "Virus"))

table(dat_r$treatment)
table(dat_r$treatment_2)
table(dat_r$treatment_1)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Cumulative queen survival in percent ####
#mean survival time accross the treatments
# drawing nice survival curves: https://rpkgs.datanovia.com/survminer/index.html


tapply(dat_r$day_of_death[dat_r$y2_survival_status==1], dat_r$treatment[dat_r$y2_survival_status==1], mean)
plot(dat_r$day_of_death ~ dat_r$treatment)
plot(survfit(Surv(dat_r$day_of_death, dat_r$y2_survival_status)~1)) #plots first basic plot with 95 confidence intervalls for all the data
fit=survfit(Surv(dat_r$day_of_death, dat_r$y2_survival_status)~dat_r$treatment) 
plot(fit)
summary(fit)
ggsurvplot(fit, data = dat_r, pval = TRUE,
           xlab = "Time in days",
           ylim = c(0.5,1),
           legend.labs = c("Control", "Mixed", "Neonic", "Virus"),
           ggtheme = theme_bw())

#calculate potential statistical difference between the groups
mod <- survdiff(Surv(day_of_death, y2_survival_status) ~ treatment, data=dat_r)
mod # no statistical difference between the four groups. 

require(dplyr)
group_by(dat_r, treatment_1, treatment_2) %>%
  summarise(
    count = n(),
    count_survivers = sum(y2_survival_status == 0),
    count_death = sum(y2_survival_status == 1)
  )


#### Number of Queens producing adults - Is there a sig. difference in the number of queens producing adults between the treatment? ####

dat_r$adult_production <- as.factor(dat_r$adult_production)
plot(dat_r$adult_production ~ dat_r$treatment)
chisq.test(dat_r$adult_production, dat_r$treatment_1)

#generalised logistic mixed model!!!!!!!!!!! (ev ordinal family, Bernulli ja / nein 1/0
mod <- glm(dat_r$adult_production ~ dat_r$treatment_1*dat_r$treatment_2, family = binomial)
summary(mod)
Anova(mod)

# test model assumtions
par(mfrow=c(2,2))
leveneTest(residuals(mod) ~ dat_r$treatment_1*dat_r$treatment_2) #non-significant = ok #homogenity of variance (plot(mod,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #does not need to be fulfilled in a binomial glm
bptest(mod) #non significant --> assumtions of heteroskedasticity not violated ok 






#### Bodymass analyses ####

### workers ###
boxplot(dat_r$weight_adults ~ dat_r$treatment)
#stat summary:
require(dplyr)

group_by(dat_r, tr) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(weight_adults)), 
    mean = mean(weight_adults, na.rm = TRUE),
    sd = sd(weight_adults, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$weight_adults)$p,
    levene = leveneTest(dat_r$weight_adults, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(weight_adults, probs = (0), na.rm = TRUE),
    q_25 = quantile(weight_adults, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(weight_adults, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(weight_adults, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(weight_adults, probs = (1), na.rm = TRUE),
  )


#model
mod <- lm(weight_adults ~ treatment_1*treatment_2, data = dat_r) # is equal to the interaction model. 
summary(mod)
Anova(mod)
baseline <- lm(weight_adults ~ 1, data = dat_r)
neonic_M <- update(baseline, .~. + treatment_1)
virus_M <- update(neonic_M, .~. + treatment_2)
interaction_M <- update(virus_M, .~. + treatment_1:treatment_2)
anova(baseline, neonic_M, virus_M, interaction_M)
Anova(mod)

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1*treatment_2, data = dat_r, adjust = "tukey")
TukeyHSD(aov(mod))
marginal = lsmeans(mod, ~ treatment_1*treatment_2, data = dat_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# test model assumtions
compareqqnorm(mod)
par(mfrow=c(2,2))
leveneTest(mod) #non-significant = ok #homogenity of variance (plot(mod,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #non significant --> assumption of normally distributed residuals is ok ##plot(mod, 2)#normality
bptest(mod) #non significant --> assumtions of heteroskedasticity not violated ok 
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

#Effect size
control_neonic <- subset(dat_r, dat_r$treatment == "Control" | dat_r$treatment == "Neonic")
control_virus <- subset(dat_r, dat_r$treatment == "Control" | dat_r$treatment == "Virus")
control_combined <- subset(dat_r, dat_r$treatment == "Control" | dat_r$treatment == "Interaction")
neonic_virus <- subset(dat_r, dat_r$treatment == "Neonic" | dat_r$treatment == "Virus")
neonic_combined <- subset(dat_r, dat_r$treatment == "Neonic" | dat_r$treatment == "Interaction")
virus_combined <- subset(dat_r, dat_r$treatment == "Virus" | dat_r$treatment == "Interaction")

control_neonic$treatment <- droplevels(control_neonic)$treatment
control_virus$treatment <- droplevels(control_virus)$treatment
control_combined$treatment <- droplevels(control_combined)$treatment
neonic_virus$treatment <- droplevels(neonic_virus)$treatment
neonic_combined$treatment <- droplevels(neonic_combined)$treatment
virus_combined$treatment <- droplevels(virus_combined)$treatment

cohen.d(control_neonic$weight_adults, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$weight_adults, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$weight_adults, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$weight_adults, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$weight_adults, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$weight_adults, virus_combined$treatment, na.rm = TRUE)



### queens ### 
boxplot(dat_r$weight_queens ~ dat_r$treatment)
#stat summary:
require(dplyr)
group_by(dat_r, tr) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(weight_queens)), 
    mean = mean(weight_queens, na.rm = TRUE),
    sign = "±",
    sd = sd(weight_queens, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$weight_queens)$p,
    levene = leveneTest(dat_r$weight_queens, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(weight_queens, probs = (0), na.rm = TRUE),
    q_25 = quantile(weight_queens, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(weight_queens, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(weight_queens, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(weight_queens, probs = (1), na.rm = TRUE),
  ) %>% as.data.frame()

#model
mod <- lm(weight_queens ~ treatment_1*treatment_2, data = dat_r)
summary(mod)
Anova(mod)
baseline <- lm(weight_queens ~ 1, data = dat_r)
neonic_M <- update(baseline, .~. + treatment_1)
virus_M <- update(neonic_M, .~. + treatment_2)
interaction_M <- update(virus_M, .~. + treatment_1:treatment_2)
anova(baseline, neonic_M, virus_M, interaction_M)

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1*treatment_2, data = dat_r, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1*treatment_2, data = dat_r)
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




#### Colony size after year 1 ####
#Y1 only for treatment 1 because 2 not yet established --> y2 for all 4 treatments
# no need for random factor beceause there are no repeated measures and conditions the same for all apart from the treatment
# In the mansucript we will only adutls and brood (combined eggs larva & pupae )

#overview
dat_r$y1_brood <- dat_r$y1_eggs + dat_r$y1_larva + dat_r$y1_pupae

boxplot(dat_r$y1_adults ~ dat_r$treatment_1)
boxplot(dat_r$y1_brood ~ dat_r$treatment_1)
boxplot(dat_r$y1_eggs ~ dat_r$treatment_1)
boxplot(dat_r$y1_larva ~ dat_r$treatment_1)
boxplot(dat_r$y1_pupae ~ dat_r$treatment_1)

### adults ### 
# summary
require(dplyr)
group_by(dat_r, treatment_1) %>%
  summarise(
    count = n(),
    surviving = sum(y1_survival_status == 0),
    dead = sum(y1_survival_status == 1),
    mean = mean(y1_adults, na.rm = TRUE),
    sd = sd(y1_adults, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$y1_adults)$p,
    levene = leveneTest(dat_r$y1_adults, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(y1_adults, probs = (0), na.rm = TRUE),
    q_25 = quantile(y1_adults, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(y1_adults, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(y1_adults, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(y1_adults, probs = (1), na.rm = TRUE),
  ) %>% as.data.frame()

#data not normally distributed resp modelassumptions are violated with linear models --> non parametric test
kruskal.test(dat$y1_adults ~ dat$treatment_1, na.action = na.omit)

### brood ###
require(dplyr)
group_by(dat_r, treatment_1) %>%
  summarise(
    count = n(),
    surviving = sum(y1_survival_status == 0),
    dead = sum(y1_survival_status == 1),
    mean = mean(y1_brood, na.rm = TRUE),
    sd = sd(y1_brood, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$y1_brood)$p,
    levene = leveneTest(dat_r$y1_brood, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(y1_brood, probs = (0), na.rm = TRUE),
    q_25 = quantile(y1_brood, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(y1_brood, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(y1_brood, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(y1_brood, probs = (1), na.rm = TRUE),
  ) %>% as.data.frame()

kruskal.test(dat_r$y1_brood ~ dat_r$treatment_1, na.action = na.omit)

#### Colony after year II - adults and brood ####
#overview
dat_r$y2_brood <- dat_r$y2_eggs + dat_r$y2_larva + dat_r$y2_pupae

boxplot(dat_r$y2_adults ~ dat_r$treatment)
boxplot(dat_r$y2_brood ~ dat_r$treatment)


### adults y2 ### 
# summary
require(dplyr)
group_by(dat_r, tr) %>%
  summarise(
    count = n(),
    surviving = sum(y2_survival_status == 0),
    dead = sum(y2_survival_status == 1),
    mean____ = mean(y2_adults, na.rm = TRUE),
    sd = sd(y2_adults, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$y2_adults)$p,
    levene = leveneTest(dat_r$y2_adults, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(y2_adults, probs = (0), na.rm = TRUE),
    q_25 = quantile(y2_adults, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(y2_adults, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(y2_adults, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(y2_adults, probs = (1), na.rm = TRUE),
  ) %>% as.data.frame()


#model
boxplot(dat_r$y2_adults ~ dat_r$treatment)
mod <- lm(y2_adults ~ treatment_1*treatment_2, data = dat_r)
summary(mod)
baseline <- lm(y2_adults ~ 1, data = dat_r)
neonic_M <- update(baseline, .~. + treatment_1)
virus_M <- update(neonic_M, .~. + treatment_2)
interaction_M <- update(virus_M, .~. + treatment_1:treatment_2)
anova(baseline, neonic_M, virus_M, interaction_M)
Anova(mod)

coefficients(mod)

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
#homogenity of variance # plot(mod, 1)
leveneTest(mod) #non-significant = ok 
#normality  # plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> ok
#heteroskedastity 
bptest(mod) #non significant --> ok


#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey") #or TukeyHSD(aov(mod))
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

non_NA <- subset(dat_r, dat_r$y2_adults != "NA")
non_NA$y2_adults

#Pairwise Cohens D
cohen.d(control_neonic$y2_adults, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$y2_adults, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$y2_adults, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$y2_adults, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$y2_adults, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$y2_adults, virus_combined$treatment, na.rm = TRUE)






### Ideas and stuff to plot interactions 
#Model.1 <- lm(y2_adults~treatment_1+treatment_2, dat_r)
#Model.2 <- lm(y2_adults~treatment_1*treatment_2, dat_r)
#library(stargazer)
#stargazer(Model.1, Model.2,type="text", 
#          column.labels = c("Main Effects", "Interaction"), 
#          intercept.bottom = FALSE, 
#          single.row=FALSE,     
#          notes.append = FALSE, 
#          header=FALSE) 
#FacetPlot1 = ggplot(dat_r, aes(x=treatment_1, y=y2_adults, fill = treatment_2)) + geom_boxplot() + facet_grid(~treatment_2)
#FacetPlot1
#par(mfrow=c(1,1))
#interaction.plot(x.factor = non_NA$treatment_1, trace.factor = non_NA$treatment_2, response = non_NA$y2_adults, 
#                 xlab = "Neonicotinoid Treatment", ylab = "Predicted mean number of workers", 
#                 trace.label = "Virus Treatment", main = "Interaction Plot - Worker number (weeek 64)")
#plot_model(mod, type = "int", 
#           title = "Predicted values of the number of workers (week 64)", 
#           axis.title = c("Thiamethoxam treatment", "Number of workers"), 
#           legend.title = "Virus treatment")



### brood ###
require(dplyr)
group_by(dat_r, tr) %>%
  summarise(
    count = n(),
    surviving = sum(y2_survival_status == 0),
    dead = sum(y2_survival_status == 1),
    mean = mean(y2_brood, na.rm = TRUE),
    sd = sd(y2_brood, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$y2_brood)$p,
    levene = leveneTest(dat_r$y2_brood, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(y2_brood, probs = (0), na.rm = TRUE),
    q_25 = quantile(y2_brood, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(y2_brood, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(y2_brood, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(y2_brood, probs = (1), na.rm = TRUE),
  ) %>% as.data.frame()
boxplot(dat_r$y2_brood ~ dat_r$treatment)

dat_r$y2_brood

mod1 <- lm(dat_r$y2_brood ~ treatment_1*treatment_2, data = dat_r)
summary(mod1)
Anova(mod1)
mod1$coefficients
#homogenity of variance # plot(mod, 1)
leveneTest(mod1) #significant = bad 
#normality  # plot(mod, 2)
aov_residuals <- residuals(object = mod1)
shapiro.test(x = aov_residuals) # non significant = ok
#heteroskedastity 
bptest(mod1) #significant --> bad
# the use of linear models seems not appropriate as homogeneity of variances are violated and/or 
# if the date are logtransformed the normality of the residuals can not be assumed

#ART MODEL
dat_noNA <- subset(dat_r, dat_r$y2_brood != "NA") #model can not process NA's!
m = art(y2_brood ~ treatment_1 * treatment_2, data=dat_noNA) # uses a linear mixed model
anova(m)
names(m)





#pairwise comparisons: 
# An good option for cross-factor pairwise comparisons is to use either nonparametric Mann-Whitney U tests or Wilcoxon signed-rank tests on the original data.
# Use the former when subjects were in only one of the conditions being compared (i.e., between-subjects), and use the latter when subjects were in both of the conditions being compared (i.e., within-subjects).
# For example, if the interaction between X1 and X2 is statistically significant from the art and anova calls above, and if X1 now, for simplicity, has levels a and b, and X2 has levels d and e, then we can do:

# cross-factor pairwise comparisons using Mann-Whitney U tests. Note each comparison assumes subjects were in only one of the compared conditions.

uncorrected_pairwise_p_values <- c(
  c_v = wilcox.test(dat_r[dat_r$treatment == "Control",]$y2_brood,   dat_r[dat_r$treatment == "Virus",]$y2_brood)$p.value, 
  c_n = wilcox.test(dat_r[dat_r$treatment == "Control",]$y2_brood,   dat_r[dat_r$treatment == "Neonic",]$y2_brood)$p.value,
  c_m = wilcox.test(dat_r[dat_r$treatment == "Control",]$y2_brood,   dat_r[dat_r$treatment == "Interaction",]$y2_brood)$p.value,
  v_n = wilcox.test(dat_r[dat_r$treatment == "Virus",]$y2_brood,     dat_r[dat_r$treatment == "Neonic",]$y2_brood)$p.value,
  v_m = wilcox.test(dat_r[dat_r$treatment == "Virus",]$y2_brood,     dat_r[dat_r$treatment == "Interaction",]$y2_brood)$p.value,
  n_m = wilcox.test(dat_r[dat_r$treatment == "Neonic",]$y2_brood,   dat_r[dat_r$treatment == "Interaction",]$y2_brood)$p.value
)

uncorrected_pairwise_p_values

# correct for multiple comparisons using Holm's sequential Bonferroni procedure (method="holm") or the bonferroni correction (Holm 1979)
p.adjust(uncorrected_pairwise_p_values, method="bonferroni")


#Pairwise Cohens D
#Effect size
control_neonic <- subset(dat_r, dat_r$treatment == "Control" | dat_r$treatment == "Neonic")
control_virus <- subset(dat_r, dat_r$treatment == "Control" | dat_r$treatment == "Virus")
control_combined <- subset(dat_r, dat_r$treatment == "Control" | dat_r$treatment == "Interaction")
neonic_virus <- subset(dat_r, dat_r$treatment == "Neonic" | dat_r$treatment == "Virus")
neonic_combined <- subset(dat_r, dat_r$treatment == "Neonic" | dat_r$treatment == "Interaction")
virus_combined <- subset(dat_r, dat_r$treatment == "Virus" | dat_r$treatment == "Interaction")
control_neonic$treatment <- droplevels(control_neonic)$treatment
control_virus$treatment <- droplevels(control_virus)$treatment
control_combined$treatment <- droplevels(control_combined)$treatment
neonic_virus$treatment <- droplevels(neonic_virus)$treatment
neonic_combined$treatment <- droplevels(neonic_combined)$treatment
virus_combined$treatment <- droplevels(virus_combined)$treatment

cohen.d(control_neonic$y2_brood, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$y2_brood, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$y2_brood, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$y2_brood, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$y2_brood, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$y2_brood, virus_combined$treatment, na.rm = TRUE)









#### create dataframe to compare the castes using models ####
dat_r$tw <- dat_r$tw_thiamethoxam_workers_dryconcentration
dat_r$cw <- dat_r$cw_clothianidin_workers_dryconcentration
dat_r$tq <- dat_r$tq_thiamethoxam_queens_dryconcentration
dat_r$cq <- dat_r$cq_clothianidin_queens_dryconcentration

#create virus titer per mg bodyweight to compare them between castes and species. 
dat_r$log10_viral_copies_per_mg_workers <- log10(dat_r$virus_copies_per_worker/(dat_r$weight_adults/20))
dat_r$log10_viral_copies_per_mg_queens <- log10(dat_r$virus_copies_per_queen/dat_r$weight_queens)
dat_r$infection_status_queens <- revalue(dat_r$pos_neg_queens, c("neg" = 0, "pos" = 1))
dat_r$infection_status_workers <- revalue(dat_r$pos_neg_workers, c("neg" = 0, "pos" = 1))

#virus titer dead individuals of colony H9V
copies_per_individual <- 1.34E+07
copies_per_mg <- copies_per_individual/(mean(dat_r$weight_adults, na.rm = TRUE)/20)
1.4212794*10^7

plot(as.factor(dat_r$infection_status_queens)~dat_r$treatment)
plot(as.factor(dat_r$infection_status_workers)~dat_r$treatment)

data <- data.frame(sample =c(1:40, 1:40), 
                         identity = rep(dat_r$identity, times = 2), 
                         treatment = rep(dat_r$treatment, times = 2), 
                         treatment_1 = rep(dat_r$treatment_1, times = 2), 
                         treatment_2 = rep(dat_r$treatment_2, times = 2), 
                         caste = c(rep("w", times = length(dat_r$identity)), rep("q", times = length(dat_r$identity))),
                         thiamethoxam = c(dat_r$tw, dat_r$tq), 
                         clothianidin = c(dat_r$cw, dat_r$cq), 
                         virus_titre = c(dat_r$log10_viral_copies_per_mg_workers, dat_r$log10_viral_copies_per_mg_queens)
)
head(data)



#### virus titers ####

# Virus detection treshold
vdt_workers_uncorrected <- 3.03E+02
vdt_queens_uncorrected <- 1.48E+04
vdt_workers <- log10(vdt_workers_uncorrected/(mean(dat_r$weight_adults, na.rm = TRUE)/20))
vdt_queens <- log10(vdt_queens_uncorrected/mean(dat_r$weight_queens, na.rm = TRUE))

# Virus detection treshold new (derived form the max instead of the average of the negative samples)
vdt_workers_uncorrected_new <- 1.00E+03
vdt_queens_uncorrected_new <- 4.07E+04
vdt_workers_new <- log10(vdt_workers_uncorrected_new/(mean(dat_r$weight_adults, na.rm = TRUE)/20))
vdt_queens_new <- log10(vdt_queens_uncorrected_new/mean(dat_r$weight_queens, na.rm = TRUE))

virus_load_dead_workers <- log10(1.34E+07/(mean(dat_r$weight_adults, na.rm = TRUE)/20))
max(dat_r$log10_viral_copies_per_mg_queens, na.rm = TRUE)
max(dat_r$log10_viral_copies_per_mg_workers, na.rm = TRUE)
max(dat_r$log10_viral_copies_per_mg_workers, na.rm = TRUE)

10^virus_load_dead_workers/10^max(dat_r$log10_viral_copies_per_mg_workers, na.rm = TRUE)


# look at data 
#queens
group_by(dat_r, tr) %>%
  summarise(
    count = n(),
    N_Q = sum(!is.na(log10_viral_copies_per_mg_queens)),
    positive_q = sum(!is.na(log10_viral_copies_per_mg_queens)) - sum(pos_neg_queens=="neg", na.rm = TRUE),
    negative_q = sum(pos_neg_queens=="neg", na.rm = TRUE),
    mean_titer_q = mean(log10_viral_copies_per_mg_queens, na.rm = TRUE),
    sd_titer_q = sd(log10_viral_copies_per_mg_queens, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$log10_viral_copies_per_mg_queens)$p,
    levene = leveneTest(dat_r$log10_viral_copies_per_mg_queens, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(log10_viral_copies_per_mg_queens, probs = (0), na.rm = TRUE),
    q_25 = quantile(log10_viral_copies_per_mg_queens, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(log10_viral_copies_per_mg_queens, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(log10_viral_copies_per_mg_queens, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(log10_viral_copies_per_mg_queens, probs = (1), na.rm = TRUE),
  )

#only positive values for summary stats ?
virus_q_pos <- subset(dat_r, dat_r$pos_neg_queens == "pos")
group_by(virus_q_pos, tr) %>%
  summarise(
    count = n(),
    N_Q = sum(!is.na(log10_viral_copies_per_mg_queens)),
    positive_q = sum(!is.na(log10_viral_copies_per_mg_queens)) - sum(pos_neg_queens=="neg", na.rm = TRUE),
    negative_q = sum(pos_neg_queens=="neg", na.rm = TRUE),
    mean_titer_q = mean(log10_viral_copies_per_mg_queens, na.rm = TRUE),
    sd_titer_q = sd(log10_viral_copies_per_mg_queens, na.rm = TRUE),
    q_0 = quantile(log10_viral_copies_per_mg_queens, probs = (0), na.rm = TRUE),
    q_25 = quantile(log10_viral_copies_per_mg_queens, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(log10_viral_copies_per_mg_queens, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(log10_viral_copies_per_mg_queens, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(log10_viral_copies_per_mg_queens, probs = (1), na.rm = TRUE),
  )





#workers
group_by(dat_r, tr) %>%
  summarise(
    count = n(),
    N_Q = sum(!is.na(log10_viral_copies_per_mg_workers)),
    positive_w = sum(!is.na(log10_viral_copies_per_mg_workers)) - sum(pos_neg_workers=="neg", na.rm = TRUE),
    negative_w = sum(pos_neg_workers=="neg", na.rm = TRUE),
    mean_titer_w = mean(log10_viral_copies_per_mg_workers, na.rm = TRUE),
    sd_titer_w = sd(log10_viral_copies_per_mg_workers, na.rm = TRUE),
    shapiro = shapiro.test(dat_r$log10_viral_copies_per_mg_workers)$p,
    levene = leveneTest(dat_r$log10_viral_copies_per_mg_workers, dat_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(log10_viral_copies_per_mg_workers, probs = (0), na.rm = TRUE),
    q_25 = quantile(log10_viral_copies_per_mg_workers, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(log10_viral_copies_per_mg_workers, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(log10_viral_copies_per_mg_workers, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(log10_viral_copies_per_mg_workers, probs = (1), na.rm = TRUE),
  )


plot(dat_r$log10_viral_copies_per_mg_queens, dat_r$log10_viral_copies_per_mg_workers)

boxplot(dat_r$log10_viral_copies_per_mg_queens ~ dat_r$treatment)
boxplot(dat_r$log10_viral_copies_per_mg_workers ~ dat_r$treatment)
plot(dat_r$pos_neg_queens ~ dat_r$treatment)
plot(dat_r$pos_neg_workers ~ dat_r$treatment)

### positive queens or workers ###
mod <- glm(infection_status_queens ~ treatment_1*treatment_2, data=dat_r, family="binomial")
Anova(mod)
mod <- glm(infection_status_workers ~ treatment_1*treatment_2, data=dat_r, family="binomial")
Anova(mod)

#Analyse castes separately
queens <- subset(data, caste == "q")
workers <- subset(data, caste == "w")


### workers ###
boxplot(log10(workers$virus_titre) ~ workers$treatment)
boxplot(workers$virus_titre ~ workers$treatment)
abline(vdt_workers, 0, lty = 2)
abline(vdt_workers_new, 0, lty = 2)

mod <- lm(log10(virus_titre) ~ treatment_1*treatment_2, data = workers)
baseline <- lm(log10(virus_titre) ~ 1, data = workers)
treatment_1_M <- update(baseline, .~. + treatment_1)
treatment_2_M <- update(treatment_1_M, .~. + treatment_2)
interaction_M <- update(treatment_2_M, .~. + treatment_1:treatment_2 )
anova(baseline, treatment_1_M, treatment_2_M, interaction_M)
Anova(mod)
mod <- interaction_M

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey") #or TukeyHSD(aov(mod))
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

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
#homogenity of variance # plot(mod, 1)
leveneTest(mod) #non-significant = ok 
#normality  # plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> ok
#heteroskedastity 
bptest(mod) #non significant --> ok

#Pairwise Cohens D
#Effect size
levels(workers$treatment)
control_neonic <- subset(workers, workers$treatment == "Control" | workers$treatment == "Neonic")
control_virus <- subset(workers, workers$treatment == "Control" | workers$treatment == "Virus")
control_combined <- subset(workers, workers$treatment == "Control" | workers$treatment == "Interaction")
neonic_virus <- subset(workers, workers$treatment == "Neonic" | workers$treatment == "Virus")
neonic_combined <- subset(workers, workers$treatment == "Neonic" | workers$treatment == "Interaction")
virus_combined <- subset(workers, workers$treatment == "Virus" | workers$treatment == "Interaction")
control_neonic$treatment <- droplevels(control_neonic)$treatment
control_virus$treatment <- droplevels(control_virus)$treatment
control_combined$treatment <- droplevels(control_combined)$treatment
neonic_virus$treatment <- droplevels(neonic_virus)$treatment
neonic_combined$treatment <- droplevels(neonic_combined)$treatment
virus_combined$treatment <- droplevels(virus_combined)$treatment

cohen.d(control_neonic$virus_titre, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$virus_titre, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$virus_titre, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$virus_titre, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$virus_titre, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$virus_titre, virus_combined$treatment, na.rm = TRUE)











### queens ###
boxplot(queens$virus_titre ~ queens$treatment)
abline(vdt_queens, 0, lty = 2)
abline(vdt_queens_new, 0, lty = 2)

mod <- lm(virus_titre ~ treatment_1*treatment_2, data = queens)
baseline <- lm(virus_titre ~ 1, data = queens)
treatment_1_M <- update(baseline, .~. + treatment_1)
treatment_2_M <- update(treatment_1_M, .~. + treatment_2)
interaction_M <- update(treatment_2_M, .~. + treatment_1:treatment_2 )
anova(baseline, treatment_1_M, treatment_2_M, interaction_M)
Anova(mod)
Anova(interaction_M)

#pairwise comparison
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey") #or TukeyHSD(aov(mod))
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = dat_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

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
#homogenity of variance # plot(mod, 1)
leveneTest(mod) #non-significant = ok 
#normality  # plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> ok
#heteroskedastity 
bptest(mod) #non significant --> ok


#Pairwise Cohens D
#Effect size
levels(queens$treatment)
control_neonic <- subset(queens, queens$treatment == "Control" | queens$treatment == "Neonic")
control_virus <- subset(queens, queens$treatment == "Control" | queens$treatment == "Virus")
control_combined <- subset(queens, queens$treatment == "Control" | queens$treatment == "Interaction")
neonic_virus <- subset(queens, queens$treatment == "Neonic" | queens$treatment == "Virus")
neonic_combined <- subset(queens, queens$treatment == "Neonic" | queens$treatment == "Interaction")
virus_combined <- subset(queens, queens$treatment == "Virus" | queens$treatment == "Interaction")
control_neonic$treatment <- droplevels(control_neonic)$treatment
control_virus$treatment <- droplevels(control_virus)$treatment
control_combined$treatment <- droplevels(control_combined)$treatment
neonic_virus$treatment <- droplevels(neonic_virus)$treatment
neonic_combined$treatment <- droplevels(neonic_combined)$treatment
virus_combined$treatment <- droplevels(virus_combined)$treatment

cohen.d(control_neonic$virus_titre, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$virus_titre, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$virus_titre, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$virus_titre, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$virus_titre, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$virus_titre, virus_combined$treatment, na.rm = TRUE)




#caste not really affecting virus titer 
mod <- lmer(virus_titre ~ caste + (1|identity) + (1|treatment), data = data, REML = FALSE)
Anova(mod)
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
#homogenity of variance # plot(mod, 1)
leveneTest(residuals(mod) ~ data$caste[!is.na(data$virus_titre)]) #non-significant = ok 
boxplot(residuals(mod) ~ data$caste[!is.na(data$virus_titre)])
#normality  # plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> ok











#### neonicotinoid residues ####

#Overview
par(mfrow = c(1,2))
boxplot(dat_r$tq ~ dat_r$treatment)
boxplot(dat_r$tw ~ dat_r$treatment)
boxplot(dat_r$cq ~ dat_r$treatment)
boxplot(dat_r$cw ~ dat_r$treatment)

mean(dat_r$tq[dat_r$treatment == "Neonic"], na.rm = TRUE)
mean(dat_r$cq[dat_r$treatment == "Neonic"], na.rm = TRUE)
mean(dat_r$tq[dat_r$treatment == "Interaction"], na.rm = TRUE)
mean(dat_r$cq[dat_r$treatment == "Interaction"], na.rm = TRUE)
mean(dat_r$tq[dat_r$treatment == "Neonic"& dat_r$death_date == "31.10.17"], na.rm = TRUE)
mean(dat_r$cq[dat_r$treatment == "Neonic"& dat_r$death_date == "31.10.17"], na.rm = TRUE)
mean(dat_r$tq[dat_r$treatment == "Interaction" & dat_r$death_date == "31.10.17"], na.rm = TRUE)
mean(dat_r$cq[dat_r$treatment == "Interaction" & dat_r$death_date == "31.10.17"], na.rm = TRUE)

#Summary stats workers only
head(workers)
table(workers$treatment)
workers$tr <- workers$treatment
workers$tr <- relevel(workers$tr, "Neonic") 
workers$tr <- relevel(workers$tr, "Control")
workers$tr<- revalue(workers$tr, c("Interaction" = "Combined"))
table(workers$tr)

#workers thiamethoxam
DF <- group_by(workers, tr) %>%
  summarise(
    count = n(),
    N = sum(!is.na(thiamethoxam)),
    Pos_count = sum(thiamethoxam > 0, na.rm = TRUE),
    mean = mean(thiamethoxam, na.rm = TRUE),
    sd = sd(thiamethoxam, na.rm = TRUE),
    shapiro = shapiro.test(workers$thiamethoxam)$p,
    levene = leveneTest(workers$thiamethoxam, workers$treatment)$`Pr(>F)`[1],
    q_0 = quantile(thiamethoxam, probs = (0), na.rm = TRUE),
    q_25 = quantile(thiamethoxam, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(thiamethoxam, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(thiamethoxam, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(thiamethoxam, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)

#workers Clothianidin
DF <- group_by(workers, tr) %>%
  summarise(
    count = n(),
    N = sum(!is.na(clothianidin)),
    Pos_count = sum(clothianidin > 0, na.rm = TRUE),
    mean = mean(clothianidin, na.rm = TRUE),
    sd = sd(clothianidin, na.rm = TRUE),
    shapiro = shapiro.test(workers$clothianidin)$p,
    levene = leveneTest(workers$clothianidin, workers$treatment)$`Pr(>F)`[1],
    q_0 = quantile(clothianidin, probs = (0), na.rm = TRUE),
    q_25 = quantile(clothianidin, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(clothianidin, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(clothianidin, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(clothianidin, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)


#summary stats queens only

table(queens$treatment)
queens$tr <- queens$treatment
queens$tr <- relevel(queens$tr, "Neonic") 
queens$tr <- relevel(queens$tr, "Control")
queens$tr<- revalue(queens$tr, c("Interaction" = "Combined"))
table(queens$tr)



#queens thiamethoxam
DF <- group_by(queens, tr) %>%
  summarise(
    count = n(),
    N = sum(!is.na(thiamethoxam)),
    P_count = sum(thiamethoxam > 0, na.rm = TRUE),
    mean = mean(thiamethoxam, na.rm = TRUE),
    sd = sd(thiamethoxam, na.rm = TRUE),
    shapiro = shapiro.test(queens$thiamethoxam)$p,
    levene = leveneTest(queens$thiamethoxam, queens$treatment)$`Pr(>F)`[1],
    q_0 = quantile(thiamethoxam, probs = (0), na.rm = TRUE),
    q_25 = quantile(thiamethoxam, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(thiamethoxam, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(thiamethoxam, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(thiamethoxam, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)

#queens Clothianidin
DF <- group_by(queens, tr) %>%
  summarise(
    count = n(),
    N = sum(!is.na(clothianidin)),
    Pos_count = sum(clothianidin > 0, na.rm = TRUE),
    mean = mean(clothianidin, na.rm = TRUE),
    sd = sd(clothianidin, na.rm = TRUE),
    shapiro = shapiro.test(queens$clothianidin)$p,
    levene = leveneTest(queens$clothianidin, queens$treatment)$`Pr(>F)`[1],
    q_0 = quantile(clothianidin, probs = (0), na.rm = TRUE),
    q_25 = quantile(clothianidin, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(clothianidin, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(clothianidin, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(clothianidin, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)




head(dat_r)
boxplot(data$clothianidin ~ data$treatment*data$caste)
boxplot(data$thiamethoxam ~ data$treatment*data$caste)

require(dplyr)
group_by(caste_reduced, caste,treatment) %>%
  summarise(
    count = n(),
    N = sum(!is.na(thiamethoxam)),
    T_mean = mean(thiamethoxam, na.rm = TRUE), 
    T_sd = sd(thiamethoxam, na.rm = TRUE),
    N = sum(!is.na(clothianidin)), 
    C_mean = mean(clothianidin, na.rm = TRUE), 
    C_sd = sd(clothianidin, na.rm = TRUE)
  )



### Thiamethoxam ###
# all worker sample and all but one queen of the virus treatment are zero. Thus, the two levels are dropped for the analysis.

caste_reduced <- subset(data, treatment_1 != "Control")
caste_reduced$treatment <- droplevels(caste_reduced)$treatment
caste_reduced$treatment_1 <- droplevels(caste_reduced)$treatment
boxplot(thiamethoxam ~ treatment_2*caste, data = caste_reduced)
table(data_noNA$treatment)

data_noNA <- subset(caste_reduced, caste_reduced$thiamethoxam != "NA") #model can not process NA's!
m = art(thiamethoxam ~ caste * treatment_2 + (1|identity), data=data_noNA) # uses a linear mixed model
anova(m)

#pairwise comparisons: 
# An good option for cross-factor pairwise comparisons is to use either nonparametric Mann-Whitney U tests or Wilcoxon signed-rank tests on the original data.
# Use the former when subjects were in only one of the conditions being compared (i.e., between-subjects), 
# and use the latter when subjects were in both of the conditions being compared (i.e., within-subjects).

table(data_noNA$treatment_2)

library(reshape2)
df2 <- dcast(data_noNA, identity ~ treatment_2 + caste, value.var="thiamethoxam") # make wide-format table
uncorrected_pairwise_p_values <- c(
  cq_vs_cw = wilcox.test(df2$Control_q, df2$Control_w, paired=TRUE)$p.value,
  cq_vs_vq = wilcox.test(df2$Control_q, df2$Virus_q, paired=FALSE)$p.value,
  cq_vs_vw = wilcox.test(df2$Control_q, df2$Virus_w, paired=FALSE)$p.value,
  cw_vs_vq = wilcox.test(df2$Control_w, df2$Virus_q, paired=FALSE)$p.value,
  cw_vs_vw = wilcox.test(df2$Control_w, df2$Virus_w, paired=FALSE)$p.value,
  vq_vs_vw = wilcox.test(df2$Virus_q, df2$Virus_w, paired=TRUE)$p.value)

# correct for multiple comparisons using Holm's sequential Bonferroni procedure (Holm 1979)
p.adjust(uncorrected_pairwise_p_values, method="holm")

#Pairwise Cohens D
#Effect size
head(data_noNA)
wNeonic_wCombined <- subset(data_noNA, data_noNA$caste =="w") #1
wNeonic_qNeonic <- subset(data_noNA, data_noNA$treatment =="Neonic") #2
wNeonic_qCombined1 <- subset(data_noNA, data_noNA$caste == "w" & data_noNA$treatment =="Neonic")
wNeonic_qCombined2 <- subset(data_noNA, data_noNA$caste == "q" & data_noNA$treatment == "Interaction")
wNeonic_qCombined <- rbind(wNeonic_qCombined1, wNeonic_qCombined2) #3
wCombined_qNeonic1 <- subset(data_noNA, data_noNA$caste == "w" & data_noNA$treatment == "Interaction")
wCombined_qNeonic2 <- subset(data_noNA, data_noNA$caste == "q" & data_noNA$treatment == "Neonic")
wCombined_qNeonic <- rbind(wCombined_qNeonic1, wCombined_qNeonic2) #4
wCombined_qCombined <- subset(data_noNA, data_noNA$treatment =="Interaction") #5
qNeonic_qCombined <- subset(data_noNA, data_noNA$caste == "q") #6

cohen.d(wNeonic_wCombined$thiamethoxam, wNeonic_wCombined$treatment, na.rm = TRUE)
cohen.d(wNeonic_qNeonic$thiamethoxam, wNeonic_qNeonic$caste, na.rm = TRUE)
cohen.d(wNeonic_qCombined$thiamethoxam, wNeonic_qCombined$caste, na.rm = TRUE)
cohen.d(wCombined_qNeonic$thiamethoxam, wCombined_qNeonic$caste, na.rm = TRUE)
cohen.d(wCombined_qCombined$thiamethoxam, wCombined_qCombined$caste, na.rm = TRUE)
cohen.d(qNeonic_qCombined$thiamethoxam, qNeonic_qCombined$treatment, da.rm = TRUE)










### clothianidin ###
boxplot(clothianidin ~ treatment_2*caste, data = caste_reduced)
data_noNA <- subset(caste_reduced, caste_reduced$clothianidin != "NA") #model can not process NA's!
m = art(clothianidin ~ caste * treatment_2 + (1|identity), data=data_noNA) # uses a linear mixed model
anova(m)

df2 <- dcast(data_noNA, identity ~ treatment_2 + caste, value.var="clothianidin")
uncorrected_pairwise_p_values <- c(
  cq_vs_cw = wilcox.test(df2$Control_q, df2$Control_w, paired=TRUE)$p.value,
  cq_vs_vq = wilcox.test(df2$Control_q, df2$Virus_q, paired=FALSE)$p.value,
  cq_vs_vw = wilcox.test(df2$Control_q, df2$Virus_w, paired=FALSE)$p.value,
  cw_vs_vq = wilcox.test(df2$Control_w, df2$Virus_q, paired=FALSE)$p.value,
  cw_vs_vw = wilcox.test(df2$Control_w, df2$Virus_w, paired=FALSE)$p.value,
  vq_vs_vw = wilcox.test(df2$Virus_q, df2$Virus_w, paired=TRUE)$p.value)

# correct for multiple comparisons using Holm's sequential Bonferroni procedure (Holm 1979)
p.adjust(uncorrected_pairwise_p_values, method="holm")

#Pairwise Cohens D
#Effect size
cohen.d(wNeonic_wCombined$clothianidin, wNeonic_wCombined$treatment, na.rm = TRUE)
cohen.d(wNeonic_qNeonic$clothianidin, wNeonic_qNeonic$caste, na.rm = TRUE)
cohen.d(wNeonic_qCombined$clothianidin, wNeonic_qCombined$caste, na.rm = TRUE)
cohen.d(wCombined_qNeonic$clothianidin, wCombined_qNeonic$caste, na.rm = TRUE)
cohen.d(wCombined_qCombined$clothianidin, wCombined_qCombined$caste, na.rm = TRUE)
cohen.d(qNeonic_qCombined$clothianidin, qNeonic_qCombined$treatment, da.rm = TRUE)






#### analyse clothianidin in queens (or workers) vs amount of workers within high treatment ! ####
#not in manuscript
plot(dat_r$cq[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
mod <- lm(dat_r$cq[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
abline(mod)
summary(mod)
plot(dat_r$tq[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
mod <- lm(dat_r$tq[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
abline(mod)
summary(mod)
plot(dat_r$cw[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
mod <- lm(dat_r$cw[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
abline(mod)
summary(mod)
plot(dat_r$tw[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
mod <- lm(dat_r$tw[dat_r$treatment_1 == "Neonicotinoid"] ~ dat_r$y2_adults[dat_r$treatment_1 == "Neonicotinoid"])
abline(mod)
summary(mod)
plot(mod)



#### neonic ratio ####

caste_ratio <- subset(caste_reduced, clothianidin > 0 & thiamethoxam > 0)
caste_ratio$ratio <- caste_ratio$clothianidin/caste_ratio$thiamethoxam
table(caste_ratio$caste, caste_ratio$treatment)
boxplot(caste_ratio$ratio ~ caste_ratio$caste)
caste_ratio$treatment <- droplevels(caste_ratio)$treatment
caste_ratio$identity <- droplevels(caste_ratio)$identity

#modeling
baseline <- lmer(log10(ratio) ~ (1|identity), data = caste_ratio, REML = FALSE)
treatment_M <- lmer(log10(ratio) ~ (1|treatment)+ (1|identity), data = caste_ratio, REML = FALSE)
mod <- lmer(log10(ratio) ~ caste + (1|treatment) + (1|identity), data = caste_ratio, REML = FALSE)
anova(baseline, treatment_M ,mod)
Anova(mod)

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
#homogenity of variance # plot(mod, 1)
leveneTest(residuals(mod) ~ caste_ratio$caste) #non-significant = ok 
#normality  # plot(mod, 2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # non significant --> ok



#### Behavioural assays ####

#preparation of dataframe
head(mov)
str(mov)

boxplot(mov$number_of_stops~mov$treatment)
boxplot(mov$active_time~mov$treatment)
boxplot(mov$overall_movement~mov$treatment)
boxplot(mov$average_speed~mov$treatment)
boxplot(mov$initial_speed~mov$treatment)

mov$treatment


mov_r <- subset(mov, treatment_1 != "low") #removed  because it does not really add valuable information to the study
mov_r$treatment_1 <- droplevels(mov_r)$treatment_1
mov_r$treatment <- droplevels(mov_r)$treatment
mov_r$id <-  droplevels(mov_r)$id
mov_r$treatment <- as.factor(mov_r$treatment)
mov_r$treatment <- revalue(mov_r$treatment, c("1"="Control", "2"="Virus","5" = "Neonic","6" = "Interaction")) 
#mov_r$treatment <- factor(mov_r$treatment, levels(mov_r$treatment)[c(1,3,2,4)]) # we did not revalue here
#mov_r$treatment <- revalue(mov_r$treatment, c("control-control"="Control", "control-virus"="Virus","high-control" = "Neonicotinoid","high-virus" = "Mixed" ))
mov_r$treatment_1 <- revalue(mov_r$treatment_1, c("control" = "Control", "high" = "Neonicotinoid"))
mov_r$treatment_2 <- revalue(mov_r$treatment_2, c("control" = "Control", "virus" = "Virus"))
mov_r$time <- times(mov_r$time)
mov_r$time_minutes <- 60*hours(mov_r$time) + minutes(mov_r$time)
table(mov_r$treatment)
table(mov_r$treatment_2)
table(mov_r$treatment_1)
table(mov_r$id)

mov_r$tr <- mov_r$treatment
mov_r$tr <- relevel(mov_r$tr, "Neonic") 
mov_r$tr <- relevel(mov_r$tr, "Control")
mov_r$tr<- revalue(mov_r$tr, c("Interaction" = "Combined"))
table(mov_r$tr)

# Define time inactive!
mov_r$time_inactive <- 120-mov_r$active_time
mov_r$time_inactive <- round(mov_r$time_inactive, digits = 0)
hist(mov_r$time_inactive)
boxplot(mov_r$time_inactive ~ mov_r$treatment)


#### overall movement ####
par(mfrow = c(1,1))
boxplot(mov_r$overall_movement ~ mov_r$datum) #day might have and effect and will be included in the model as random factor
plot(mov_r$overall_movement ~ mov_r$time_minutes)
abline(lm(mov_r$overall_movement ~ mov_r$time_minutes))              # time of day not relevant
hist(mov_r$overall_movement)
boxplot(mov_r$overall_movement ~ mov_r$treatment)

head(mov_r)


# summary stats
DF <- group_by(mov_r, tr) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(overall_movement)), 
    mean = mean(overall_movement, na.rm = TRUE),
    sd = sd(overall_movement, na.rm = TRUE),
    shapiro = shapiro.test(mov_r$overall_movement)$p,
    levene = leveneTest(mov_r$overall_movement, mov_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(overall_movement, probs = (0), na.rm = TRUE),
    q_25 = quantile(overall_movement, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(overall_movement, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(overall_movement, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(overall_movement, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)

mean(mov_r$overall_movement[mov_r$treatment=="Control"])

baseline <- lmer(overall_movement ~ (1|datum), data=mov_r, REML=FALSE)
colony_M <- lmer(overall_movement ~ (1|id) + (1|datum), data=mov_r, REML=FALSE)
treatment_1_M <- lmer(overall_movement ~ treatment_1 + (1|id) + (1|datum), data=mov_r, REML=FALSE)
treatment_2_M <- lmer(overall_movement ~ treatment_1 + treatment_2 + (1|id) + (1|datum), data=mov_r, REML=FALSE)
interaction_M <- lmer(overall_movement ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|datum), data=mov_r, REML=FALSE)
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)
lrtest(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)
mod <- interaction_M
Anova(mod)

summary(mod)
estimates <- fixef(mod)
estimates

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1*treatment_2, data = mov_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# test model assumtions
compareqqnorm(mod)
leveneTest(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2) # homogenity of variance ok
boxplot(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # normality of residues is ok
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$datum[,1])
qqline(ranef(mod)$datum[,1])
boxplot(resid(mod)~mov_r$datum)
par(mfrow=c(1,1))

#pairwise comparisons Effect size 
control_neonic <- subset(mov_r, mov_r$treatment == "Control" | mov_r$treatment == "Neonic")
control_virus <- subset(mov_r, mov_r$treatment == "Control" | mov_r$treatment == "Virus")
control_combined <- subset(mov_r, mov_r$treatment == "Control" | mov_r$treatment == "Interaction")
neonic_virus <- subset(mov_r, mov_r$treatment == "Neonic" | mov_r$treatment == "Virus")
neonic_combined <- subset(mov_r, mov_r$treatment == "Neonic" | mov_r$treatment == "Interaction")
virus_combined <- subset(mov_r, mov_r$treatment == "Virus" | mov_r$treatment == "Interaction")
control_neonic$treatment <- droplevels(control_neonic)$treatment
control_virus$treatment <- droplevels(control_virus)$treatment
control_combined$treatment <- droplevels(control_combined)$treatment
neonic_virus$treatment <- droplevels(neonic_virus)$treatment
neonic_combined$treatment <- droplevels(neonic_combined)$treatment
virus_combined$treatment <- droplevels(virus_combined)$treatment

cohen.d(control_neonic$overall_movement, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$overall_movement, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$overall_movement, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$overall_movement, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$overall_movement, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$overall_movement, virus_combined$treatment, na.rm = TRUE)










#### active_time - time_inactive ####
mov_r$time_inactive <- 120-mov_r$active_time
mov_r$time_inactive <- round(mov_r$time_inactive, digits = 0)
hist(mov_r$time_inactive)
boxplot(mov_r$time_inactive ~ mov_r$treatment)

# summary stats
DF <- group_by(mov_r, tr) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(time_inactive)), 
    mean = mean(time_inactive, na.rm = TRUE),
    sd = sd(time_inactive, na.rm = TRUE),
    shapiro = shapiro.test(mov_r$time_inactive)$p,
    levene = leveneTest(mov_r$time_inactive, mov_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(time_inactive, probs = (0), na.rm = TRUE),
    q_25 = quantile(time_inactive, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(time_inactive, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(time_inactive, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(time_inactive, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)

mean(mov_r$time_inactive[mov_r$treatment=="Control"])



boxplot(mov_r$time_inactive ~ mov_r$datum) #day might have and effect and will be included in the model as random factor
plot(mov_r$time_inactive ~ mov_r$time_minutes)
abline(lm(mov_r$time_inactive ~ mov_r$time_minutes)) # time of day seems not relevant and will be neglected in the model

baseline <- glmer(time_inactive ~ (1|datum), data=mov_r, family = "poisson")
colony_M <- glmer(time_inactive ~ (1|id) + (1|datum), data=mov_r, family = "poisson")
treatment_1_M <- glmer(time_inactive ~ treatment_1 + (1|id) + (1|datum), data=mov_r, family = "poisson")
treatment_2_M <- glmer(time_inactive ~ treatment_1 + treatment_2 + (1|id) + (1|datum), data=mov_r, family = "poisson")
interaction_M <- glmer(time_inactive ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|datum), data=mov_r, family = "poisson")
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)
lrtest(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)
mod <- interaction_M
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

# test model assumtions
compareqqnorm(mod)
leveneTest(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2) # homogenity of variance ok
boxplot(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # normality of residues is ok
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$datum[,1])
qqline(ranef(mod)$datum[,1])
boxplot(resid(mod)~mov_r$datum)
par(mfrow=c(1,1))

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1*treatment_2, data = mov_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# Pairwise effect size differences
cohen.d(control_neonic$time_inactive, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$time_inactive, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$time_inactive, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$time_inactive, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$time_inactive, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$time_inactive, virus_combined$treatment, na.rm = TRUE)






#### average speed ####
boxplot(mov_r$average_speed ~ mov_r$datum) #day might have and effect and will be included in the model as random factor
plot(mov_r$average_speed ~ mov_r$time_minutes)
abline(lm(mov_r$average_speed ~ mov_r$time_minutes))              # time of day not relevant
hist(mov_r$average_speed)
boxplot(mov_r$average_speed ~ mov_r$treatment)
stripchart(average_speed ~ treatment, data = mov_r, 
           vertical = TRUE, at = c(1:4) ,method = "jitter", 
           pch = c(1), col= "grey78",
           add = TRUE, cex= 0.8)

# summary stats
DF <- group_by(mov_r, tr) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(average_speed)), 
    mean = mean(average_speed, na.rm = TRUE),
    sd = sd(average_speed, na.rm = TRUE),
    shapiro = shapiro.test(mov_r$average_speed)$p,
    levene = leveneTest(mov_r$average_speed, mov_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(average_speed, probs = (0), na.rm = TRUE),
    q_25 = quantile(average_speed, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(average_speed, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(average_speed, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(average_speed, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)
mean(mov_r$average_speed[mov_r$treatment=="Control"])

baseline <- lmer(average_speed ~ (1|datum), data=mov_r, REML=FALSE)
colony_M <- lmer(average_speed ~ (1|id) + (1|datum), data=mov_r, REML=FALSE)
treatment_1_M <- lmer(average_speed ~ treatment_1 + (1|id) + (1|datum), data=mov_r, REML=FALSE)
treatment_2_M <- lmer(average_speed ~ treatment_1 + treatment_2 + (1|id) + (1|datum), data=mov_r, REML=FALSE)
interaction_M <- lmer(average_speed ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|datum), data=mov_r, REML=FALSE)
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)
lrtest(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)
mod <- interaction_M
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = mov_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# test model assumtions
compareqqnorm(mod)
leveneTest(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2) # homogenity of variance ok
boxplot(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # normality of residues is ok
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$datum[,1])
qqline(ranef(mod)$datum[,1])
boxplot(resid(mod)~mov_r$datum)
plot(mod)
par(mfrow=c(1,1))

# Pairwise effect size differences
cohen.d(control_neonic$average_speed, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$average_speed, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$average_speed, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$average_speed, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$average_speed, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$average_speed, virus_combined$treatment, na.rm = TRUE)







#### initial speed ####
boxplot(mov_r$initial_speed ~ mov_r$datum) #day might have and effect and will be included in the model as random factor
plot(mov_r$initial_speed ~ mov_r$time_minutes)
abline(lm(mov_r$initial_speed ~ mov_r$time_minutes))              # time of day not relevant
hist(mov_r$initial_speed)
boxplot(mov_r$initial_speed ~ mov_r$treatment)
stripchart(initial_speed ~ treatment, data = mov_r, 
           vertical = TRUE, at = c(1:4) ,method = "jitter", 
           pch = c(1), col= "grey78",
           add = TRUE, cex= 0.8)

# summary stats
DF <- group_by(mov_r, tr) %>%
  summarise(
    count = n(),
    N_W = sum(!is.na(initial_speed)), 
    mean = mean(initial_speed, na.rm = TRUE),
    sd = sd(initial_speed, na.rm = TRUE),
    shapiro = shapiro.test(mov_r$initial_speed)$p,
    levene = leveneTest(mov_r$initial_speed, mov_r$treatment)$`Pr(>F)`[1],
    q_0 = quantile(initial_speed, probs = (0), na.rm = TRUE),
    q_25 = quantile(initial_speed, probs = (0.25), na.rm = TRUE),
    q_50 = quantile(initial_speed, probs = (0.5), na.rm = TRUE),
    q_75 = quantile(initial_speed, probs = (0.75), na.rm = TRUE),
    q_100 = quantile(initial_speed, probs = (1), na.rm = TRUE),
  )
as.data.frame(DF)
mean(mov_r$initial_speed[mov_r$treatment=="Control"])


baseline <- glmer(initial_speed ~ (1|datum), data=mov_r, family = "poisson")
colony_M <- glmer(initial_speed ~ (1|id) + (1|datum), data=mov_r, family = "poisson")
treatment_1_M <- glmer(initial_speed ~ treatment_1 + (1|id) + (1|datum), data=mov_r, family = "poisson")
treatment_2_M <- glmer(initial_speed ~ treatment_1 + treatment_2 + (1|id) + (1|datum), data=mov_r, family = "poisson")
interaction_M <- glmer(initial_speed ~ treatment_1 + treatment_2 + treatment_1:treatment_2 + (1|id) + (1|datum), data=mov_r, family = "poisson")
anova(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)
lrtest(baseline, colony_M, treatment_1_M, treatment_2_M, interaction_M)

mod <- interaction_M
Anova(mod)
summary(mod)
estimates <- fixef(mod)
estimates

# test model assumtions
compareqqnorm(mod)
leveneTest(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2) # homogenity of variance ok
boxplot(residuals(mod) ~ mov_r$treatment_1*mov_r$treatment_2)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # normality of residues is not fullfilled but also not expected from a glmer
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
qqnorm(ranef(mod)$id[,1]) # qq of random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$datum[,1])
qqline(ranef(mod)$datum[,1])
boxplot(resid(mod)~mov_r$datum)
plot(mod)
par(mfrow=c(1,1))

#pairwise differences
lsmeans(mod, pairwise ~ treatment_1:treatment_2, adjust = "tukey")
marginal = lsmeans(mod, ~ treatment_1:treatment_2, data = mov_r)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
CLD

# Pairwise effect size differences
cohen.d(control_neonic$initial_speed, control_neonic$treatment, na.rm = TRUE)
cohen.d(control_virus$initial_speed, control_virus$treatment, na.rm = TRUE)
cohen.d(control_combined$initial_speed, control_combined$treatment, na.rm = TRUE)
cohen.d(neonic_virus$initial_speed, neonic_virus$treatment, na.rm = TRUE)
cohen.d(neonic_combined$initial_speed, neonic_combined$treatment, na.rm = TRUE)
cohen.d(virus_combined$initial_speed, virus_combined$treatment, na.rm = TRUE)











#### graphs for the manuscript ####
par(mar=c(4.1, 4.1, 3.1, 1.1),
    oma = c(0,0,0,0)+0.5,
    cex.lab=1.5, 
    cex.axis=1.5, 
    cex.main= 1.5,
    mfrow=c(1,2))
label <- c("Control", "Virus", "Neonicotinoid", "Combined")

#### Graph Bodymass #### 
#Bodymass exported as 1200 X 880  (later ev. 10x8 inch pdfs)
boxplot(dat_r$weight_adults ~ dat_r$treatment, main = "Workers", cex.main = 1.5 ,xlab = "Treatments", ylab="Mass [mg]", ylim=c(14.5, 25), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 24.5, labels = c("a", "ab", "b", "b"), font = 2, cex = 1.5)
mtext('a', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(dat_r$weight_queens ~ dat_r$treatment, main = "Queens", xlab = "Treatments", ylab="Mass [mg]", ylim=c(22, 47), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 45.6, labels = c("a", "a", "a", "a"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)

#### graph colonysize Y1 + Y2 #### exportet 1000X880
# adults + brood
par(mfrow = c(2,2))
label_y1 <- c("Control", "Neonicotinoid")

boxplot(dat_r$y1_adults ~ dat_r$treatment_1, main = "Week 13", cex.main = 1.5, xlab = "Treatments", ylab = "Number of workers", ylim=c(0, 26), las=1, xaxt = "n")
axis(1, at=1:2, labels = label_y1, cex.axis = 1.15)
text(c(1:2), 25, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('a', side=3, line=1.5, at=0.2, cex = 2, font = 2)
boxplot(dat_r$y2_adults ~ dat_r$treatment, main = "Week 64", cex.main = 1.5, xlab = "Treatments", ylab = "Number of workers", ylim=c(80, 420), las=1, xaxt = "n", yaxt = "n")
axis(2, cex.axis = 1.25, las = 1)
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 410, labels = c("a", "b", "b", "b"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(dat_r$y1_brood ~ dat_r$treatment_1, main = "Week 13", xlab = "Treatments", ylab = "Amount of brood", ylim=c(0, 75), las=1, xaxt = "n")
axis(1, at=1:2, labels = label_y1, cex.axis = 1.15)
text(c(1:2), 65, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('c', side=3, line=1.5, at=0.2, cex = 2, font = 2)
boxplot(dat_r$y2_brood ~ dat_r$treatment, main = "Week 64", xlab = "Treatments", ylab = "Amount of brood", ylim=c(0, 220), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 210, labels = c("a", "ab", "b", "b"), font = 2, cex = 1.5)
mtext('d', side=3, line=1.5, at=0, cex = 2, font = 2)

#### graph virus titers ####
# export as 1200 x 880 to clipboard
boxplot(workers$virus_titre ~ workers$treatment, main = "Workers", xlab = "Treatments", ylab="Genomic ABPV copies / mg tissue [log]", ylim=c(1, 6), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 5.8, labels = c("a", "b", "a", "b"), font = 2, cex = 1.5)
abline(vdt_workers_new, 0, lty = 2)
mtext('a', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(queens$virus_titre ~ queens$treatment, main = "Queens", xlab = "Treatments", ylab="Genomic ABPV copies / mg tissue [log]", ylim=c(1.5, 5), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 4.9, labels = c("a", "b", "ab", "b"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)
abline(vdt_queens_new, 0, lty = 2)

#### graph neonic residues ####
label_2 <- c("Neonicotinoid", "Combined")
boxplot(thiamethoxam ~ caste*treatment_2, data = caste_reduced, main = "", xlab = "Treatments", ylab = "ng thiamethoxam / g dryweight",
        col  = c("white","grey"), ylim=c(0, 2.7), las=1, xaxt = "n", at = c(0.8, 1.8, 3.2, 4.2))
axis(1, at = c(1.3,3.7), labels=label_2, cex.axis = 1.15)
text(c(0.8, 1.8, 3.2, 4.2), 2.6, labels = c("a", "a", "a", "a"), font = 2, cex = 1.5)
mtext('a', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(clothianidin ~ caste*treatment_2, data = caste_reduced, main = "", xlab = "Treatments", ylab = "ng clothianidin / g dryweight",
        col  = c("white","grey"), ylim=c(0, 4), las=1, xaxt = "n", at = c(0.8, 1.8, 3.2, 4.2))
axis(1, at = c(1.3,3.7), labels=label_2, cex.axis = 1.15)
text(c(0.8, 1.8, 3.2, 4.2), 3.8, labels = c("a", "b", "a", "b"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)


#### graph behavioural assay ####
par(mfrow = c(2,2))
boxplot(mov_r$overall_movement ~ mov_r$treatment, main = "", xlab = "Treatments", ylab="Distance [~cm]", ylim=c(0, 480), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 470, labels = c("a", "ab", "a", "b"), font = 2, cex = 1.5)
mtext('a', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(mov_r$time_inactive ~ mov_r$treatment, main = "", xlab = "Treatments", ylab="Time inactive [s]", ylim=c(0, 140), las=1, xaxt = "n", yaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
axis(2, at=c(0, 20, 40, 60, 80, 100, 120), labels = c(0, 20, 40, 60, 80, 100, 120), las = 1)
text(c(1:4), 130, labels = c("a", "a", "a", "a"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(mov_r$average_speed ~ mov_r$treatment, main = "", xlab = "Treatments", ylab="Average speed [~cm/s]", ylim=c(0, 4.8), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 4.5, labels = c("ab", "bc", "a", "c"), font = 2, cex = 1.5)
mtext('c', side=3, line=1.5, at=0, cex = 2, font = 2)
boxplot(mov_r$initial_speed ~ mov_r$treatment, main = "", xlab = "Treatments", ylab="Initial speed [~cm/s]", ylim=c(0, 75), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 71, labels = c("a", "b", "a", "b"), font = 2, cex = 1.5)
mtext('d', side=3, line=1.5, at=0, cex = 2, font = 2)



#### NEW GRAPHS ####
par(mar=c(4.1, 4.1, 3.1, 1.1),
    oma = c(0,0,0,0)+0.5,
    cex.lab=1.5, 
    cex.axis=1.5, 
    cex.main= 1.5,
    mfrow=c(1,2))
label <- c("Control", "Virus", "Neonicotinoid", "Combined")
label_2 <- c("Neonicotinoid", "Combined")

par(mfrow = c(3,2))

boxplot(dat_r$weight_adults ~ dat_r$treatment, main = "Worker bodymass", cex.main = 1.5 ,xlab = "Treatments", ylab="Mass [mg]", ylim=c(14.5, 25), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 24.5, labels = c("a", "ab", "b", "b"), font = 2, cex = 1.5)
mtext('a', side=3, line=1.5, at=0, cex = 2, font = 2)

boxplot(dat_r$weight_queens ~ dat_r$treatment, main = "Queen bodymass", xlab = "Treatments", ylab="Mass [mg]", ylim=c(22, 47), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 45.6, labels = c("a", "a", "a", "a"), font = 2, cex = 1.5)
mtext('b', side=3, line=1.5, at=0, cex = 2, font = 2)

boxplot(workers$virus_titre ~ workers$treatment, main = "Virus titres - Workers", xlab = "Treatments", ylab="Genomic ABPV copies / mg tissue [log]", ylim=c(1, 6), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 5.8, labels = c("a", "b", "a", "b"), font = 2, cex = 1.5)
abline(vdt_workers_new, 0, lty = 2)
mtext('c', side=3, line=1.5, at=0, cex = 2, font = 2)

boxplot(queens$virus_titre ~ queens$treatment, main = "Virus titres - Queens", xlab = "Treatments", ylab="Genomic ABPV copies / mg tissue [log]", ylim=c(1.5, 5), las=1, xaxt = "n")
axis(1, at=1:4, labels = label, cex.axis = 1.15)
text(c(1:4), 4.9, labels = c("a", "b", "ab", "b"), font = 2, cex = 1.5)
mtext('d', side=3, line=1.5, at=0, cex = 2, font = 2)
abline(vdt_queens_new, 0, lty = 2)

boxplot(thiamethoxam ~ caste*treatment_2, data = caste_reduced, main = "Thiamethoxam", xlab = "Treatments", ylab = "ng thiamethoxam / g dryweight",
        col  = c("white","grey"), ylim=c(0, 2.7), las=1, xaxt = "n", at = c(0.7, 1.7, 3.3, 4.3))
axis(1, at = c(1.2,3.8), labels=label_2, cex.axis = 1.15)
text(c(0.7, 1.7, 3.3, 4.3), 2.6, labels = c("a", "a", "a", "a"), font = 2, cex = 1.5)
mtext('e', side=3, line=1.5, at=-0.35, cex = 2, font = 2)
legend(0.1, 2.4, legend=c("Queens", "Workers"), col=c("white", "grey"), fill = c("white", "grey"), box.lty = 0, cex= 1.15)

boxplot(clothianidin ~ caste*treatment_2, data = caste_reduced, main = "Clothianidin", xlab = "Treatments", ylab = "ng clothianidin / g dryweight",
        col  = c("white","grey"), ylim=c(0, 4), las=1, xaxt = "n", at = c(0.7, 1.7, 3.3, 4.3))
axis(1, at = c(1.2,3.8), labels=label_2, cex.axis = 1.15)
text(c(0.7, 1.7, 3.3, 4.3), 3.8, labels = c("a", "b", "a", "b"), font = 2, cex = 1.5)
mtext('f', side=3, line=1.5, at=-0.35, cex = 2, font = 2)
legend(3.7, 3.5, legend=c("Queens", "Workers"), col=c("white", "grey"), fill = c("white", "grey"), box.lty = 0, cex = 1.15)

# saved as 1000x1320 jpeg
