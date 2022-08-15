### ### ### ### ### ### ### ### ### ### ### ### ### ###
### HB-VIS - Honey bee viruses in Spiders   ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### load libraries ####
suppressPackageStartupMessages(library(survival))# contains kalan meier plot function
suppressPackageStartupMessages(library(plyr)) #contains revalue
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(chron))
suppressPackageStartupMessages(library(car))#contains levene's test
suppressPackageStartupMessages(library(dunn.test)) #contains the dunn test
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

#### define functions used later ####

### campell2x2 optimized chi test chi-test ###
# compare to https://measuringu.com/what-test/   &   https://measuringu.com/ab-cal/
campbell2x2.test <- function(t) {
  min_exp_val <- min(colSums(t))*min(rowSums(t))/sum(t)
  
  if (min_exp_val < 1) {
    # in Campbell's naming: Fisher–Irwin test by Irwin’s rule
    result <- fisher.test(t, alternative = "two.sided")
    result$method <- paste("Optimal 2x2 test according to Campbell(2007) recommendation\n\n",
                           paste("Minimal expected cell count: ", round(min_exp_val, 3), "\n\n", sep = ""),
                           paste("Performing", result$method)
    )
    return(result)
  } else {
    #  'N − 1' Pearson's Chi-squared test
    n1chisq.test <- function(t) {
      chisqtst <- chisq.test(t, correct = FALSE)
      N <- sum(chisqtst$observed)
      chisqtst$statistic = ((N-1)/N) * chisqtst$statistic
      chisqtst$p.value <- 1 - pchisq(chisqtst$statistic, chisqtst$parameter)
      chisqtst$method <- paste("'N-1'", chisqtst$method)
      return(chisqtst)
    }
    result <- n1chisq.test(t)
    result$method <- paste("Optimal 2x2 test according to Campbell(2007) recommendation\n\n",
                           paste("Minimal expected cell count: ", round(min_exp_val, 3), "\n\n", sep = ""),
                           paste("Performing", result$method)
    )
    return(result)
  }
}

n1chisq.test <- function(t) {
  chisqtst <- chisq.test(t, correct = FALSE)
  N <- sum(chisqtst$observed)
  chisqtst$statistic = ((N-1)/N) * chisqtst$statistic
  chisqtst$p.value <- 1 - pchisq(chisqtst$statistic, chisqtst$parameter)
  chisqtst$method <- paste("'N-1'", chisqtst$method)
  return(chisqtst)
}



#### prerequisites ####
# setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV") #office computer
setwd("/Users/gismo/Desktop/R/spiders") # homeoffice mac
dat <- read.table("viruses_in_spiders.txt", header = TRUE)


head(dat)
str(dat)
dat$treatment <- revalue(dat$treatment, c("control"="Control", "virus"="Virus"))
table(dat$treatment)
levels(dat$sample_id)
table(dat$sample_id)


# Virus detection treshold 
vdt_ABPV <- log10(7.45E+03) # corresponds to the treshold for 1 body part, but it needs to be doubled to stand for the full body? check later in ghraphs!
vdt_ABPV_full <- log10(2*7.45E+03) # 
vdt_DWV <- log10(2.40E+04)
vdt_DWV_full <- log10(2*2.40E+04)

# preparation of titres used in analyses: substract detection treshold and apply 0 to negative samples
dat$P_ABPV_titer_cor <- dat$P_ABPV_titer - 7.45E+03
dat$O_ABPV_titer_cor <- dat$O_ABPV_titer - 7.45E+03
dat$P_DWV_titer_cor <- dat$P_DWV_titer - 2.40E+04
dat$O_DWV_titer_cor <- dat$O_DWV_titer - 2.40E+04

dat$P_ABPV_titer_cor[dat$P_ABPV_PN == "neg"] <- 0 
dat$O_ABPV_titer_cor[dat$O_ABPV_PN == "neg"] <- 0
dat$P_DWV_titer_cor[dat$P_DWV_PN == "neg"] <- 0 
dat$O_DWV_titer_cor[dat$O_DWV_PN == "neg"] <- 0

#logtransform plus 1 to avoid log10(0) = -inf
dat$ABPV_pro <- log10(dat$P_ABPV_titer_cor +1)
dat$ABPV_opi <- log10(dat$O_ABPV_titer_cor +1)
dat$DWV_pro <- log10(dat$P_DWV_titer_cor +1)
dat$DWV_opi <- log10(dat$O_DWV_titer_cor +1)
# titres of the samples when combining the pro- and opisthosoma
dat$ABPV_full <- log10(dat$P_ABPV_titer_cor + dat$O_ABPV_titer_cor + 1)
dat$DWV_full <- log10(dat$P_DWV_titer_cor + dat$O_DWV_titer_cor +1)
# titre as a sum of both virus titres for the entire body
dat$virus_full <- log10(dat$P_ABPV_titer_cor + dat$O_ABPV_titer_cor+dat$P_DWV_titer_cor + dat$O_DWV_titer_cor +1)

#### descriptive statistics ####
#### ABPV ####
# prosoma
group_by(dat, treatment) %>%
  summarise(
    count = n(),
    pos = sum(P_ABPV_PN =="pos"),
    neg = sum(P_ABPV_PN == "neg"),
    p_percent = sum(P_ABPV_PN =="pos")/n(),
    mean_tit = mean(ABPV_pro[P_ABPV_PN == "pos"]), 
    sd_tit = sd(ABPV_pro[P_ABPV_PN == "pos"]),
    q_0 = quantile(ABPV_pro[P_ABPV_PN == "pos"], probs = (0), na.rm = TRUE),
    q_25 = quantile(ABPV_pro[P_ABPV_PN == "pos"], probs = (0.25), na.rm = TRUE),
    q_50 = quantile(ABPV_pro[P_ABPV_PN == "pos"], probs = (0.5), na.rm = TRUE),
    q_75 = quantile(ABPV_pro[P_ABPV_PN == "pos"], probs = (0.75), na.rm = TRUE),
    q_100 = quantile(ABPV_pro[P_ABPV_PN == "pos"], probs = (1), na.rm = TRUE),
    HLI = sum(ABPV_pro >= 7)
  )

# opisthosoma
group_by(dat, treatment) %>%
  summarise(
    count = n(),
    pos = sum(O_ABPV_PN =="pos"),
    neg = sum(O_ABPV_PN == "neg"),
    p_percent = sum(O_ABPV_PN =="pos")/n(),
    mean_tit = mean(ABPV_opi[O_ABPV_PN == "pos"]), 
    sd_tit = sd(ABPV_opi[O_ABPV_PN == "pos"]),
    q_0 = quantile(ABPV_opi[O_ABPV_PN == "pos"], probs = (0), na.rm = TRUE),
    q_25 = quantile(ABPV_opi[O_ABPV_PN == "pos"], probs = (0.25), na.rm = TRUE),
    q_50 = quantile(ABPV_opi[O_ABPV_PN == "pos"], probs = (0.5), na.rm = TRUE),
    q_75 = quantile(ABPV_opi[O_ABPV_PN == "pos"], probs = (0.75), na.rm = TRUE),
    q_100 = quantile(ABPV_opi[O_ABPV_PN == "pos"], probs = (1), na.rm = TRUE),
    HLI = sum(ABPV_opi >= 7)
  )

#### DWV ####
# prosoma
group_by(dat, treatment) %>%
  summarise(
    count = n(),
    pos = sum(P_DWV_PN =="pos"),
    neg = sum(P_DWV_PN == "neg"),
    p_percent = sum(P_DWV_PN =="pos")/n(),
    mean_tit = mean(DWV_pro[P_DWV_PN == "pos"]), 
    sd_tit = sd(DWV_pro[P_DWV_PN == "pos"]),
    q_0 = quantile(DWV_pro[P_DWV_PN == "pos"], probs = (0), na.rm = TRUE),
    q_25 = quantile(DWV_pro[P_DWV_PN == "pos"], probs = (0.25), na.rm = TRUE),
    q_50 = quantile(DWV_pro[P_DWV_PN == "pos"], probs = (0.5), na.rm = TRUE),
    q_75 = quantile(DWV_pro[P_DWV_PN == "pos"], probs = (0.75), na.rm = TRUE),
    q_100 = quantile(DWV_pro[P_DWV_PN == "pos"], probs = (1), na.rm = TRUE),
    HLI = sum(DWV_pro >= 7)
  )
# opisthosoma
group_by(dat, treatment) %>%
  summarise(
    count = n(),
    pos = sum(O_DWV_PN =="pos"),
    neg = sum(O_DWV_PN == "neg"),
    pos_per = sum(O_DWV_PN =="pos")/n(),
    mean_tit = mean(DWV_opi[O_DWV_PN == "pos"]), 
    sd_tit = sd(DWV_opi[O_DWV_PN == "pos"]),
    q_0 = quantile(DWV_opi[O_DWV_PN == "pos"], probs = (0), na.rm = TRUE),
    q_25 = quantile(DWV_opi[O_DWV_PN == "pos"], probs = (0.25), na.rm = TRUE),
    q_50 = quantile(DWV_opi[O_DWV_PN == "pos"], probs = (0.5), na.rm = TRUE),
    q_75 = quantile(DWV_opi[O_DWV_PN == "pos"], probs = (0.75), na.rm = TRUE),
    q_100 = quantile(DWV_opi[O_DWV_PN == "pos"], probs = (1), na.rm = TRUE),
    HLI = sum(DWV_opi >= 7),
    min = min(DWV_opi[O_DWV_PN == "pos"]),
    max = max(DWV_opi[O_DWV_PN == "pos"])
  )

# Full titres
group_by(dat) %>% 
  summarise(
    ABPV_min = min(ABPV_full), 
    ABPV_max = max(ABPV_full), 
    DWV_min = min(DWV_full[DWV_full!=0], na.rm = TRUE), 
    DWV_max = max(DWV_full, na.rm = TRUE) 
  )


#### Tests ####
#### comparison of treatments with regards of virus titres

#ABPV
plot(dat$ABPV_full ~ dat$treatment, ylim = c(0,11))
abline(vdt_ABPV_full, 0, lty = 2)
shapiro.test(dat$ABPV_full) # significant --> not normally distributed
leveneTest(dat$ABPV_full, dat$treatment) # non significant --> nullhypothesis of equal variance can be assumed
wilcox.test(dat$ABPV_full ~ dat$treatment)

#DWV
plot(dat$DWV_full ~ dat$treatment, ylim = c(0,11))
abline(vdt_DWV_full, 0, lty = 2)
shapiro.test(dat$DWV_full) # significant --> not normally distributed
leveneTest(dat$DWV_full, dat$treatment) # non significant --> nullhypothesis of equal variance can be assumed
wilcox.test(dat$DWV_full ~ dat$treatment)

#### High level infection (HLI)
control_total <- sum(dat$treatment=="Control") #count(dat, treatment)[1,2] #total number of samples in control group
treatment_total <- sum(dat$treatment=="Virus") #count(dat, treatment)[2,2] #total number of samples in treatment group
control_hli <- sum(dat$treatment == "Control" & (dat$ABPV_pro >= 7 | dat$ABPV_opi >= 7 | dat$DWV_pro >= 7 | dat$DWV_opi >= 7)) #number of cases (spiders with at least one body part having a high level infection with at least one of the two viruses) in control group
treatment_hli <- sum(dat$treatment == "Virus" & (dat$ABPV_pro >= 7 | dat$ABPV_opi >= 7 | dat$DWV_pro >= 7 | dat$DWV_opi >= 7)) #number of cases (high level infections) in treatment group
values <- c(control_hli,
            treatment_hli,
            control_total - control_hli,
            treatment_total - treatment_hli) #create the vector with the values
chi_matrix <- matrix(values,  nrow = 2) # create a matrix from the vector
chi_matrix
# alternatively create a more informative datatable
treatment = c("control", "virus treatment", "control", "virus treatment")
hli = c("yes", "yes", "no", "no")
N = values #c(4, 7, 1, 12)
data <- data.frame(treatment, hli, N)
#create crosstabulation
highlevelinfections <- xtabs(N ~ treatment + hli, data = data)
highlevelinfections
#all the same just for testing
#campbell2x2.test(matrix(c(1, 12, 4, 7), nrow = 2))
campbell2x2.test(highlevelinfections)
#campbell2x2.test(chi_matrix)
#campbell2x2.test(matrix(values, nrow = 2))
#Output for the one tailed campell test: 
value <- as.numeric(campbell2x2.test(highlevelinfections)$p.value/2)
paste("p-value for one tailed campell test = ", value)


####
# comparison of pro and opisthosoma

boxplot(dat$ABPV_pro, dat$ABPV_opi)
boxplot(dat$DWV_pro, dat$DWV_opi)
boxplot(dat$ABPV_pro+dat$DWV_pro, dat$DWV_opi+dat$ABPV_opi)

# ABPV
shapiro.test(c(dat$ABPV_pro, dat$ABPV_opi)) #significant --> not normally distributed
bodypart = c(rep("prosoma", nrow(dat)), rep("opistosoma", nrow(dat)))
ABPV_titer = c(dat$ABPV_pro, dat$ABPV_opi)
levene_dat <- data.frame(bodypart, ABPV_titer)
leveneTest(levene_dat$ABPV_titer, levene_dat$bodypart)# non significant --> nullhypothesis of equal variance can be assumed
test <- wilcox.test(dat$ABPV_pro, dat$ABPV_opi, paired = TRUE)
test
Z_score <-  qnorm(test$p.value/2)
Z_score

#DWV
shapiro.test(c(dat$DWV_pro, dat$DWV_opi)) #significant --> not normally distributed
DWV_titer = c(dat$DWV_pro, dat$DWV_opi)
levene_dat <- data.frame(bodypart, DWV_titer)
leveneTest(levene_dat$DWV_titer, levene_dat$bodypart)# significant --> variance not equal. but wilcox test is ok. 
test <- wilcox.test(dat$DWV_pro, dat$DWV_opi, paired = TRUE)
test
Z_score <-  qnorm(test$p.value/2)
Z_score

#### 
# correlation of pro and opisthosoma
# normally use pearson correlation if data is normally distributed else use Kendall tau or Spearman rho as non-parametric correlation coefficients
library("ggpubr")
ggscatter(dat, x = "ABPV_pro", y = "ABPV_opi", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "ABPV titer prosoma", ylab = "ABPV titer opisthosoma")

cor.test(dat$ABPV_pro, dat$ABPV_opi, method = "kendall")

ggscatter(dat, x = "DWV_pro", y = "DWV_opi", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "DWV titer prosoma", ylab = "DWV titer opisthosoma")

cor.test(dat$DWV_pro, dat$DWV_opi, method = "kendall")


# cor(x, y, method = c("pearson", "kendall", "spearman")) # calculate correlation coefficient 
# cor.test(x, y, method=c("pearson", "kendall", "spearman")) # test between paired samples returns both correlation coefficient and significance level



coco <- subset(dat, dat$sex == "female")
table(coco$treatment)
control_total <- 5 #total number of samples in control group
treatment_total <- 16 #total number of samples in treatment group
control_cocoons <- 5 #number of cases in control group
treatment_cocoons <- 4 #number of cases in treatment group

values <- c(control_cocoons,
            treatment_cocoons,
            control_total - control_cocoons,
            treatment_total - treatment_cocoons) #create the vector with the values
chi_matrix <- matrix(values,  nrow = 2) # create a matrix from the vector
chi_matrix
# alternatively create a more informative datatable
treatment = c("control", "virus treatment", "control", "virus treatment")
cocoons = c("yes", "yes", "no", "no")
N = values #c(5, 4, 0, 12)
data <- data.frame(treatment, cocoons, N)
#create crosstabulation
cocoon_tab <- xtabs(N ~ treatment + cocoons, data = data)
cocoon_tab

#all the same just for testing
#campbell2x2.test(matrix(c(xxx), nrow = 2))
campbell2x2.test(cocoon_tab)
campbell2x2.test(chi_matrix)
campbell2x2.test(matrix(values, nrow = 2))





#### proportion of samples positive for DWV in the two groups 
control_total <- 5 #total number of samples in control group
treatment_total <- 19 #total number of samples in treatment group
cases_control <- 3 #number of cases in control group
cases_treatment <- 17 #number of cases in treatment group
values <- c(cases_control,
            cases_treatment,
            control_total - cases_control,
            treatment_total - cases_treatment) #create the vector with the values
treatment = c("control", "virus treatment", "control", "virus treatment")
cases = c("yes", "yes", "no", "no")
N = values
data <- data.frame(treatment, cases, N)
tab <- xtabs(N ~ treatment + cases, data = data) #create crosstabulation
tab
test <- n1chisq.test (tab)
test$p.value/2
test$statistic
test$parameter


#### Cocoons dependant on virus titres ####
# only with females. Males have cocoon == "NA" and are thus neglected --> no need for subset

#### overall virus 

#old version
shapiro.test(dat$virus_full) #non significant --> normally distributed
leveneTest(dat$virus_full, dat$cocoon) # non significant --> nullhypothesis of equal variance can be assumed
mod <- lm(dat$virus_full ~ dat$cocoon)
#mod <- glm(dat$cocoon ~ dat$virus_full, family = "binomial")
summary(mod)
compareqqnorm(mod)
par(mfrow=c(2,2))
#leveneTest(mod) #non-significant = ok #homogenity of variance (plot(mod,1))
leveneTest(dat$virus_full, dat$cocoon)
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

# new version accounting for differences in the food source
mod <- lmer(virus_full ~ cocoon + (1|treatment), data = dat, REML = FALSE)
Anova(mod)
mod <- lm(virus_full ~ cocoon + treatment, data = dat)
summary(mod)



#ABPV full
plot(dat$ABPV_full ~ dat$cocoon)
abline(vdt_ABPV_full, 0, lty = 2)
shapiro.test(dat$ABPV_full) # significant --> not normally distributed
leveneTest(dat$ABPV_full, dat$treatment) # non significant --> nullhypothesis of equal variance can be assumed
wilcox.test(dat$ABPV_full ~ dat$cocoon)


#DWV full
mod<- lm(dat$DWV_full ~ dat$cocoon)
summary(mod)
compareqqnorm(mod)
par(mfrow=c(2,2))
#leveneTest(mod) #non-significant = ok #homogenity of variance (plot(mod,1))
leveneTest(dat$DWV_full, dat$cocoon) #non significant --> ok
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



#### models controlling for food ####
# are differences still significant if we exclude/include food as explanatory or random variable?# 


#### test just in the treated spiders (excluding the controls)
# cocoon ~ titres, controls excluded -------------------------------------------------------------------------------------------
# cocoons dependant on virus titres with controls excluded to control for differences in food
no_control <- subset(dat, dat$treatment == "Virus")
table(no_control$treatment)
boxplot(no_control$virus_full ~ no_control$cocoon)
boxplot(no_control$DWV_full ~ no_control$cocoon)
boxplot(no_control$ABPV_full ~ no_control$cocoon)

shapiro.test(no_control$virus_full) #non significant --> normally distributed
leveneTest(no_control$virus_full, no_control$cocoon) # non significant --> nullhypothesis of equal variance can be assumed
mod <- lm(no_control$virus_full ~ no_control$cocoon)
summary(mod)
mod1 <- glm(no_control$cocoon ~ no_control$virus_full, family = "binomial")
summary(mod1)
mod <- lm(no_control$DWV_full ~ no_control$cocoon)
#same trends as before but not significant due to small sample size? 


#### test again with mixed effect models that include treatment as random factor or additional explanatory variable?
#### thats the version that got into the manuscripts, replacing the older simple 
#But an alternative version to do the analysis would be using generalized models with treatment as random factor 

# full virus
mod <- lmer(virus_full ~ cocoon + (1|treatment), data = dat, REML = FALSE)
Anova(mod)
# would work. but the model fit is singular... overfitted? Probably not. Discussion with Nathalie: should be ok to use this model --> go with this solution!
boxplot(dat$virus_full ~ dat$cocoon)
compareqqnorm(mod)
par(mfrow=c(2,2))
leveneTest(dat$virus_full, dat$cocoon) #non significant --> ok
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #non significant --> assumption of normally distributed residuals is ok ##plot(mod, 2)#normality
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

# DWV
mod <- lmer(DWV_full ~ cocoon + (1|treatment), data = dat, REML = FALSE)
Anova(mod)

boxplot(dat$DWV_full ~ dat$cocoon)
compareqqnorm(mod)
leveneTest(dat$DWV_full, dat$cocoon) #non significant --> ok
shapiro.test(x = residuals(object = mod)) #non significant --> assumption of normally distributed residuals is ok ##plot(mod, 2)#normality
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))

# ABPV
mod <- lmer(ABPV_full ~ cocoon + (1|treatment), data = dat, REML = FALSE)
Anova(mod)

boxplot(dat$ABPV_full ~ dat$cocoon)
compareqqnorm(mod)
leveneTest(dat$ABPV_full, dat$cocoon) #non significant --> ok
shapiro.test(x = residuals(object = mod))
par(mfrow=c(2,2))
plot(mod)
scatter.smooth(fitted(mod),resid(mod)); abline(h=0, lty=2)  # residuals vs. fitted
title("Tukey-Anscombe Plot")
qqnorm(resid(mod), main="normal QQ-plot, residuals") 
qqline(resid(mod))  # qq of residuals
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))  # homogeneity of variance
par(mfrow=c(1,1))








# maybe it would be right to change the way of analysing the data > to this point we always looked at differences between the spiders that
# built a cocoon and spiders that did not. But we can turn it round and to see if virus titers explain cocoon yes or no?
# repeat model with treatment as a second explanatory variable to make sure what we record is not only due to differences in food
new_model <- glmer(cocoon ~ virus_full + (1|treatment), family = "binomial", data = coco)
summary(new_model)
new_model_2 <- glm(cocoon ~ virus_full + treatment, family = "binomial", data = coco)
summary(new_model_2)
# non siginificant... Why is there a significant difference in virus load between cocooners and non cocooners while virus load does not
# help to predict cocoon or not? 






#### graphs for the manuscript ####

par(mar=c(4.1, 4.1, 3.1, 1.1),
    oma = c(0,0,0,0)+0.5,
    cex.lab=1.5, 
    cex.axis=1.5, 
    cex.main= 1.5,
    mfrow=c(2,2))

label <- c("Prosoma", "Opisthosoma")
 
boxplot(dat$ABPV_full ~ dat$treatment, main ="", xlab = "Treatments", ylab = "Genomic ABPV copies [log]", las=1, ylim = c(4,11))
abline(vdt_ABPV_full, 0, lty = 2)
text(c(1,2), 10.8, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('A', side=3, line=1, at=0.15, cex = 2, font = 2)
boxplot(dat$DWV_full ~ dat$treatment, main ="", xlab = "Treatments", ylab = "Genomic DWV copies [log]", las=1, ylim = c(0,10.5))
abline(vdt_DWV_full, 0, lty = 2)
text(c(1,2), 10.2, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('B', side=3, line=1, at=0.15, cex = 2, font = 2)
boxplot(dat$ABPV_pro, dat$ABPV_opi, main = "", xlab = "Bodyparts", ylab = "Genomic ABPV copies [log]", las=1, ylim = c(0,11.2))
axis(1, at=(1:2), labels=label)
abline(vdt_ABPV, 0, lty = 2)
text(c(1,2), 10.9, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('C', side=3, line=1, at=0.15, cex = 2, font = 2)
boxplot(dat$DWV_pro, dat$DWV_opi, main = "", xlab = "Bodyparts", ylab = "Genomic DWV copies [log]", las=1, ylim = c(0,10.7))
axis(1, at=(1:2), labels=label)
abline(vdt_ABPV, 0, lty = 2)
text(c(1,2), 10.4, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('D', side=3, line=1, at=0.15, cex = 2, font = 2)

# new version

boxplot(dat$ABPV_full ~ dat$cocoon, main ="", xlab = "Cocoon", ylab = "Genomic ABPV copies [log]", las=1, ylim = c(4,11)) #
abline(vdt_ABPV_full, 0, lty = 2)
text(c(1,2), 10.8, labels = c("a", "a"), font = 2, cex = 1.5)
mtext('C', side=3, line=1, at=0.15, cex = 2, font = 2)
boxplot(dat$DWV_full ~ dat$cocoon, main ="", xlab = "Cocoon", ylab = "Genomic DWV copies [log]", las=1, ylim = c(0,10.5))
abline(vdt_DWV_full, 0, lty = 2)
text(c(1,2), 10.2, labels = c("a", "b"), font = 2, cex = 1.5)
mtext('D', side=3, line=1, at=0.15, cex = 2, font = 2)








#### Randomization test ####
# Aim: test whether the results (cocoon-spiders have lower virus loads compared are independent of food source (treatment: cricket vs. honey bee)
# H0: there is no difference in virus load between cocoon and no cocoon spiders, the only relevant aspect is which treatment the spiders were in (i.e. what type of food they received)
# In the observed data set certain spiders were fed with crickets (controls) instead of honey bees)
# The randomization will consists in randomly drawing the same number of spiders from the total and saying that in these spiders were fed with crickets (irrespective of truth), and in all others served food was honey bees (irrespective of truth)
# If the food source is the only thing that plays a role, we should observe similar difference in the virus load between spiders with and without cocoon in the randomized and the real data

# If however the observed difference between cocoon and no cocoon spiders is bigger than in the random data, it is 

# use the dataset coco which does not contain males!

#define variables
number_of_cricket_spiders <- aggregate(sample_id~treatment, FUN=length, data= dat[which(dat$sex=="female"),])
observed_diff <- mean(dat$virus_full[which(dat$treatment=="Virus"&dat$sex=="female")])-mean(dat$virus_full[which(dat$treatment=="Control")])
observed_diff <- mean(coco$virus_full[which(coco$cocoon=="no"&coco$sex=="female")])-mean(coco$virus_full[which(coco$cocoon=="yes")])

head(coco)

random_data <- NULL
#randomisation loop
for(i in 1:1000){
  rand_nb <- NULL
  rand_diff <- NULL
  subset_data <- coco
  nb_control <- 5
  spider_list <- unique(coco$sample_id)
  hypothethic_controls <- sample(spider_list, size=nb_control, replace = F)
  #create  ####now create first_source_discovered_RAND column according to hypothetical status
  rand_subset <- subset_data #first put the original data set
  rand_subset$treatment_rand <- "virus" #create a new column with the newly assigned random treatment
  rand_subset[which(rand_subset$sample_id%in%hypothethic_controls), "treatment_rand"] <- "control"
  
  #concatenate and store
  random_data  <- rbind(random_data, data.frame(RAND=i, rand_subset))
}


#get the summary statistics from the new random data for each of the data sets

mean_data <- aggregate(virus_full ~ RAND + treatment_rand, FUN = mean, data=random_data)
#mean_data <- aggregate(virus_full ~ RAND + cocoon, FUN = mean, data=random_data)
wide_means <-  reshape(data = mean_data, 
                            idvar=c("RAND"),
                            v.names = "virus_full",
                            timevar="treatment_rand",
                            direction = "wide"
)
wide_means$diff <- wide_means$virus_full.virus - wide_means$virus_full.control
hist(wide_means$diff)

# at this point we abandoned the idea of the randomisation. 