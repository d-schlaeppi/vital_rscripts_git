### ### ### ### ### ### ### ### ### ###
### ANV - Ants, Neonics and Viruses ###
### ### ### ### ### ### ### ### ### ###

#### prerequisites ####

#read tables
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV")
# Neonic and virus table 
ANV <- read.table("20190226_ANV.txt", header = TRUE)

### ev include some of the count data from anv full into the firs table? 
# table with full colony developement
dat <- read.table("ANVlong_full.txt", header = TRUE)


#load libraries
library(lme4)
library(plyr)
library(survival) # contains kalan meier plot function
install.packages("ggpubr")
library("ggpubr")

#Citing R Studio
RStudio.Version()




#### Modifying the table for Analyses ####
head(ANV)
str(ANV)

ANV$treatment_1 <- relevel(ANV$treatment_1, ref ="low")
ANV$treatment_1 <- relevel(ANV$treatment_1, ref ="control")
ANV$treatment<- as.factor(ANV$treatment)
ANV$TREATMENT_LARS <- as.factor(ANV$TREATMENT_LARS)
ANV$death_date <- as.Date(ANV$death_date, "%d.%m.%Y")


#### Creating Subsets and subtables ####

survivors <- subset(ANV, id!= "C9C" & id!="C12C" & id!="C15V" & id!= "C16V" & id!= "L3C" & id!="L10V" & id!="L14C" & id!= "H4C" & id!="H11V" & id!="H18V" & id!="H19C")
analyzed <- subset(ANV, neonic_analyses == "y" & virus_analyses == "y")
analyzed_w <- subset(analyzed, sample=="worker")
analyzed_q <- subset(analyzed, sample=="queen")

analyzed_for_neonics <- subset(ANV, neonic_analyses == "y")
analyzed_for_virus <- subset(ANV, virus_analyses =="y")
analyzed_for_neonics_w <- subset(ANV, neonic_analyses == "y" & sample == "worker")
analyzed_for_virus_w <- subset(ANV, virus_analyses =="y" & sample == "worker")
analyzed_for_neonics_q <- subset(ANV, neonic_analyses == "y" & sample == "queen")
analyzed_for_virus_q <- subset(ANV, virus_analyses =="y"  & sample == "queen")

View(data) #look at dataframe in separate window


analyzed_w <- subset(ANV, sample == "worker")
analyzed_q <- subset(ANV, sample == "queen")

viruspositive <- subset(analyzed, virus_qualitative == "y")
viruspositive_w <- subset(viruspositive, sample == "worker")
viruspositive_q <- subset(viruspositive, sample == "queen")


#### First look at data ####
table(analyzed_for_virus_w$treatment)
table(analyzed_for_virus_q$treatment)
table(analyzed_for_neonics_w$treatment)
table(analyzed_for_neonics_q$treatment)
table(analyzed_w$treatment)
table(analyzed_q$treatment)

table(viruspositive_w$treatment)
table(viruspositive_q$treatment)

View(viruspositive_w)
boxplot(viruspositive_w$virus_titer ~ viruspositive_w$treatment)
boxplot(viruspositive_q$virus_titer ~ viruspositive_q$treatment)
viruspositive_q_reduced <- subset(viruspositive_q, id != "C19C")
boxplot(viruspositive_q_reduced$virus_titer ~ viruspositive_q_reduced$treatment)

plot(viruspositive_w$virus_titer ~ viruspositive_w$conc_T_dry)
plot(viruspositive_w$virus_titer ~ viruspositive_w$conc_C_dry)
ggscatter(viruspositive_w, x = "virus_titer", y = "conc_T_dry", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "virustiter", ylab = "thiamethoxam concentration")

ggscatter(viruspositive_w, x = "virus_titer", y = "conc_C_dry", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "virustiter", ylab = "clothianidin concentration")

plot(viruspositive_q$virus_titer ~ viruspositive_q$conc_T_dry)
ggscatter(viruspositive_q, x = "virus_titer", y = "conc_T_dry", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "virustiter", ylab = "thiamethoxam concentration")

plot(viruspositive_q$virus_titer ~ viruspositive_q$conc_C_dry)
plot(viruspositive_q_reduced$virus_titer ~ viruspositive_q_reduced$conc_T_dry)
ggscatter(viruspositive_q_reduced, x = "virus_titer", y = "conc_T_dry", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "virustiter", ylab = "thiamethoxam concentration")
plot(viruspositive_q_reduced$virus_titer ~ viruspositive_q_reduced$conc_C_dry)


#maybe a correlation over all samples (workers and queen?)
analyzed_reduced <- subset(analyzed, id != "C19C")
ggscatter(analyzed, x = "conc_T_dry", y = "virus_titer", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "thiamethoxam concentration", ylab = "virustiter", 
          title = "thiamethoxam")
ggscatter(analyzed, x = "conc_C_dry", y = "virus_titer", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "clothianidin concentration", ylab = "virustiter", 
          title = "clothianidin")
ggscatter(analyzed_reduced, x = "conc_T_dry", y = "virus_titer", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "thiamethoxam concentration", ylab = "virustiter", 
          title = "thiamethoxam")

ggscatter(analyzed_reduced, x = "conc_C_dry", y = "virus_titer", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "clothianidin concentration", ylab = "virustiter", 
          title = "clothianidin")

ggscatter(analyzed, x = "conc_T_dry", y = "conc_C_dry", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "thiamethoxam concentration", ylab = "clothianidin concentration")

ggscatter(analyzed, x = "ctd_q", y = "ctd_w", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "thiamethoxam concentration in queens", ylab = "thiam. concentration in workers")





#number of adults
boxplot(analyzed$nr_adults ~ analyzed$treatment, main = "final colony size", xlab = "treatments", ylab = "# workers")

ggscatter(analyzed_w, x = "conc_T_dry", y = "nr_adults", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "thiamethoxam concentration in workers", ylab = "final number of adults", title = "thiamethoxam in workers and final colony size")
ggscatter(analyzed_q, x = "conc_T_dry", y = "nr_adults", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "thiamethoxam concentration in queens", ylab = "final number of adults", title = "thiamethoxam in queens and final colony size")

ggscatter(analyzed_w, x = "conc_T_dry", y = "nr_adults", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "clothianidin concentration in workers", ylab = "final number of adults", title = "clothianidin in workers and final colony size")
ggscatter(analyzed_q, x = "conc_T_dry", y = "nr_adults", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "clothianidin concentration in queens", ylab = "final number of adults", title = "clothianidin in queens and final colony size")

thiam_workers
thiam_queens



boxplot(viruspositive_w$virus_titer~viruspositive_w$treatment, main = "workers", xlab = "treatments", ylab = "virus titer")

boxplot(viruspositive_q_reduced$virus_titer~viruspositive_q_reduced$treatment, main = "queen", xlab = "treatments", ylab = "virus titer")


#### Cumulative survival of queens in % ####
#kaplan meier oder so in %
survival <- read.table("ANV_queen-survival.txt", header = TRUE)
levels(survival$treatment_1)
survival$treatment_1 <- relevel(survival$treatment_1, ref ="low")
survival$treatment_1 <- relevel(survival$treatment_1, ref ="control")

survival$treatment <- factor(survival$treatment, levels(survival$treatment)[c(1,2,5,6,3,4)])
levels(survival$treatment)


#Data ANV only
ANV <- subset(survival, treatment_2=="control")
table(survival$treatment)


#mean survival time accross the treatments
tapply(survival$deathday_2[survival$sy2==1], survival$treatment[survival$sy2==1], mean)
plot(survfit(Surv(survival$deathday_2, survival$sy2)~1)) #plots first basic plot with 95 confidence intervalls for all the data
fit=survfit(Surv(survival$deathday_2, survival$sy2)~survival$treatment)
plot(fit) 

#plot cumulative queen survival
plot(fit, col=c(1:6), xlab="days since queens were placed in chambers", ylab="cumulative queen survival (%) ", lty = c(1:6), lwd = 2, yaxt="n", ylim=c(0.4, 1), cex.axis = 1.5, cex.lab = 1.5)
axis(2, at = c(0, 0.20, 0.40, 0.60, 0.80, 1), labels=c("0","20","40","60","80","100") , cex.axis = 1.5)
legend(12, 0.71, legend = c("0.0 µg thiamethoxam/L; no virus", "0.0 µg thiamethoxam/L; virus", "4.5 µg thiamethoxam/L; no virus", "4.5 µg thiamethoxam/L; virus", "30 µg thiamethoxam/L; no virus", "30 µg thiamethoxam/L; virus"),  col=(1:6), lty= c(1:6) , lwd=2, cex = 1.5, bty = "n")

#calculate potential statistical difference between the groups
mod <- survdiff(Surv(deathday_2, sy2) ~ treatment, data=survival)
summary(mod)
mod



