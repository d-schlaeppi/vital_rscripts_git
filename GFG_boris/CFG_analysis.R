rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###  GFG (Group level flugus experiment) analysis ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Read me ####
# Analysis script for boris data files based on a script written by Natahlie Stroeymeyt and modified by Rachel Brown
# So far a Baustelle for a first glimpse at the data solely based on the videos manually annotated by Ellie Smith 
# Ellie annotated videos from day 5-7 typically 20 min per video
# So far only the manually annotated videos are included in the analysis, rest will hopefully done with DLC2Action

#### prerequisites ####
# working directory and data
setwd("/Users/gismo/Documents/GitHub/vital_rscripts_git/GFG_boris/") # mac daniel
events <- read.csv(file="GFG_flugus_all.csv", header =T, sep=",") # boris file
key <- read.csv(file="key.csv", header =T, sep=",") # key to un-blind the GFG data (annotations in boris were done with the blinded video names)

# packages and functions
{
library(curl)
library(usethis)
library(devtools)
library(performance)
library(stringr)
library(lme4)
library(carData)
library(car)
library(pscl)
library(ggplot2)
library(multcomp)
library(glmmTMB)
library(ggbeeswarm)
library(numDeriv)
library(e1071)
library(dplyr)
library(gridExtra)
#library(glmmADMB)
standard_error <- function(x) sd(x,na.rm=T) / sqrt(length(x))
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...)) #allows expand grid function to pass a data frame
}

test_norm <- function(resids){ # adopted from Adriano
  print("Testing normality")
  print(paste0("Number of data points: ",length(resids)))
  print("")
  print("Fewer than 300 data points so performing Shapiro-Wilk's test")
  print(shapiro.test(resids))
  print("If p-value is less than 0.05 (significant), residuals are NOT normally distributed")
  print("")
  print("More than 300 data points so using the skewness and kurtosis approach")
  print("Skewness should be between -3 and +3 (best around zero")
  print(skewness(resids))
  print("")
  print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
  print(kurtosis(resids))
}
### prepare dataframes and creating additional variables for analyses ###
{
behavior_list <- c("G","A","P","F","R","B", "S", "N")
names(behavior_list) <- c("self_grooming", "allo_grooming", "poison_spraying", "aggression", "trophallaxis", "brood_tending", "alert_state", "escape_or_death")
ant_list <- c("BLUE","YELLOW", "WHITE", "RED", "GREEN", "PINK")
#removing spaces in events: first part removes leading spaces and then adds in _ instead of spaces
events$Duration..s. <- trimws(events$Duration..s., which=c("left"))
events$Start..s. <- trimws(events$Start..s., which=c("left"))
events$Stop..s.<- trimws(events$Stop..s., which=c("left"))
events <- as.data.frame(apply(events, 2, function(x)
  str_replace_all(string=x, pattern=" ", repl="_")))
events <- events[, -which(names(events) == "Media.file")]
# adding the information from the key df to unblind the data and changing a few variable names matching previous scripts of other users
names(events)[names(events) == "Observation.id"] <- "video_id"
names(events)[which(names(events)=="Duration..s.")] <- "event_duration"
events <- merge(events, key, by = "video_id", all.x = TRUE)
names(events)[names(events) == "rep"] <- "replicate"
events$treatment_full <- substr(events$treatment_id, nchar(events$treatment_id) - 1, nchar(events$treatment_id))
treatment_levels_flupy <- c("control", "low", "high")
events$t_flupy <- factor(events$t_flupy, levels = treatment_levels_flupy, ordered = TRUE)
treatment_levels_fungus <- c("sham", "fungus")
events$t_fungus <- factor(events$t_fungus, levels = treatment_levels_fungus, ordered = TRUE)
events$ant_ID <- events$Subject
events$trial <- as.character(events$video_id)
#removing the . in variable names and changing it to _ + making variables lower case
original_names <- names(events)
new_names <- gsub("\\.$", "", gsub("_$", "", gsub("^\\.", "_", tolower(gsub("\\.+", "_", original_names)))))
names(events) <- new_names
events <- within (events, petri_dish <- paste(replicate,treatment_full,sep="_"))
events$n <- 1 #adding column to count each event
events$ant_status <- "NESTMATE" 
events[ which (events$ant_id=="BLUE")  , "ant_status"	] <- "FOCAL"
events$behavior_explicit <- names(behavior_list)[match(events$behavior,behavior_list)] #replace behaviour names with more explicit behaviour names
events <- events[events$discard_ignore != "yes", ] #Discard and igenore all replicats and trials that were lost so far... needs to modified in the original dataframe first
trial_list <- unique(events[c("video_id", "treatment_full", "t_fungus", "t_flupy", "petri_dish")])

# adding grooming modifiers into allo grooming
events[which(events$behavior == "A"),"behavior_explicit"] <- events[which(events$behavior == "A"),"modifiers"]
unique(events$behavior_explicit)

# seconds not characters:
events <- events %>% mutate(start_s = as.numeric(start_s), stop_s = as.numeric(stop_s), event_duration = as.numeric(event_duration))

ant_bloc_behaviour <- expand.grid (
  behaviour = na.omit(unique(events$behavior_explicit)),
  recording_duration = 20*60
)

}

###  create behavior table counting the number of behaviors of each type for each ant ###
{
behaviour_table <- NULL
for (repl in 1:nrow(trial_list)){ 
  ants <- data.frame(ant_id=ant_list[ 1:6 ])
  repl_behaviours <- expand.grid.df(   data.frame( trial_list[repl,] )
                                       ,
                                       ants
                                       ,
                                       ant_bloc_behaviour
  )
  behaviour_table <- rbind(behaviour_table , repl_behaviours)
}

### count number of behaviours of each type and fill behaviour_table
events_Nb <- events[c(names(behaviour_table)[which(names(behaviour_table)%in%names(events))],"ant_status","behavior_explicit", "n")]
names(events_Nb)[which(names(events_Nb)=="behavior_explicit")] <- "behaviour"
events_Nb <- aggregate(na.rm=T,na.action="na.pass",n ~ . , FUN=sum, data= events_Nb)

### do the same for event duration
events_duration <- events[c(names(behaviour_table)[which(names(behaviour_table)%in%names(events))],"ant_status","behavior_explicit", "event_duration")]
names(events_duration)[which(names(events_duration)=="behavior_explicit")] <- "behaviour"
events_duration$event_duration <- as.numeric(events_duration$event_duration)
events_duration <- aggregate(na.rm=T,na.action="na.pass",event_duration ~ . , FUN=sum, data= events_duration)

### merge into behaviour table
behaviour_table <- merge(behaviour_table,events_Nb,all.x=T)
behaviour_table <- merge(behaviour_table,events_duration,all.x=T)

### fill in missing cases
behaviour_table  [ which(is.na(behaviour_table$n) ) ,"n"] <- 0
behaviour_table  [ which(is.na(behaviour_table$event_duration)&behaviour_table$behaviour!="twitch" ) ,"event_duration"] <- 0

### complement ant status column
behaviour_table[ which(behaviour_table$ant_id=="BLUE")   ,	"ant_status"] <- "FOCAL"
behaviour_table[ which(behaviour_table$ant_id!="BLUE")   ,	"ant_status"] <- "NESTMATE"

behaviour_table$performs <- as.numeric(behaviour_table$n>0) # add column for whether behaviour was performed or not
behaviour_table$full_ant_ID <- with (behaviour_table,paste(petri_dish,ant_id,sep="_")) #add a full ant ID column
behaviour_table$replicate <- substr(behaviour_table$petri_dish, 1, 2)
}










### ### ### ### ### ### ### ### ### ###
####             GRAPHS            ####
### ### ### ### ### ### ### ### ### ###


#### Graph: Number of self grooming events ####

self_groom <- subset(behaviour_table, behaviour=="self_grooming")
self_groom$recording_duration <- self_groom$recording_duration / 60
self_groom$frequency <- self_groom$n / self_groom$recording_duration

mean_status1 <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=standard_error, data=self_groom)

mean_status1$se <- error_status$frequency
colnames(mean_status1) <- c("Treatment_Fungus", "Treatment_Flupy", "ant_status", "Mean", "Standard_error")


x <- ggplot(mean_status1, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - Standard_error, ymax = Mean + Standard_error), 
                position = position_dodge(0.9), width = 0.2, size = 0.7) +
  facet_wrap(~ ant_status) +
  labs(x = "Fungus Treatment", y = "Mean number self-groom", 
       title = "Mean number of self-grooming events per ant per minute") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "sans"), 
        axis.title = element_text(size = 15),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "italic"),
        plot.title.position = "plot"
  )


### just flupy

mean_status2 <- aggregate( frequency ~  t_flupy + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( frequency ~  t_flupy + ant_status, FUN=standard_error, data=self_groom)

mean_status2$se <- error_status$frequency
colnames(mean_status2) <- c("Treatment_Flupy", "ant_status", "Mean", "standard_error")

y <- ggplot(mean_status2, aes(fill = Treatment_Flupy, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ ant_status) +
  labs(x = "Treatment", y = "Mean number self-groom", 
       title = "Treatment only") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "italic")
  )


### just fungus

mean_status4 <- aggregate( frequency ~  t_fungus + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( frequency ~  t_fungus + ant_status, FUN=standard_error, data=self_groom)
mean_status4$se <- error_status$frequency
colnames(mean_status4) <- c("Treatment_Fungus", "ant_status", "Mean", "standard_error")
mean_status4$Treatment_fungus <- factor(mean_status4$Treatment_fungus,
                                 levels=c("sham", "fungus"))

z <- ggplot(mean_status4, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Fungus)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ ant_status) +
  labs(x = "Treatment", y = "Mean number self-groom", 
       title = "Treatment only") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "italic")
  )


# plot all three plots toghether
hlay <- rbind(c(1,1,1,1),
              c(1,1,1,1),
              c(2,2,3,3),
              c(2,2,3,3))

grid.arrange(x,y,z, layout_matrix=hlay)


#### Stats: Number of self grooming events ####
# predictions without looking at graphs: 
# fungus should have higher number of self grooming events than sham
# with increasing flupy concentration the number of self grooming events get lower

self_grooming_Nb_focal <- glmer(n ~ t_fungus * t_flupy + (1|replicate), family="poisson", data=self_groom[which(self_groom$ant_status=="FOCAL"),])
test_norm(residuals(self_grooming_Nb_focal))
Anova(self_grooming_Nb_focal, type="III")
qqnorm(residuals(self_grooming_Nb_focal))
qqline(residuals(self_grooming_Nb_focal))
hist(residuals(self_grooming_Nb_focal))

# posthoc and compact letter display
cld(lsmeans(self_grooming_Nb_focal, pairwise~t_fungus*t_flupy, adjust="tukey"))

# consider if the full comparison is required - the focals were never exposed to the pesticide so I do not expect any interaction with the flupy treatment regarding self grooming


# following the same for nestmate
self_grooming_Nb_nestmate <- glmer(n ~ t_fungus * t_flupy + (1|replicate), family="poisson", data=self_groom[which(self_groom$ant_status=="NESTMATE"),])
Anova(self_grooming_Nb_nestmate)
test_norm(residuals(self_grooming_Nb_nestmate)) # residuals not normally distributed but skewness and kurtosis are within acceptable range... might need some modifications
Anova(self_grooming_Nb_nestmate, type="III")
qqnorm(residuals(self_grooming_Nb_nestmate))
qqline(residuals(self_grooming_Nb_nestmate))
hist(residuals(self_grooming_Nb_nestmate))

cld(lsmeans(self_grooming_Nb_nestmate, pairwise~t_fungus*t_flupy, adjust="tukey"))


### Take away ###
# sham fungus not difference regarding self grooming... 
# in focals the pattern is quite weird but as they were never exposed to pesticides it is probably not relevant and only sham fungus needs to be compared
# in nestmates self grooming is lower and there might be a trend that with increasing flupy concentration the number of self grooming events get lower
# fungus high as well as sham high have lower amounts of self grooming than fungus low and fungus control.


#### Graph: Duration of self grooming ####

self_groom <- subset(behaviour_table, behaviour=="self_grooming")
mean_status6 <- aggregate( event_duration ~ t_flupy + t_fungus + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( event_duration ~ t_flupy + t_fungus + ant_status, FUN=standard_error, data=self_groom)

mean_status6$se <- error_status$event_duration
colnames(mean_status6) <- c("Treatment_Fungus", "Treatment_Flupy", "ant_status", "Mean", "Standard_error")

a <- ggplot(mean_status6, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - Standard_error, ymax = Mean + Standard_error), 
                position = position_dodge(0.9), width = 0.2, size = 0.7) +
  facet_wrap(~ ant_status) +
  labs(x = "Fungus Treatment", y = "Mean number self-groom", 
       title = "Mean time spent self-grooming in 20 min (s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "sans"), 
        axis.title = element_text(size = 15),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "italic"),
        plot.title.position = "plot"
  )
# print(a)


### just flupy
mean_status7 <- aggregate( event_duration ~  t_flupy + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( event_duration ~  t_flupy + ant_status, FUN=standard_error, data=self_groom)
mean_status7$se <- error_status$event_duration
colnames(mean_status7) <- c("Treatment_Flupy", "ant_status", "Mean", "standard_error")

b <- ggplot(mean_status7, aes(fill = Treatment_Flupy, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ ant_status) +
  labs(x = "Treatment", y = "Duration of self grooming (s)", 
       title = "Treatment only") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "italic")
  )
# print(b)



### just fungus
mean_status8 <- aggregate( event_duration ~  t_fungus + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( event_duration ~  t_fungus + ant_status, FUN=standard_error, data=self_groom)
mean_status8$se <- error_status$event_duration
colnames(mean_status8) <- c("Treatment_Fungus", "ant_status", "Mean", "standard_error")
mean_status8$Treatment_Fungus <- factor(mean_status8$Treatment_Fungus,
                                        levels=c("sham", "fungus"))

c <- ggplot(mean_status8, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Fungus)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ ant_status) +
  labs(x = "Treatment", y = "Duration (s)", 
       title = "Fungus only") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "italic")
  )
#print(c)

hlay <- rbind(c(1,1,1,1),
              c(1,1,1,1),
              c(2,2,3,3),
              c(2,2,3,3))


grid.arrange(a,b,c, layout_matrix=hlay)



#### Stats: Duration of self grooming events ####
# predictions without looking at graphs: 
# fungus should have higher duration of self grooming events than sham
# with increasing flupy concentration the duration of self grooming events is expected to be lower

self_grooming_duration_focal <- glmer(round(event_duration) ~ t_fungus * t_flupy + (1|replicate), family="poisson", data=self_groom[which(self_groom$ant_status=="FOCAL"),])
test_norm(residuals(self_grooming_duration_focal))
Anova(self_grooming_duration_focal, type="III")
qqnorm(residuals(self_grooming_duration_focal))
qqline(residuals(self_grooming_duration_focal))
hist(residuals(self_grooming_duration_focal))

# posthoc and compact letter display
cld(lsmeans(self_grooming_duration_focal, pairwise~t_fungus*t_flupy, adjust="tukey"))
# weird pattern: controls fungus and sham the same but opposing trends in combination with flupy: reduced selfgrooming in sham but upregulated selfgrooming in fungus treatment
# consider if the full comparison is required - the focals were never exposed to the pesticide so I do not expect any interaction with the flupy treatment regarding self grooming
# when ignoring the pesticide treatment for the focals: sefl grooming is upregulated when exposed to fungus

# Nestmates
self_grooming_duration_nestmate <- glmer(round(event_duration) ~ t_fungus * t_flupy + (1|replicate), family = "poisson", data=self_groom[which(self_groom$ant_status=="NESTMATE"),])
Anova(self_grooming_duration_nestmate) # Anova(self_grooming_duration_nestmate, type="III")
test_norm(residuals(self_grooming_duration_nestmate))
qqnorm(residuals(self_grooming_duration_nestmate))
qqline(residuals(self_grooming_duration_nestmate))
hist(residuals(self_grooming_duration_nestmate))
cld(lsmeans(self_grooming_duration_nestmate, pairwise~t_fungus*t_flupy, adjust="tukey"))

zero_poisson_model <- glmmTMB(round(event_duration) ~ t_fungus * t_flupy + (1|replicate), ziformula=~1, family = poisson, data = self_groom)
summary(zero_poisson_model)
test_norm(residuals(zero_poisson_model)) # residuals not normally distributed & skewness and kurtosis are not within acceptable range
Anova(zero_poisson_model, type="III")
qqnorm(residuals(zero_poisson_model))
qqline(residuals(zero_poisson_model))
hist(residuals(zero_poisson_model))
cld(lsmeans(zero_poisson_model, pairwise~t_fungus*t_flupy, adjust="tukey"))


### try out lukes box cox transformation: 
library(MASS)
b<-boxcox(lm(self_groom$event_duration+1 ~ 1))
lambda <- b$x[which.max(b$y)]
self_groom$boxcox_duration<-((self_groom$event_duration+1)^lambda-1)/lambda
self_grooming_duration_cox <- lmer(boxcox_duration ~ t_fungus * t_flupy + (1|replicate),  data=self_groom[which(self_groom$ant_status=="NESTMATE"),])
#did not really work either... :-(

# not yet sure about the stats  - needs more time for the nest mate model... might need to be a zero inflated model but I am not sure about it... 
# either way the pattern is not as clear..



#### Graph: Allogrooming by nestmates ####

allo <- subset(behaviour_table, ant_status == "NESTMATE" & 
                 (behaviour == "grooming_other_nestmates" | behaviour == "grooming_focal_ant"))
allo$recording_duration <- allo$recording_duration / 60
allo$frequency <- allo$n / allo$recording_duration


mean_status_a <- aggregate( frequency ~ t_flupy + t_fungus + ant_status + behaviour, FUN=mean, data=allo)
error_status <- aggregate( frequency ~ t_flupy + t_fungus + ant_status + behaviour, FUN=standard_error, data=allo)
mean_status_a$se <- error_status$frequency
colnames(mean_status_a) <- c("Treatment_Flupy", "Treatment_Fungus", "ant_status", "Behaviour", "Mean", "standard_error")
mean_status_a$Behaviour <- factor(mean_status_a$Behaviour,
                                levels=c("grooming_focal_ant", "grooming_other_nestmates"))


a <- ggplot(mean_status_a, aes(fill = Treatment_Flupy, y = Mean, x = Treatment_Fungus)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ Behaviour, labeller = labeller(Behaviour = 
                                                c("grooming_other_nestmates" = "Grooming other nestmates",
                                                  "grooming_focal_ant" = "Grooming focal ant"))) +
  labs(x = "Fungus Treatment", y = "Mean number events", 
       title = "Allo-grooming performed by nestmates") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"), 
        axis.title.x = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 20)
  )
# print(a)

mean_status_b <- aggregate( event_duration ~ t_flupy + t_fungus + ant_status + behaviour, FUN=mean, data=allo)
error_status <- aggregate( event_duration ~ t_flupy + t_fungus + ant_status + behaviour, FUN=standard_error, data=allo)
mean_status_b$se <- error_status$event_duration
colnames(mean_status_b) <- c("Treatment_Flupy", "Treatment_Fungus", "ant_status", "Behaviour", "Mean", "standard_error")
mean_status_b$Behaviour <- factor(mean_status_b$Behaviour,
                                levels=c("grooming_focal_ant", "grooming_other_nestmates"))

b <- ggplot(mean_status_b, aes(fill = Treatment_Flupy, y = Mean, x = Treatment_Fungus)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ Behaviour, labeller = labeller(Behaviour = 
                                                c("grooming_other_nestmates" = "Grooming other nestmates",
                                                  "grooming_focal_ant" = "Grooming focal ant"))) +
  labs(x = "Fungus Treatment", y = "Mean duration (s)", 
       title = "") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"), 
        axis.title.x = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 20)
  )
# print(b)
hlay <- rbind(c(1,1,1),
              c(1,1,1),
              c(2,2,2),
              c(2,2,2))
grid.arrange(a,b, layout_matrix=hlay)




#### Stats: Number of allogrooming events performed by nestmates ####
# predictions without looking at graphs: 
# fungus should have higher number of allo grooming events than sham
# with increasing flupy concentration the number of allo grooming events get lower

head(allo)
hist(allo$n)

allo_grooming_Nb_focal <- glmer(n ~ t_fungus * t_flupy + (1|replicate), family="poisson", data=allo[which(allo$behaviour=="grooming_focal_ant"),])
test_norm(residuals(allo_grooming_Nb_focal))
Anova(allo_grooming_Nb_focal, type="III")
qqnorm(residuals(allo_grooming_Nb_focal))
qqline(residuals(allo_grooming_Nb_focal))
hist(residuals(allo_grooming_Nb_focal))

hist(allo$n)

# posthoc and compact letter display
cld(lsmeans(allo_grooming_Nb_focal, pairwise~t_fungus*t_flupy, adjust="tukey"))
zero_poisson_model <- glmmTMB(n ~ t_fungus * t_flupy + (1|replicate), ziformula=~1, family = poisson, data = allo[which(allo$behaviour=="grooming_focal_ant"),])
summary(zero_poisson_model)
test_norm(residuals(zero_poisson_model)) # residuals not normally distributed & skewness and kurtosis are not within acceptable range
Anova(zero_poisson_model, type="III")
qqnorm(residuals(zero_poisson_model))
qqline(residuals(zero_poisson_model))
hist(residuals(zero_poisson_model))
cld(lsmeans(zero_poisson_model, pairwise~t_fungus*t_flupy, adjust="tukey"))


# consider if the full comparison is required - the focals were never exposed to the pesticide so I do not expect any interaction with the flupy treatment regarding self grooming


# following the same for nestmate
self_grooming_Nb_nestmate <- glmer(n ~ t_fungus * t_flupy + (1|replicate), family="poisson", data=self_groom[which(self_groom$ant_status=="NESTMATE"),])
Anova(self_grooming_Nb_nestmate)
test_norm(residuals(self_grooming_Nb_nestmate)) # residuals not normally distributed but skewness and kurtosis are within acceptable range... might need some modifications
Anova(self_grooming_Nb_nestmate, type="III")
qqnorm(residuals(self_grooming_Nb_nestmate))
qqline(residuals(self_grooming_Nb_nestmate))
hist(residuals(self_grooming_Nb_nestmate))

cld(lsmeans(self_grooming_Nb_nestmate, pairwise~t_fungus*t_flupy, adjust="tukey"))









#### Graph: Poison spraying ####

poison <- subset(behaviour_table, behaviour=="poison_spraying")
poison$recording_duration <- poison$recording_duration / 60
poison$frequency <- poison$n / poison$recording_duration

mean_status_p <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=mean, data=poison)
error_status <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=standard_error, data=poison)
mean_status_p$se <- error_status$frequency
colnames(mean_status_p) <- c("Treatment_Fungus", "Treatment_Flupy", "ant_status", "Mean", "Standard_error")

p <- ggplot(mean_status_p, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - Standard_error, ymax = Mean + Standard_error), 
                position = position_dodge(0.9), width = 0.2, size = 0.7) +
  facet_wrap(~ ant_status) +
  labs(x = "Fungus Treatment", y = "Mean number poison spraying", 
       title = "Mean number of poison spraying events per ant per minute") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "sans"), 
        axis.title = element_text(size = 15),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "italic"),
        plot.title.position = "plot"
  )


# just flupy
mean_status_p2 <- aggregate( frequency ~  t_flupy + ant_status, FUN=mean, data=poison)
error_status <- aggregate( frequency ~  t_flupy + ant_status, FUN=standard_error, data=poison)
mean_status_p2$se <- error_status$frequency
colnames(mean_status_p2) <- c("Treatment_Flupy", "ant_status", "Mean", "standard_error")

p2 <- ggplot(mean_status_p2, aes(fill = Treatment_Flupy, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ ant_status) +
  labs(x = "Treatment", y = "Mean number poison spraying", 
       title = "Flupy treatment only") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "italic")
  )

# fungus only
mean_status_p3 <- aggregate( frequency ~  t_fungus + ant_status, FUN=mean, data=poison)
error_status <- aggregate( frequency ~  t_fungus + ant_status, FUN=standard_error, data=poison)
mean_status_p3$se <- error_status$frequency
colnames(mean_status_p3) <- c("Treatment_Fungus", "ant_status", "Mean", "standard_error")

p3 <- ggplot(mean_status_p3, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Fungus)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Mean - standard_error, ymax = Mean + standard_error), 
                width = 0.6, position = position_dodge(0.9), size = 1) +
  facet_wrap(~ ant_status) +
  labs(x = "Treatment", y = "Mean number poison spraying", 
       title = "Fungus treatment only") +
  theme_bw() +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "italic")
  )


# plot all three plots toghether
hlay <- rbind(c(1,1,1,1),
              c(1,1,1,1),
              c(2,2,3,3),
              c(2,2,3,3))

grid.arrange(p,p2,p3, layout_matrix=hlay)

#### Graph: Aggression ####

aggro <- subset(behaviour_table, behaviour=="aggression")
aggro$recording_duration <- aggro$recording_duration / 60
aggro$frequency <- aggro$n / aggro$recording_duration

mean_status_a1 <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=mean, data=aggro)
error_status <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=standard_error, data=aggro)
mean_status_a1$se <- error_status$frequency
colnames(mean_status_a1) <- c("Treatment_Fungus", "Treatment_Flupy", "ant_status", "Mean", "Standard_error")

a1 <- ggplot(mean_status_a1, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - Standard_error, ymax = Mean + Standard_error), 
                position = position_dodge(0.9), width = 0.2, size = 0.7) +
  facet_wrap(~ ant_status) +
  labs(x = "Fungus Treatment", y = "Mean number aggression", 
       title = "Mean number of aggression events per ant per minute") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "sans"), 
        axis.title = element_text(size = 15),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "italic"),
        plot.title.position = "plot"
  )

print(a1)

# seems to be only interesting for nestmates... there might be something with nest mates beeing agressive towards pathogen! but only in controls and low but not in high! 
# this is interesting explore further. 

#### Graph: Trophallaxis ####

trophy <- subset(behaviour_table, behaviour=="trophallaxis")
trophy$recording_duration <- trophy$recording_duration / 60
trophy$frequency <- trophy$n / trophy$recording_duration

mean_status_t1 <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=mean, data=trophy)
error_status <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=standard_error, data=trophy)
mean_status_t1$se <- error_status$frequency
colnames(mean_status_t1) <- c("Treatment_Fungus", "Treatment_Flupy", "ant_status", "Mean", "Standard_error")

t1 <- ggplot(mean_status_t1, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - Standard_error, ymax = Mean + Standard_error), 
                position = position_dodge(0.9), width = 0.2, size = 0.7) +
  facet_wrap(~ ant_status) +
  labs(x = "Fungus Treatment", y = "Mean number trophallaxis events", 
       title = "Mean number of trophallaxis events per ant per minute") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "sans"), 
        axis.title = element_text(size = 15),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "italic"),
        plot.title.position = "plot"
  )

print(t1)

# there might be something  - a slight trend that fungus exposure reduces trophallaxis among nestmates but probably no influence of pesticide... 

#### brood_tending ####

brood <- subset(behaviour_table, behaviour=="brood_tending")
brood$recording_duration <- brood$recording_duration / 60
brood$frequency <- brood$n / brood$recording_duration

mean_status_b1 <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=mean, data=brood)
error_status <- aggregate( frequency ~ t_flupy + t_fungus + ant_status, FUN=standard_error, data=brood)
mean_status_b1$se <- error_status$frequency
colnames(mean_status_b1) <- c("Treatment_Fungus", "Treatment_Flupy", "ant_status", "Mean", "Standard_error")

b1 <- ggplot(mean_status_b1, aes(fill = Treatment_Fungus, y = Mean, x = Treatment_Flupy)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - Standard_error, ymax = Mean + Standard_error), 
                position = position_dodge(0.9), width = 0.2, size = 0.7) +
  facet_wrap(~ ant_status) +
  labs(x = "Fungus Treatment", y = "Mean number brood tending events", 
       title = "Mean number of brood tending events per ant per minute") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "sans"), 
        axis.title = element_text(size = 15),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "italic"),
        plot.title.position = "plot"
  )

print(b1)