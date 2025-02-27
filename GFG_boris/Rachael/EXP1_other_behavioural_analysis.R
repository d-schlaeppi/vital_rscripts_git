rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###  GFG (Group level flugus experiment) analysis ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Read me ####
# Analysis script for boris data files based on a script written by Natahlie Stroeymeyt and modified by Rachel Brown
# So far a Baustelle for a first glimpse at the data solely based on the videos manually annotated by Ellie Smith

#### prerequisites ####
rm(list=ls())
getwd()
# working directory
setwd("/Users/gismo/Documents/GitHub/vital_rscripts_git/GFG_boris/") # mac daniel
events <- read.csv(file="GFG_flugus_all.csv", header =T, sep=",") # load boris file
key <- read.csv(file="key.csv", header =T, sep=",")# load the key to un-blind the GFG data asannotations in boris were done with the blinded video names)
#events_rachael <- read.csv(file="Rachael/events.csv", header =T, sep=",")

str(events)

?read.csv()

# install.packages("pscl")
# install.packages("performance")
# install.packages("ggbeeswarm")
# install.packages("R2admb")
# install.packages("glmmADMB",
#            	repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                    	getOption("repos")),
#            	type="source")
# install.packages("glmmTMB")
# install.packages("carData")
# install.packages("usethis")
# install.packages("devtools")
# install.packages('TMB', type = 'source')
library(curl)
library(usethis)
library(devtools)
library(performance)
library(stringr)
library(lme4)
library(carData)
library(car)
library(pscl)
#library(glmmADMB)
library(ggplot2)
library(multcomp)
library(glmmTMB)
library(ggbeeswarm)
library(numDeriv)
library(e1071)




behavior_list <- c("G","A","P","F","R","B", "S", "N")
names(behavior_list) <- c("self_grooming", "allo_grooming", "poison_spraying", "aggression", "trophallaxis", "brood_tending", "alert_state", "escape_or_death")
#behavior_list <- c("T","S","N","G","A","B")
#names(behavior_list) <- c("twitch","alarm","escape_or_death","self_grooming","allo_grooming","brood_tending")

#allows expand grid function to pass a data frame
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
ant_list <- paste("ANT",c(1:6),sep="_")

##removing spaces in events
##first part removes leading spaces
#events$Replicate <- trimws(events$Replicate, which=c("left"))
events$Duration..s. <- trimws(events$Duration..s., which=c("left"))
events$Start..s. <- trimws(events$Start..s., which=c("left"))
events$Stop..s.<- trimws(events$Stop..s., which=c("left"))
#adds in _ in spaces
events <- as.data.frame(apply(events, 2, function(x)
  str_replace_all(string=x, pattern=" ", repl="_")))





##adding unique identifiers for each trial
events$Trial <- as.character(events$Trial)
events$group_size <- 1+as.numeric(sapply(strsplit(events$Trial, ",_"), getElement, 1))
events$treatment <- sapply(strsplit(events$Trial, ",_"), getElement, 2)
events$ant_ID <- events$Subject
events <- within (events, Petri_dish <- paste(Replicate,group_size,treatment,sep="_"))

# replicate time
events$time_of_year <- "WINTER"
events[which(grepl("\\/08\\/",events$Observation_date)),"time_of_year"] <- "SUMMER"

no_events <- events[ which(   events$ant_ID  == "NO_EVENTS"  	)  ,   ]
discards  <- events[ which(   events$ant_ID  == "DISCARD"  	)  ,   ]
events	<- events[ which( !  events$ant_ID   %in% c("DISCARD","NO_EVENTS")  	)  ,   ]



###get non-discarded replicate list
replicate_list <- rbind(unique(events[c("Replicate","treatment","group_size","Petri_dish","time_of_year")]),unique(no_events[c("Replicate","treatment","group_size","Petri_dish","time_of_year")]))

#ant status
events$ant_status <- "NESTMATE"
events[ which (events$ant_ID=="ANT_1")  , "ant_status"	] <- "FOCAL"
# unique(events$ant_status)



#replace behaviour names with more explicit behaviour names
events$Behavior_explicit <- names(behavior_list)[match(events$Behavior,behavior_list)]


##adding twitch modifiers into separate columns
#twitch size, repetition, proximity to brood and nestmates
events$Behavior <- as.character(events$Behavior)
events$Modifiers <- as.character(events$Modifiers)
events[which(events$Behavior == "T"),"twitch_size"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 1)
events[which(events$Behavior == "T"),"twitch_size"] <- sapply(strsplit(events[which(events$Behavior == "T"),"twitch_size"], "_"), getElement, 1)
events [ which(events$twitch_size == "Small") , "twitch_size" ] <- "small"
unique(events$twitch_size)

events[which(events$Behavior == "T"),"twitch_repetition"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 2)
events[ which(events$twitch_repetition == "single_twitch")  , "twitch_repetition" ] <- "single"
events[ which(events$twitch_repetition == "part_of_a_repetitive_set")  , "twitch_repetition" ] <- "repetitive"
unique(events$twitch_repetition)

events[which(events$Behavior == "T"),"proximity_to_nestmates"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 3)
events[ which(events$proximity_to_nestmates == "isolated_twitch/alone") , "proximity_to_nestmates" ] <- "alone"
events[ which(events$proximity_to_nestmates == "in_contact_with_nestmate") , "proximity_to_nestmates" ] <- "contact"
events[ which(events$proximity_to_nestmates == "in_contact_with_multiple_nestmates") , "proximity_to_nestmates" ] <- "multiple_contact"
unique(events$proximity_to_nestmates)

events[which(events$Behavior == "T"),"proximity_to_brood"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 4 )
unique(events$proximity_to_brood)

##adding grooming modifiers into grooming
events[which(events$Behavior == "A"),"Behavior_explicit"] <- events[which(events$Behavior == "A"),"Modifiers"]
unique(events$Behavior_explicit)

#adding column to count each event
events$N <- 1

###rename event duration column
names(events)[which(names(events)=="Duration_.s.")] <- "event_duration"

#each ant gets a bloc of all possible combinations of twitch modifiers
ant_bloc_twitches <- expand.grid (  
  twitch_size = na.omit(unique(events$twitch_size))
  ,
  twitch_repetition =  na.omit(unique (events$twitch_repetition))
  ,
  proximity_to_nestmates =  na.omit(unique (events$proximity_to_nestmates))
  ,
  proximity_to_brood =  na.omit(unique (events$proximity_to_brood))
  ,
  recording_duration = 20*60
  
)
ant_bloc_behaviour <- expand.grid (  
  behaviour = na.omit(unique(events$Behavior_explicit))
  ,
  recording_duration = 20*60
)

###create a twitch table counting the number of twitches of each type for each ant
###and a behaviour table counting the number of behaviours of each type for each ant
twitch_table	<- NULL
behaviour_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list[ 1:   replicate_list[repl,"group_size"] ])
  repl_twitches <- expand.grid.df(   data.frame( replicate_list[repl,] )
                                     ,
                                     ants
                                     ,
                                     ant_bloc_twitches
  )
  repl_behaviours <- expand.grid.df(   data.frame( replicate_list[repl,] )
                                       ,
                                       ants
                                       ,
                                       ant_bloc_behaviour
  )
  
  twitch_table	<- rbind(twitch_table , repl_twitches)
  behaviour_table <- rbind(behaviour_table , repl_behaviours)
}

###count number of twitches of each type and fill twitch_table
events_Nb <- events[which(events$Behavior=="T"),c(names(twitch_table)[which(names(twitch_table)%in%names(events))],"ant_status","Behavior", "N")]
events_Nb <- aggregate(N ~ . , FUN=sum, data= events_Nb)
twitch_table <- merge(twitch_table,events_Nb,all.x=T)
###fill in missing cases
twitch_table  [ is.na(twitch_table$N)  ,"N"] <- 0
twitch_table[ which(twitch_table$ant_ID=="ANT_1")   ,	"ant_status"] <- "FOCAL"
twitch_table[ which(twitch_table$ant_ID!="ANT_1")   ,	"ant_status"] <- "NESTMATE"
###add column for whether behaviour was performed or not
twitch_table$performs <- as.numeric(twitch_table$N>0)

###count number of behaviours of each type and fill behaviour_table
events_Nb <- events[c(names(behaviour_table)[which(names(behaviour_table)%in%names(events))],"ant_status","Behavior_explicit", "N")]
names(events_Nb)[which(names(events_Nb)=="Behavior_explicit")] <- "behaviour"
events_Nb <- aggregate(na.rm=T,na.action="na.pass",N ~ . , FUN=sum, data= events_Nb)

###do the same for event duration
events_duration <- events[c(names(behaviour_table)[which(names(behaviour_table)%in%names(events))],"Behavior_explicit", "event_duration")]
names(events_duration)[which(names(events_duration)=="Behavior_explicit")] <- "behaviour"
events_duration$event_duration <- as.numeric(events_duration$event_duration)
events_duration <- aggregate(na.rm=T,na.action="na.pass",event_duration ~ . , FUN=sum, data= events_duration)
###merge into behaviour table
behaviour_table <- merge(behaviour_table,events_Nb,all.x=T)
behaviour_table <- merge(behaviour_table,events_duration,all.x=T)
###fill in missing cases
behaviour_table  [ which(is.na(behaviour_table$N) ) ,"N"] <- 0
behaviour_table  [ which(is.na(behaviour_table$event_duration)&behaviour_table$behaviour!="twitch" ) ,"event_duration"] <- 0
behaviour_table  [ which(behaviour_table$behaviour=="twitch" ) ,"event_duration"] <- NA
###complement ant status column
behaviour_table[ which(behaviour_table$ant_ID=="ANT_1")   ,	"ant_status"] <- "FOCAL"
behaviour_table[ which(behaviour_table$ant_ID!="ANT_1")   ,	"ant_status"] <- "NESTMATE"
###add column for whether behaviour was performed or not
behaviour_table$performs <- as.numeric(behaviour_table$N>0)

######add a full ant ID column
twitch_table$full_ant_ID	<- with (twitch_table,paste(Replicate,treatment,group_size,ant_ID,sep="_"))
behaviour_table$full_ant_ID <- with (behaviour_table,paste(Replicate,treatment,group_size,ant_ID,sep="_"))

##deaths
deaths <- subset(events, Comment_start == "DEAD")
dead_ants <- with (deaths, paste(Replicate,treatment,group_size,ant_ID,sep="_"))
###remove dead ants from twitch_table and behaviour_table
twitch_table <- twitch_table [ which( ! twitch_table$full_ant_ID %in% dead_ants) ,   ]
behaviour_table <- behaviour_table [ which( ! behaviour_table$full_ant_ID %in% dead_ants) ,   ]
###correct the group size
actual_group_size <- aggregate (ant_ID ~ Replicate + treatment + group_size, function(x)length(unique(x)),data=twitch_table)
names(actual_group_size)[which(names(actual_group_size)=="ant_ID")] <- "actual_group_size"
twitch_table <- merge(twitch_table,actual_group_size,all.x=T)
behaviour_table <- merge(behaviour_table,actual_group_size,all.x=T)

###escapes
escapes <- subset(events, Comment_start == "ESCAPE")
escapes$full_ant_ID <-  with (escapes, paste(Replicate,treatment,group_size,ant_ID,sep="_"))
twitch_table[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] <- twitch_table[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] - as.numeric(escapes[match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  )[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  )))],"event_duration"])
behaviour_table[which(!is.na(match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] <- behaviour_table[which(!is.na(match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] - as.numeric(escapes[match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  )[which(!is.na(match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  )))],"event_duration"])

behaviour_table$interac_size_treat <- with (behaviour_table, interac_size_treat <- paste(group_size,treatment,sep="_"))




#remove summer focals

behaviour_table$Replicate <- as.numeric(behaviour_table$Replicate)
behaviour_table$rep_set <- "1"
behaviour_table [ which (behaviour_table$Replicate > 16 ) , "rep_set" ] <- "2"


###remove summer focals
rep_set_one <- subset(behaviour_table, rep_set=="1")
rep_set_two <- subset(behaviour_table, rep_set=="2")

summer_nestmates <- subset(rep_set_two, ant_status=="NESTMATE")

behaviour_table <- rbind(rep_set_one, summer_nestmates)





#########################################
###### GRAPHS FOR OTHER BEHAVIOURS ######
#########################################

#install.packages("gridExtra")
library(gridExtra)



### SELF GROOMING NB EVENTS ###


self_groom <- subset(behaviour_table, behaviour=="self_grooming")

self_groom <- subset(behaviour_table, behaviour=="self_grooming")
self_groom$recording_duration <- self_groom$recording_duration / 60
self_groom$frequency <- self_groom$N / self_groom$recording_duration


mean_status1 <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=standard_error, data=self_groom)

mean_status1$se <- error_status$frequency
colnames(mean_status1) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
mean_status1$Group_size <- factor(mean_status1$Group_size,
                                 levels=c("1", "2", "6", "11"))

mean_status1[which(mean_status1$Treatment=="C"), "Treatment"] <- "Control"
mean_status1[which(mean_status1$Treatment=="S"), "Treatment"] <- "Sham"
mean_status1[which(mean_status1$Treatment=="P"), "Treatment"] <- "Pathogen"


x <- ggplot(mean_status1, aes(fill=Treatment, y=Mean, x=Group_size)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~ant_status) +
  labs(x="Group size", y="Mean number self-groom", title="Mean number of self-grooming events per ant per minute")+
  #ylim(0,0.3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=15),
        legend.title=element_blank(), legend.text=element_text(size=16), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), legend.position = "bottom", plot.title =element_text(size=18, face="italic"), plot.title.position="panel" )


## just treatment

mean_status2 <- aggregate( frequency ~  treatment + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( frequency ~  treatment + ant_status, FUN=standard_error, data=self_groom)

mean_status2$se <- error_status$frequency
colnames(mean_status2) <- c("Treatment", "ant_status", "Mean", "standard_error")


mean_status2[which(mean_status2$Treatment=="C"), "Treatment"] <- "Control"
mean_status2[which(mean_status2$Treatment=="S"), "Treatment"] <- "Sham"
mean_status2[which(mean_status2$Treatment=="P"), "Treatment"] <- "Pathogen"

y <- ggplot(mean_status2, aes(fill=Treatment, y=Mean, x=Treatment)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~ant_status) +
  labs(x="Treatment", y="Mean number self-groom", title="Treatment only")+
  #ylim(0,0.3)+
  theme_bw()+
  guides(fill="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=16),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), plot.title = element_text(size=16, face="italic") )


## just group size

mean_status4 <- aggregate( frequency ~  group_size + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( frequency ~  group_size + ant_status, FUN=standard_error, data=self_groom)

mean_status4$se <- error_status$frequency
colnames(mean_status4) <- c("Group_size", "ant_status", "Mean", "standard_error")
mean_status4$Group_size <- factor(mean_status4$Group_size,
                                 levels=c("1", "2", "6", "11"))


z <- ggplot(mean_status4, aes(fill=Group_size, y=Mean, x=Group_size)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~ant_status) +
  labs(x="Group size", y="Mean number self-groom", title="Group size only")+
  #ylim(0,0.3)+
  theme_bw()+
  guides(fill="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto", colour = "black"), axis.title = element_text(size=16),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), plot.title = element_text(size=16, face="italic") )


hlay <- rbind(c(1,1,1,1),
              c(1,1,1,1),
              c(2,2,3,3),
              c(2,2,3,3))


grid.arrange(x,y,z, layout_matrix=hlay)





### SELF GROOMING DURATION ###


self_groom <- subset(behaviour_table, behaviour=="self_grooming")

self_groom <- subset(behaviour_table, behaviour=="self_grooming")
# self_groom$recording_duration <- self_groom$recording_duration / 60
self_groom$event_duration <- self_groom$event_duration / self_groom$recording_duration


mean_status6 <- aggregate( event_duration ~ group_size + treatment + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( event_duration ~ group_size + treatment + ant_status, FUN=standard_error, data=self_groom)

mean_status6$se <- error_status$event_duration
colnames(mean_status6) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
mean_status6$Group_size <- factor(mean_status6$Group_size,
                                 levels=c("1", "2", "6", "11"))

mean_status6[which(mean_status6$Treatment=="C"), "Treatment"] <- "Control"
mean_status6[which(mean_status6$Treatment=="S"), "Treatment"] <- "Sham"
mean_status6[which(mean_status6$Treatment=="P"), "Treatment"] <- "Pathogen"


a <- ggplot(mean_status6, aes(fill=Treatment, y=Mean, x=Group_size)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~ant_status) +
  labs(x="Group size", y="Mean duration (s)", title="Mean duration of self-grooming events per ant per minute (s)")+
  #ylim(0,0.3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=15),
        legend.title=element_blank(), legend.text=element_text(size=16), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), legend.position = "bottom", plot.title =element_text(size=18, face="italic"), plot.title.position="panel" )


## just treatment

mean_status7 <- aggregate( event_duration ~  treatment + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( event_duration ~  treatment + ant_status, FUN=standard_error, data=self_groom)

mean_status7$se <- error_status$event_duration
colnames(mean_status7) <- c("Treatment", "ant_status", "Mean", "standard_error")


mean_status7[which(mean_status7$Treatment=="C"), "Treatment"] <- "Control"
mean_status7[which(mean_status7$Treatment=="S"), "Treatment"] <- "Sham"
mean_status7[which(mean_status7$Treatment=="P"), "Treatment"] <- "Pathogen"

b <- ggplot(mean_status7, aes(fill=Treatment, y=Mean, x=Treatment)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~ant_status) +
  labs(x="Treatment", y="Mean duration (s)", title="Treatment only")+
  #ylim(0,0.3)+
  theme_bw()+
  guides(fill="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=16),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), plot.title = element_text(size=16, face="italic") )


## just group size

mean_status8 <- aggregate( event_duration ~  group_size + ant_status, FUN=mean, data=self_groom)
error_status <- aggregate( event_duration ~  group_size + ant_status, FUN=standard_error, data=self_groom)

mean_status8$se <- error_status$event_duration
colnames(mean_status8) <- c("Group_size", "ant_status", "Mean", "standard_error")
mean_status8$Group_size <- factor(mean_status8$Group_size,
                                 levels=c("1", "2", "6", "11"))


c <- ggplot(mean_status8, aes(fill=Group_size, y=Mean, x=Group_size)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~ant_status) +
  labs(x="Group size", y="Mean duration (s)", title="Group size only")+
  #ylim(0,0.3)+
  theme_bw()+
  guides(fill="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto", colour = "black"), axis.title = element_text(size=16),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), plot.title = element_text(size=16, face="italic") )


hlay <- rbind(c(1,1,1,1),
              c(1,1,1,1),
              c(2,2,3,3),
              c(2,2,3,3))


grid.arrange(a,b,c, layout_matrix=hlay)




### NESTMATE ALLOGROOMING PLOTS ###

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

allo <- subset(behaviour_table, behaviour=="grooming_other_nestmates" | behaviour=="grooming_focal_ant")



allo$recording_duration <- allo$recording_duration / 60
allo$frequency <- allo$N / allo$recording_duration


mean_status <- aggregate( frequency ~ group_size + treatment + ant_status + behaviour, FUN=mean, data=allo)
error_status <- aggregate( frequency ~ group_size + treatment + ant_status + behaviour, FUN=standard_error, data=allo)

mean_status$se <- error_status$frequency
colnames(mean_status) <- c("Group_size", "Treatment", "ant_status", "Behaviour", "Mean", "standard_error")
mean_status$Group_size <- factor(mean_status$Group_size,
                                 levels=c("1", "2", "6", "11"))

mean_status <- subset(mean_status, ant_status=="NESTMATE")
mean_status[which(mean_status$Treatment=="P"), "Treatment"] <- "Pathogen"
mean_status[which(mean_status$Treatment=="C"), "Treatment"] <- "Control"
mean_status[which(mean_status$Treatment=="S"), "Treatment"] <- "Sham"

mean_status$Behaviour <- factor(mean_status$Behaviour,
                                levels=c("grooming_focal_ant", "grooming_other_nestmates"))

a <- ggplot(mean_status, aes(fill=Treatment, y=Mean, x=Group_size)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~Behaviour, labeller = labeller(Behaviour = 
                                               c("grooming_other_nestmates" = "Grooming other nestmates",
                                                 "grooming_focal_ant" = "Grooming focal ant"))) +
  labs(x="Group size", y="Mean Nb events", title="Allo-grooming performed by nestmates")+
  #ylim(0,0.3)+
  theme_bw()+
  guides(fill="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title.x = element_blank(),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), plot.title=element_text(size=20) )



# legend <- get_legend(a)




mean_status <- aggregate( events_duration ~ group_size + treatment + ant_status + behaviour, FUN=mean, data=allo)
error_status <- aggregate( events_duration ~ group_size + treatment + ant_status + behaviour, FUN=standard_error, data=allo)

mean_status$se <- error_status$events_duration
colnames(mean_status) <- c("Group_size", "Treatment", "ant_status", "Behaviour", "Mean", "standard_error")
mean_status$Group_size <- factor(mean_status$Group_size,
                                 levels=c("1", "2", "6", "11"))

mean_status <- subset(mean_status, ant_status=="NESTMATE")
mean_status[which(mean_status$Treatment=="P"), "Treatment"] <- "Pathogen"
mean_status[which(mean_status$Treatment=="C"), "Treatment"] <- "Control"
mean_status[which(mean_status$Treatment=="S"), "Treatment"] <- "Sham"

mean_status$Behaviour <- factor(mean_status$Behaviour,
                                levels=c("grooming_focal_ant", "grooming_other_nestmates"))

b <- ggplot(mean_status, aes(fill=Treatment, y=Mean, x=Group_size)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5) +
  facet_wrap(~Behaviour, labeller = labeller(Behaviour = 
                                               c("grooming_other_nestmates" = "Grooming other nestmates",
                                                 "grooming_focal_ant" = "Grooming focal ant"))) +
  labs(x="Group size", y="Mean Duration (s)")+
  #ylim(0,0.3)+
  theme_bw()+
  #guides(fill="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), legend.position = "bottom" )


hlay <- rbind(c(1,1,1),
              c(1,1,1),
              c(2,2,2),
              c(2,2,2))


grid.arrange(a,b, layout_matrix=hlay)








###########################################
##### STATISTICS FOR OTHER BEHAVIOURS #####
###########################################






behaviour_table$interac_size_treat <- with (behaviour_table, interac_size_treat <- paste(group_size,treatment,sep="_"))
behaviour_table$group_size <- factor(behaviour_table$group_size)



#########################
##### SELF GROOMING #####
#########################



### FOCAL SELF GROOMING NB EVENTS ###

### STRUGGLING TO GET BELOW NORMALITY THRESHOLD ###
## TRY WITH POISSON ##



self_grooming_Nb <- glmer(log(N+1/sqrt(2)) ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate),
                         family="gaussian", data=self_groom[which(self_groom$ant_status=="FOCAL"),])


test_norm(residuals(self_grooming_Nb))


Anova(self_grooming_Nb, type="III")

qqnorm(residuals(self_grooming_Nb))
qqline(residuals(self_grooming_Nb))
hist(residuals(self_grooming_Nb))


self_grooming_Nb <- glmer(log(N+1/sqrt(2)) ~ treatment+group_size + offset(log(recording_duration)) + (1|Replicate),
                          family="gaussian", data=self_groom[which(self_groom$ant_status=="FOCAL"),])

test_norm(residuals(self_grooming_Nb))


Anova(self_grooming_Nb)


p <- summary(glht(self_grooming_Nb, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
print(cld(p))
g <- summary(glht(self_grooming_Nb, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
print(cld(g))
## these results match the graph so maybe its okay, double check




### FOCAL SELF GROOMING DURATION ###

### CANNOT GET THIS TO BE NORMALLY DISTRIBUTED, CANNOT USE POISSON AS CONTINUOUS DATE ###


self_grooming_d <- glmer(log10(event_duration+1/sqrt(2)) ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate),
                         family="gaussian", data=self_groom[which(self_groom$ant_status=="FOCAL"),])



test_norm(residuals(self_grooming_d))

## notes: could not get data to be normal (shapiro wilk test) so looked at histogram and chose poisson dist
## looked at the qqplot of a poisson and it fits well
## cannot be poisson because it isnt count data


Anova(self_grooming_d, type="III")

qqnorm(residuals(self_grooming_d))
qqline(residuals(self_grooming_d))
hist(residuals(self_grooming_d))




#non significant interaction
self_grooming_d <- glmer(log(event_duration + 1/sqrt(2)) ~ treatment+group_size + offset(log(recording_duration)) + (1|Replicate),
                         family="gaussian", data=self_groom[which(self_groom$ant_status=="FOCAL"),])

test_norm(residuals(self_grooming_d))

Anova(self_grooming_d)

p <- summary(glht(self_grooming_d, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))

p <- summary(glht(self_grooming_d, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
print(cld(p))




### NESTMATE SELFGROOMING NB EVENTS ###



self_grooming_Nb <- glmer(sqrt(N) ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                          + (1|time_of_year), family="gaussian", data=self_groom[which(self_groom$ant_status=="NESTMATE"),])


test_norm(residuals(self_grooming_Nb))

Anova(self_grooming_Nb, type="III")

##significant interaction
self_grooming_Nb <- glmer(sqrt(N) ~ treatment+group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                          + (1|time_of_year), family="gaussian", data=self_groom[which(self_groom$ant_status=="NESTMATE"),])


Anova(self_grooming_Nb)


p <- summary(glht(self_grooming_Nb, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
print(cld(p))




### NESTMATE SELFGROOMING DURATION ###


self_grooming_d <- glmer(sqrt(event_duration) ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate),
                         family="gaussian", data=self_groom[which(self_groom$ant_status=="NESTMATE"),])


test_norm(residuals(self_grooming_d))

Anova(self_grooming_d, type="III")

##significant interaction
self_grooming_d <- glmer(sqrt(event_duration) ~ interac_size_treat + offset(log(recording_duration)) + (1|Replicate),
                         family="gaussian", data=self_groom[which(self_groom$ant_status=="NESTMATE"),])



p <- summary(glht(self_grooming_d, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
print(cld(p))
#matches graphs :)





#########################
##### ALLO GROOMING #####
#########################


allo_focal <- subset(behaviour_table, behaviour=="grooming_focal_ant")


### ALLO GROOMING FOCAL ANTS NB EVENTS ###


allo_grooming_f_Nb <- glmer(sqrt(N) ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                            + (1|time_of_year), family="gaussian", data=allo_focal[which(allo_focal$ant_status=="NESTMATE"),])


test_norm(residuals(allo_grooming_f_Nb))

qqnorm(residuals(allo_grooming_f_Nb))
qqline(residuals(allo_grooming_f_Nb))
hist(residuals(allo_grooming_f_Nb))



Anova(allo_grooming_f_Nb, type="III")


##significant interaction
allo_grooming_f_Nb <- glmer(sqrt(N) ~ interac_size_treat + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                            + (1|time_of_year), family="gaussian", data=allo_focal[which(allo_focal$ant_status=="NESTMATE"),])


p <- summary(glht(allo_grooming_f_Nb, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
print(cld(p))




### ALLO GROOMING FOCAL ANTS DURATION ###

allo_grooming_d <- glmer(log(event_duration+1/sqrt(2)) ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate),
                         family="gaussian", data=allo_focal[which(allo_focal$ant_status=="NESTMATE"),])


test_norm(residuals(allo_grooming_d))

Anova(allo_grooming_d, type="III")


##significant interaction
allo_grooming_d <- glmer(log(event_duration+1/sqrt(2)) ~ interac_size_treat + offset(log(recording_duration)) + (1|Replicate),
                         family="gaussian", data=allo_focal[which(allo_focal$ant_status=="NESTMATE"),])

p <- summary(glht(allo_grooming_d, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
print(cld(p))





### ALLO GROOMING NESTMATES NB EVENTS ###
### VERY HIGH KURTOSIS

allo_nestmates <- subset(behaviour_table, behaviour=="grooming_other_nestmates")



allo_grooming_n_Nb <- glmer(N^0.1 ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                            + (1|time_of_year), family="gaussian", data=allo_nestmates[which(allo_nestmates$ant_status=="NESTMATE"),])



test_norm(residuals(allo_grooming_n_Nb))

# qqnorm(residuals(allo_grooming_n_Nb))
# qqline(residuals(allo_grooming_n_Nb))
# hist(residuals(allo_grooming_n_Nb))



Anova(allo_grooming_n_Nb, type="III")

#interaction non signifcant

allo_grooming_n_Nb <- glmer(N^0.1 ~ treatment+group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                            + (1|time_of_year), family="gaussian", data=allo_nestmates[which(allo_nestmates$ant_status=="NESTMATE"),])



test_norm(residuals(allo_grooming_n_Nb))

Anova(allo_grooming_n_Nb, type="II")

#no significance bad fitting model




### ALLO GROOMING NESTMATES DURATION ###

allo_grooming_n_Nb <- glmer(event_duration^0.1 ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                            + (1|time_of_year), family="gaussian", data=allo_nestmates[which(allo_nestmates$ant_status=="NESTMATE"),])



test_norm(residuals(allo_grooming_n_Nb))

Anova(allo_grooming_n_Nb, type="III") #non sig interaction


allo_grooming_n_Nb <- glmer(event_duration^0.1 ~ treatment+group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) 
                            + (1|time_of_year), family="gaussian", data=allo_nestmates[which(allo_nestmates$ant_status=="NESTMATE"),])


Anova(allo_grooming_n_Nb)



#no significance bad fitting model












