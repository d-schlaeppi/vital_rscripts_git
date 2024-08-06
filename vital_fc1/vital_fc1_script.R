# rm(list = ls())
rm(list = setdiff(ls(), "first_time_use_working_directory"))

### ### ### ### ### ### ### ### ### ### ### ###
### ### ###  VITAL - FC1 Analyses   ### ### ###
### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####
#' This script was written by Nathalie Stroexmeyt and Daniel Schläppi
#' It contains the analysis of the binary food choice experient conducted prior to the vital tracking experiment
#' Of 16 colonies 2 feeding events were recorded (one week apart) at each the ants were presented with a virus or a control food source
#' the side of the where the virus was placed has been chosen pseudo randomly to avoid learning effects from 1th to 2nd feeding
#' The first data table contains the information on the number of ants present at each of the food sources at any given time 
#' the second table contains annotations of ant feeding durations. For the first feeding event we covered 1h until after the second food source has been discovered
#' for the second feeding event only the first 5 ants feeding on each food source were recorded.


#### prerequisites ####

# libraries
library(pacman)
pacman::p_load(lubridate, plotrix, scales, car, lme4, Hmisc, dplyr, tidyverse, blmeco, lmtest, lsmeans, emmeans, multcompView, multcomp, viridis, crayon)

#functions
sem <- function(x) {sd(x,na.rm=T)/sqrt(lenght(na.omit(x)))} # function to get standard error of means


# choose directory 
if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") {
  library(tcltk)
  setwd(tk_choose.dir(default = "~/", caption = "Select Working Directory"))
  first_time_use_working_directory <- getwd()
  setwd(first_time_use_working_directory)
  cat(blue(getwd()))
} else {setwd(first_time_use_working_directory)
  cat(blue(getwd()))}

# load data
dat_numbers <- read.csv("vital_fc1_data_numbers.csv", header = TRUE)
dat_duration <- read.csv("vital_fc1_data_feedingdurations.csv", header = TRUE)


#### Analysis of the number of ant present on each of the food sources at a given time since discovery ####

#### data preparation ####

# in the loaded data block corresponds to the eight days of video recordings -> needs to be renamed to video_session
# a new variable "block" gets created grouping the 4x4 colonies that were recorded as a group (twice)
# 1-4 instead of 1-8
dat_numbers <- rename(dat_numbers, video_session = block)
block_corrected <- NULL
for (i in 1:nrow(dat_numbers)) {
  block_old           <- dat_numbers[i, "video_session"]
  block_new       <- if(block_old<5) {block_old} else {block_old-4}
  block_corrected <- rbind(block_corrected, data.frame(block_new))
}
dat_numbers$block <- block_corrected$block_new



# create empty variables
dynamic_dat_same_t0_Nb      <- NULL
dynamic_dat_shifted_t0_Nb   <- NULL
dynamic_dat_same_t0_Diff      <- NULL
dynamic_dat_shifted_t0_Diff   <- NULL

# create time point vecotr (13 time points from t00, t05...t60)
# Time_points_to_include <- 1:13 # full hour 
Time_points_to_include <- 1:5  # first 20 minutes (was significant in the past)


# fill up the created data frames

for (i in 1:nrow(dat_numbers)){ # i <- 1
  colony_id                          <- dat_numbers[i,"colony_id"]
  feeding_session                    <- dat_numbers[i,"feeding_session"]
  block                              <- dat_numbers[i,"block"]
  virus_position                     <- dat_numbers[i,"position_virus_corrected"]
  healthy_position                   <- c("left","right")[which(c("left","right")!=virus_position)]
  discovery_time_virus               <- period_to_seconds(hms(dat_numbers[i,paste("time_first_ant_",virus_position,sep="")]))
  discovery_time_healthy             <- period_to_seconds(hms(dat_numbers[i,paste("time_first_ant_",healthy_position,sep="")]))
  
  Nb_ant_virus_same_t0               <- as.numeric(dat_numbers[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(virus_position,1,1)  ,sep="_")])[Time_points_to_include]
  Nb_ant_healthy_same_t0             <- as.numeric(dat_numbers[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(healthy_position,1,1),sep="_")])[Time_points_to_include]
  
  if (discovery_time_virus<discovery_time_healthy){
    first_source_discovered          <- "virus"
    Nb_ant_virus_shifted_t0          <- Nb_ant_virus_same_t0
    Nb_ant_healthy_shifted_t0         <- as.numeric(dat_numbers[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
  }else{
    first_source_discovered          <- "healthy"
    Nb_ant_virus_shifted_t0          <- as.numeric(dat_numbers[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
    Nb_ant_healthy_shifted_t0        <- Nb_ant_healthy_same_t0
  }
  
  time_since_discovery               <- period_to_seconds(hms(dat_numbers[i,c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60")]))[Time_points_to_include]  - period_to_seconds(hms(dat_numbers[i,"time_first_ant"]))
  
  dynamic_dat_same_t0_Nb                <- rbind(dynamic_dat_same_t0_Nb, data.frame (colony_id,
                                                                                     feeding_session,
                                                                                     block,
                                                                                     virus_position,
                                                                                     discovery_time_virus,
                                                                                     discovery_time_healthy,
                                                                                     first_source_discovered,
                                                                                     time=rep(time_since_discovery,2),
                                                                                     status = rep(c("virus","healthy"),each=length(time_since_discovery)),
                                                                                     Nb_ants=c(Nb_ant_virus_same_t0,Nb_ant_healthy_same_t0),
                                                                                     stringsAsFactors = F))
  
  dynamic_dat_shifted_t0_Nb              <- rbind(dynamic_dat_shifted_t0_Nb, data.frame (colony_id,
                                                                                         feeding_session,
                                                                                         block,
                                                                                         virus_position,
                                                                                         discovery_time_virus,
                                                                                         discovery_time_healthy,
                                                                                         first_source_discovered,
                                                                                         time=rep(time_since_discovery,2),
                                                                                         status = rep(c("virus","healthy"),each=length(time_since_discovery)),
                                                                                         Nb_ants=c(Nb_ant_virus_shifted_t0,Nb_ant_healthy_shifted_t0),
                                                                                         stringsAsFactors = F))
  
  dynamic_dat_same_t0_Diff                <- rbind(dynamic_dat_same_t0_Diff, data.frame (colony_id,
                                                                                         feeding_session,
                                                                                         block,
                                                                                         virus_position,
                                                                                         discovery_time_virus,
                                                                                         discovery_time_healthy,
                                                                                         first_source_discovered,
                                                                                         time=time_since_discovery,
                                                                                         status = "virus_minus_healthy",
                                                                                         Delta_Nb_ants=Nb_ant_virus_same_t0-Nb_ant_healthy_same_t0,
                                                                                         stringsAsFactors = F))
  
  dynamic_dat_shifted_t0_Diff              <- rbind(dynamic_dat_shifted_t0_Diff, data.frame (colony_id,
                                                                                             feeding_session,
                                                                                             block,
                                                                                             virus_position,
                                                                                             discovery_time_virus,
                                                                                             discovery_time_healthy,
                                                                                             first_source_discovered,
                                                                                             time=time_since_discovery,
                                                                                             status = "virus_minus_healthy",
                                                                                             Delta_Nb_ants=Nb_ant_virus_shifted_t0-Nb_ant_healthy_shifted_t0,
                                                                                             stringsAsFactors = F))
  
}


#### the analysis ####

for (time_origin in c("same","shifted")){
  ###get dynamic data
  dynamic_dat_Nb      <- get(paste("dynamic_dat",time_origin,"t0","Nb",sep="_")) 
  dynamic_dat_Diff    <- get(paste("dynamic_dat",time_origin,"t0","Diff",sep="_")) 
  
  ###calculate means and SE create dynamic_summary_dat_Nb objects
  dynamic_summary_dat_Nb    <- data.frame(as.matrix(aggregate(Nb_ants   ~feeding_session+status+time,function(x)cbind(mean(x),std.error(x)),data=dynamic_dat_Nb)),stringsAsFactors = F)
  names(dynamic_summary_dat_Nb)[which(grepl(".1",names(dynamic_summary_dat_Nb)))] <- gsub(".1","_Mean",names(dynamic_summary_dat_Nb)[which(grepl(".1",names(dynamic_summary_dat_Nb)))])
  names(dynamic_summary_dat_Nb)[which(grepl(".2",names(dynamic_summary_dat_Nb)))] <- gsub(".2","_SE",names(dynamic_summary_dat_Nb)[which(grepl(".2",names(dynamic_summary_dat_Nb)))])
  dynamic_summary_dat_Nb$Nb_ants_Mean <- as.numeric( dynamic_summary_dat_Nb$Nb_ants_Mean)
  dynamic_summary_dat_Nb$Nb_ants_SE   <- as.numeric( dynamic_summary_dat_Nb$Nb_ants_SE)
  dynamic_summary_dat_Nb <- dynamic_summary_dat_Nb[order(dynamic_summary_dat_Nb$time),]
  
  ###calculate means and SE create dynamic_summary_dat_Diff objects
  dynamic_summary_dat_Diff    <- data.frame(as.matrix(aggregate(Delta_Nb_ants   ~feeding_session+status+time,function(x)cbind(mean(x),std.error(x)),data=dynamic_dat_Diff)),stringsAsFactors = F)
  names(dynamic_summary_dat_Diff)[which(grepl(".1",names(dynamic_summary_dat_Diff)))] <- gsub(".1","_Mean",names(dynamic_summary_dat_Diff)[which(grepl(".1",names(dynamic_summary_dat_Diff)))])
  names(dynamic_summary_dat_Diff)[which(grepl(".2",names(dynamic_summary_dat_Diff)))] <- gsub(".2","_SE",names(dynamic_summary_dat_Diff)[which(grepl(".2",names(dynamic_summary_dat_Diff)))])
  dynamic_summary_dat_Diff$Delta_Nb_ants_Mean <- as.numeric( dynamic_summary_dat_Diff$Delta_Nb_ants_Mean)
  dynamic_summary_dat_Diff$Delta_Nb_ants_SE   <- as.numeric( dynamic_summary_dat_Diff$Delta_Nb_ants_SE)
  dynamic_summary_dat_Diff <- dynamic_summary_dat_Diff[order(dynamic_summary_dat_Diff$time),]
  
  ### calculate one mean number of ants for each source, each colony and each feeding session
  summary_dat_Nb      <- aggregate(Nb_ants   ~colony_id+feeding_session+block+virus_position+discovery_time_virus+discovery_time_healthy+first_source_discovered+status,FUN=mean,data=dynamic_dat_Nb)
  ### stats
  print( paste("Statistics -",time_origin,"T0"))
  model <- lmer(Nb_ants ~ status*feeding_session + (1|colony_id) + (1|first_source_discovered), data=summary_dat_Nb )
  print(Anova(model))
  print(shapiro.test(residuals(model)))
  ### plot  TO IMPROVE: SHOULD BE MEANS AND STANDARD ERRORS RATHR THAN BOXPLOT
  boxplot(Nb_ants~status+feeding_session,data=summary_dat_Nb)
  
  ### NOTE: AS FEEDING SESSION IS HIGHLY NON-SIGNIFICANT, WE DON'T NEED TO DISTINGUISH BETWEEN THE TWO IN THE ANALYSIS OF RANDOMISATION TESTS
  
  ###Plot Dynamic data For Each Session
  for (session in c(1,2)){
    ymin <- min( dynamic_summary_dat_Nb$Nb_ants_Mean-dynamic_summary_dat_Nb$Nb_ants_SE )
    ymax <- max( dynamic_summary_dat_Nb$Nb_ants_Mean+dynamic_summary_dat_Nb$Nb_ants_SE )
    
    plot(Nb_ants_Mean ~ time, pch=16, col="red",
         data=dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),],
         ylim=c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin)),bty="l", 
         main =paste("Session",session,"-", time_origin ,"t0"))
    polygon( x= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"time"]
                   ,
                   rev(dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"time"]) ) , 
             y= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_Mean"] - dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_SE"]  
                   ,
                   rev( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_Mean"] + dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="virus"),"Nb_ants_SE"])),
             col=alpha("red",0.5),border=NA
    )
    points(Nb_ants_Mean ~ time, pch=16, col="blue",
           data=dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),])
    polygon(  x= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"time"]
                    ,
                    rev(dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"time"]) ) , 
              y= c( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_Mean"] - dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_SE"]  
                    ,
                    rev( dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_Mean"] + dynamic_summary_dat_Nb[which(dynamic_summary_dat_Nb$feeding_session==session & dynamic_summary_dat_Nb$status=="healthy"),"Nb_ants_SE"])),
              col=alpha("blue",0.5),border=NA
    )
    legend ("topright",legend=c("virus","healthy"),col=c("red","blue"),pch=16,bty="n")
    
  }
  
  ### Randomisation test for each session
  ### Aim: test whether the results are an artifact from virus food tending to be discovered earlier than healthy food
  ### H0: there is no difference in ant behaviour towards either virus or healthy food, the only relevant point for recruitment is which food was discovered first
  ### In the observed dataset the virus was discovered first a certain number of time
  ### The randomisation will consist in randomly drawing the same number of colonies from the total and saying that in these colonies the first source had virus (irrespective of truth), and in all others the first source was healthy irrespective of truth)
  ### If observed data is an artifact, we should observe a similar difference in number of ants between the two sources in randomised and real data
  
  ### Pre-calculations: nb of colonies that found the virus first in each session, and observed difference in number of ants across all colonies and all times
  nb_discoveries <- aggregate(colony_id~first_source_discovered+feeding_session ,FUN=length, data= summary_dat_Nb[which(summary_dat_Nb$status=="healthy"),])
  observed_Diff    <- aggregate(Delta_Nb_ants   ~ status,FUN=mean,data=dynamic_dat_Diff)
  
  random_data_dynamic_Nb <- NULL
  random_data_Diff       <- NULL 
  
  for (i in 1 : 1000){ ### randomisation loop
    rand_Nb   <- NULL
    rand_Diff <- NULL
    ### randomise separately for each session
    for (session in c(1,2)){
      
      nb_virus       <- nb_discoveries[which(nb_discoveries$first_source_discovered=="virus"&nb_discoveries$feeding_session==session),"colony_id"]
      
      subset_Nb   <- dynamic_dat_Nb[which(dynamic_dat_Nb$feeding_session==session),]
      subset_Diff <- dynamic_dat_Diff[which(dynamic_dat_Diff$feeding_session==session),]
      colony_list <- unique(subset_Nb$colony_id)
      
      hypothetic_virus_first    <- sample(colony_list,size=nb_virus,replace = F)
      ####now create first_source_discovered_RAND column according to hypothetical status
      rand_subset_Nb                                                                                           <- subset_Nb             ; rand_subset_Diff                                                                                             <- subset_Diff  
      rand_subset_Nb$first_source_discovered_RAND                                                              <- "healthy"             ; rand_subset_Diff$first_source_discovered_RAND                                                                <- "healthy" ;
      rand_subset_Nb[which(rand_subset_Nb$colony_id%in%hypothetic_virus_first),"first_source_discovered_RAND"] <- "virus"               ; rand_subset_Diff[which(rand_subset_Diff$colony_id%in%hypothetic_virus_first),"first_source_discovered_RAND"] <- "virus"
      
      #### for Nb: copy status column into a new column status_RAND, and switch the values for those colonies who have been assigned a different first source discovered than in the observed data
      rand_subset_Nb$status_RAND  <- rand_subset_Nb$status                                                                           <- rand_subset_Nb$status 
      rand_subset_Nb[which(rand_subset_Nb$first_source_discovered_RAND!=rand_subset_Nb$first_source_discovered &  rand_subset_Nb$status=="virus"),"status_RAND"] <- "healthy"
      rand_subset_Nb[which(rand_subset_Nb$first_source_discovered_RAND!=rand_subset_Nb$first_source_discovered &  rand_subset_Nb$status=="healthy"),"status_RAND"] <- "virus"
      
      #### for Diff: multiply the rows of colonies who have been assigned a different first source discovered than in the observed data by minus 1
      rand_subset_Diff[which(rand_subset_Diff$first_source_discovered_RAND!=rand_subset_Diff$first_source_discovered),"Delta_Nb_ants"] <- (-1)*rand_subset_Diff[which(rand_subset_Diff$first_source_discovered_RAND!=rand_subset_Diff$first_source_discovered),"Delta_Nb_ants"] 
      
      #### concatenate and store
      rand_Nb   <- rbind(rand_Nb  ,data.frame(RAND=i, feeding_session=session, rand_subset_Nb))
      rand_Diff <- rbind(rand_Diff,data.frame(RAND=i, feeding_session=session, rand_subset_Diff))
      
       
    }
    
    ###then we use aggregate to get the mean number of ants at each time point in the hypothetical scenario, across all colonies and both feeding events
    mean_rand_dat_Nb <- aggregate(Nb_ants ~ status_RAND + time, FUN=mean, data=rand_Nb)
    #### and we use aggregate to calculate the mean delta nb ant across all times,  all colonies and both feeding events in the hypothetical scenario
    mean_rand_dat_Diff <- aggregate(Delta_Nb_ants ~  1, FUN=mean, data=rand_Diff)
    ###and we concatenate
    random_data_dynamic_Nb <- rbind(random_data_dynamic_Nb,data.frame(RAND=i, feeding_session=session, mean_rand_dat_Nb))
    random_data_Diff       <- rbind(random_data_Diff      ,data.frame(RAND=i, feeding_session=session, mean_rand_dat_Diff))
    
    
  }
  
  #### Plot expected vs observed, Delta ants
  xmin <- min (c(random_data_Diff$Delta_Nb_ants),observed_Diff[,"Delta_Nb_ants"])
  xmax <- max (c(random_data_Diff$Delta_Nb_ants),observed_Diff[,"Delta_Nb_ants"])
  
  hist(random_data_Diff$Delta_Nb_ants,col=alpha("grey",0.5),xlim=c(xmin,xmax),main=paste("Observed vs Expected,",capitalize(time_origin),"T0"), xlab=expression(paste(Delta, " number of ants")))
  arrows(x0=observed_Diff[,"Delta_Nb_ants"],
         y0=0,
         y1=-10,
         code=1,
         col="black",
         lwd=4,
         length=0.1)
  
  #### P value
  if (observed_Diff[,"Delta_Nb_ants"]>median(random_data_Diff$Delta_Nb_ants)){
    pval <- 2*length(which(random_data_Diff$Delta_Nb_ants>=observed_Diff[,"Delta_Nb_ants"]))/length(random_data_Diff$Delta_Nb_ants)
    mtext(paste("p=",pval),3, cex = 1.3)
  }else{
    pval <- 2*length(which(random_data_Diff$Delta_Nb_ants<=observed_Diff[,"Delta_Nb_ants"]))/length(random_data_Diff$Delta_Nb_ants)
    mtext(paste("p=",pval),3, cex = 1.3)
    
  }
  
}




#### update here with the parts usefull from reconstruciton and learning ####

# not used anymore:     text(x=1, y = 120, label = paste("p=",pval),3, cex = 1.3)



#### Analysis of feeding duration ####

# get duration in seconds: 
dat_duration <- dat_duration %>%
  mutate(feeding_duration_seconds = as.numeric(as.difftime(feeding_duration, format="%H:%M:%S", units="secs")))


# get correct position of virus and data from other numbers data frame
# create new variable food source saying for each of the feeding events whether it was no a virus of control food source
for (i in 1:nrow(dat_duration)) { # i <- 1
  index <- which(dat_numbers$colony_id == dat_duration$colony_id[i] &
                   dat_numbers$feeding_session == dat_duration$feeding_session[i]) 
  dat_duration$position_virus[i] <- dat_numbers$position_virus_corrected[index]
  dat_duration$treatment[i]      <- dat_numbers$treatment[index]
  if (dat_duration$feeding_side[i] == "r") { # get the side at which the feeding event was observed
    mapped_feeding_side <- "right"
  } else if (dat_duration$feeding_side[i] == "l") {
    mapped_feeding_side <- "left"
  } else {
    mapped_feeding_side <- NA  
  }
  if (mapped_feeding_side == dat_numbers$position_virus_corrected[index]) { #define whether feeding event occured on a virus or control food source
    dat_duration$food_source[i] <- "virus"
  } else {
    dat_duration$food_source[i] <- "control"
  }
}


# mean duration of feeding  events based on the two food sources just for the first 5 events per food source
# subset based on event_id (for each event a number saying this was feeding nr X on this food source for left or right)
dat_duration <- dat_duration %>% # extract the number from event_id for new variable
  mutate(event_number = as.numeric(gsub("^[^0-9]*(\\d+).*", "\\1", event_id))) %>% 
  filter(!is.na(event_number))  
dat_duration_sub <- subset(dat_duration, event_number <= 5 )

dat_duration_sub %>%
  group_by(food_source) %>%
  summarize(
    mean_feeding_duration = mean(feeding_duration_seconds, na.rm = TRUE),
    sd_feeding_duration = sd(feeding_duration_seconds, na.rm = TRUE)
  )

# It was worth investigating but I already had a feeling when watching the videos that feeding duration does not necessarily reflect feeding volume very well... 
# Some ants are full after 30 sec others take 5 min to get full, might depend more on other variables e.g. like consistency of food and such...

ggplot(dat_duration_sub, aes(x = food_source, y = feeding_duration_seconds, fill = food_source)) +
  geom_boxplot() +
  labs(
    title = "Feeding Duration by Food Source (first 5 events only)",
    x = "Food Source",
    y = "Feeding Duration (seconds)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("virus" = "blue", "control" = "green"))

lmm <- lmer(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block) + (1|feeding_session), data = dat_duration_sub)
summary(lmm)
Anova(lmm)
# might need another model, but for now it does the job --> consider changing to glmer with poisson







#### Analysis of completely annotated first feeding session ####
dat_duration_first <- subset(dat_duration, feeding_session == 1)

ggplot(dat_duration_first, aes(x = food_source, y = feeding_duration_seconds)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, color = "black", alpha = 0.1) +
  labs(title = "Feeding Duration by Food Source",
       x = "Food Source",
       y = "Feeding Duration (seconds)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


mod <- glmer(feeding_duration_seconds ~ food_source + (1|colony_id) + (1|block), data = dat_duration_first, family = "poisson")
summary(mod)
Anova(mod)
### does it for now, but check if model assumptions need to be met for this thing and if not consider alternative versions 
mod <- lmer(sqrt(feeding_duration_seconds) ~ food_source + (1|colony_id) + (1|block), data = dat_duration_first)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow = c(1,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) # not significant would be good 



# boxcox transformation: 
prelim_model <- lm(feeding_duration_seconds ~ food_source, data = dat_duration_first)

bc <- boxcox(prelim_model, lambda = seq(-2, 2, by = 0.1)) # determine the optimal lambda to transform based on max log likelihood
lambda <- bc$x[which.max(bc$y)]
if (lambda == 0) {
  dat_duration_first$transformed_duration <- log(dat_duration_first$feeding_duration_seconds + 1)
} else {
  dat_duration_first$transformed_duration <- (dat_duration_first$feeding_duration_seconds^lambda - 1) / lambda
}

#full model
mod <- lmer(transformed_duration ~ food_source + (1 | colony_id) + (1 | block), data = dat_duration_first)
# Summary of the model
summary(mod)
# Anova to check the significance of fixed effects
Anova(mod)
# Check residuals for normality
par(mfrow = c(1, 1))
aov_residuals <- residuals(object = mod)
qqnorm(aov_residuals)
qqline(aov_residuals, col = "red")
shapiro.test(aov_residuals)

# still not ok... check Nathalies or Adrianos scipts for what they do because it is very easy for norm distribution not to be ok if the dataset is very large


### colony level summary  
dat_summary <- dat_duration_first %>%
  group_by(food_source, colony_id) %>%
  summarise(
    num_feeding_events = n(),
    total_feeding_duration = sum(feeding_duration_seconds, na.rm = TRUE)
  )
dat_summary <-  as.data.frame(dat_summary)

# number of feeding events per colony per food source 
mod <- lm(num_feeding_events ~ food_source, data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod); par(mfrow = c(1,1))
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #not significant good 

ggplot(dat_summary, aes(x = food_source, y = num_feeding_events)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "Number of Feeding Events by Food Source",
       x = "Food Source",
       y = "Number of Feeding Events") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# summed up feeding duration per colony and food source 
mod <- lm(log10(total_feeding_duration) ~ food_source, data = dat_summary)
summary(mod)
Anova(mod)
compareqqnorm(mod)
aov_residuals <- residuals(object = mod)
shapiro.test(x = aov_residuals) #not significant good

ggplot(dat_summary, aes(x = food_source, y = total_feeding_duration)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "Number of Feeding Events by Food Source",
       x = "Food Source",
       y = "Number of Feeding Events") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




# then plot a time curve showing the exploitation of the two food sources
# is first discovery the driving factor? maybe needs another randomisation as above... 



