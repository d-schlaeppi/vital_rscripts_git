### ### ### ### ### ### ### ### ### ###
### ### ###  FC1 Analyses   ### ### ###
### ### ### ### ### ### ### ### ### ###
# Nathalies script #


#### prerequisites ####
library(lubridate)
library(plotrix)
library(scales)
library(car)
library(lme4)
library(Hmisc)
#setwd("~/Dropbox/SeniorLectureship_Bristol_Professorship_Fribourg/Students_postdocs/Post-Docs/Daniel Schlaeppi/FoodChoice_Experiment1") # homeoffice mac
setwd("/Users/gismo/Desktop/R/vital") # Daniel



#### data preparation ####
rm(list=ls())
dat <- read.table("Food_choice_1_data_2021_DS.txt", header = TRUE,stringsAsFactors = F)

# create empty variables (future data frames?)
dynamic_dat_same_t0_Nb      <- NULL
dynamic_dat_shifted_t0_Nb   <- NULL
dynamic_dat_same_t0_Diff      <- NULL
dynamic_dat_shifted_t0_Diff   <- NULL

# create time point vecotr (13 time points from t00, t05...t60)
Time_points_to_include <- 1:13
Time_points_to_include <- 1:5


# fill up the created 
# reminder for myself: dat[row,column]

for (i in 1:nrow(dat)){
  colony_id                          <- dat[i,"colony_id"]
  feeding_session                    <- dat[i,"feeding_session"]
  block                              <- dat[i,"block"]
  virus_position                     <- dat[i,"pos_corrected"]
  healthy_position                   <- c("left","right")[which(c("left","right")!=virus_position)]
  discovery_time_virus               <- period_to_seconds(hms(dat[i,paste("time_first_ant_",virus_position,sep="")]))
  discovery_time_healthy             <- period_to_seconds(hms(dat[i,paste("time_first_ant_",healthy_position,sep="")]))
  
  Nb_ant_virus_same_t0               <- as.numeric(dat[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(virus_position,1,1)  ,sep="_")])[Time_points_to_include]
  Nb_ant_healthy_same_t0             <- as.numeric(dat[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(healthy_position,1,1),sep="_")])[Time_points_to_include]
  
  if (discovery_time_virus<discovery_time_healthy){
    first_source_discovered          <- "virus"
    Nb_ant_virus_shifted_t0          <- Nb_ant_virus_same_t0
    Nb_ant_healthy_shifted_t0         <- as.numeric(dat[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
  }else{
    first_source_discovered          <- "healthy"
    Nb_ant_virus_shifted_t0          <- as.numeric(dat[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])[Time_points_to_include]
    Nb_ant_healthy_shifted_t0        <- Nb_ant_healthy_same_t0
  }
  
  time_since_discovery               <- period_to_seconds(hms(dat[i,c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60")]))[Time_points_to_include]                           - period_to_seconds(hms(dat[i,"time_first_ant"]))
  
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

#not used anymore:     text(x=1, y = 120, label = paste("p=",pval),3, cex = 1.3)




