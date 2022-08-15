rm(list=ls())


library(lubridate)
library(plotrix)
library(scales)
library(car)
# setwd("~/Dropbox/SeniorLectureship_Bristol_Professorship_Fribourg/Students_postdocs/Post-Docs/Daniel Schlaeppi/FoodChoice_Experiment1") # Nathalie
setwd("/Users/gismo/Desktop/R/vital") # Daniel
dat <- read.table("food_choice_1_data_2021_DS.txt", header = TRUE,stringsAsFactors = F)

summary_dat              <- NULL
dynamic_dat_same_t0      <- NULL
dynamic_dat_shifted_t0   <- NULL

for (i in 1:nrow(dat)){
  colony_id                          <- dat[i,"colony_id"]
  feeding_session                    <- dat[i,"feeding_session"]
  block                              <- dat[i,"block"]
  virus_position                     <- dat[i,"pos_corrected"]
  healthy_position                   <- c("left","right")[which(c("left","right")!=virus_position)]
  discovery_time_virus               <- period_to_seconds(hms(dat[i,paste("time_first_ant_",virus_position,sep="")]))
  discovery_time_healthy             <- period_to_seconds(hms(dat[i,paste("time_first_ant_",healthy_position,sep="")]))
  
  Nb_ant_virus_same_t0               <- as.numeric(dat[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(virus_position,1,1)  ,sep="_")])
  Nb_ant_healthy_same_t0             <- as.numeric(dat[i,paste(c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60"),substr(healthy_position,1,1),sep="_")])

  if (discovery_time_virus<discovery_time_healthy){
    first_source_discovered          <- "virus"
    Nb_ant_virus_shifted_t0          <- Nb_ant_virus_same_t0
    Nb_ant_healthy_shifted_t0         <- as.numeric(dat[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])
  }else{
    first_source_discovered          <- "healthy"
    Nb_ant_virus_shifted_t0          <- as.numeric(dat[i,c("rl00","rl05","rl10","rl15","rl20","rl25","rl30","rl35","rl40","rl45","rl50","rl55","rl60")])
    Nb_ant_healthy_shifted_t0        <- Nb_ant_healthy_same_t0
  }
  
  time_since_discovery               <- period_to_seconds(hms(dat[i,c("t00","t05","t10","t15","t20","t25","t30","t35","t40","t45","t50","t55","t60")]))        - period_to_seconds(hms(dat[i,"time_first_ant"]))
  
  dynamic_dat_same_t0                <- rbind(dynamic_dat_same_t0, data.frame (colony_id,
                                                                       feeding_session,
                                                                       block,
                                                                       virus_position,
                                                                       discovery_time_virus,
                                                                       discovery_time_healthy,
                                                                       first_source_discovered,
                                                                       time=rep(time_since_discovery, each = 2),
                                                                       status = rep(c("virus","healthy"),each=length(time)),
                                                                       Nb_ant_same_t0=c(Nb_ant_virus_same_t0,Nb_ant_healthy_same_t0),
                                                                       stringsAsFactors = F))
  
  dynamic_dat_shifted_t0              <- rbind(dynamic_dat_shifted_t0, data.frame (colony_id,
                                                                               feeding_session,
                                                                               block,
                                                                               virus_position,
                                                                               discovery_time_virus,
                                                                               discovery_time_healthy,
                                                                               first_source_discovered,
                                                                               time=rep(time_since_discovery, each = 2),
                                                                               status = rep(c("virus","healthy"),each=length(time)),
                                                                               Nb_ant_shifted_t0=c(Nb_ant_virus_shifted_t0,Nb_ant_healthy_shifted_t0),
                                                                               stringsAsFactors = F))
  
}




###calculate means and create summary_dat objects
summary_dat_same_t0    <- aggregate(Nb_ant_same_t0   ~colony_id+feeding_session+block+virus_position+discovery_time_virus+discovery_time_healthy+first_source_discovered+status,FUN=mean,data=dynamic_dat_same_t0)
summary_dat_shifted_t0 <- aggregate(Nb_ant_shifted_t0~colony_id+feeding_session+block+virus_position+discovery_time_virus+discovery_time_healthy+first_source_discovered+status,FUN=mean,data=dynamic_dat_shifted_t0)

###calculate means and SE create dynamic_summary_dat objects
dynamic_summary_dat_same_t0    <- data.frame(as.matrix(aggregate(Nb_ant_same_t0   ~feeding_session+status+time,function(x)cbind(mean(x),std.error(x)),data=dynamic_dat_same_t0)),stringsAsFactors = F)
names(dynamic_summary_dat_same_t0)[which(grepl(".1",names(dynamic_summary_dat_same_t0)))] <- gsub(".1","_Mean",names(dynamic_summary_dat_same_t0)[which(grepl(".1",names(dynamic_summary_dat_same_t0)))])
names(dynamic_summary_dat_same_t0)[which(grepl(".2",names(dynamic_summary_dat_same_t0)))] <- gsub(".2","_SE",names(dynamic_summary_dat_same_t0)[which(grepl(".2",names(dynamic_summary_dat_same_t0)))])

dynamic_summary_dat_shifted_t0    <- data.frame(as.matrix(aggregate(Nb_ant_shifted_t0   ~feeding_session+status+time,function(x)cbind(mean(x),std.error(x)),data=dynamic_dat_shifted_t0)),stringsAsFactors = F)
names(dynamic_summary_dat_shifted_t0)[which(grepl(".1",names(dynamic_summary_dat_shifted_t0)))] <- gsub(".1","_Mean",names(dynamic_summary_dat_shifted_t0)[which(grepl(".1",names(dynamic_summary_dat_shifted_t0)))])
names(dynamic_summary_dat_shifted_t0)[which(grepl(".2",names(dynamic_summary_dat_shifted_t0)))] <- gsub(".2","_SE",names(dynamic_summary_dat_shifted_t0)[which(grepl(".2",names(dynamic_summary_dat_shifted_t0)))])

###Plot - Mean number of ants
boxplot(Nb_ant_same_t0~status+feeding_session,data=summary_dat_same_t0)
boxplot(Nb_ant_shifted_t0~status+feeding_session,data=summary_dat_shifted_t0)

###Plot - Dynamic data - Feeding session 1
dynamic_summary_dat_same_t0$Nb_ant_same_t0_Mean <- as.numeric( dynamic_summary_dat_same_t0$Nb_ant_same_t0_Mean)
dynamic_summary_dat_same_t0$Nb_ant_same_t0_SE   <- as.numeric( dynamic_summary_dat_same_t0$Nb_ant_same_t0_SE)
dynamic_summary_dat_same_t0 <- dynamic_summary_dat_same_t0[order(dynamic_summary_dat_same_t0$time),]

ymin <- min( dynamic_summary_dat_same_t0$Nb_ant_same_t0_Mean-dynamic_summary_dat_same_t0$Nb_ant_same_t0_SE )
ymax <- max( dynamic_summary_dat_same_t0$Nb_ant_same_t0_Mean+dynamic_summary_dat_same_t0$Nb_ant_same_t0_SE )


virus_first_session    <- dynamic_summary_dat_same_t0[which(dynamic_summary_dat_same_t0$status=="virus"  &dynamic_summary_dat_same_t0$feeding_session==1,),]
healthy_first_session  <- dynamic_summary_dat_same_t0[which(dynamic_summary_dat_same_t0$status=="healthy"&dynamic_summary_dat_same_t0$feeding_session==1,),]
virus_second_session   <- dynamic_summary_dat_same_t0[which(dynamic_summary_dat_same_t0$status=="virus"  &dynamic_summary_dat_same_t0$feeding_session==2,),]
healthy_second_session <- dynamic_summary_dat_same_t0[which(dynamic_summary_dat_same_t0$status=="healthy"&dynamic_summary_dat_same_t0$feeding_session==2,),]

plot(Nb_ant_same_t0_Mean ~ time, pch=16, col="red",
     data=virus_first_session,
     ylim=c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin)),bty="l", main ="First session - same t0")
polygon( x= c( virus_first_session$time,rev(virus_first_session$time) ) , 
         y= c( virus_first_session$Nb_ant_same_t0_Mean - virus_first_session$Nb_ant_same_t0_SE  , rev( virus_first_session$Nb_ant_same_t0_Mean + virus_first_session$Nb_ant_same_t0_SE)),
         col=alpha("red",0.5),border=NA
         )
points(Nb_ant_same_t0_Mean ~ time, pch=16, col="blue",
       data=healthy_first_session)
polygon( x= c( healthy_first_session$time,rev(healthy_first_session$time) ) , 
         y= c( healthy_first_session$Nb_ant_same_t0_Mean - healthy_first_session$Nb_ant_same_t0_SE  , rev( healthy_first_session$Nb_ant_same_t0_Mean + healthy_first_session$Nb_ant_same_t0_SE)),
         col=alpha("blue",0.5),border=NA
)
legend ("topright",legend=c("virus","healthy"),col=c("red","blue"),pch=16,bty="n")


plot(Nb_ant_same_t0_Mean ~ time, pch=16, col="red",
     data=virus_second_session,
     ylim=c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin)),bty="l", main ="Second session - same t0")
polygon( x= c( virus_second_session$time,rev(virus_second_session$time) ) , 
         y= c( virus_second_session$Nb_ant_same_t0_Mean - virus_second_session$Nb_ant_same_t0_SE  , rev( virus_second_session$Nb_ant_same_t0_Mean + virus_second_session$Nb_ant_same_t0_SE)),
         col=alpha("red",0.5),border=NA
)
points(Nb_ant_same_t0_Mean ~ time, pch=16, col="blue",
       data=healthy_second_session)
polygon( x= c( healthy_second_session$time,rev(healthy_second_session$time) ) , 
         y= c( healthy_second_session$Nb_ant_same_t0_Mean - healthy_second_session$Nb_ant_same_t0_SE  , rev( healthy_second_session$Nb_ant_same_t0_Mean + healthy_second_session$Nb_ant_same_t0_SE)),
         col=alpha("blue",0.5),border=NA
)
legend ("topright",legend=c("virus","healthy"),col=c("red","blue"),pch=16,bty="n")


###Plot - Dynamic data - Feeding session 1 - Shifted t0
dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_Mean <- as.numeric( dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_Mean)
dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_SE   <- as.numeric( dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_SE)
dynamic_summary_dat_shifted_t0 <- dynamic_summary_dat_shifted_t0[order(dynamic_summary_dat_shifted_t0$time),]

ymin <- min( dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_Mean-dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_SE )
ymax <- max( dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_Mean+dynamic_summary_dat_shifted_t0$Nb_ant_shifted_t0_SE )


virus_first_session    <- dynamic_summary_dat_shifted_t0[which(dynamic_summary_dat_shifted_t0$status=="virus"  &dynamic_summary_dat_shifted_t0$feeding_session==1,),]
healthy_first_session  <- dynamic_summary_dat_shifted_t0[which(dynamic_summary_dat_shifted_t0$status=="healthy"&dynamic_summary_dat_shifted_t0$feeding_session==1,),]
virus_second_session   <- dynamic_summary_dat_shifted_t0[which(dynamic_summary_dat_shifted_t0$status=="virus"  &dynamic_summary_dat_shifted_t0$feeding_session==2,),]
healthy_second_session <- dynamic_summary_dat_shifted_t0[which(dynamic_summary_dat_shifted_t0$status=="healthy"&dynamic_summary_dat_shifted_t0$feeding_session==2,),]

plot(Nb_ant_shifted_t0_Mean ~ time, pch=16, col="red",
     data=virus_first_session,
     ylim=c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin)),bty="l", main ="First session - shifted t0")
polygon( x= c( virus_first_session$time,rev(virus_first_session$time) ) , 
         y= c( virus_first_session$Nb_ant_shifted_t0_Mean - virus_first_session$Nb_ant_shifted_t0_SE  , rev( virus_first_session$Nb_ant_shifted_t0_Mean + virus_first_session$Nb_ant_shifted_t0_SE)),
         col=alpha("red",0.5),border=NA
)
points(Nb_ant_shifted_t0_Mean ~ time, pch=16, col="blue",
       data=healthy_first_session)
polygon( x= c( healthy_first_session$time,rev(healthy_first_session$time) ) , 
         y= c( healthy_first_session$Nb_ant_shifted_t0_Mean - healthy_first_session$Nb_ant_shifted_t0_SE  , rev( healthy_first_session$Nb_ant_shifted_t0_Mean + healthy_first_session$Nb_ant_shifted_t0_SE)),
         col=alpha("blue",0.5),border=NA
)
legend ("topright",legend=c("virus","healthy"),col=c("red","blue"),pch=16,bty="n")


plot(Nb_ant_shifted_t0_Mean ~ time, pch=16, col="red",
     data=virus_second_session,
     ylim=c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin)),bty="l", main ="Second session - shifted t0")
polygon( x= c( virus_second_session$time,rev(virus_second_session$time) ) , 
         y= c( virus_second_session$Nb_ant_shifted_t0_Mean - virus_second_session$Nb_ant_shifted_t0_SE  , rev( virus_second_session$Nb_ant_shifted_t0_Mean + virus_second_session$Nb_ant_shifted_t0_SE)),
         col=alpha("red",0.5),border=NA
)
points(Nb_ant_shifted_t0_Mean ~ time, pch=16, col="blue",
       data=healthy_second_session)
polygon( x= c( healthy_second_session$time,rev(healthy_second_session$time) ) , 
         y= c( healthy_second_session$Nb_ant_shifted_t0_Mean - healthy_second_session$Nb_ant_shifted_t0_SE  , rev( healthy_second_session$Nb_ant_shifted_t0_Mean + healthy_second_session$Nb_ant_shifted_t0_SE)),
         col=alpha("blue",0.5),border=NA
)
legend ("topright",legend=c("virus","healthy"),col=c("red","blue"),pch=16,bty="n")


###statistics 
model <- glmer(    Nb_ant_same_t0 ~ status*feeding_session + (1| block/colony_id)  , family=poisson(link = "log"), data=dynamic_dat_same_t0)
Anova(model)
