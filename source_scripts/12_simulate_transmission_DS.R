
### ### ### ### ### ### ### ### ### ### ### ###
#### 12_simulate_transmission.R #####
### ### ### ### ### ### ### ### ### ### ### ###

### Sources C++ function simulate_transmission from source folder

### Takes an interaction list and simulates the transmission of that agent from a list of originally contaminated workers to the rest of the colony

### Created by Nathalie Stroeymeyt
### Adaptations by Adriano Wanderlingh & Linda Sartoris
### Adjusted to the needs of Daniel Schläppi

#### Data preparation ####

cat(red(paste("### ### ### Starting transmission simulations for ", interaction_type, "### ### ###")), "\n")

to_keep_ori <- to_keep

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(code_path,"/simulate_transmission.cpp",sep=""))
loop_start_time <- Sys.time()

### ### ### ### ### ### ### ### ### ### ###
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

#### get input file list
input_path         <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists",sep="")
setwd(input_path)  
# input_files        <- paste(input_path,"/",list.files(recursive=T,pattern="interactions"),sep="")  
input_files <- list.files(path = input_path, pattern = "interactions", recursive = TRUE, full.names = TRUE)



### arguments 
N_SIM  <- 500 # N_SIM  <- 2  
#if (!grepl("survival",data_path)){ # DS never has survival
  seed_files <- c("treated_workers.txt","random_workers.txt","foragers.txt","nurses.txt")
# }else{
#   seed_files <- c("treated_workers.txt")
# }

to_keep <- c(ls(),"to_keep","seed_file","outputfolder","interac_folders","interac_folder","interac_list","summary_collective","summary_individual","interac", "loop_start_time")


for (seed_file in seed_files ){ # seed_file <- "treated_workers.txt"
  print(paste("### ### Seeds =",gsub("\\.txt","",seed_file), "### ###"))
  
  if (seed_file=="treated_workers.txt"){
    #if (!grepl("survival",data_path)){
      outputfolder    <- paste(data_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="")
      interac_folders <- "observed"
    # }else{
    #   outputfolder    <- paste(data_path,"/transmission_simulations/post_treatment/experimentally_exposed_seeds",sep="")
    #   interac_folders <- "observed"
    # }
  }else{
    outputfolder    <- paste(data_path,"/transmission_simulations/random_vs_observed/",gsub("\\.txt","",seed_file),sep="")
    interac_folders <- c("observed",1:100)
  }
  
  if (!file.exists(outputfolder)){dir.create(outputfolder,recursive = T)}  
  
  for (interac_folder in interac_folders){ # interac_folder <- "observed"
    if (interac_folder!="observed"){interac_folder <- paste("random_",paste(rep(0,3-nchar(interac_folder)),collapse=""),interac_folder,sep="")} 
    if (!file.exists(paste(outputfolder,"/individual_simulation_results_",interac_folder,".txt",sep=""))){
      summary_collective <- NULL
      summary_individual <- NULL
      interac_list <- input_files[grepl(interac_folder, input_files)]
      if (seed_file!="treated_workers.txt"){interac_list <- interac_list[grepl("PreTreatment",interac_list)]}
      for (interac in interac_list){ # interac  <- interac_list[6] 
        ### get and display colony info
        root_name          <- file_name <- basename(interac)
        colony             <- unlist(strsplit(root_name, split="_"))[1]
        treatment          <- unlist(strsplit(root_name, split="_"))[2]
        period_detail      <- sub("(.*)Treatment.*", "\\1", unlist(strsplit(root_name,split="_"))[3]) # pre post?
        if (grepl("Pre",root_name)){period="before"}else{period="after"}
        cat(blue(paste0("starting ", colony, " ",  period_detail, "Treatment - ",  tolower(interaction_type)), " - ", interac_folder, "\n"))
            
        ### read interactions
        interaction_table <- read_table_with_auto_delim(interac) # DS using function to read file to work both with space and tab separated tables (observed and randomized might have different format)
        interaction_table <- interaction_table[order(interaction_table$Starttime),]
        interaction_table <- interaction_table[which(!duplicated(interaction_table)),] # removes duplicates should there be any
        
        ### DS: additional checks because the trophallaxis interaction_table is missing some columns that were not created due to my simplified design with only 1 3h block to analyze - technically not required but easier with subsequent scripts
        if (!"time_hours" %in% colnames(interaction_table)) { # Check if 'interaction_table' contains the column 'time_hours'
          interaction_table$time_hours <- ifelse(period == "before", -24, ifelse(period == "after", 0, NA)) # Create and fill 'time_hours' based on the value of 'period'
        }
        if (!"time_of_day" %in% colnames(interaction_table)) { # Check if 'interaction_table' contains the column 'time_of_day'
          interaction_table$time_of_day <- 11  # Create and fill 'time_of_day' with the value 11
        }
        
        ### get appropriate task_group list 
        colony_task_group  <- task_groups[which(task_groups$colony==colony),] # contains both task group types i.e. percentage and facet net based
        queenid            <- as.character(colony_task_group[which(colony_task_group$task_group_prop=="queen"),"tag"]) #call specific queen tag instead of fixed 665. it has to be a character to work with igraph, technically could be just set to to "1" if appropriate queen tags were used during tagging
      
        ### read tag list to define list of live ants: get tag list for colony and filter out alive ones. 
        tag <- read.tag(tag_list, colony) # function defined in main analysis
        #if (!grepl("survival",data_path)){ # for DS never survival so  simplified here... 
        alive      <- tag[which(tag$final_status=="alive"),"tag"]
        # }else{
        #   alive      <- tag[which(as.numeric(tag$death)==0|as.numeric(tag$death)>=max(interaction_table$Stopframe,na.rm=T)),"tag"]
        # }
        #
        
        ### make sure interaction table only contains ants that were alive the end
        all_interacting_ants <- unique(c(interaction_table$Tag1, interaction_table$Tag2))
        dead_ants_among_actors <-  all_interacting_ants[!all_interacting_ants %in% alive] # get all ants that die before the end of the experiment but occur in the interaction list prior to that - will be excluded form the analysis
        if(length(dead_ants_among_actors) > 0) {
          cat(red(paste("Warning:")), "\n")
          cat("In", basename(interac),"the following ants are present in the interaction table but are not alive according to tag_list:","\n",
            paste(dead_ants_among_actors, collapse = ", "),"   -->     ","The interaction list gets subsetted to remove interactions with those ants ", "\n")
          interaction_table <- interaction_table[which(interaction_table$Tag1%in%alive & interaction_table$Tag2%in%alive),] # subset interaction table
        }

        ### read seeds
        seeds              <- read.table(paste(seed_path,seed_file,sep=""),header=T,stringsAsFactors = F)
        seeds              <- seeds[which(seeds$colony==colony),"tag"]
        seeds              <- seeds[which(seeds%in%alive)]
        
        ### read time aggregation info for both current file and PostTreatment file if relevant (i.e. if reorder=T), and get start time from that
        # DS: I never have survival and as linda do not have reorder either
        #if (!grepl("survival",data_path)){
          splitfile          <- split_list[grepl(sub("(.*Treatment).*", "\\1", root_name),split_list)]  
          splitinfo          <- read.table(splitfile,header=T,stringsAsFactors = F)
          # LS: I don't have "reorder" (would always be FALSE anyway) so following statement was commented out
          # if (reorder){
          #   splitfile_Post   <- split_list[grepl(gsub("Pre","Post",root_name),split_list)]
          #   splitinfo_Post   <- read.table(splitfile_Post,header=T,stringsAsFactors = F)
          # } 
          time_start         <- min(splitinfo$time,na.rm=T)
        # }else{
        #   time_start         <- min(interaction_table$Starttime,na.rm=T)
        # }
        # 
        ####for simulations to be comparable between the Pre-Treatment and Post-treatment periods, they should start at the same time of day
        ####so in that case, Pre-Treatment interactions should be reordered
        # LS: I don't have "reorder" (would always be FALSE anyway) so following statement was commented out
        # if (reorder & period=="before") {
        #   ####start time must be the same time of day as start time of after.
        #   time_start <- min(splitinfo_Post$time,na.rm=T) - (24*3600)
        #   ####get the corresponding row number in the interaction table
        #   start_index <- min(which(interaction_table$Starttime>=time_start))
        #   ####split the interaction table into two halves
        #   interactions_unchanged <- interaction_table[start_index:nrow(interaction_table),]
        #   interactions_to_change <- interaction_table[1:(start_index-1),]
        #   ####modify frame number and time of interactions_to_change
        #   interactions_to_change[c("Starttime","Stoptime")] <- interactions_to_change[c("Starttime","Stoptime")] + 24 * 3600
        #   frame_offset <- (interactions_unchanged[nrow(interactions_unchanged),"Startframe"]+floor((interactions_to_change[1,"Starttime"]-interactions_unchanged[nrow(interactions_unchanged),"Starttime"])*2))-interactions_to_change[1,"Startframe"]
        #   interactions_to_change[c("Startframe","Stopframe")] <- interactions_to_change[c("Startframe","Stopframe")] + frame_offset
        #   interaction_table <- rbind(interactions_unchanged,interactions_to_change)
        # }
          
        ### Create antlist object
        antlist                                                 <- data.frame(tag=alive,status=0,load=0,stringsAsFactors = F)
        antlist[which(antlist$tag%in%seeds),c("status","load")] <- 1
        
        ### Change formats to match formats expected by the C++ program
        antlist$tag                                             <- as.integer(antlist$tag )
        antlist$status                                          <- as.integer(antlist$status )
        antlist$load                                            <- as.numeric(antlist$load )
        
        interaction_table$Tag1                                  <- as.integer(interaction_table$Tag1 )
        interaction_table$Tag2                                  <- as.integer(interaction_table$Tag2 )
        if ("Startframe" %in% colnames(interaction_table) && "Stopframe" %in% colnames(interaction_table)) { # fix because interaction table has slightly different names for classic or trophallactic interactions
          interaction_table$Startframe                          <- as.integer(interaction_table$Startframe)
          interaction_table$Stopframe                           <- as.integer(interaction_table$Stopframe)
        } else {
          interaction_table$Startframe                          <- as.integer(interaction_table$frame_start)
          interaction_table$Stopframe                           <- as.integer(interaction_table$frame_stop)
        }
        interaction_table$Starttime                             <- as.numeric(interaction_table$Starttime )
        interaction_table$Stoptime                              <- as.numeric(interaction_table$Stoptime )
        
        # create vector of names matching the interac_folder format
        interac_folders_named <- character(length(interac_folders))
        for (i in seq_along(interac_folders)) {
          if (interac_folders[i] == "observed") {interac_folders_named[i] <- interac_folders[i]}
          if (grepl("^\\d+$", interac_folders[i])) { # If it is numeric, apply sprintf to format it as "random_###"
            interac_folders_named[i] <- sprintf("random_%03d", as.numeric(interac_folders[i]))
          } else {interac_folders_named[i] <- interac_folders[i]}} # If it is not "observed" or numeric, keep the element as is
        
        
        #### Perform simulations ####
        cat(yellow(paste("Performing simulations...")), "\n")
        simulations <- NULL
        for (i in 1:N_SIM){ # i=1
          simulations <- rbind(simulations, 
                               data.frame(sim_number=i,
                                          simulate_transmission(i_table=interaction_table,
                                                                ant_list=antlist,
                                                                t0=time_start,
                                                                seed=i*which(seed_file==seed_files)*which(interac_folder==interac_folders_named)*which(interac==interac_list),
                                                                frame_rate=FRAME_RATE)))
          # Sys.sleep(1)
          # gc()
        }
        
        #### Summarise simulations ####
        
        #### Step 1: modify columns ####
        simulations$relative_contamination_time                      <- simulations$relative_contamination_time/3600 # express contamination time in hours instead of seconds
        simulations                                                  <- simulations[,!names(simulations)%in%c("absolute_contamination_time","infectious")] # remove unncessary columns
        simulations["contaminated"]                                  <- 1          # add a column containing information on whether the ant was contaminated or not during the simulation
        simulations["status"]                                        <- "untreated"
        simulations[which(simulations$tag%in%seeds),"status"]        <- "treated"
        
        
        #### Step 2: add individuals that did not get contaminated during the simulations ####
        all_individuals <- expand.grid(sim_number=c(1:N_SIM),tag=alive)
        simulations     <- merge(all_individuals,simulations,all.x=T,all.y=T)
        if (nrow( simulations[which(is.na(simulations$relative_contamination_time)),])>0){
          simulations[which(is.na(simulations$relative_contamination_time)),"status"]     <- "untreated"
          simulations[which(is.na(simulations$relative_contamination_time)),c("relative_contamination_time","contaminated_by","final_load","contaminated")]     <- rep(c(24,-1,0,0) ,each=nrow(simulations[which(is.na(simulations$relative_contamination_time)),]))
        }
        simulations     <- simulations[order(simulations$sim_number,simulations$relative_contamination_time),]
        
        #### Step3: add further individual-level information ####
        # Add infection rank
        ranks <- aggregate(relative_contamination_time~sim_number,function(x)rank(x,ties.method="min"),data=simulations[which(simulations$status!="treated"),])
        for (sim_number in unique(ranks$sim_number)){ #sim_number <- 1
          simulations[which(simulations$sim_number==sim_number&simulations$status=="treated"),"rank"] <- 0
          simulations[which(simulations$sim_number==sim_number&simulations$status!="treated"),"rank"] <- as.numeric(ranks[which(ranks$sim_number==sim_number),c(2:ncol(ranks))])
        }
        
        
        # Add high load/low load
        simulations["high_level"] <- as.numeric(simulations$final_load>high_threshold)
        simulations["low_level"] <- as.numeric(simulations$final_load>0&simulations$final_load<=high_threshold)
        
        
        #### Step4: summarise simulations - colony-level data ####
        
        #### Step 4.1: Prevalence, Mean load, Prop. high level, Prop. low level, Mean load ####
        colony_level         <- aggregate(na.rm=T,na.action="na.pass",cbind(contaminated,final_load,high_level,low_level)~sim_number,FUN=mean,data=simulations[which(simulations$status!="treated"),])
        names(colony_level)  <- c("sim_number","Prevalence","Mean_load","Prop_high_level","Prop_low_level")
        
        #### Step 4.2: Load skewness ####
        skewness             <- aggregate(na.rm=T,na.action="na.pass",final_load~sim_number,FUN=skewness,data=simulations[which(simulations$status!="treated"),]) # LS: only works if object "skewness" does not exist yet, when re-running this line during testing an error will occur -> rm(ls = skewness)
        names(skewness)      <- c("sim_number","Load_skewness")
        colony_level         <- merge(colony_level,skewness)
        
        # Step 4.3: Queen load
        queen                <- simulations[which(simulations$tag==queenid),c("sim_number","final_load")]
        names(queen)         <- c("sim_number","Queen_load")
        colony_level         <- merge(colony_level,queen)
        
        #### Step 4.4: Transmission rate ####
        colony_level["logistic_r"] <- NA
        for (sim in unique(simulations$sim_number)){  # sim <- 1 # LS: warnings but only that elements produced in try() do not exist in some cases
          # define table that will be used in the non linear fit
          subset <- simulations[which(simulations$sim_number==sim),]
          # remove from subset the ants that were not infected during the simulation and were added later
          subset <- subset[which(subset$contaminated!=0),]
          # sort by time
          subset <- subset[order(subset$relative_contamination_time,subset$status),]
          # get population size
          pop_size <- sum(subset[1,c("nb_susceptible","nb_contaminated")])
          #get relevant data
          spread                 <- subset[c("relative_contamination_time","rank")]
          spread["nb_seeds"]     <- nrow(subset[which(subset$status=="treated"),])
          spread["nb_contaminated"] <- spread$rank+spread$nb_seeds
          spread <- spread[which(!duplicated(spread[c("relative_contamination_time","nb_contaminated")])),]
          spread["proportion_contaminated"] <- spread$nb_contaminated/pop_size
          
          # try logistic fit
          y <- spread$proportion_contaminated
          x <- spread$relative_contamination_time
          P_zero <- spread[which(spread$rank==0),"proportion_contaminated"]
          # K <- 1
          if (exists("fit_logistic")){rm(list=c("fit_logistic"))}
          try(fit_logistic <- nls(
            y ~ (K*P_zero*exp(r*x))
            /
              (K+(P_zero*(-1+exp(r*x))))
            ,
            start=list(r=1,K=1)
            #start=list(r=1)
            ,
            control=list(maxiter = 1000)
          )
          ,silent=T)
          if (exists("fit_logistic")){
            r <- summary(fit_logistic)$coefficients["r","Estimate"]
            K <- summary(fit_logistic)$coefficients["K","Estimate"]
            colony_level[which(colony_level$sim_number==sim),"logistic_r"] <- r
          }
          rm(list=c("fit_logistic","r","K"))
        }
        
        #### Step 4.5: Take average over all simulations ####
        colony_level         <- colMeans(colony_level[names(colony_level)!="sim_number"],na.rm=T)
        
        #### Step 4.6: Store ####
        summary_collective   <- rbind(summary_collective, cbind(data.frame(colony=colony,
                                                                           colony_size=info[which(info$colony==colony),"colony_size"],
                                                                           treatment=info[which(info$colony==colony),"treatment"],
                                                                           period=period,
                                                                           period_detail=period_detail,
                                                                           time_hours=interaction_table[1,"time_hours"],
                                                                           time_of_day=interaction_table[1,"time_of_day"],
                                                                           stringsAsFactors = F),
                                                                           t(colony_level))) 
        
        #### Step5: summarize simulations - individual-level data ####
        
        #### Step 5.1: Final_load, Prob. contamination, Prob. high level, Prob. low level ####
        individual_level         <- aggregate(na.rm=T,na.action="na.pass", cbind(final_load,contaminated,high_level,low_level) ~ tag,FUN=mean,data=simulations)
        names(individual_level)  <- c("tag","simulated_load","probability_of_transmission","probability_high_level","probability_low_level")
        
        #### Step 5.2: Transmission latency and transmission rank ####
        individual_level[c("transmission_latency","transmission_rank")] <- NA
        for (ant in alive){
          if (all(simulations[which(simulations$tag==ant),"contaminated"]==1)){
            individual_level[which(individual_level$tag==ant),"transmission_latency"] <- mean(simulations[which(simulations$tag==ant),"relative_contamination_time"])
            individual_level[which(individual_level$tag==ant),"transmission_rank"]    <- mean(simulations[which(simulations$tag==ant),"rank"])
            
          }else{
            # transmission latency
            model                                                                     <- coxph(Surv(relative_contamination_time,contaminated)~1,data=simulations[which(simulations$tag==ant),])
            mean_data                                                                 <- summary(survfit(model),rmean="common")$table
            individual_level[which(individual_level$tag==ant),"transmission_latency"] <-  mean_data["rmean"] # LS: column name changed to "rmean"
            # transmission rank
            model                                                                     <- coxph(Surv(rank,contaminated)~1,data=simulations[which(simulations$tag==ant),])
            mean_data                                                                 <- summary(survfit(model),rmean="common")$table
            individual_level[which(individual_level$tag==ant),"transmission_rank"]    <-  mean_data["rmean"] # LS: column name changed to "rmean"
            
          }
        }
        #### Step 5.3: Store ####
        individual_level$antid          <- as.character(interaction(colony,individual_level$tag))
        individual_level$status         <- "untreated";individual_level[which(individual_level$tag%in%seeds),"status"] <- "treated"
        individual_level["colony"]      <- colony
        individual_level["colony_size"] <- info[which(info$colony==colony),"colony_size"]
        individual_level["treatment"]   <- info[which(info$colony==colony),"treatment"]
        individual_level["period"]      <- period
        individual_level["period_detail"]      <- period_detail     
        individual_level["time_hours"]  <- interaction_table[1,"time_hours"]
        individual_level["time_of_day"] <- interaction_table[1,"time_of_day"]
        
        individual_level       <- individual_level[c("colony","colony_size","treatment","tag","antid","status","period","period_detail",
                                                     "time_hours","time_of_day","simulated_load","probability_of_transmission","probability_high_level",
                                                     "probability_low_level","transmission_latency","transmission_rank")] 
        summary_individual     <- rbind(summary_individual,individual_level)
        clean()
      }
      
      #ordering based on time hours. 
      summary_collective <- summary_collective[order(summary_collective$colony,summary_collective$time_hours),]
      summary_individual <- summary_individual[order(summary_individual$colony,summary_individual$tag,summary_individual$time_hours),]
      
      if (grepl("random_vs_observed",outputfolder)){
        summary_collective["randy"] <- interac_folder
        summary_individual["randy"] <- interac_folder
      }
      
      write.table(summary_collective,paste(outputfolder,"/collective_simulation_results_",interac_folder,".txt",sep=""),col.names=T,row.names=F,quote=F,append=F)
      write.table(summary_individual,paste(outputfolder,"/individual_simulation_results_",interac_folder,".txt",sep=""),col.names=T,row.names=F,quote=F,append=F)
      
    }
  }
  
  if (grepl("random_vs_observed",outputfolder )){
    setwd(outputfolder)
    file_list <- list.files(pattern="individual")
    all_indiv_results <- NULL
    for (file in file_list){
      temp <- read.table(file, header=T,stringsAsFactors = F)
      if (grepl("random",unique(temp$randy))){
        temp$randy <- "random"
      }
      all_indiv_results <- rbind(all_indiv_results,temp)
    }
    all_indiv_results <- aggregate(na.rm=T,na.action="na.pass",cbind(simulated_load,probability_of_transmission,probability_high_level,probability_low_level,transmission_latency,transmission_rank)~.,FUN=mean,data=all_indiv_results)
    all_indiv_results$treatment <- all_indiv_results$randy 
    write.table(all_indiv_results,file=paste(outputfolder,"/summarised_individual_results.dat",sep=""),append=F,quote=F,row.names=F,col.names=T)
    
  }
  
  
  
}

loop_end_time <- Sys.time()
print(paste("loop took ", as.numeric(difftime(loop_end_time, loop_start_time, units = "mins")), " minutes to complete"))
to_keep <- to_keep_ori
