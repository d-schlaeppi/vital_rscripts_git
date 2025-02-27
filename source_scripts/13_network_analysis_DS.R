### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### NETWORK ANALYSIS - 13_network_analysis_DS.R ####
### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Read me ####

# Takes an interaction list as an input, builds a network, and analyses its properties

# Created by Nathalie Stroeymeyt
# Modified by Adriano Wanderlingh to work with FORT formicidae Tracking data.
# Modified by Nathalie Stroeymeyt to include number of events in addition to duration
# with adaptations by Linda Sartoris and then adjusted to the needs of Daniel Schläppi's data


#### TODO's 
#' Change to facet net task group allocation
#' Update colony metadata so that food types can be added and we can get distance to treated but split by food type (control, virus or pseudo virus)
#' 
#' 
#' RECENT UPDATES: 
#' change communtiy detection so it is fixed at three communities...



### ### ### ### ### ### ### ### ### ### ### ###
to_keep_ori <- to_keep
### ### ### ### ### ### ### ### ### ### ### ###


options(digits = 16) ; options(digits.secs = 6) ; options("scipen" = 10) #set options to control how numerical values are printed and handled in the R environment
edge_weights <- c("number", "duration") # DS: script now loops over both types, to undo change here back to one of the two, remove loop start plus the end bracket and change edge_weight back to edge_weights


### output file: if it already exists remove or relable previous version - step only needed if the script is run multiple times - maybe check if there are other output files that need to be deleted in this instance... ?
output_path <-  paste(data_path,"/processed_data/individual_behaviour/post_treatment/",sep="")
output_filename_old <- paste(output_path,"interactions_with_treated.txt",sep="")
if (file.exists(output_filename_old)){ 
  system_date <- Sys.Date()
  new_file_name <- paste(output_path, "interactions_with_treated_old_", Sys.Date(), ".txt", sep="")
  # file.remove(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))
  file.rename(output_filename_old, new_file_name)}

### UPDATE SO WE GET FACET NET MODULARITY! 
# facet_net_modularity <- paste(data_path, "")
 
if (file.exists()) # if facet net table exists read it... 

#### get input file list
if (!grepl("survival",data_path)){ 
  input_path           <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists",sep="")
  setwd(input_path)  
  input_folders        <- list.dirs(recursive=T,path="PreTreatment",full.names=F)
  # input_folders        <- input_folders[which(input_folders!="")]
  input_folders <- "observed" # delete again and activate line above 
} # else{                 # DS never has survival . 
#   input_path           <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists",sep="")
#   setwd(input_path)  
#   input_folders        <- list.dirs(recursive=T,path="PostTreatment",full.names=F)
#   input_folders        <- input_folders[which(input_folders!="")]
# }

queen_community_summary <- NULL
to_keep <- c(ls(),"to_keep","input_folder","network_files","options","option","summary_collective","summary_individual","outputfolder","network_file","queenid", "edge_weights", "edge_weight")
for (input_folder in input_folders){ # input_folder <- input_folders[1]
  cat(yellow(paste0(input_folder, "                                            ")), "\n")
  setwd(input_path)
  network_files <- list.files(path=paste("PreTreatment/",input_folder,sep=""),full.names=T)
  
  if (input_folder=="observed"){
    network_files <- c(network_files,list.files(path=paste("PostTreatment/",input_folder,sep=""),full.names=T))
    if(grepl("main",data_path)){
      options <- c("all_workers","untreated_only")
    }else{options <- c("all_workers")}
  }else{options <- c("all_workers")}
  
  for (option in options){ # option <- options[1]
    cat(blue(paste0(option, "                                            "), "\n"))
    for (edge_weight in edge_weights) { # edge_weight <- edge_weights[1] # #DS: for loop to loop over edge weights inserted
    outputfolder <- paste(data_path,"/processed_data/network_properties_edge_weights_",edge_weight,sep="")
    
    summary_collective <- NULL
    summary_individual <- NULL
    for (network_file in network_files){ # network_file <- network_files[28]    network_file <- network_files[1]
      cat("\r",network_file)
      
      ### get file metadata
      root_name          <- gsub("_interactions.txt","",unlist(strsplit(network_file,split="/"))[grepl("interactions",unlist(strsplit(network_file,split="/")))]) # LS: replace grepl("colony", ...) with grepl("interactions")
      components         <- unlist(strsplit(root_name,split="_"))
      colony             <- unlist(strsplit(root_name,split="_"))[1]
      treatment          <- unlist(strsplit(root_name,split="_"))[2] 
      colony_size        <- info[which(info$colony==colony),"colony_size"]
      if (!all(!grepl("PreTreatment",components))) {period <- "pre"} else {period <- "post"} 
      
      # here a section specific to Lindas design was deleted

      # if(!grepl("survival",data_path)){ # DS never set to survival 
        time_hours         <- as.numeric(gsub("TH","",components[which(grepl("TH",components))]))
        time_of_day        <- as.numeric(gsub("TD","",components[which(grepl("TD",components))]))
      # } 
      # here some more of lindas code was not needed
    
      ### get appropriate task_group list, treated list and tag
      colony_treated     <- treated[which(treated$colony==colony),"tag"] #AW  # DS might need double checking to see if this is the best source for treated ants. 
      colony_task_group  <- task_groups[which(task_groups$colony==colony),]
      queenid            <- as.character(colony_task_group[which(colony_task_group$task_group_FACETNET_0.5=="queen"),"tag"]) # call specific queen as a character to work with igraph | could also just be set to 1 if the queen has always been tagged with tag 0x000

      ### read interactions & remove dead ants from interactions list
      interactions       <- read.table(network_file,header=T,stringsAsFactors = F)
      tag <- read.tag(tag_list, colony)
      alive <- tag$tag
      interactions <- subset(interactions, Tag1 %in% alive)
      interactions <- subset(interactions, Tag2 %in% alive)

      ### if untreated only, reduce tag , colony_task_group, and interactions
      if (option=="untreated_only"){
        colony_task_group <- colony_task_group[which(!colony_task_group$tag%in%colony_treated),]
        tag               <- tag[which(!tag$tag%in%colony_treated),]
        interactions      <- interactions[which(!((interactions$Tag1%in%colony_treated)|(interactions$Tag2%in%colony_treated))),]}
      
      actors             <- data.frame(name=as.character(tag$tag))
      
      ### add a column containing interaction duration in min
      interactions["duration_min"] <- (interactions$Stoptime - interactions$Starttime + (1/FRAME_RATE)) /60 # duration in minutes (one frame = 0.125 second)
      interactions$N               <- 1 # each interaction is given a count of 1
      
      ### add a column containing the status of tag 1 and the status of tag2
      interactions[c("status_Tag1","status_Tag2")] <- "untreated"
      interactions[which(interactions$Tag1%in%colony_treated),"status_Tag1"] <- "treated"
      interactions[which(interactions$Tag2%in%colony_treated),"status_Tag2"] <- "treated"
      
      ### use this information to calculate, for each worker, the cumulative duration of interaction with treated workers
      if (input_folder=="observed"&option=="all_workers"){
        
        ### AW: To overcome the "Error in aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...) : no rows to aggregate", create empty df if no rows to aggregate on
        aggregated1 <- tryCatch(
          { aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Tag1+status_Tag2,FUN=sum,data=interactions[which(interactions$status_Tag2=="treated"),]) },
          error = function(e) {data.frame(tag= integer(),partner_status=character(), duration_min=numeric(), number_contacts=numeric())} # Return an empty data frame in case of an error
        )
        #aggregated1                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag1+status_Tag2,FUN=sum,data=interactions[which(interactions$status_Tag2=="treated"),])
        names(aggregated1)          <- c("tag","partner_status","duration_min","number_contacts")
        aggregated2 <- tryCatch(
          { aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Tag2+status_Tag1,FUN=sum,data=interactions[which(interactions$status_Tag1=="treated"),]) },
          error = function(e) {data.frame(tag= integer(),partner_status=character(), duration_min=numeric(), number_contacts=numeric())} # Return an empty data frame in case of an error
        )
        #aggregated2                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag2+status_Tag1,FUN=sum,data=interactions[which(interactions$status_Tag1=="treated"),])
        names(aggregated2)          <- c("tag","partner_status","duration_min","number_contacts")
        aggregated                  <- rbind(aggregated1,aggregated2)
        aggregated <- tryCatch(
          { aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,number_contacts)~tag+partner_status,FUN=sum,data=aggregated) },
          error = function(e) {data.frame(tag= integer(),partner_status=character(), duration_min=numeric(),number_contacts=numeric())} # Return an empty data frame in case of an error
        )
        #aggregated                  <- aggregate(na.rm=T,na.action="na.pass",duration_min~tag+partner_status,FUN=sum,data=aggregated)
        interactions_with_treated   <- merge(data.frame(tag=tag[which(tag$tag%in%alive),"tag"],stringsAsFactors = F),aggregated[c("tag","duration_min","number_contacts")],all.x=T)
        interactions_with_treated[is.na(interactions_with_treated$duration_min),"duration_min"] <- 0
        interactions_with_treated[is.na(interactions_with_treated$number_contacts),"number_contacts"] <- 0
        names(interactions_with_treated) <- c("tag","duration_of_contact_with_treated_min","number_of_contact_with_treated")
        interactions_with_treated["colony"] <- colony
        interactions_with_treated["time_hours"] <- time_hours
        
        ### write results
        if (grepl("main",data_path)){
          behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
          if (!"duration_of_contact_with_treated_min"%in%names(behav)){
            behav[c("duration_of_contact_with_treated_min")] <- NA
            behav[c("number_of_contact_with_treated")] <- NA
          }
          behav[match(as.character(interaction(interactions_with_treated$colony,interactions_with_treated$tag,interactions_with_treated$time_hours)),as.character(interaction(behav$colony,behav$tag,behav$time_hours))),c("duration_of_contact_with_treated_min")]  <- interactions_with_treated$duration_of_contact_with_treated_min
          behav[match(as.character(interaction(interactions_with_treated$colony,interactions_with_treated$tag,interactions_with_treated$time_hours)),as.character(interaction(behav$colony,behav$tag,behav$time_hours))),c("number_of_contact_with_treated")]        <- interactions_with_treated$number_of_contact_with_treated
          options(digits=3)
          write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
          options(digits=16)
        }
        ### interactions with treated: post-treatment
        if (period=="post"){
          outputfoldy <- paste(data_path,"/processed_data/individual_behaviour/post_treatment",sep="")
          if(!file.exists(outputfoldy)){dir.create(outputfoldy,recursive=T)}
          int_with_treated <- data.frame(colony_id = colony, colony_size=colony_size,treatment=treatment,period=period, time_of_day=time_of_day,
                                         interactions_with_treated,stringsAsFactors = F)
          if (!file.exists(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))){
            write.table(int_with_treated,file=paste(outputfoldy,"/interactions_with_treated.txt",sep=""),col.names=T,row.names=F,quote=F,append=F)
          }else{
            write.table(int_with_treated,file=paste(outputfoldy,"/interactions_with_treated.txt",sep=""),col.names=F,row.names=F,quote=F,append=T)
          }
        }
      }
      
      ### build NETWORK
      if (!grepl("survival",data_path)){
        net <- graph_from_data_frame(interactions[c("Tag1","Tag2")],directed=F,vertices=actors)
        ### add edge weights
        if (edge_weight=="number"){
          E(net)$weight <- interactions[,"N"]
        }else if (edge_weight=="duration"){
          E(net)$weight <- interactions[,"duration_min"]        
         }

        ### simplify graph (merge all edges involving the same pair of ants into a single one whose weight = sum of these weights)
        net <- simplify(net,remove.multiple=TRUE,remove.loop=TRUE,edge.attr.comb="sum")
        
        ### get all degress without removing any node
        degrees_all         <- degree(net,mode="all")
        
        ### remove unconnected nodes
        unconnected <- actors[degree(net)==0,]
        net <- net - as.character(unconnected)
        ### get actor list from net
        actors <- vertex_attr(net,"name")
        
      #### Part 1: individual network properties ####
        ### prepare table
        tag["status"] <- "untreated"; tag[which(tag$tag%in%colony_treated),"status"] <- "treated"
        individual <- data.frame(randy=input_folder,
                                 colony=colony, 
                                 colony_size=colony_size,
                                 treatment=treatment,
                                 tag=tag$tag,
                                 age=tag$age,
                                 status=tag$group,
                                 period=period,
                                 time_hours=time_hours,
                                 time_of_day=time_of_day,
                                 degree=NA,
                                 aggregated_distance_to_queen=NA,
                                 mean_aggregated_distance_to_treated=NA, 
                                 mean_aggregated_distance_to_f1v_treated=NA, # adjust/add a section below to get distance to treated ants with access to specific food source (virus/pseudovirus)
                                 mean_aggregated_distance_to_f2c_treated=NA, # treated ants with control food
                                 community_id = NA,
                                 same_community_as_queen=NA)
        ## degree
         individual[match(names(degrees_all),individual$tag),"degree"] <- degrees_all
        
        ## skip queen in grooming & trophallactic interactions interactions
        if (!grepl("grooming|trophallaxis", input_path)) {
          # communities             <- cluster_louvain(net, weights = E(net)$weight) # updated with the below lines so one can fix the number of communities to two. 
          best_partition <- net %>% cluster_fast_greedy(weights = E(net)$weight) %>% cut_at(no = 2) # no = x defines the number of communities. 
          # community_membership    <- communities$membership
          community_membership <- best_partition
          
          ## same community as queen
          queen_comm <- community_membership[which(V(net)$name==queenid)]
          community_membership <- community_membership==queen_comm
          individual[match(V(net)$name,individual$tag),"same_community_as_queen"] <- community_membership
          ## path length to queen
          if (queenid%in%actors){
            path_length_to_queen <- t(distances(net,v=actors,to=queenid,weights=1/E(net)$weight))
            individual[match(colnames(path_length_to_queen),individual$tag),"aggregated_distance_to_queen"] <- as.numeric(path_length_to_queen )
          }
        }
        
        #### Mean path length to treated; aggregated_network
        if(option!="untreated_only"){
          
          ### ADD LINES TO DUPLICATE DISTANCE TO TREATED BUT SO THAT IT BECOMES DISTANCE TO TREATED THAT WERE FED ON SPECIFIC FOOD TYPE 
          path_length_to_treated                             <- as.data.frame(as.matrix(shortest.paths(net,v=actors,to=as.character(colony_treated)[as.character(colony_treated)%in%V(net)$name],weights=1/E(net)$weight)))
          ### here modify
          path_length_to_treated["mean_distance_to_treated"] <- NA
          path_length_to_treated$mean_distance_to_treated    <- as.numeric(rowMeans(path_length_to_treated,na.rm=T))
          individual[match(rownames(path_length_to_treated),individual$tag),"mean_aggregated_distance_to_treated"] <- path_length_to_treated[,"mean_distance_to_treated"]
        }
         
        ### Add data to main data table
        summary_individual <- rbind(summary_individual,individual)
        
        
        #### Part 2: collective network properties ####
        
        ## remove outliers (e.g. ants that have only 1 or 2 connections)
        # outlier_p <- 0
        # outlier_removed <- c()
        # while(outlier_p<0.05){
        #   outlier_p <- grubbs.test(sort(degree(net))[1:min(30,length(degree(net)))])$p.value ###do the test on lowest 30 values first as test does not work as well when there are too many nodes
        #   if (outlier_p<0.05){
        #     outlier_removed <- c(outlier_removed,actors[which.min(degree(net))])
        #     net <- net - as.character(actors[which.min(degree(net))])
        #   }
        # }
        outlier_removed <- actors[which(strength(net)<(5/60))]  # this probably replaced the above outlier removal part... 
        net <- net - as.character(outlier_removed)
        
        ### update actor list
        actors <- get.vertex.attribute(net,"name")
        
        ## Assortativity - Age
        ## if age experiment, get colony ages
        # if (grepl("age",data_path)){
        #   colony_ages <- ages [which(ages$colony==colony),]
        #   ####set queen age to NA as this would bias the result (the queen is the oldest individual and interacts mostly with the young nurses)
        #   colony_ages[which(colony_ages$tag==queenid),"age"] <- NA
        #   ####order the age according to the order of the network's vertices
        #   ordered_ages        <- colony_ages[match(V(net)$name,as.character(colony_ages$tag)),"age"]
        #   #### calculate age assortativity
        #   age_assortativity <- assortativity(net-V(net)$name[is.na(ordered_ages)],types1=ordered_ages[!is.na(ordered_ages)],directed=F)
        # }else{
          age_assortativity <- NA
        # }
        
        ## Assortativity - Task
        ordered_task_groups <- colony_task_group[match(actors,as.character(colony_task_group$tag)),"task_group_FACETNET_0.5"]     #### done with the facet net task group allocation
        ordered_task_groups[is.na(ordered_task_groups)] <- "nurse"
        ordered_task_groups <- as.numeric(as.factor(ordered_task_groups))
        task_assortativity  <- assortativity_nominal(net,types=ordered_task_groups,directed=F)  ### 
        ## Clustering
        clustering <- mean(transitivity(net,type="barrat",weights=E(net)$weight,isolates = c("NaN")),na.rm=T)
        ## Degree mean and max
        degrees         <- degree(net,mode="all")
        degree_mean     <- mean(degrees,na.rm=T)
        degree_maximum  <- max(degrees,na.rm=T)
        ## Density
        density  <- igraph::edge_density(net)
        ## Diameter
        diameter <- igraph::diameter(net,directed=F,unconnected=TRUE,weights=(1/E(net)$weight)) ### here use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
        ## Efficiency
        net_dist                    <- shortest.paths(net, weights=1/E(net)$weight, mode="all") ## again use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
        net_dist[net_dist==0]       <- NA ##remove distances to self
        efficiency                  <- 1/net_dist ##transform each distance into an efficiency
        efficiency <- (1/((vcount(net)*(vcount(net)-1))))*(sum(efficiency,na.rm=TRUE))
        ## Modularity
        best_partition             <-  net %>% cluster_fast_greedy(weights = E(net)$weight) %>% cut_at(no = 2)
        community_membership    <- as.numeric(best_partition)
        names(community_membership) <- V(net)$name
        individual <- individual %>% mutate(community_id = community_membership[as.character(tag)])
        #
        nr_communities <- length(unique(community_membership))
        modularity              <- modularity(net,community_membership,weights=E(net)$weight)
        modularity              <- modularity(net,best_partition,weights=E(net)$weight)
        nr_community_members      <- paste("[",paste(table(best_partition),collapse=", "),"]",sep="")
        
        modularity_facetnet <- NA
        ####CONTINUE HERE!!!
        # modularity_facetnet <- UPDATE SO WE CAN GET MODULARITY CALCULATED VIA FACET NET... 
        
        ### Add to data
        summary_collective <- rbind(summary_collective,data.frame(randy=input_folder,
                                                                  colony=colony,
                                                                  colony_size=colony_size,
                                                                  treatment=treatment,
                                                                  period=period,
                                                                  time_hours=time_hours,
                                                                  time_of_day=time_of_day,
                                                                  age_assortativity=age_assortativity,
                                                                  task_assortativity=task_assortativity,
                                                                  clustering=clustering,
                                                                  degree_mean=degree_mean,
                                                                  degree_maximum=degree_maximum,
                                                                  density=density,
                                                                  diameter=diameter,
                                                                  efficiency=efficiency,
                                                                  modularity=modularity,
                                                                  modularity_facetnet=modularity_facetnet,
                                                                  nr_communities=nr_communities,
                                                                  nr_community_members=nr_community_members,
                                                                  nb_unconnected=length(unconnected),
                                                                  nb_outliers_removed=length(outlier_removed),
                                                                  stringsAsFactors = F))
        
       }
      clean()
    #}# TEMPORARY EXCLUSION 
    }
    
    #print progress AW
    cat("\n End of network_files processing >> writing")
    
    #### write #####
    if (!grepl("survival",data_path)){
      if (input_folder=="observed"){
        ### Main experiment: write pre_vs_post_treatment data
        if (!grepl("age",data_path)){
          outputfolder2 <- paste(outputfolder,"pre_vs_post_treatment",option,sep="/")
          if(!file.exists(outputfolder2)){dir.create(outputfolder2,recursive=T)}
          write.table(summary_collective[,names(summary_collective)!="randy"],file=paste(outputfolder2,"/colony_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          write.table(summary_individual[,names(summary_individual)!="randy"],file=paste(outputfolder2,"/individual_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        }
        ### All workers: write pre_treatment data into random_vs_observed folder
        if (option=="all_workers"){
          outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
          if(!file.exists(outputfolder3)){dir.create(outputfolder3,recursive=T)}
          write.table(summary_collective[which(summary_collective$period=="pre"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          write.table(summary_individual[which(summary_individual$period=="pre"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)

          outputfolder4 <- paste(outputfolder,"post_treatment",sep="/")
          if(!file.exists(outputfolder4)){dir.create(outputfolder4,recursive=T)}
          write.table(summary_collective[which(summary_collective$period=="post"),],file=paste(outputfolder4,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) 
          write.table(summary_individual[which(summary_individual$period=="post"),],file=paste(outputfolder4,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) 
        }
        ### Main experiment, All workers: add pre_treatment node_properties information to pre_treatment behaviour file
        if (!grepl("age",data_path)&option=="all_workers"){
          pre_treatment_behav_file <- paste(data_path,"/processed_data/individual_behaviour/pre_treatment/network_position_vs_time_outside.dat",sep="")
          pre_treatment_behav      <- read.table(pre_treatment_behav_file,header=T,stringsAsFactors = F)
          names(summary_individual)[which(names(summary_individual)=="aggregated_distance_to_queen")] <- paste("aggregated_distance_to_queen_edge_weights_",edge_weight,sep="")
          pre_treatment_behav      <- merge(pre_treatment_behav,summary_individual[which(summary_individual$period=="pre"),c("colony","tag","time_hours","degree",paste("aggregated_distance_to_queen_edge_weights_",edge_weight,sep=""))],all.x=T,all.y=T) 
          pre_treatment_behav      <- pre_treatment_behav[order(pre_treatment_behav$colony,pre_treatment_behav$tag,pre_treatment_behav$time_hours),]
          write.table(pre_treatment_behav, file=pre_treatment_behav_file,col.names=T,row.names=F,quote=F,append=F) # LS: this file adds columns to "network_position_vs_time_outside.dat" which already has the information of period_detail in it
        }
      }else{
        outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
        if(!file.exists(outputfolder3)){dir.create(outputfolder3,recursive=T)}
        write.table(summary_collective[which(summary_collective$period=="pre"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) #
        write.table(summary_individual[which(summary_individual$period=="pre"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) #
        
      }
    }
    
    ### Get characteristics of queen community vs. other communities (worker age, prop. of foragers)
    if (!grepl("survival",data_path)&option=="all_workers"){ 
      print(" End of writing >> Get characteristics of queen community vs. other communities")
      summary_individual_pre <- read.table(paste(outputfolder,"/random_vs_observed/node_properties_",input_folder,".txt",sep=""),header = T,stringsAsFactors = F) 
      summary_individual_pre <- summary_individual_pre[which(summary_individual_pre$period=="pre"),] 
      
      # ### if necessary: add age
      # if (grepl("age",data_path)){
      #   summary_individual_pre <- merge(summary_individual_pre[,which(names(summary_individual_pre)!="age")],ages,all.x=T,all.y=F)
      # }else{
        summary_individual_pre$age <- NA
      # }
      
      ### add task_group
      summary_individual_pre <- merge(summary_individual_pre,task_groups,all.x=T,all.y=F)
      ### remove queen
      summary_individual_pre <- summary_individual_pre[which(summary_individual_pre$tag!=queenid),]
      
      ### 1. calculate mean proportion of foragers depending on with vs without queen
      summary_individual_pre["forager"] <- 0
      summary_individual_pre[which(summary_individual_pre$task_group_FACETNET_0.5=="forager"),"forager"] <- 1
      ## skip queen in grooming an trophallactic interactions
      if (!grepl("grooming|trophallaxis",input_path)) { 
      prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+time_hours+same_community_as_queen,FUN=mean,data=summary_individual_pre) 
      prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+same_community_as_queen,FUN=mean,data=prop_foragers) 
      names(prop_foragers)[names(prop_foragers)=="same_community_as_queen"] <- "in_queen_comm";names(prop_foragers)[names(prop_foragers)=="forager"] <- "proportion_of_foragers"
      }else{
        prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+time_hours,FUN=mean,data=summary_individual_pre) 
        prop_foragers$in_queen_comm <- NA
        names(prop_foragers)[names(prop_foragers)=="forager"] <- "proportion_of_foragers"
      }
      ### 2. calculate mean age of workers depending on with vs without queen
      # if (grepl("age",data_path)){
      #   mean_age <- aggregate(na.rm=T,na.action="na.pass",age~colony+randy+colony_size+treatment+period+time_hours+same_community_as_queen,FUN=mean,data=summary_individual_pre) 
      #   mean_age <- aggregate(na.rm=T,na.action="na.pass",age~colony+randy+colony_size+treatment+period+same_community_as_queen,FUN=mean,data=mean_age) 
      #   names(mean_age)[names(mean_age)=="same_community_as_queen"] <- "in_queen_comm"
      #   prop_foragers <- merge(prop_foragers,mean_age,all.x=T)
      # }else{
        prop_foragers$age <- NA
      # }
      
      
      prop_foragers[which(prop_foragers$in_queen_comm=="FALSE"),"in_queen_comm"] <- "not_with_queen"
      prop_foragers[which(prop_foragers$in_queen_comm=="TRUE"),"in_queen_comm"] <- "with_queen"
      if (grepl("random",input_folder)){
        prop_foragers["randy"] <- "random"
      }
      queen_community_summary <- rbind(queen_community_summary,prop_foragers)
     }
    } # DS one bracket added here for edge weight (before the script was run twice, now it is loops over both options: duration and number. )
  } 
}


if (!grepl("survival",data_path)){
  queen_community_summary <- aggregate(na.rm=T,na.action="na.pass",cbind(proportion_of_foragers,age)~.,FUN=mean,data=queen_community_summary)
  queen_community_summary <- queen_community_summary[order(queen_community_summary$randy,queen_community_summary$colony),]
  queen_community_summary$treatment <- queen_community_summary$randy
  if (!file.exists(paste(outputfolder,"/random_vs_observed",sep=""))){dir.create(paste(data_path,"/processed_data/network_properties/random_vs_observed",sep=""),recursive=T)}
  write.table(queen_community_summary,file=paste(outputfolder,"/random_vs_observed/queen_community.dat",sep=""),append=F,quote=F,row.names=F,col.names=T)
}
to_keep <- to_keep_ori



