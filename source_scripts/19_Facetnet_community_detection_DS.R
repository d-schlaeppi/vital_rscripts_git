
### ### ### ### ### ### ### ### ### ### ### ###
#### Facetnet_community_detection.R ####
### ### ### ### ### ### ### ### ### ### ### ###

#### Read me ####
# Created by TO Richardson and Adriano Wanderlingh & adjusted by DS to the needs of Daniel Schläppi
# Takes an interaction list as an input, builds a network to provide a continuous measure of workers' social maturity
# and from social maturity creates the nurse and worker community. 
# Recent update: Paralellization implemented so that multiple cores are used.

# requires the following:
# https://c4science.ch/source/facet_unil/
# TO Richardson, T Kay, R Braunschweig, OA Journeau, M Rüegg, ... Ant behavioral maturation is mediated by a stochastic transition between two fundamental states. Current Biology 31, 1-8



### Second part still needs to be updated to my's needs

#### Prerequisites ####
# Set up directories and parameters

# FACENET community partition parameters
m                   <- 2          ## how many communities do we want to instruct FACETNET to search for...? Nurses and Foragers so 2?
alpha               <- 0.5        ## used for modulating the `memory` - how much the community structure @ time t influences that at t+1 - only matters when we ask facetnet to pay attention to this...
t_step              <- 0          ## not important
N_ITERATIONS        <- 50        ## number of iterations to find best modularity: the modularity of the found solutions vary quite a bit --> repeat multiple times & select the highest-modularity solution 
UPDATE_DATA         <- FALSE      ## updating of tables - e.g. adding modularity based on facet_net, might need some fixing as well...
DISPLAY_RESULTS_I   <- FALSE       ## optional plotting of results in loop #
DISPLAY_RESULTS_II  <- TRUE       ## optional plotting of results after loop # might need some fixing... 

# set directories to match original version of script ### should already be defined as we are sourcing this from the main_analysis script anyways. 
# DATADIR <- paste("/media",usr, hd, "vital/fc2",sep="/") 
# SCRIPTDIR <- paste("/home",usr,"Documents/vital_rscripts_git",sep="/") # does not (yet) work on mac - just keep working on linux

# to match with the code below
data_path <- paste0(DATADIR, "/vital_experiment/main_experiment")
# code_path  <- paste0(SCRIPTDIR, "/source_scripts") # already exists due to running main script - if things run it can be deleted
# source(paste(code_path,"/functions_and_parameters_DS.R",sep=""))  # so far I have not seen the use of this functions. Clean() already defined due to running main analysis

FACETNET_DIR <- paste(DATADIR, "facetNet", sep = "/")  ## the Python script is here - can be anywhere, but this must point to it...

# interaction lists (observed)
input_path           <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists",sep="")
setwd(input_path)  
input_folders        <- list.dirs(recursive=T,path="PreTreatment",full.names=F)
input_folders        <- input_folders[which(input_folders!="")]

## load the file containing the % time each ant spent outside
TaskStats_File <- paste(data_path, "processed_data/individual_behaviour/pre_treatment/network_position_vs_time_outside.dat", sep="/")
TaskStats_all      <- read.table(TaskStats_File , header=T, stringsAsFactors = F)

# where scores are saved
WORKDIR      <- paste(data_path,"Soft_community_scores_duration",sep="/")
to_keep <- c(ls(),"to_keep","network_files","network_file","output_folders","output_folder", "num_cores", "process_iteration")


# define function to run computation of community partition in parallel on multiple cores 
num_cores <- detectCores() - 2  # number of cores to use
process_iteration <- function(ITER) {
  FACETNET_REP_OUTPUT_DIR_M <- paste(FACETNET_REP_folder, paste(Cassette, ",m=", m, ",Iteration=", ITER, sep=""), sep="/") # Create new subdirectory for community results with m
  if (!file.exists(FACETNET_REP_OUTPUT_DIR_M)) {dir.create(FACETNET_REP_OUTPUT_DIR_M)}
  
  # Check if the outputs already exist and if not proceed with calculation
  if (!file.exists(paste(FACETNET_REP_OUTPUT_DIR_M, "soft_comm_step_alpha0.5_nw0.csv", sep="/"))) {
    FACETNET_INPUT_brackers <- paste0("'", FACETNET_INPUT, "'")
    FACETNET_REP_OUTPUT_DIR_M_brackers <- paste0("'", FACETNET_REP_OUTPUT_DIR_M, "'")
    command <- paste("python3", paste(FACETNET_DIR, "facetnet_step.py", sep = "/"),
                     FACETNET_INPUT_brackers, alpha, m,
                     FACETNET_REP_OUTPUT_DIR_M_brackers, t_step, sep=" ")
    OutPut <- system(command, intern=TRUE) # run command and capture output
    
    # error handling
    if ("TRUE" %in% grepl("rror", OutPut)) {
      print("ERROR in OutPut")
      print(OutPut)
      return(NULL)  # move on to next... 
    } else { # get modularity 
      Item <- grep("modularity", OutPut)
      Modularity <- as.numeric(gsub(" ", "", gsub("\\(", "", gsub("\\)", "", strsplit(OutPut[Item], "=")[[1]][2]))))
    }
    # Stack modularity for each m
    Modules <- data.frame(colony, treatment, period, ITER, alpha, MODULARITY=Modularity)
    # Check if the modularity for m communities has already been calculated and recorded
    if (file.exists(Module_File)) {
      Modules_precomputed <- read.table(Module_File, header=TRUE)
      if (!ITER %in% Modules_precomputed$ITER) {write.table(Modules, file=Module_File, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)}
    } else {write.table(Modules, file=Module_File, row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)}
  }
}








if (RUN_19_SECOND_SUBSECTION_ONLY != TRUE){
  #### ifelse() here insert something so that inpüt folder only includes observed when running first part only.
  # input_folders        <- input_folders[which(input_folders!="random")]
for (input_folder in input_folders){ # input_folder <- input_folders[2]
  #print(input_folder)
  cat(blue("\n input folder = ",input_folder, "\n"))
  setwd(input_path)
  network_files <- list.files(path=paste("PreTreatment/",input_folder,sep=""),full.names=T)
  
  # add a section to filter out any files that have already been processed completely
  completed_files <- list.files(path=paste(WORKDIR, input_folder,sep="/"),full.names=T)
  completed_files_with_iteration_100 <- completed_files[sapply(completed_files, function(file) {   # Filter completed files by checking if a subfolder with "Iteration=100" exists
    subfolders <- list.dirs(path = file, recursive = FALSE, full.names = TRUE)
    any(grepl("Iteration=100", subfolders))})]
  completed_files_base <- gsub("_FACETNET_iters$", "", basename(completed_files_with_iteration_100))
  network_files_base <- basename(network_files); network_files_base <- gsub("_interactions.txt$", "", network_files_base)
  network_files <- network_files[!network_files_base %in% completed_files_base]
  
  output_folder <- file.path(WORKDIR, input_folder)
  if (!file.exists(output_folder)){dir.create(output_folder, recursive = TRUE)}
  
  ## file in which queen community is not the same as the one where ants spend least time outside  ###these two lines were adriano specific and can probably be deleted.
  # network_files <- network_files[!grepl("colony07SP_pathogen.small_PreTreatment",network_files)] 
  
  pb <- progress_bar$new(
    format = "Progress: :current/:total [:bar] :percent ETA: :eta",
    total = length(network_files),  # Use the length of network_files for total
    clear = FALSE,
    width = 80
  )
  
  if (length(network_files) == 0) {cat("No files left to process --> All files have been processed before (or ERROR?)  \n")}
  
  for (network_file in network_files){ # network_file <- network_files[1]
    
    ### re-open the task_groups (where results are saved) as it gets updated at the end of the loop (the loaded version has to be the newest one)
    task_groups    <- read.table(paste(data_path,"original_data/task_groups.txt",sep="/"),header=T,stringsAsFactors = F)
    # Create new columns if they don't exist in 'task_groups'
    if(!"task_group_FACETNET_0.5" %in% names(task_groups)) {task_groups$task_group_FACETNET_0.5 <- NA}
    if(!"Forager_score" %in% names(task_groups)) {task_groups$Forager_score <- NA}
    
    
    ### get file metadata
    root_name          <- gsub("_interactions.txt","",unlist(strsplit(network_file,split="/"))[grepl("interactions",unlist(strsplit(network_file,split="/")))]) # LS: replace grepl("colony", ...) with grepl("interactions")
    Cassette           <- paste(root_name,"_interactions", sep="")
    #components         <- unlist(strsplit(root_name,split="_"))
    colony             <- unlist(strsplit(root_name,split="_"))[1]
    treatment          <- unlist(strsplit(root_name,split="_"))[2]
    #colony_size        <- info[which(info$colony==colony),"colony_size"]
    period             <- "Pre"
    
    # print(basename(root_name)) #"\r",
    
    #MAKE A FOLDER PER COLONY
    FACETNET_REP_folder <- file.path(output_folder, paste(root_name,"_FACETNET_iters",sep="")); if (!file.exists(FACETNET_REP_folder)){dir.create(FACETNET_REP_folder)}
    Module_File        <- paste(FACETNET_REP_folder, paste(Cassette,"Modularities.txt",sep="_"), sep="/")
    
    ### get appropriate task_group list, treated list and tag
    tag <- read.table(paste0(data_path, "/original_data/tag_files/", colony, "_", treatment, ".txt") , header=T,stringsAsFactors = F)
    alive <- tag$tag #AW
    colony_task_group  <- task_groups[which(task_groups$colony==colony), c("tag","task_group_prop")];   colony_task_group <- subset(colony_task_group, tag %in% alive)
    queenid            <- as.character(colony_task_group[which(colony_task_group$task_group=="queen"),"tag"])
    TaskStats          <- TaskStats_all[which(TaskStats_all$colony==colony),]
    TaskStats          <- subset(TaskStats, tag %in% alive)
    
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
    ##### PART 1: re-format the aggregated interactions list for facetnet soft community labelling    #####
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    
    if (!file.exists("FACETNET_INPUT")){
      # read interactions & #remove dead ants from interactions list
      interactions       <- read.table(network_file,header=T,stringsAsFactors = F)
      interactions <- subset(interactions, Tag1 %in% alive)
      interactions <- subset(interactions, Tag2 %in% alive)
      ## time-aggregated edge list
      #EdgeList              <- aggregate(Box ~ Tag1 + Tag2, FUN=length, data=interactions)  ; colnames(EdgeList)[ncol(EdgeList)] <- "Count"       ## pairwise contact count
      EdgeList              <- aggregate(duration ~ Tag1 + Tag2, FUN=sum, data=interactions)  ; colnames(EdgeList)[ncol(EdgeList)] <- "duration"       ## pairwise contact duration 
      EdgeList$duration <- round(EdgeList$duration)
      
      ## DEFINE INPUT FILE - needs to be formatted for facetnet python to understand ---????? What do adriano and tom mean here? hopefully it is just right as is... 
      FACETNET_INPUT <- paste(FACETNET_REP_folder,"/",Cassette,"interactions_Time-aggregated_network,FACETNET_format.txt",sep="")
      
      ## Export the aggregated edge list
      write.table(EdgeList, file=`FACETNET_INPUT`, row.names=F, col.names = F, quote=F)
    }
    
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
    #### PART 2: Repeatedly apply facetnet to generate a community partition                           ####
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
    
    mclapply(1:N_ITERATIONS, process_iteration, mc.cores=num_cores)
    
    # here the old for (ITER in 1:N_ITERATIONS)  {} loop was removed. 
   
    
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    #### PART 3: Find the top-modularity solution & assign biological labels to both communities       ####
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    
    if (input_folder=="observed") {
      # reduce these files to the current colony, treatment,period
      TaskStats_ColTreatPer <- TaskStats[which(TaskStats$colony==colony & TaskStats$treatment==treatment & tolower(TaskStats$period)==tolower(period) ), c("tag","time_hours","prop_time_outside")]
      
      ## **KEY STAGE** average the % time outside across the 3 hour pre-introduction bins ## IS THIS RIGHT???
      
      TaskStats_ColTreatPer_MEAN <- aggregate(prop_time_outside ~ tag, FUN=mean, na.rm=TRUE, na.action=NULL, TaskStats_ColTreatPer)
      if (nrow(TaskStats_ColTreatPer_MEAN) != nrow(colony_task_group)) {print("ERROR : the number of tags in TaskStats_ColTreatPer_MEAN is different from that in colony_task_group")}
      
      if (file.exists(Module_File)){
        Modules <- read.table(Module_File, header = T)
        ITER    <- Modules$ITER [which.max(Modules$MODULARITY)]
        
        ## load community scores for the top solution 
        FACETNET_REP_OUTPUT_DIR_M <- paste(FACETNET_REP_folder, paste(Cassette,",m=",m,",Iteration=",ITER,sep=""), sep="/")
        ## load the community scores
        Ant_Modules <- read.csv( file = paste(FACETNET_REP_OUTPUT_DIR_M, "soft_comm_step_alpha0.5_nw0.csv", sep="/")); colnames(Ant_Modules)[match("id",colnames(Ant_Modules))] <- "tag"
        
        ## convert continuous -> binary communities
        Ant_Modules$cluster_binary <- apply( Ant_Modules[,c("cluster_0","cluster_1")], 1, function(x) {c(0,1)[which.max(x)]} ) 
        
        ## add the caste information 
        Combined <- merge(x = TaskStats_ColTreatPer_MEAN, y = colony_task_group, by="tag", all =TRUE)
        ## add the scores for the 2 as-yet unlabelled communities to the df containing % time outside
        Combined <- merge(x = Combined,                   y = Ant_Modules,    by="tag", all =TRUE)
        
        ##  Label the foragers using average time spent outside (& maybe also the distance to the queen)
        ClusterMeans <- aggregate(cbind(prop_time_outside) ~ cluster_binary, FUN=mean, na.rm=T, Combined)
        
        ## which community has the lower mean % time outside? --> THE NURSES
        Nurses_by_prop_time_outside <-  ClusterMeans$cluster_binary [which.min(ClusterMeans$prop_time_outside)]
        # ## which community contains the queen? ---> THE NURSES
        # Nurses_by_queen_affiliation <-  Combined$cluster_binary [ which(Combined$task_group=="queen") ]
        # 
        # COMMENTED OUT AS colony07SP_pathogen.small_PreTreatment has disagreement between time spent with queen and time outside
        # ##  check agreement between the 2 diagnosis methods - if they don't agree, we will need to discuss
        # if (Nurses_by_prop_time_outside != Nurses_by_queen_affiliation){
        #   print("PROBLEM: community identification by % time outside DISAGREES with that by queen affiliation")}
        # else{  
        ## the nurse community is the one that spends the least time outside and that includes the queen (assuming the 2 methods agree...)
        NURSE_COMMUNITY   <- ClusterMeans$cluster_binary [ which.min(ClusterMeans$prop_time_outside)]
        FORAGER_COMMUNITY <- ClusterMeans$cluster_binary [ which.max(ClusterMeans$prop_time_outside)]
        
        ## alter the column names in Ant_Modules to reflect task group labels - a continuous score in the range 0-1 for both modules, SUMMING TO 1 ACROSS MODULES!!
        colnames(Combined) [match( paste("cluster_",FORAGER_COMMUNITY,sep=""), colnames(Combined))] <- "Forager_score"
        colnames(Combined) [match( paste("cluster_",NURSE_COMMUNITY,sep=""),   colnames(Combined))] <- "Nurse_score"
        
        ## add categorical versions of the FACETNET scores - just to compare with the original threshold-based task_group:
        Combined$cluster_binary <- NULL
        Combined$task_group_FACETNET_0.5 <- NA
        Combined$task_group_FACETNET_0.5 [which(Combined$Forager_score > Combined$Nurse_score)] <- "forager"
        Combined$task_group_FACETNET_0.5 [which(Combined$Forager_score < Combined$Nurse_score)] <- "nurse"
        Combined$task_group_FACETNET_0.5 [which(Combined$tag == queenid)] <- "queen"
        #}
        
        if(DISPLAY_RESULTS_I){
          ## tally the agreement between the original and new (facetnet) task labels
          Agreement <- table(Combined$task_group_prop, Combined$task_group_FACETNET_0.5)
          rownames(Agreement) <- paste("orig", rownames(Agreement),sep="_")
          colnames(Agreement) <- paste("facet", colnames(Agreement), sep="_")
          print(Agreement)
        }
        
        # Save output to the task_groups file
        Combined$colony <- colony
        
        # Merge by specified columns
        merged_data <- merge(task_groups, Combined[ ,which(names(Combined) %in%  c("colony", "tag", "Forager_score", "task_group_FACETNET_0.5"))],
                             by=c("colony","tag"), all.x=TRUE, suffixes=c("", ".update"))
        
        # Update FACETNET results based on the specific colony
        # !is.na(merged_data$task_group_FACETNET_0.5.update) checks if there is an updated value available for the task_group_FACETNET_0.5 column
        merged_data$task_group_FACETNET_0.5 <- ifelse(merged_data$colony == colony & !is.na(merged_data$task_group_FACETNET_0.5.update), merged_data$task_group_FACETNET_0.5.update, merged_data$task_group_FACETNET_0.5)
        merged_data$Forager_score <- ifelse(merged_data$colony == colony & !is.na(merged_data$Forager_score.update), merged_data$Forager_score.update, merged_data$Forager_score)
        merged_data$Forager_score <- round(merged_data$Forager_score,3)
        merged_data <- merged_data[order(merged_data$colony),]
        
        write.table(merged_data[, names(task_groups)],
                    file=paste(data_path,"original_data/task_groups.txt",sep="/"),
                    append = F, col.names = T, row.names = F, quote = F, sep = "\t")
      }
    }
    pb$tick()
  }
  clean()
 }
}










### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

if (RUN_19_FIRST_SUBSECTION_ONLY != TRUE) {

#### DISPLAY RESULTS  ####

if(DISPLAY_RESULTS_II){

  ###  fetch best modularity values and save output  ###
  # # Define directory and patterns
  folder_pattern <- "FACETNET_iters"
  file_pattern <- "interactions_Modularities.txt"
  data <- NULL
  
  # Get list of folders matching the pattern
  randys <- list.dirs(WORKDIR, full.names = TRUE, recursive = FALSE)
  for(randy in randys){ # randy <- randys[1]
    folders <- list.dirs(randy, full.names = TRUE, recursive = FALSE)
    folders <- folders[grep(folder_pattern, folders)]
    
    # Read and rbind files together
    data_list <- lapply(folders, function(folder) {
      file_path <- list.files(folder, pattern = file_pattern, full.names = TRUE)
      if(length(file_path) > 0) {
        return(read.table(file_path, header = TRUE, stringsAsFactors = FALSE))
      }
      return(NULL)
    })
    
    data_randy <- do.call(rbind, data_list)
    data_randy$randy <- basename(randy)
    
    data_randy <- data_randy %>%
      group_by(colony) %>%
      filter(MODULARITY == max(MODULARITY)) %>%
      ungroup()
    
    data_randy$randy <- basename(randy)
    
    data <- rbind(data, data_randy)
  }
  #save best modularity output
  write.table(data, file= paste0(WORKDIR,"_modularity.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
  
   
  ### ### ###      PLOT results of OBSERVED ONLY     ### ### ### 
  # task_groups_A    <- read.table("/media/cf19810/Seagate Portable Drive/Lasius-Bristol_pathogen_experiment/main_experiment/original_data/task_groups.txt",header=T,stringsAsFactors = F)
  task_group_file <- paste(data_path, "original_data/task_groups.txt", sep = "/" )
  task_groups_A    <- read.table(task_group_file ,header=T,stringsAsFactors = F)
  task_groups_A$treatment      <- as.factor(task_groups_A$treatment) 
  
  # Generate plots of Social Maturity
  ggplot(task_groups_A, aes(x = Forager_score, colour = treatment)) + 
    #geom_density(alpha = 0.6,size=1.5, adjust=1/1.2) + # 'adjust' changes the smoothing
    geom_line(aes(color=treatment,group = colony), stat="density", linewidth=1, alpha=0.2, adjust=1/1.2) +
    geom_line(aes(color=treatment), stat="density", linewidth=2, alpha=1, adjust=1/1.2) +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed",colour="grey20") + 
    geom_vline(aes(xintercept = 1/4), linetype = "dashed",colour="grey20") + 
    #facet_wrap(~ colony, scales = "free") +
    theme_minimal() +
    xlab("Social maturity (duration contacts)")
  
  
  ### Check when best modularity value is achieved ###
  data_obs <- data[which(data$randy=="observed"),]
  
  # Identify and rank the top modularity values for each colony
  threshold <- 0.95
  
  data_ranked <- data_obs %>%
    group_by(colony) %>%
    mutate(max_modularity = max(MODULARITY),
           is_top = ifelse(MODULARITY >= max_modularity * threshold, 1, 0)) %>%
    filter(is_top == 1) %>%
    mutate(rank = dense_rank(desc(MODULARITY)))
  
  # Count occurrences of each ITER for each rank
  iter_counts <- data_ranked %>%
    group_by(ITER, rank) %>%
    tally()
  
  # Dynamically generate breaks and labels for ranks
  max_rank <- max(iter_counts$rank)
  all_labels <- c("Best", "Second Best", "Third", "Fourth", "Fifth")  # Extend this list if more ranks are expected
  breaks_ranks <- 1:max_rank
  labels_ranks <- all_labels[1:max_rank]
  
  # Check if any NA labels and remove them
  valid_indices <- !is.na(labels_ranks)
  breaks_ranks <- breaks_ranks[valid_indices]
  labels_ranks <- labels_ranks[valid_indices]
  
  iter_counts <- iter_counts[which(iter_counts$rank<=5),]
  
  # Plot
  ggplot(iter_counts, aes(x = factor(ITER), y = n, fill = factor(rank))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colorRampPalette(c("lightblue", "blue"))(max(iter_counts$rank)),
                      name = "Rank",
                      breaks = breaks_ranks,
                      labels = labels_ranks) +
    labs(title = "Distribution of ITER for Top Modularity Ranks",
         x = "ITER",
         y = "Count") +
    theme_minimal()
}

### Not really sure what the last plot is supposed to show? If we run the script as above, then there is only one value per colony? and I do not see ŵhy iteration matters?
### This probably only matters if we are looking at multiple time bins? I only look at one three hour block while usually it is multiple 3h bins. --> confirm if this is true when doing Flugus



# maybe the following needs to be run with once the randomised or simulated things are ones are run as well? so we have observed or simulated in addition to observed? 
# #### Merge cluster output ####
# # Set the path to the directory
# path <- paste0(data_path,"/Soft_community_scores_duration")
# 
# # List all the text files in the directory
# files <- list.files(path, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
# 
# # Read and combine all the files
# combined_data <- do.call(rbind, lapply(files, read.table, header=TRUE, stringsAsFactors=FALSE))
# 
# combined_data$JOBTYPE <- NULL
# combined_data$ITER <- NULL
# combined_data$alpha <- NULL
# # Save the combined data to a new file
# output_path <- "/media/cf19810/Seagate Portable Drive/Lasius-Bristol_pathogen_experiment/main_experiment/Soft_community_scores_duration_1ITER_0.5alpha_full.txt"
# write.table(combined_data, file=output_path, row.names=FALSE, sep="\t", quote=FALSE)


### Continue adjusting here once the randomisations have been run?!

#### Update Data ####
# add Facet net info to relevant files

if(UPDATE_DATA){
  combined_data <- read.table( "/media/cf19810/Seagate Portable Drive/Lasius-Bristol_pathogen_experiment/main_experiment/Soft_community_scores_duration_1ITER_0.5alpha_full.txt", header=TRUE, stringsAsFactors=FALSE)
  combined_data$colony <- paste0("colony",combined_data$colony)
  combined_data$treatment <- paste(combined_data$treatment,combined_data$colony_size,sep=".")
  combined_data$colony_size <- NULL
  combined_data$period <- tolower(gsub("Treatment", "", combined_data$period))
  combined_data$time_hours <- as.numeric(gsub("TH","",combined_data$bin)); combined_data$bin <- NULL
  combined_data$time_of_day <- as.numeric(gsub("TD","",combined_data$TD)); combined_data$TD <- NULL
  combined_data$modularity_FacetNet <- combined_data$MODULARITY; combined_data$MODULARITY <- NULL
  combined_data$randy <- combined_data$ObsRand; combined_data$ObsRand <- NULL
  
  # combined_data_file <- paste0(data_path,"/Soft_community_scores_duration_1ITER_0.5alpha_full.txt")
  # combined_data <- read.table(combined_data_file, header=TRUE, stringsAsFactors=FALSE)
  # combined_data$colony <- paste0("colony",combined_data$colony)
  # combined_data$treatment <- paste(combined_data$treatment,combined_data$colony_size,sep=".")
  # combined_data$colony_size <- NULL
  # combined_data$period <- tolower(gsub("Treatment", "", combined_data$period))
  # combined_data$time_hours <- as.numeric(gsub("TH","",combined_data$bin)); combined_data$bin <- NULL
  # combined_data$time_of_day <- as.numeric(gsub("TD","",combined_data$TD)); combined_data$TD <- NULL
  # combined_data$modularity_FacetNet <- combined_data$MODULARITY; combined_data$MODULARITY <- NULL
  # combined_data$randy <- combined_data$ObsRand; combined_data$ObsRand <- NULL
  
  
  outputfolder <- "/media/cf19810/Seagate Portable Drive/Lasius-Bristol_pathogen_experiment/main_experiment/processed_data/network_properties_edge_weights_duration"
  options <- c("all_workers","untreated_only")
  
  for (input_folder in unique(combined_data$randy)){
    
    if (input_folder=="observed"){
      for (option in options) {
        ####Main experiment: write pre_vs_post_treatment data
        outputfolder2 <- paste(outputfolder,"pre_vs_post_treatment",option,sep="/")
        
        summary_collective <- combined_data[which(combined_data$randy=="observed"),]
        
        #append columns for the new efficiency measures #AW mod for task efficiency
        summary_colony_existing <- read.table(paste(outputfolder2,"/colony_data.txt",sep=""),header = T,stringsAsFactors = F)
        common_names <- names(summary_colony_existing)[names(summary_colony_existing) %in% names(summary_collective)]
        summary_colony_existing <- merge(summary_colony_existing, summary_collective, by = common_names, all.x = TRUE)
        
        write.table(summary_colony_existing[,names(summary_colony_existing)!="randy"],file=paste(outputfolder2,"/colony_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        #write.table(summary_individual[,names(summary_individual)!="randy"],file=paste(outputfolder2,"/individual_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        
        
        ####All workers: write pre_treatment data into random_vs_observed folder
        if (option=="all_workers"){
          outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
          
          #append columns for the new efficiency measures #AW mod for task efficiency
          summary_colony_existing <- read.table(paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),header = T,stringsAsFactors = F)
          common_names <- names(summary_colony_existing)[names(summary_colony_existing) %in% names(summary_collective)]
          summary_colony_existing <- merge(summary_colony_existing, summary_collective, by = common_names, all.x = TRUE)
          
          write.table(summary_colony_existing[which(summary_colony_existing$period=="pre"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          #write.table(summary_individual[which(summary_individual$period=="pre"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          
          outputfolder4 <- paste(outputfolder,"post_treatment",sep="/")
          
          #append columns for the new efficiency measures #AW mod for task efficiency
          summary_colony_existing <- read.table(paste(outputfolder4,"/network_properties_",input_folder,".txt",sep=""),header = T,stringsAsFactors = F)
          common_names <- names(summary_colony_existing)[names(summary_colony_existing) %in% names(summary_collective)]
          summary_colony_existing <- merge(summary_colony_existing, summary_collective, by = common_names, all.x = TRUE)
          
          write.table(summary_colony_existing[which(summary_colony_existing$period=="post"),],file=paste(outputfolder4,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          #write.table(summary_individual[which(summary_individual$period=="post"),],file=paste(outputfolder4,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        }
        
      }
      
      #for all not observed
    }else{
      summary_collective <- combined_data[which(combined_data$randy==input_folder),]
      
      outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
      
      #append columns for the new efficiency measures #AW mod for task efficiency
      summary_colony_existing <- read.table(paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),header = T,stringsAsFactors = F)
      common_names <- names(summary_colony_existing)[names(summary_colony_existing) %in% names(summary_collective)]
      summary_colony_existing <- merge(summary_colony_existing, summary_collective, by = common_names, all.x = TRUE)
      
      write.table(summary_colony_existing[which(summary_colony_existing$period=="pre"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
      #write.table(summary_individual[which(summary_individual$period=="pre"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
    }
  }
}
}  
