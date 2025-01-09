### ### ### ### ### ### ### ### ### ### 
#### 14_summarise_interactions.R #####
### ### ### ### ### ### ### ### ### ### 

#### Read me ####
# Takes an interaction list as an input, and calculates:
#                     - the total duration of interactions of each ant to each of four categories of individuals (queen, nurses, untreated foragers and treated foragers)
#                     - the overall distribution of interactions according to task groups and ages

# Created by Nathalie Stroeymeyt
# Modified by Adriano Wanderlingh to work with FORT formicidae Tracking data. Mods tagged with the comment "AW". script wide mods cited here below.
# Modified by Nathalie Stroeymeyt to include number of events in addition to duration
# with adaptations by Linda Sartoris
# Adjusted to specific needs of the vital experiment by Daniel Schläppi... 

# Script wide mods AW
# - replaced before/after with pre/post
# when ran for the grooming interactions it should only be run for the "observed" folders

# TO DO's
#' Line 50: Check what is needed for trophallactic interactions... potentially the same as main so that we can look at pro and post treatmnet
#' See if and how to implement the trophallactic interaction in this script instead of the grooming interactions which will probably be ignored because of the framerate...
#' Check what to use for task groups... probably Facet Net!!!! update accordingly! 
#' The Script would benefit from a check to see which files have already been processed and which files did not, to only run what is required...

### Why only for PreTreatment at the beginning?

### ### ### ### ### ### ### ### ### ### ### ###
to_keep_ori <- to_keep
### ### ### ### ### ### ### ### ### ### ### ###

options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10); options(readr.show_progress = FALSE)


#### Start ####

### get input file list
input_path           <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists",sep="")
setwd(input_path)  
input_folders        <- list.dirs(recursive=T,path="PreTreatment",full.names=F); input_folders        <- input_folders[which(input_folders!="")]

outputfolder_1 <- paste(data_path,"/processed_data/collective_behaviour/random_vs_observed",sep="")
if(!file.exists(outputfolder_1)){dir.create(outputfolder_1,recursive=T)}

output_file_1  <- paste(outputfolder_1, "processed_interaction_list.txt", sep = "/")
output_file_2  <- paste(outputfolder_1, "interactions.dat"              , sep = "/")
if (file.exists(output_file_1)) {
  cat(red("\n !!! WARNING: !!! \n"))
  cat(yellow(basename(output_file_1)))
  cat(blue(" already exists and newly processed interactions will be appended."))
  flush.console()
  # user_response <- menu(c("Yes", "No"), title = "Do you want to delete the file (please type 1 or 2)?")
  # if (user_response == 1) {
  #   message("The file will be deleted before proceeding.")
  #   file.remove(output_file_1)
  # } else if (user_response == 2) {
  #   message("The file will not be deleted. Processed interactions will be appended.")
  # } else {
  #   message("Invalid response.")
  # }
  repeat{
  user_response <- readline(prompt = " Do you want to delete the file? (y/n): ")
  if(tolower(user_response) == "y") {message("the file will be deleted before proceeding")
                                     file.remove(output_file_1) ; break
                                     }
  else if (tolower(user_response) == "n") {message("The file will not be deleted. Processed interactions will be appended.") ; break}
  else {message("Invalid response. Please enter 'y' or 'n'.")}}
}

summary_dol <- NULL
to_keep <- c(ls(),"to_keep","input_folder","network_file","network_files","summary_interactions","summary_interactions_grooming","summary_pairs","all_interactions", "pb")  

cat(green("\nSummarising interactions \n"))

for (input_folder in input_folders){ # input_folder <- "random"
  cat(blue("\n input folder = ",input_folder, "\n"))
  setwd(input_path)
  network_files <- list.files(path=paste("PreTreatment/",input_folder,sep=""),full.names=T) # get all PreTreatment files into the list
  if (input_folder=="observed"&grepl("main",data_path)){                                    # if it is observed and main exp. we also get the PostTreatment files 
    network_files <- c(network_files,list.files(path=paste("PostTreatment/",input_folder,sep=""),full.names=T))
  }
  summary_interactions <- NULL
  summary_interactions_grooming <- NULL
  summary_interactions_trophy <- NULL # WILL BE NEEDED ONCE WE IMPLEMENT TROPHALLACTIC INTERACTIONS
  summary_pairs        <- NULL
  # all_interactions     <- NULL
  
  pb <- progress_bar$new(
    format = "Progress: :current/:total [:bar] :percent ETA: :eta",
    total = length(network_files),  # Use the length of network_files for total
    clear = FALSE,
    width = 80
  )
  
  for (network_file in network_files){ # network_file <- network_files[2]
    
    #### collect information on network file and get interactions ####
    # cat(yellow("\r",network_file, )) replaced by progress bar below

    ### get file metadata
    root_name          <- gsub("_interactions.txt","",unlist(strsplit(network_file,split="/"))[grepl("interactions",unlist(strsplit(network_file,split="/")))]) 
    components         <- unlist(strsplit(root_name,split="_"))
    colony             <- unlist(strsplit(root_name,split="_"))[1] 
    treatment          <- info[which(info$colony==colony),"treatment"]
    colony_size        <- info[which(info$colony==colony),"colony_size"]
    period             <- ifelse(any(grepl("PreTreatment", components)), "pre", "post") # if (!all(!grepl("PreTreatment",components))){period <- "pre"}else{period <- "post"}
    time_hours         <- as.numeric(gsub("TH","",components[which(grepl("TH",components))]))
    time_of_day        <- as.numeric(gsub("TD","",components[which(grepl("TD",components))]))
    # randy_nr           <- ifelse(input_folder == "observed", "NA", sub("(\\d{3})\\.txt", "\\1", components[8])) # 
      
    # DS: code specific to Linda and Adriano cleaned out

    ### get appropriate task_group list, treated list and tag
    colony_treated     <- treated[which(treated$colony==colony),"tag"]
    colony_task_group  <- task_groups[which(task_groups$colony==colony),]
    queenid            <- as.character(colony_task_group[which(colony_task_group$task_group_prop =="queen"),"tag"]) # has to be a character to work with igraph
    tag <- read.tag(tag_list, colony) # DS deleted some old stuff commented out.
    tag[which(tag$age==0),"age"]   <- NA # unknown ages are coded as 0 in the tag file - replaced with NA instead
    
    ### read interactions 
    interactions       <- read.table(network_file,header=T,stringsAsFactors = F)
    foragers <- colony_task_group %>% filter(task_group_FACETNET_0.5 == "forager") %>% pull(tag)
    nurses <- colony_task_group %>% filter(task_group_FACETNET_0.5 == "nurse") %>% pull(tag)
    
    interactions[c("status_Tag1","status_Tag2")] <- NA
  
    interactions <- interactions %>%
      filter(Tag1 %in% tag$tag & Tag2 %in% tag$tag) %>%  # Filter interactions based on valid tags
      mutate(duration_min = (Stoptime - Starttime + (1/FRAME_RATE)) / 60, N = 1) %>% # Calculate duration in minutes and add N column
      mutate(status_Tag1 = case_when(
             Tag1 %in% foragers ~ "forager",
             Tag1 %in% nurses   ~ "nurse",
             Tag1 == queenid    ~ "queen",
             TRUE               ~ status_Tag1),
           status_Tag2 = case_when(
             Tag2 %in% foragers ~ "forager",
             Tag2 %in% nurses   ~ "nurse",
             Tag2 == queenid    ~ "queen",
             TRUE               ~ status_Tag2)) %>% as.data.frame()
    
    # actor and receiver definition for directional behaviors such as allogrooming was removed here
    if(grepl("grooming",input_path)){ cat(red("WARNING: \n get Adriano's code to define actor and receiver and insert it here"))}
    
    if (period=="pre"){
      interactions_data <- interactions %>%
        dplyr::select(sort(names(.))) %>%
        mutate(randy = input_folder, file = network_file, colony_size = colony_size)
      
      if (!file.exists(output_file_1)) {
        invisible(write_csv(interactions_data, output_file_1))
        } else {
            invisible(write_csv(interactions_data, output_file_1, append = TRUE))
        }
      } 
    
    #### 2. continue calculations for pre vs post ####
      interactions <- interactions %>% mutate(
        status_Tag1 = ifelse(Tag1 %in% colony_treated, "treated", status_Tag1),
        status_Tag2 = ifelse(Tag2 %in% colony_treated, "treated", status_Tag2))
    
    # use this information to calculate, for each worker, the accumulated duration of interaction with treated workers
    aggregated <- interactions %>% 
      group_by(Tag1, status_Tag2) %>%
      summarise(duration_min = sum(duration_min, na.rm = TRUE), number_contacts = sum(N, na.rm = TRUE), .groups = "drop") %>%
      rename(tag = Tag1, partner_status = status_Tag2) %>%
      bind_rows(interactions %>%group_by(Tag2, status_Tag1) %>%
          summarise(duration_min = sum(duration_min, na.rm = TRUE), number_contacts = sum(N, na.rm = TRUE), .groups = "drop") %>%
          rename(tag = Tag2, partner_status = status_Tag1)) %>%
      group_by(tag, partner_status) %>%
      summarise(duration_min = sum(duration_min, na.rm = TRUE), number_contacts = sum(number_contacts, na.rm = TRUE), .groups = "drop") %>% 
      arrange(partner_status) %>% 
      as.data.frame()
    
    full_table <- tag %>% # filter(final_status == "alive") %>% filter not required as all are alive per definition - see tables to match stroeymeyt where the file got created
      dplyr::select(tag) %>% expand_grid(partner_status = c("queen", "nurse", "forager", "treated")) %>%
      left_join(aggregated %>% # join with aggregated interactions
                  dplyr::select(tag, partner_status, duration_min, number_contacts), 
                by = c("tag", "partner_status")) %>%
      mutate(  # replacing of  NA values with 0 in duration_min and number_contacts
        duration_min = if_else(is.na(duration_min), 0, duration_min),
        number_contacts = if_else(is.na(number_contacts), 0, number_contacts)
      ) %>%
      left_join(tag %>% dplyr::select(tag, group_facetNet) %>%  # add status (caste/task) - information stored in tag file
                  rename(status = group_facetNet), by = "tag") %>%
      left_join(colony_task_group %>% dplyr::select(tag, task_group_FACETNET_0.5) %>% # add task group from colony_task_groop or the tag file
                  rename(task_group = task_group_FACETNET_0.5), by = "tag") %>%
      mutate(task_group = if_else(tag %in% colony_treated, "treated", task_group)) %>% # for treated ants we overwrite task_group nurse or forager with treated
      dplyr::select(tag, task_group, status, partner_status, duration_min, number_contacts) %>% arrange(tag,partner_status) %>% as.data.frame() # Select and reorder columns # DS: age not needed
    
    summary_interactions <- bind_rows(summary_interactions, data.frame(randy = input_folder,colony = colony,colony_size = colony_size,treatment = treatment,period = period,time_hours = time_hours,time_of_day = time_of_day,full_table,stringsAsFactors = FALSE))
    
    # HERE, A MASSIVE SECION ON GROOMIN HAS BEEN CLEANED OUT FOR NOW 
    if (grepl("grooming|trophallaxis", input_path)) { cat(red("WARNING: \n See Adriano's code and reinsert some of it here so it works with trophallacis and/or grooming"))}
  
    clean()
    pb$tick()
  }
  
  #### if folder = observed, use the summary_interactions table to compute inter-caste contacts ####
  
  if (grepl("main",data_path)&input_folder=="observed"){ 
    # filter so that interactions with queens and treated ants are no longer in the dataframe then calculate duration and number of inter caste interactions for each ant
    summary_interactions <- summary_interactions %>%
      filter(!partner_status %in% c("treated", "queen"),
             !task_group %in% c("treated", "queen")) %>%
      mutate(within_vs_between = partner_status == task_group) %>%
      group_by(colony, period, time_hours, time_of_day, tag, task_group, status, partner_status, within_vs_between) %>%
      summarise(
        duration_min = sum(duration_min, na.rm = TRUE),
        number_contacts = sum(number_contacts, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(!within_vs_between) %>%
      dplyr::select(-partner_status, -within_vs_between) %>%
      rename(
        inter_caste_contact_duration = duration_min,
        inter_caste_contact_number = number_contacts
      ) %>% arrange(desc(status)) %>%  as.data.frame()
    
    # add information to individual behaviour file 
    behav <- read.table(paste(data_path, "/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt", sep=""), header = TRUE, stringsAsFactors = FALSE) %>%
      dplyr::select(-c(inter_caste_contact_duration, inter_caste_contact_number, 
                       grep("duration_grooming", names(.), value = TRUE))) %>%
      left_join(summary_interactions %>%
                  dplyr::select(colony, tag, time_hours, inter_caste_contact_duration, inter_caste_contact_number),
                by = c("colony", "tag", "time_hours")) %>% mutate(time_hours = as.integer(time_hours)) %>% arrange(desc(period)) %>% as.data.frame()
    
    
    if (grepl("trophallaxis",data_path)&input_folder=="observed"){
      cat(red(" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n THIS SECTION NEEDS AN UPDATE TO RUN CORRECTLY WITH THE TROPHALLACTIC INTERACTIONS"))} 
    if(grepl("grooming",input_path)){
      behav <- merge(behav,summary_interactions_grooming[c("colony","tag","time_hours",names(summary_interactions_grooming)[which(grepl("duration",names(summary_interactions_grooming)))])],all.x=T,all.y=F,sort=F)
      behav <- merge(behav,summary_interactions_grooming[c("colony","tag","time_hours",names(summary_interactions_grooming)[which(grepl("number",names(summary_interactions_grooming)))])],all.x=T,all.y=F,sort=F)}
    
    behav <- behav[order(behav$colony,behav$tag,behav$time_hours),]
    write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
  }
  

  cat(blue("\n Final Step: Processing the massive file containing all interactions combined... \n"))
  ### Use all_interactions to obtain information about age-based DoL pre treatment
  all_interactions <- fread(output_file_1)
  all_interactions <- all_interactions[period == "pre"] # filter to only include pre-treatment
  ### sum interactions for each pair of ants using data.table's aggregation
  all_interactions <- all_interactions[, .(duration_min = sum(duration_min, na.rm = TRUE),
                                           N = sum(N, na.rm = TRUE)),
                                       by = .(randy, colony, colony_size, period, Tag1, Tag2, 
                                              status_Tag1, status_Tag2, treatment)]
  
  ### Calculate intra_caste_over_inter_caste_WW_contact_duration
  all_interactions[, same_caste := status_Tag1 == status_Tag2]
  # the above file can be massive and might exceed the memory of your computer.. on linux  you can use swap space so that in case you run out of RAM it will
  # use swap space from the disk to have additional memory... much slower, but it works. On linux you can manually set up how much swap space is available --> consider increasing it
  
  same_caste_interactions <- all_interactions[status_Tag1 != "queen" & status_Tag2 != "queen" & same_caste == TRUE,
                                              .(duration_min_within = sum(duration_min, na.rm = TRUE),
                                                number_contacts_within = sum(N, na.rm = TRUE)),
                                              by = .(randy, colony, colony_size, treatment, period)]
  inter_caste_interactions <- all_interactions[status_Tag1 != "queen" & status_Tag2 != "queen" & same_caste == FALSE,
                                               .(duration_min_between = sum(duration_min, na.rm = TRUE),
                                                 number_contacts_between = sum(N, na.rm = TRUE)),
                                               by = .(randy, colony, colony_size, treatment, period)]
  inter_intra_caste_interactions <- merge(same_caste_interactions, inter_caste_interactions, by = c("randy", "colony", "colony_size", "treatment", "period"), all = TRUE) # merge same caste and inter caste interactions
  inter_intra_caste_interactions[, intra_caste_over_inter_caste_WW_contact_duration := duration_min_within / duration_min_between] # calculate ratios
  inter_intra_caste_interactions[, intra_caste_over_inter_caste_WW_contact_number := number_contacts_within / number_contacts_between]
  dol <- inter_intra_caste_interactions[, !c("duration_min_within", "duration_min_between", "number_contacts_within", "number_contacts_between")]
  
  ### calculate queen contact with nurses vs. workers (for regular interactions only but not grooming or traphallaxis)
  if (!grepl("grooming|trophallaxis",input_path)) {
    queen_interactions <- all_interactions[status_Tag1 == "queen" | status_Tag2 == "queen"] # Subset queen interactions
    queen_interactions[status_Tag1 == "queen", `:=`(partner = Tag2, partner_status = status_Tag2)] # Assign 'partner' and 'partner_status' columns
    queen_interactions[status_Tag2 == "queen", `:=`(partner = Tag1, partner_status = status_Tag1)]
    # aggregate interactions with nurses and foragers and merge
    interaction_with_nurses <- queen_interactions[partner_status == "nurse", 
                                                  .(duration_min_with_nurses = sum(duration_min, na.rm = TRUE), 
                                                    number_contacts_with_nurses = sum(N, na.rm = TRUE)), 
                                                  by = .(randy, colony, colony_size, period, treatment)]
    interaction_with_forager <- queen_interactions[partner_status == "forager", 
                                                   .(duration_min_with_foragers = sum(duration_min, na.rm = TRUE), 
                                                     number_contacts_with_foragers = sum(N, na.rm = TRUE)), 
                                                   by = .(randy, colony, colony_size, period, treatment)]
    queen_interac <- merge(interaction_with_nurses, interaction_with_forager, by = c("randy", "colony", "colony_size", "period", "treatment"), all = TRUE)
    # calculate ratios and merge with dol
    queen_interac[, `:=`(
      QNurse_over_QForager_contact_duration = duration_min_with_nurses / duration_min_with_foragers,
      QNurse_over_QForager_contact_number = number_contacts_with_nurses / number_contacts_with_foragers)]
    dol <- merge(dol, queen_interac[, .(randy, colony, period, QNurse_over_QForager_contact_duration, QNurse_over_QForager_contact_number)], 
                 by = c("randy", "colony", "period"), all = TRUE)}
  
  # DS: Here, a large chunk of code from AW's age experiment was cut out.
  summary_dol <- rbind(summary_dol,dol)  
}

if(!file.exists(output_file_2)){write.table(summary_dol,file=output_file_2,col.names=T,row.names=F,append=F,quote=F)}
to_keep <- to_keep_ori


#### TO DO
# See in the subsequent scripts what is needed from here and what is not needed. 
# e.g. if age based division of is probably not needed so it could be commented out and in that case we probably also do not need the all_ineractions file which is massive 
# and in this case the whole loop might not need to run over all the randomised interactions??? 



