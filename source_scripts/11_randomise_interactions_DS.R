
### ### ### ### ### ### ### ### ### ### ### ###
#### 11_randomise_interactions.R #####
### ### ### ### ### ### ### ### ### ### ### ###


### Sources C++ function randomised_edges 

#### Takes an interaction list as an input and returns a the same interaction interaction list in which the interaction partners have been randomised
#### Follows the 'Randomized edges (RE)' algorithm laid out by Holme and Saramäki 2012 (Physics Reports)

### Created by Nathalie Stroeymeyt 
### Adjusted for D.Schlaeppi

### ### ### ### ### ### ### ### ### ### ### ###
to_keep_ori <- to_keep

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(code_path,"/randomise_edges.cpp",sep=""))



# Get all interaction lists for which randomized versions should be created:
# in theory only needed for the pre treatment but should be fine with the post treatment as well... 
timings <- c("PreTreatment", "PostTreatment")# Specify timings
interact_list <- list()

for (timing in timings) {
  subdirs <- c("binned_interaction_lists", "full_interaction_lists")
  for (subdir in subdirs) {
    observed_path <- paste(data_path, "intermediary_analysis_steps", subdir, timing, "observed", sep = "/")
    files <- list.files(observed_path, full.names = TRUE)
    interact_list <- append(interact_list, files)
  }
}

to_keep <- c(ls(),"to_keep","i","interac")

for (i in 1:100){ ### perform 100 randomization                  
  print(paste("Performing randomisations",i,"out of 100..."))
  for (interac in interact_list){
    # interac <- interact_list[[1]]
    folder_name <- dirname(interac) 
    file_name <- basename(interac)
    outputfolder       <- gsub("observed", paste("random",sep=""),folder_name)
    #
    if (!file.exists(outputfolder)){dir.create(outputfolder)}
    new_file_name <- paste0(sub("\\.txt$", "", file_name),"_random_",paste(rep(0,3-nchar(i)),collapse=""),i, ".txt")
    outfile <- paste(outputfolder, new_file_name,sep="/")
    
    if (!file.exists(outfile)){
      interactions                 <- read.table(interac,header=T,stringsAsFactors=F, sep = "\t")
      randomised_partners          <- randomise_edges(interactions[c("Tag1","Tag2","Starttime","Stoptime")])
      randomised_interactions      <- interactions
      randomised_interactions$Tag1 <- randomised_partners$Tag1;randomised_interactions$Tag2 <- randomised_partners$Tag2;
      write.table(randomised_interactions,file=outfile,col.names=T,row.names=F,quote=F,append=F)
    }
    clean();
  }
}

to_keep <- to_keep_ori

# to do next: 
# Organize the trophallactic interactions:
# load the the original file and create subsets for the pre 3h period and the post 3 h period. 
