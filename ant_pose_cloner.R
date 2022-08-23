rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1. ANT-POSE-CLONER  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####
# If ants are used in two separate tracking systems it is useful if the ant-pose, size and orientation are identical for both tracking files
# (e.g. tracking of a colony in one tracking system and then a few workers get sampled for treatment that is also tracked and then the workers get returned to their colony)
# especially if the treatment tracking is only short and there are no pictures saved yet that can be used for manual orientation
# (when doing short tracking make sure to adjust the frequency at which pictures are taken - e.g every 10-15 minutes or so for a 2h tracking so you have enough to manually orient ants)
# This script will copy the ant orientation/pose and size from the main tracking and apply it the second tracking file with a scaling factor based on tag size (as not both tracking systems are set up identically)


#### prerequisites ####

# load libraries
library(FortMyrmidon) #R bindings
#library(Rcpp)
#library(circular)
#library(R.utils)

# set working directory
directory <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(directory)



# once we have the post-processed date rewrite the script below to do all colonies in one go: 

for (i in 1:nrow(data_collection)) {
  main_file_name  <- paste0("colony_nr and rest of the main file name", "with colony nr set as i")
  secondary_file_name <- paste0("colony_nr and rest of the corresponding feeding file name", "with colony nr set as i")
  #then include here the code below 
}
  

# select the two related files you want to work with (both need to have ants already created) 
main_file_name <- "vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined_metaID.myrmidon" # main tracking file
#secondary_file_name <- "vital_fc2_esterhase_c02_feeding_DS_AntsCreated_metaID.myrmidon" # treatment tracking file
secondary_file_name <- "vital_fc2_esterhase_c02_feeding_DS_AntsCreated.myrmidon" # treatment tracking file


# file containing the oriented ants with capsules
source_file <- paste0(directory, main_file_name)
destination_files <- list(
  paste0(directory, secondary_file_name)
)

# files containing the ants that need automatic orienting and capsules based on the source file

#### Extraction of information on AntPose and Capsules from the source file ####

data_list <- list(source_file)
source_metadata <- NULL
capsule_list<- list()
ant_lenght_px <- NULL

# define name of output file
output_name <- file.path(paste0(directory, "Mean_ant_length_colonies.txt")) 

for (file in data_list) {
  experiment_name <- unlist(strsplit(file,split="/"))[length(unlist(strsplit(file,split="/")))]
  source_data <- fmExperimentOpen(file)
  source_ants <- source_data$ants
  capsule_names <- source_data$antShapeTypeNames
  for (ant in source_ants) {
    # extract ant length and capsuöes
    ant_length_px <- mean(fmQueryComputeMeasurementFor(source_data,antID=ant$ID)$length_px)
    capsules <- ant$capsules
    for (caps in 1:length(capsules)){
      capsule_name  <- capsule_names[[capsules[[caps]]$type]]
      capsule_coord <- capsules[[caps]]$capsule
      capsule_info <- data.frame(experiment = experiment_name,
                                 antID      = ant$ID,
                                 c1_ratio_x = capsule_coord$c1[1]/ant_length_px,
                                 c1_ratio_y = capsule_coord$c1[2]/ant_length_px,
                                 c2_ratio_x = capsule_coord$c2[1]/ant_length_px,
                                 c2_ratio_y = capsule_coord$c2[2]/ant_length_px,
                                 r1_ratio   = capsule_coord$r1[1]/ant_length_px,
                                 r2_ratio   = capsule_coord$r2[1]/ant_length_px)
      if (!capsule_name %in%names(capsule_list)){ ###if this is the first time we encounter this capsule, add it to capsule list...
        capsule_list <- c(capsule_list,list(capsule_info)) 
        if(length(names(capsule_list))==0){
          names(capsule_list) <- capsule_name
        }else{
          names(capsule_list)[length(capsule_list)] <- capsule_name
        }
      }else{###otherwise, add a line to the existing dataframe within capsule_list
        capsule_list[[capsule_name]] <- rbind(capsule_list[[capsule_name]] , capsule_info)
      }
    }
    for (id in ant$identifications){
      source_metadata <- rbind(source_metadata, data.frame(experiment       = experiment_name,
                                                              antID            = ant$ID,
                                                              tagIDdecimal     = id$tagValue,
                                                              angle            = id$antAngle,
                                                              x_tag_coord      = id$antPosition[1], 
                                                              y_tag_coord      = id$antPosition[2],
                                                              x_ant_coord      = id$antPosition[1]*cos(-id$antAngle) - id$antPosition[2]*sin(-id$antAngle),
                                                              y_ant_coord      = id$antPosition[1]*sin(-id$antAngle) + id$antPosition[2]*cos(-id$antAngle),
                                                              length_px        = ant_length_px,
                                                              stringsAsFactors = F))
    }
  }
 }

# remove the queen
interquartile_range <- quantile(source_metadata$length_px,probs=c(0.25,0.75), na.rm =TRUE)
outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
### apply outlier exclusion to source_metadata
source_metadata <- source_metadata[which(source_metadata$length_px>=outlier_bounds[1]&source_metadata$length_px<=outlier_bounds[2]),]  
### apply outlier exclusion to capsule list
for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(source_metadata $experiment, source_metadata $antID)),]
}
for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
}

mean_ant_size_without_queen_source_file <- mean(source_metadata$length_px)

ANT.LENGTH <- NULL





#### Write capsule data to each manually oriented file ####
for (destination_file in destination_files) {
  # open tracking data which need new capsule
  tracking_data <- fmExperimentOpen(destination_file) 
  source_data <- fmExperimentOpen(source_file)
  ants <- tracking_data$ants
  for (ant in ants){
    ### extract ant length and capsules
    ant_length_px <- mean_ant_size_without_queen_source_file
    ANT.LENGTH <- rbind(ANT.LENGTH,data.frame(
      length_px        = ant_length_px,
      stringsAsFactors = F))
  }
  #create dataset with mean values x colony
  ant_length_colony <- data.frame(ant.length=mean(ANT.LENGTH$length_px,na.rm=T), colony=destination_file)    
  for (caps in 1:length(capsule_list)){
    tracking_data$createAntShapeType(names(capsule_list)[caps])
  }
  for (i in 1:length(ants)){
    #use mean size of each manually oriented file that needs the capsule
    ant_length_px <- mean_ant_size_without_queen_source_file
    ants[[i]]$clearCapsules()
    ##assign capule numbers that match the order of the looped capsule names positions
    capsule_number <- 0
    for (capsule_name in unlist(tracking_data$antShapeTypeNames)) {
      capsule_number <- capsule_number +1
      # the file information
      #MAKE SURE THERE IS CAPSULE MATCHING, TO AVOID MIXING UP SHAPE INDEXES
      capsule_ratios <- capsule_list[[capsule_name]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
      capsule_coords <- ant_length_px*capsule_ratios
      
      ants[[i]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
      
    }#IF THIS RESULTS IN A ERROR, CHECK PREVIOUS VERSION ON GIT OR COMPARE WITH DataPrep4_Clone-capule-queens-only_v082.R
  }
  
  #tracking_data$save(destination_file) 
  tracking_data$save(paste0(sub("\\..*", "", destination_file),"_CapsAutoDefined_NEW.myrmidon"))
  
  ## save
  if (file.exists(output_name)){
    write.table(ant_length_colony,file=output_name,append=T,col.names=F,row.names=F,quote=T,sep=",")
  }else{
    write.table(ant_length_colony,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
  }
  
  #close experiment
  rm(list=(c("tracking_data")))
  
}











#### ant orientation cloner ####

for each ant in the feeding file 
check for the corresponding ant in the source file and assign its values to this ant. 























