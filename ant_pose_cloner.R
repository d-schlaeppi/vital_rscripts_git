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

library(FortMyrmidon) #R bindings

# set working directory
directory <- '/media/gw20248/gismo_hd2/vital/fc2/'  # HD with extrapolated data
setwd(directory)

source_files <- list.files(directory, pattern = "m_queen_modified_c\\d{2}\\.myrmidon$")
tracking_data_type <- c("source", "dest")



for (source_file in source_files) {
    colony_name <- substr(source_file, 18, nchar(source_file)-9)
    source_data <- fmExperimentOpen(source_file)
    source_ants <- source_data$ants
    capsule_names <- source_data$antShapeTypeNames
    dest_data <- fmExperimentOpen(paste0("f_AntsCreated_", substr(source_file, 18, nchar(source_file))))
    dest_ants <- dest_data$ants  
    tag_value_vector <- data.frame(tag_values = numeric(0)) # create vector of all ants from the source file 
    tag_values <- NULL
    for (i in source_ants) {
      tryCatch({
        if (is.numeric(i$identifications[[1]]$tagValue)) {
          tag_values <- i$identifications[[1]]$tagValue
          tag_value_vector <- rbind(tag_value_vector, data.frame(tag_values))
        }}, error = function(e) {
          # ignore the error and continue to the next iteration
          NULL
        })
      }
    # for each ant in the destination file copy the size, orientation and capsule (and optionally all metadata) form the source ant. 
    # here get the information on the capsules from the source file
    for (ant in source_ants){
      if (ant$identifications[[1]]$tagValue==0) {next} # skip the queen
      # extract ant length and capsules
      ant_length_px <- mean(fmQueryComputeMeasurementFor(source_data,antID=ant$ID)$length_px)
      capsules      <- ant$capsules
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
                                   r2_ratio   = capsule_coord$r2[1]/ant_length_px
        )
        if (!capsule_name %in%names(capsule_list)){ # if this is the first time we encounter this capsule, add it to capsule list...
          capsule_list <- c(capsule_list,list(capsule_info)) 
          if(length(names(capsule_list))==0){
            names(capsule_list) <- capsule_name
          }else{
            names(capsule_list)[length(capsule_list)] <- capsule_name
          }
        }else{ #otherwise, add a line to the existing dataframe within capsule_list
          capsule_list[[capsule_name]] <- rbind(capsule_list[[capsule_name]] , capsule_info)
        }
      }
      # extract offset btewen tag centre and ant centre
      for (id in ant$identifications){
        oriented_metadata <- rbind(oriented_metadata,data.frame(experiment       = experiment_name,
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
    mean_x_ant_coord <- mean(oriented_metadata$x_ant_coord)
    mean_y_ant_coord <- 0 ##set it to zero manually because we don't expect there to be a consistent bias in deviation / #mean_y_ant_coord <- mean(oriented_metadata$y_ant_coord) ###this is expected to be 0 or near 0 (check!) because the tag should be as likely to be to the right or to the left of the ant's bilateral symmetry line
    mean_worker_length_px <- mean(oriented_metadata$length_px) # Get the average worker length from the data
    for (caps in 1:length(capsule_list)){ # Finally, get information on each capsule
      capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
    }
    ### now copy orientation and capsule on the destination ants 
    for (caps in 1:length(capsule_list)){ #add the caps from the source file above
      dest_data$createAntShapeType(names(capsule_list)[caps])
    }
    for (ant in dest_ants) {
      if(is.element(ant$identifications[[1]]$tagValue, as.matrix(tag_value_vector))) {
          id$setUserDefinedAntPose(c(source_ants[["the one that has the same tag value"]]$identifications[[1]]$antPosition[1],  # x_tag_coord
                                     source_ants[["the one that has the same tag value"]]$identifications[[1]]$antPosition[2]), # y_tag_coord
                                     source_ants[["the one that has the same tag value"]]$identifications[[1]]$antAngle) # AntAngle
      }
    # finally, for each ant, add capsules using mean_ant_length and capsule_list
      for (caps in 1:length(capsule_list)){
        capsule_ratios <- capsule_list[[caps]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
        capsule_coords <- mean_worker_length_px*capsule_ratios
        ant[[i]]$addCapsule(caps, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
      }
      
      print(paste0(source_file, " treated ant ", ant$identifications[[1]]$tagValue, " got oriented :-) "))
      
    }
    dest_data$save(paste0("f_AntsModified_", substr(source_file, 18, nchar(source_file))))  
}

  
    
source_ants[[2]]$capsules
dest_ants[[2]]$capsules

dest_ants[[2]]$capsules <- source_ants[[2]]$capsules
dest_ants[[2]]$antPosition[1]

dest_ants[[2]]$identifications[[1]]$antPosition




source_ants[[2]]$identifications
fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_px

id$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)

































