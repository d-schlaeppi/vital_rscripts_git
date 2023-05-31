rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1. ANT ORIENT EXPRESS ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###



#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

#### prerequisites ####

#load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp) # contains sourceCpp
library(circular) # contains circular function

# directory of data and myrmidon files
directory <- "/media/gw20248/gismo_hd2/vital/fc2/"

# copy the "Get_Movement_Angle.cpp" file from sharepoint into your directory mentioned above
# source C++ movement direction program
sourceCpp(paste0(directory,"Get_Movement_Angle.cpp"))


# create a list containing a source file: manually oriented myrmidon file with predefined capsules
source_data_list <- list(paste0(directory,"vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsAutoDefined.myrmidon"))


#### Step 1 - Extract information on capsules from source file ####

{#run this whole part below here to extract all the information required from the capsule source file. 
oriented_metadata <- NULL
capsule_list <- list()
for (myrmidon_file in source_data_list){
  experiment_name <- unlist(strsplit(myrmidon_file,split="/"))[length(unlist(strsplit(myrmidon_file,split="/")))]
  oriented_data <- fmExperimentOpen(myrmidon_file)
  oriented_ants <- oriented_data$ants
  capsule_names <- oriented_data$antShapeTypeNames
  for (ant in oriented_ants){
    # extract ant length and capsules
    ant_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px)
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
}
# Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen
# So we need to remove the queen from the computation by removing outliers in ant lenght measurements
interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75))
outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  # apply outlier exclusion to oriented_metadata and to capsule list
for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
}
# Once queen(s) has(have) been removed, get the mean coordinates of the offset between tag centre and ant centre
mean_x_ant_coord <- mean(oriented_metadata$x_ant_coord)
mean_y_ant_coord <- 0 ##set it to zero manually because we don't expect there to be a consistent bias in deviation / #mean_y_ant_coord <- mean(oriented_metadata$y_ant_coord) ###this is expected to be 0 or near 0 (check!) because the tag should be as likely to be to the right or to the left of the ant's bilateral symmetry line
mean_worker_length_px <- mean(oriented_metadata$length_px) # Get the average worker length from the data
for (caps in 1:length(capsule_list)){ # Finally, get information on each capsule
  capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
}
}


#### Step 2 - Copy the capsule information onto a list of files that will get oriented and the capsules assigned #### 

# create a list of files that contains all the myrmidon files to be modified
files <- list.files(directory)
files <- grep(files, pattern = 'c29', invert = FALSE, value = TRUE)
files <- grep(files, pattern = 'done', invert = FALSE, value = TRUE)

# ADJUST OUTPUT FILENAME IN LOOP ACCORDING TO YOUR OWN FILENAMING

{
not_oriented <- NULL
to_orient_manually <- NULL
for (file in files) {
  tracking_data <- fmExperimentOpen(paste0(directory, file))
  for (caps in 1:length(capsule_list)){ #add the caps from the source file above
    tracking_data$createAntShapeType(names(capsule_list)[caps])
  }
  ants <- tracking_data$ants
  # get trajectory data to extract ant orientation:
  # short tracking using all data -> using start and end time from the experiment metadata
  #from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start) # experiment start time # from <- fmTimeSinceEver()
  #to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$end  ) # experiment end time   # to   <- fmTimeForever()
  from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$end -12*3600) # longer tracking data - computation takes very long -> only use a subset of the data e.g. last 12 hours
  to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$end)
  max_gap <- fmHour(24*365)  # use a  large value to make sure you get only one trajectory per ant (larger than the time difference between from to defined above)
  positions <- fmQueryComputeAntTrajectories(tracking_data, start = from, end = to, maximumGap = max_gap, computeZones = TRUE)
  #Then compute trajectories & hard-wire ant correspondence between trajectories_summary and trajectories
  positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") # creates a ID string for each ant: ant1, ant2,...
  names(positions$trajectories)       <- positions$trajectories_summary$antID_str # and use the content of that column to rename the objects within trajectory list
  max_time_gap <- 0.5 # define a max temporal gap for which you are happy to calculate a movement angle; e.g. 0.5 s
  min_dist_moved <- 30 # define a minimum distance moved, as you don't want to use noise or small shifts in position in this calculation; e.g. 30 pix (to think about)
  for (i in 1:length(ants)){
    if(length(tracking_data$ants[[i]]$identifications) == 0) {
      print(paste(file, i, "no ant", sep = " -> "))
      not_oriented <- append(not_oriented, paste(file, i, "no ant", sep = " -> "))
      next
    }
    #skip ants that got retagged and create a table (to be saved in a file) with all the skipped ants that need manual orientation
    if (ants[[i]]$getValue("tag_reoriented", fmTimeForever())) {
      to_orient_manually <- append(to_orient_manually, paste(file, i, ants[[i]]$ID, "retagged", sep = " -> "))
      print(paste(file, i, "skipped because of retagging", sep = " -> "))
      next
    }
    if (tracking_data$ants[[i]]$identifications[[1]]$tagValue==0) {next} # skip the queen
    if (is.element(i,positions$trajectories_summary$antID)) {
      # to be fool proof, and be sure you extract the trajectory corresponding the correct ant, make sure you make use of the antID_str column!
      traj <- positions$trajectories [[   positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"antID_str"]    ]]
      traj <- cbind(traj,add_angles(traj,max_time_gap,min_dist_moved))   # feed traj to c++ program
      AntAngle <- as.numeric(- mean(circular(na.omit(traj$Tag_minus_Movement_Angle),units="radians",zero=0)))   #  get mean deviation angle between body and tag - the ant angle is equal to minus the Tag minus Movement angle output by C++ program
      x_tag_coord <- mean_x_ant_coord*cos(AntAngle) - mean_y_ant_coord*sin(AntAngle)  # now use trigonometry to calculate the pose, using AntAngle
      y_tag_coord <- mean_x_ant_coord*sin(AntAngle) + mean_y_ant_coord*cos(AntAngle)
      for (id in ants[[i]]$identifications){  # write this into ant metadata
        id$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)
      }
      # also add this to trajectories_summary
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"ant_angle"] <- AntAngle
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"x_tag_coord"] <- x_tag_coord
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"y_tag_coord"] <- y_tag_coord
      # finally, for each ant, add capsules using mean_ant_length and capsule_list
      for (caps in 1:length(capsule_list)){
        capsule_ratios <- capsule_list[[caps]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
        capsule_coords <- mean_worker_length_px*capsule_ratios
        ants[[i]]$addCapsule(caps, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
      }
    } else {
      print(paste(file, i, "FALSE", "exit before selected period", sep = " -> "))
      not_oriented <- append(not_oriented, paste(file, i, "FALSE", "exit before selected period", sep = " -> "))
    }
  }
  # tracking_data$save(paste0(directory, substr(file, 1, nchar(file)-24),'oriented.myrmidon')) ### adjust the name of the output file according to your input file names!!!
  tracking_data$save(paste0(directory, "oriented", substr(file, 12, nchar(file))))             ### adjust the name of the output file according to your input file names!!!
  not_oriented <- append(not_oriented, paste("time", Sys.time() ,sep = " : "))
}
print(not_oriented)
print(to_orient_manually)
}


# optionally run the code below to save the information on which ants to orient manually into files
# fwrite(list(not_oriented), file = paste(Sys.Date(), "",format(Sys.time(), "%H-%M-%S"), "not_oriented.txt", sep = "_"))
# fwrite(list(to_orient_manually), file = paste(Sys.Date(), "",format(Sys.time(), "%H-%M-%S"), "to_orient_manually.txt", sep = "_"))



# Go to fort-studio to check things look all right
# Next - modify queen information - 
