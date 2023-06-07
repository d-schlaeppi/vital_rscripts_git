### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### CAPSULE CLONER ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


#### READ ME ####
# This script was written by Daniel Schläppi, based on the capsule_assignment script from Adriano Wanderling and Nathalie Stroeymeyt

# This script contains a function that can be sourced to clone capsules from a source myrmidon file containing predefined capsules into a list of destination_myrmidon file.access
# The function takes as arguments 
# a list of source myrmidon files containting the capsule definitions for workers & queen "capsule_source_files_list" (one file is sufficient)
# a list of destination myrmidon files which will get the new capsules
# a metadata table containing at least the unique identifier specific and mean ant length in pixel
# Queen data orientation and size needs to be done manually after worker auto orientation (see vital main script step)
# This script assumes that you used Queen tags only 0x000 and thus the queen has antID = 1 per definition! If you used other tags for the queen adjust the way the queens are identified (see Adrianos scripts)
# capsule structure should be the same between workers and queens (i.e. both have the same capsules types defined e.g. head, body & antenna) #if not, script adjustment might be needed


# for more informations go to the vital main script 
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"

# this script was written for data for which pre-processing was completed -> files should have Queen infos  overwritten (tag size, manual orientation, manual capsules).


#### prerequisites, i.e. things required to run this function ####

# load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp) # contains sourceCpp
library(circular) # contains circular function
library(R.utils)
library(reader)
library(stringr)


directory_scripts <- "/home/gw20248/Documents/vital_rscripts_git/" # directory with the R scripts linked with github
# set directory to where you have your myrmidon files and assosciated data
directory_data <- '/media/gw20248/gismo_hd2/vital/fc2/'
setwd(directory_data)

# load the meta data dataframe 
source(paste0(directory_scripts,"vital_meta_data.R")) # will add the meta data dataframe to your environment so it can be accessed within this script
head(meta_data)

# list of myrmidon files containing the capsule definitions to be applied/cloned to the destination files (for the capsule definition you need to start with a manually oriented colony)
# save your source myrmidon files like this: filename_CapsuleDefXX_source.myrmidon 
# get the list of all the filenames (with path form directory) containting the capsule definition (listing it manually or getting it automatically from the directory)
add_directory <- function(filename) {# Function to add directory path to each filename
  paste0(directory_data, filename)
}

capsule_source_files_list <- list.files()
capsule_source_files_list <- grep(capsule_source_files_list, pattern = 'source', invert = FALSE, value = TRUE)
capsule_source_files_list <- grep(capsule_source_files_list, pattern = 'CapsuleDef', invert = FALSE, value = TRUE)
add_directory <- function(filename) {# Function to add directory path to each filename
  paste0(directory_data, filename)
}
capsule_source_files_list <- lapply(capsule_source_files_list, add_directory)

# get a list af all the filenames that shall have the capsule definitions (destination files)
# save your destination myrmidon files like this: final_uniqueIdentifier.myrmidon # e.g. final_c01.myrmidon

capsule_destination_files_list <- list.files()
capsule_destination_files_list <- grep(capsule_destination_files_list, pattern = 'final', invert = FALSE, value = TRUE)
capsule_destination_files_list <- grep(capsule_destination_files_list, pattern = 'CapsuleDef', invert = TRUE, value = TRUE)
capsule_destination_files_list <- lapply(capsule_destination_files_list, add_directory)




#### capsule cloner function ####
# function that in step 1 will extract the information about AntPose and Capsules
#               in step 2 will apply this information (extracted capsules) to new myrmidon files (destination files)
###############C CHECK if this should be done as well! # then, get the corresponding queen capsule and apply it as well 
# Step 1:
clone_capsules <- function(capsule_source_files_list, capsule_destination_files_list, meta_data) {
  for (source_file in capsule_source_files_list){
    oriented_metadata <- NULL
    queen_oriented_metadata <- NULL
    capsule_list <- list()
    queen_capsule_list <- list()
    experiment_name <- "vital"
    oriented_data <- fmExperimentOpen(source_file)
    oriented_ants <- oriented_data$ants
    capsule_names <- oriented_data$antShapeTypeNames
    queen_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=1)$length_px) 
    for (ant in oriented_ants){
      # extract ant length and capsules
      ant_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px)
      capsules      <- ant$capsules
      queen_capsules      <- oriented_ants[[1]]$capsules
      for (caps in 1:length(capsules)){ # get capsule information for workers
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
      
      for (caps in 1:length(queen_capsules)){ # same as above but this time we get capsule information for queens
        queen_capsule_name  <- capsule_names[[queen_capsules[[caps]]$type]]
        queen_capsule_coord <- queen_capsules[[caps]]$capsule
        queen_capsule_info <- data.frame(experiment = experiment_name,
                                         antID      = 1,
                                         c1_ratio_x = queen_capsule_coord$c1[1]/queen_length_px,
                                         c1_ratio_y = queen_capsule_coord$c1[2]/queen_length_px,
                                         c2_ratio_x = queen_capsule_coord$c2[1]/queen_length_px,
                                         c2_ratio_y = queen_capsule_coord$c2[2]/queen_length_px,
                                         r1_ratio   = queen_capsule_coord$r1[1]/queen_length_px,
                                         r2_ratio   = queen_capsule_coord$r2[1]/queen_length_px
        )
        if (!queen_capsule_name %in%names(queen_capsule_list)){ # if this is the first time we encounter this capsule, add it to capsule list...
          queen_capsule_list <- c(queen_capsule_list,list(queen_capsule_info)) 
          if(length(names(queen_capsule_list))==0){
            names(queen_capsule_list) <- queen_capsule_name
          }else{
            names(queen_capsule_list)[length(queen_capsule_list)] <- queen_capsule_name
          }
        }else{ #otherwise, add a line to the existing dataframe within queen_capsule_list
          queen_capsule_list[[queen_capsule_name]] <- rbind(queen_capsule_list[[queen_capsule_name]] , queen_capsule_info)
        }
      }
      # extract offset btewen tag centre and ant centre ibce for workers and once for the queen
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
      for (id in oriented_ants[[1]]$identifications){
        queen_oriented_metadata <- rbind(queen_oriented_metadata,data.frame(experiment       = experiment_name,
                                                                            antID            = 1,
                                                                            tagIDdecimal     = id$tagValue,
                                                                            angle            = id$antAngle,
                                                                            x_tag_coord      = id$antPosition[1], 
                                                                            y_tag_coord      = id$antPosition[2],
                                                                            x_ant_coord      = id$antPosition[1]*cos(-id$antAngle) - id$antPosition[2]*sin(-id$antAngle),
                                                                            y_ant_coord      = id$antPosition[1]*sin(-id$antAngle) + id$antPosition[2]*cos(-id$antAngle),
                                                                            length_px        = queen_length_px,
                                                                            stringsAsFactors = F))
      }
    }
  }
  # Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen -> queen removed by removing outliers in ant lenght measurements
  # outlier exclusion not required for queen capsules!
  interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75))
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  # apply outlier exclusion to oriented_metadata and to capsule list
  for (caps in 1:length(capsule_list)){
    capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
  }
  for (caps in 1:length(capsule_list)){ # Finally, get information on each capsule (workers)
    capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
  }
  for (caps in 1:length(queen_capsule_list)){ # Finally, get information on each capsule (queen)
    queen_capsule_list[[caps]] <- colMeans(queen_capsule_list[[caps]][,which(grepl("ratio",names(queen_capsule_list[[caps]])))])
  }
  
  # Step 2 apply capsules to destination files
  for (destination_file in capsule_destination_files_list) {
    tracking_data <-fmExperimentOpen(destination_file)
    identifier <- sub(".*_(c\\d{2}).*", "\\1", destination_file) # extract the colony identifier from the destination file
    mean_worker_length_px <- meta_data$mean_ant_lenght_px[meta_data$colony_id == identifier]
    mean_queen_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=1)$length_px)
    ants <- tracking_data$ants
    for (i in 1:length(ants)) { # skip potential errors due to missing ants (i.e. a lost tag that was miss-identified as an ant and then removed during manual post processing results in a missing ant-id)
      if(length(tracking_data$ants[[i]]$identifications) == 0) {
        print(paste(destination_file, i, "no ant", sep = " -> "))
        next
      }
      ants[[i]]$clearCapsules()  # delete/clear individuals' capsule data and capsule shapes in general (if present)
    }
    if (length(tracking_data$antShapeTypeNames)>0) {
      for (caps in 1:length(tracking_data$antShapeTypeNames)){
        tracking_data$deleteAntShapeType(caps)
      }
    }
    for (caps in 1:length(capsule_list)){ # recreate the ant shapes in the destination file by calling and applying the new names from the source files (will be the same for queens and workers)
      tracking_data$createAntShapeType(names(capsule_list)[caps])
    }
    for (i in 1:length(ants)) { # skip potential errors due to missing ants (i.e. a lost tag that was miss-identified as an ant and then removed during manual post processing results in a missing ant-id)
      if(length(tracking_data$ants[[i]]$identifications) == 0) {
        next
      }
      if (tracking_data$ants[[i]]$identifications[[1]]$tagValue==0) { # Identify the queen apply the queen capsules and define tag size 
        try(ants[[i]]$identifications[[1]]$tagSize <- 1.52, silent = TRUE) # adjust tag size for the queen spits out an error saying the property is read only but it seems to work - thus it is made silent so the loop does not interrupt. 
        capsule_number <- 0
        for (capsule_name in unlist(tracking_data$antShapeTypeNames)) {
          capsule_number <- capsule_number +1
          capsule_ratios <- queen_capsule_list[[capsule_number]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
          capsule_coords <- mean_queen_length_px*capsule_ratios
          ants[[i]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
        }
      } else {
        capsule_number <- 0
        for (capsule_name in unlist(tracking_data$antShapeTypeNames)) { # apply the workers capsules
          capsule_number <- capsule_number +1
          capsule_ratios <- capsule_list[[capsule_number]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
          capsule_coords <- mean_worker_length_px*capsule_ratios
          ants[[i]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
        }
      }
    }
  }
  # finally, save each myrmidon file before moving on to the next set of capsule definitions.
  tracking_data$save(paste0(substr(destination_file, 0  , nchar(destination_file)-9), substr(source_file, nchar(source_file)-28, nchar(source_file)-16), ".myrmidon"))
}



clone_capsules(capsule_source_files_list, capsule_destination_files_list, meta_data)


#### To do next ####
# test if queen capsule is cloned correctly
# adjust vital main script to source this function!
# rename tracking files that are ready
# define different capsule definitions
# run function to clone defined capsules to all files
# give stuff to Nathalie to adjust her script accordingly. 



