### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### CAPSULE CLONER ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


#### READ ME ####
# This script was written by Daniel Schläppi, based on the capsule_assignment script from Adriano Wanderling and Nathalie Stroeymeyt

# This script contains a function that can be sourced to clone capsules from a source myrmidon file containing predefined capsules into a list of destination_myrmidon file.access
# The function takes as arguments 
# a list of source myrmidon files containting the capsule definitions "capsule_source_files_list"
# a list of destination myrmidon files which will get the new capsules
# Note: the queens are excluded here due to the size difference and will need a separate capsul

# for more informations go to the vital main script 
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"

#### prerequisites, i.e. things required to run this function ####

# load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp) # contains sourceCpp
library(circular) # contains circular function
library(R.utils)

# set directory to where you have your myrmidon files and assosciated data 
directory <- '/media/gw20248/gismo_hd2/vital/fc2/'
setwd(directory)

# list of myrmidon files containing the capsule definitions to be applied/cloned to the destination files (for the capsule definition you need to start with a manually oriented colony)
# save your source myrmidon files like this: filename_CapsuleDefXX.myrmidon 
# get the list of all the filenames (with path form directory) containting the capsule definition (listing it manually or getting it automatically from the directory)
capsule_source_files_list <- list.files()
capsule_source_files_list <- grep(capsule_source_files_list, pattern = 'CapsuleDef', invert = FALSE, value = TRUE)
add_directory <- function(filename) {# Function to add directory path to each filename
  paste0(directory, filename)
}
capsule_source_files_list <- lapply(capsule_source_files_list, add_directory)

# get a list af all the filenames that shall have the capsule definitions (destination files)
# save your destination myrmidon files like this: final_uniqueIdentifier.myrmidon # e.g. final_c01.myrmidon

# List of remaining manually oriented files without capsules
capsule_destination_files_list <- list.files()
capsule_destination_files_list <- grep(capsule_destination_files_list, pattern = 'final', invert = FALSE, value = TRUE)
add_directory <- function(filename) {# Function to add directory path to each filename
  paste0(directory, filename)
}
capsule_destination_files_list <- lapply(capsule_destination_files_list, add_directory)


#### capsule cloner function ####
clone_capsules <- function(capsule_source_files_list, capsule_destination_files_list) {
  for (source_file in capsule_source_files_list){
    oriented_metadata <- NULL
    capsule_list <- list()
    experiment_name <- unlist(strsplit(source_file,split="/"))[length(unlist(strsplit(source_file,split="/")))]
    oriented_data <- fmExperimentOpen(source_file)
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
# Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen -> queen removed by removing outliers in ant lenght measurements
  interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75))
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  # apply outlier exclusion to oriented_metadata and to capsule list
  for (caps in 1:length(capsule_list)){
    capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
  }
  mean_x_ant_coord <- mean(oriented_metadata$x_ant_coord)   # Once queen(s) has(have) been removed, get the mean coordinates of the offset between tag centre and ant centre
  mean_y_ant_coord <- 0 ##set it to zero manually because we don't expect there to be a consistent bias in deviation / #mean_y_ant_coord <- mean(oriented_metadata$y_ant_coord) ###this is expected to be 0 or near 0 (check!) because the tag should be as likely to be to the right or to the left of the ant's bilateral symmetry line
  mean_worker_length_px <- mean(oriented_metadata$length_px) # Get the average worker length from the data
  for (caps in 1:length(capsule_list)){ # Finally, get information on each capsule
    capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
  }
  
  
  
  # Next apply the extracted capsules to all workers of all destination files
  # then, get the corresponding queen capsule and apply it as well 
  # finally, save each myrmidon file before moving on to the next set of capsule definitions.
  for (destination_file in capsule_destination_files_list) {
    tracking_data <-fmExperimentOpen(paste0(directory, destination_file))
    for (caps in 1:length(capsule_list)) {#add the caps from the source file above
      tracking_data$createAntShapeType(names(capsule_list)[caps])
    }
    ant_length_px <- 163# get mean worker length for each colony from a vector
    queen_length_px <- y# get queen length from a vector
    ants <- tracking_data$ants
    for (i in 1:length(ants)) {
      if(length(tracking_data$ants[[i]]$identifications) == 0) {
        print(paste(destination_file, i, "no ant", sep = " -> "))
        next
      }
      if (tracking_data$ants[[i]]$identifications[[1]]$tagValue==0) {next} # skip the queen  ----------- <- <- <- <- <- <- <-  here apply a different function to assign the queen? 
      ants[[i]]$clearCapsules() # clear previous capsules 
      if (length(tracking_data$antShapeTypeNames)>0) {# delete the capsule shapes IF PRESENT
        for (caps in 1:length(tracking_data$antShapeTypeNames)){
          tracking_data$deleteAntShapeType(caps)
        }
      }
      
      
      ####### continue here!!!!!!!
      for (caps in 1:length(capsule_list)){
        capsule_ratios <- capsule_list[[caps]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
        capsule_coords <- ant_length_px*capsule_ratios
        ants[[i]]$addCapsule(caps, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
      }
    }
    tracking_data$save(paste0(substr(destination_file, 0  , nchar(destination_file)-9), substr(source_file, nchar(source_file)-21, nchar(source_file))))

      }

      
    }
    

  }


  
}






       







#### queen capsule cloner ####

#### run the function on your data ####

clone_capsules(capsule_source_files_list, capsule_destination_files_list)
