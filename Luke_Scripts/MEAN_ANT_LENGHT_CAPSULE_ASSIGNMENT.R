
rm(list=ls())
gc()

########################################################################
########### CLONE CAPSULES FROM MANUAL TO MANUAL FILES #################
# CREATE THE MANUALLY ORIENTED BASE FILES
# Orient 1 large colony per tracking system used (5 total) by hand.
# pick 1 out of this 5 and, in FortStudio, create a capsule definition for a medium sized ant and replicate the shape for all of the ants of the colony.
# Copy this capsule for the remaining 4 colonies using Clone_capsule_manual_to_manual.R . The originals of these files have been stored as *.myrmidon.old

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

######load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)

### directory of data and myrmidon files
#dir_data <- '/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/'
dir_data <- '/media/ll16598/SeagateDesktopDrive/SICC_DATA/'


### tracking data folder name
#track_data_name <- 'EG_NTM_s30_MODa.0000'

### manually oriented myrmidon file name
#defined_capsule_file <- paste(dir_data,"REP1/R1SP_27-02-21_Oriented.myrmidon",sep='') #polyakov
defined_capsule_file <- paste(dir_data,"13P_IV31_110522_HDN/13P_IV31_110522_HDN_manual_orient_test_withMetaData-base.myrmidon",sep='') #prideaux


### manually oriented myrmidon file name
#no_capsule_file <- 'EG_NTM_s25_DEHb.myrmidon

# List of manually oriented files without capsules
no_capsule_list <- list(
  paste(dir_data,"3P_IV9_160222_ALL/3P_IV3_160222_ALL_manual_orient.myrmidon",sep=''), 
  paste(dir_data,"6P_IV37_220322_POL/6P_IV37_220322_POL_manual_orient.myrmidon",sep=''), 
  paste(dir_data,"8S_IV21_060422_WST/8S_IV21_060422_WST_manual_orient.myrmidon",sep=''), 
  paste(dir_data,"13P_IV31_110522_HDN/13P_IV31_110522_HDN_manual_orient.myrmidon",sep=''), 
  paste(dir_data,"14P_IV20_170522_LEA/14P_IV20_170522_LEA_manual_orient.myrmidon",sep='')) 

#create output folder
output_name <- file.path(paste0('/media/ll16598/SeagateDesktopDrive/', "Mean_ant_length_colonies1.txt")) # (saved INSIDE the Network_analysis folder)


#################################################################################################################################################################################################################
###########STEP 1 - use manually oriented data to extract important information about AntPose and Capsules
###########          IN YOUR PARTICULAR SPECIES AND EXPERIMENTAL SETTINGS  
###########IMPORTANT: FOR THIS TO WORK YOU NEED TO HAVE DEFINED THE SAME CAPSULES WITH THE SAME NAMES ACROSS ALL YOUR MANUALLY ORIENTED DATA FILES
#################################################################################################################################################################################################################
#data_list         <- list ("/home/eg15396/Documents/Data/NTM/NTM_s30_auto_orient.myrmidon") ###here list all the myrmidon files containing oriented data

##### DANIEL's ammendments: calculate mean ant length

for (element in no_capsule_list) {
  ant_measurements <- NULL
  ANT.LENGTH = NULL
  tracking_data <- fmExperimentOpen(element)
  ants <- tracking_data$ants
  for (ant in ants){
    ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_px)
    ant_length_mm <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_mm)
    ant_measurements <- rbind(ant_measurements, data.frame(length_px = ant_length_px,
                                                           length_mm = ant_length_mm,
                                                           stringsAsFactors = F))
  }
  #queen exclusion
  interquartile_range <- quantile(ant_measurements$length_px,probs=c(0.25,0.75), na.rm =TRUE)
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  ant_measurements <- ant_measurements[which(ant_measurements$length_px>=outlier_bounds[1]&ant_measurements$length_px<=outlier_bounds[2]),]
  
  
  #printing and saving
  print(element)  
  print(mean(ant_measurements$length_px, na.rm=TRUE))
  print(mean(ant_measurements$length_mm, na.rm=TRUE))
  table <- NULL
  table <- rbind(table, data.frame(mean(ant_measurements$length_px, na.rm=TRUE),
                                   mean(ant_measurements$length_mm, na.rm=TRUE),
                                   element,
                                   stringsAsFactors = F))
  if (file.exists(output_name)){
    write.table(table,file=output_name,append=T,col.names=F,row.names=F,quote=T,sep=",")
  }else{
    write.table(table,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
  }
    }#IF THIS RESULTS IN A ERROR, CHECK PREVIOUS VERSION ON GIT OR COMPARE WITH DataPrep4_Clone-capule-queens-only_v082.R
  
  
  #close experiment
  rm(list=(c("tracking_data")))
##### END DANIEL AMMEND
  

###ADRIANO capsule assign - this did not remove queen from capsule mean
data_list         <- list (defined_capsule_file)
oriented_metadata <- NULL
capsule_list <- list()
for (myrmidon_file in data_list){
  experiment_name <- unlist(strsplit(myrmidon_file,split="/"))[length(unlist(strsplit(myrmidon_file,split="/")))]
  oriented_data <- fmExperimentOpen(myrmidon_file)
  oriented_ants <- oriented_data$ants
  capsule_names <- oriented_data$antShapeTypeNames
  for (ant in oriented_ants){
    ###extract ant length and capsules
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
    
    
    ###extract offset btewen tag centre and ant centre
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


### Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen #########
### So we need to remove the queen from the computation
### One way of doing so is to find and remove outliers in the ant length measure (provided there is enough variation between queen and worker size)
  interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75), na.rm = TRUE)
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  ###apply outlier exclusion to oriented_metadata...
  oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  
  ###...and to capsule list
  for (caps in 1:length(capsule_list)){
    capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
  }
  
  for (caps in 1:length(capsule_list)){
    capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
  }
  
  ### Write capsule data in each manually oriented file
  for (no_capsule_file in no_capsule_list) {
    
    # open tracking data which need new capsule
    tracking_data <- fmExperimentOpen(no_capsule_file) 
    ants <- tracking_data$ants
    for (caps in 1:length(capsule_list)){
      tracking_data$createAntShapeType(names(capsule_list)[caps])
    }
    
    for (i in 1:length(ants)){
      #use mean size of each manually oriented file that needs the capsule
      ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ants[[i]]$ID)$length_px)
      interquartile_range <- quantile(ant_measurements$length_px,probs=c(0.25,0.75), na.rm =TRUE)
      outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
      ant_measurements <- ant_measurements[which(ant_measurements$length_px>=outlier_bounds[1]&ant_measurements$length_px<=outlier_bounds[2]),]
      ant_length_px_mean <- (mean(ant_measurements$length_px, na.rm=TRUE))
      ants[[i]]$clearCapsules()
      
      ##assign capule numbers that match the order of the looped capsule names positions
      capsule_number <- 0
      for (capsule_name in unlist(tracking_data$antShapeTypeNames)) {
        capsule_number <- capsule_number +1
        # the file information
        #MAKE SURE THERE IS CAPSULE MATCHING, TO AVOID MIXING UP SHAPE INDEXES
        capsule_ratios <- capsule_list[[capsule_name]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
        capsule_coords <- ant_length_px_mean*capsule_ratios
        
        ants[[i]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
        
      }#IF THIS RESULTS IN A ERROR, CHECK PREVIOUS VERSION ON GIT OR COMPARE WITH DataPrep4_Clone-capule-queens-only_v082.R
    }
    
    
    #tracking_data$save(no_capsule_file) 
    tracking_data$save(paste0(sub("\\..*", "", no_capsule_file),"_CapsuleDef4.myrmidon"))
    
    #close experiment
    rm(list=(c("tracking_data")))
  
}#no_capsule_file
}
#IMPORTANT:
# WOULD BE BETTER TO DO THIS STEP AUTOMATICALLY BUT IT HAS BEEN DONE MANUALLY AT THE MOMENT:
# SAVE THE FILES WITH A NEW NAME "NAME_TrackSystemName-base.myrmidon" 
# THIS IS THE BASE INPUT FOR THE FOLLOWING STEP IN auto_orientation_loop.R

##LASTLY!! Go to fort-studio to check things look all right

###end adriano capsule


