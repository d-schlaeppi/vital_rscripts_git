# CHECK THIS SCRIPT FOR MISTAKES BEFORE USE


# this script contains: 
# 1. CLONE CAPSULES FROM MANUAL TO MANUAL FILES
# 2. Extraction of in formation on AntPose and Capsules from manually oriented data

# for more information on fort-myrmidon and fort-studio see: 
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"

rm(list=ls())
gc()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1. CLONE CAPSULES FROM MANUAL TO MANUAL FILES ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Steps to take first:
# Create manually oriented base files
# Orient 1 lanrge colony per tracking system (food flow experiment Daniel - 3 tracking systems used for main tracking and 2 tracking systems used to record feeding sessions)
# used tracking systems Daniel: Maintracking - trojan, prideaux, guillam ; feeding sessions - guillam, esterhase 
# choose 1 of the manually oriented colonies to define a capsule definition in fort studio for a medium sized ant and replicate the shape for all of the ants of the colony
# apply this capsule definition for the remaining manually oriented colonies using the script below. The originals of these files have been stored as *.myrmidon.old

# manually oriented colonies: 
#     guillam_c03
#     prideaux_c02
#     trojan_c27

# question: do I also need the mean body size for the feeding tracking systems? guillam and esterhase?



#### prerequisites ####

### load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)

### set directory of data and myrmidon files
# dir_data <- '/media/gw20248/gismo_hd2/vital/fc2/'
dir_data <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(dir_data)


### myrmidon file (manually oriented) containing the capsule definition which was applied/cloned to all ants of the colony
defined_capsule_file <- paste(dir_data,"vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented_CapsuleManuallyDefined.myrmidon",sep='')

### List of remaining manually oriented files without capsules
no_capsule_list <- list(
  paste(dir_data,"vital_fc2_guillam_c03_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''), 
  paste(dir_data,"vital_fc2_trojan_c27_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''), 
  paste(dir_data,"vital_fc2_esterhase_c02_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_guillam_c12_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_guillam_c27_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep='')
  )

# define name of output file
#output_name <- file.path(paste0('/media/gw20248/gismo_hd2/vital/fc2/', "Mean_ant_length_colonies.txt")) 
output_name <- file.path(paste0(dir_data,"Mean_ant_length_colonies.txt")) 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2. Extraction of information on AntPose and Capsules from manually oriented data ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# prerequisites
# all manually oriented data files need to have the same capsule definition assigned for the following to work

data_list         <- list (defined_capsule_file) # here list all the myrmidon files containing oriented data
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
}






#### the queen ####
# Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen
# The queen needs to be removed as an outlier to get a good measure of mean worker size
# => find and remove outliers in the ant length measure (provided there is enough variation between queen and worker size)

### defining outliers
interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75), na.rm =TRUE)
outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))

### apply outlier exclusion to oriented_metadata
oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  

### apply outlier exclusion to capsule list
for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
}

for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
}



#### Write capsule data to each manually oriented file ####
for (no_capsule_file in no_capsule_list) {
  ANT.LENGTH <- NULL
  # open tracking data which need new capsule
  tracking_data <- fmExperimentOpen(no_capsule_file) 
  ants <- tracking_data$ants
  
  for (ant in ants){
    ###extract ant length and capsules
    ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_px)
    
    ANT.LENGTH <- rbind(ANT.LENGTH,data.frame(
      length_px        = ant_length_px,
      stringsAsFactors = F))
  }
  
  #create dataset with mean values x colony
  ant_length_colony <- data.frame(ant.length=mean(ANT.LENGTH$length_px,na.rm=T), colony=no_capsule_file)
  
  
  for (caps in 1:length(capsule_list)){
    tracking_data$createAntShapeType(names(capsule_list)[caps])
  }
  
  for (i in 1:length(ants)){
    #use mean size of each manually oriented file that needs the capsule
    ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ants[[i]]$ID)$length_px)
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
  
  
  #tracking_data$save(no_capsule_file) 
  tracking_data$save(paste0(sub("\\..*", "", no_capsule_file),"_CapsAutoDefined.myrmidon"))
  
  ## save
  if (file.exists(output_name)){
    write.table(ant_length_colony,file=output_name,append=T,col.names=F,row.names=F,quote=T,sep=",")
  }else{
    write.table(ant_length_colony,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
  }
  
  
  #close experiment
  rm(list=(c("tracking_data")))
  
}





 
#### Extract ant length for defined capsule file ####

for (defined_capsule_file in defined_capsule_file) {
  # open tracking data which need new capsule
  tracking_data <- fmExperimentOpen(defined_capsule_file) 
  ants <- tracking_data$ants
}
for (ant in ants){
  ###extract ant length and capsules
  ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_px)
  
  ANT.LENGTH <- rbind(ANT.LENGTH,data.frame(
    length_px        = ant_length_px,
    stringsAsFactors = F))
}

#### store mean values for each colony ####

# create dataset with mean values x colony

ant_length_colony1 <- data.frame(ant.length=mean(ANT.LENGTH$length_px,na.rm=T), colony=defined_capsule_file)

#IMPORTANT:
# WOULD BE BETTER TO DO THIS STEP AUTOMATICALLY BUT IT HAS BEEN DONE MANUALLY AT THE MOMENT:
# SAVE THE FILES WITH A NEW NAME "NAME_TrackSystemName-base.myrmidon" 
# THIS IS THE BASE INPUT FOR THE FOLLOWING STEP IN auto_orientation_loop.R

#### LAST! ####

# Go to fort-studio to check things look all right 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ####













