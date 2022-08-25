# this script contains: 
# 1. ANT-RULER to measure mean ant size in pixel and mm 

rm(list=ls())
gc()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1. ANT RULER  #### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### INFORMATION ####

# create a list with all the tracking files for which a measurement is needed and run the loop to create a text file containing mean size in pixel and mm
# for more information on fort-myrmidon and fort-studio see: 
# https://formicidae-tracker.github.io/myrmidon/latest/index.html

#### PREREQUISITES ####

# load required libraries
library(FortMyrmidon) # R bindings

# set and save working directory
### set directory of data and myrmidon files
# dir_data <- '/media/gw20248/gismo_hd2/vital/fc2/'
dir_data <- "/home/gw20248/Documents/data_copy_for_trials/"
setwd(dir_data)

# List of files for which you would like the mean ant size (they need to be manually oriented because thats how the head-tail measurement is made)
files <- list(
  paste(dir_data,"vital_fc2_guillam_c03_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''), 
  paste(dir_data,"vital_fc2_trojan_c27_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_prideaux_c02_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''), 
  paste(dir_data,"vital_fc2_esterhase_c02_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_guillam_c12_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep=''),
  paste(dir_data,"vital_fc2_guillam_c27_feeding_DS_AntsCreated_ManuallyOriented.myrmidon",sep='')
)

#### ANT-RULER ####

output_name <- file.path(paste0(dir_data,"Mean_ant_length_colonies_new.txt")) # define the name of the textfile containing your measurements

for (element in files) {
  # get tracking data
  ant_measurements <- NULL
  tracking_data <- fmExperimentOpen(element) 
  ants <- tracking_data$ants
  for (ant in ants){ # for each ant get the mean size and put it in a dataframe
    ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_px)
    ant_length_mm <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ant$ID)$length_mm)
    ant_measurements <- rbind(ant_measurements, data.frame(length_px = ant_length_px,
                                                           length_mm = ant_length_mm,
                                                           stringsAsFactors = F))
  }
  # queen exclusion from the dataframe
  interquartile_range <- quantile(ant_measurements$length_px,probs=c(0.25,0.75), na.rm =TRUE)
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  ant_measurements <- ant_measurements[which(ant_measurements$length_px>=outlier_bounds[1]&ant_measurements$length_px<=outlier_bounds[2]),]
  # printing and saving of the means from the dataframes without the queen
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
}

