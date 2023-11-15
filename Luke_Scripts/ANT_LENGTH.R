output_name <- file.path(paste0(dir_data,"Mean_ant_length_colonies.txt"))

for (element in no_capsule_list) {
  ant_measurements <- NULL
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
}
