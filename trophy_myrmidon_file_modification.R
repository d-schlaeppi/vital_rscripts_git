# trophy myrmidon file modifications to make it fit with the FCII files so they can go together into the trophallyxis inferral process

directory <- "/media/gw20248/gismo_hd2/vital/fc2/"
setwd()
data_list <- list.files()
data_list <- grep(data_list, pattern = 'trophy', invert = FALSE, value = TRUE)
data_list <- grep(data_list, pattern = 'myrmidon', invert = FALSE, value = TRUE)

#create metadata keys
for (i in data_list) {
  fort_data <- fmExperimentOpen(i)
  fort_data$setMetaDataKey(key = "meta_ID",     default_Value = 001)  # create the key variables you want as metadata in your data sets 
  fort_data$setMetaDataKey(key = "IsQueen",     default_Value = FALSE)
  fort_data$setMetaDataKey(key = "IsTreated",   default_Value = FALSE)
  fort_data$setMetaDataKey(key = "IsAlive",     default_Value = TRUE)
  fort_data$setMetaDataKey(key = "treatment",   default_Value = "NA") # treated ants will get control or virus
  fort_data$setMetaDataKey(key = "glass_beads", default_Value = "NA") # treated ants will get yellow or blue
  fort_data$setMetaDataKey(key = "comment",     default_Value = "NA")
  fort_data$setMetaDataKey(key = "tag_reoriented", default_Value = FALSE)
  fort_data$spaces[[1]]$createZone(name = "nest") # create zones to be defined manually in the fort files 
  fort_data$spaces[[1]]$createZone(name = "arena")
  fort_data$spaces[[1]]$createZone(name = "water_left")
  fort_data$spaces[[1]]$createZone(name = "sugar_right")
  for (y in 1:length(fort_data$ants)) {
    fort_data$ants[[y]]$setValue(key="meta_ID", value = c(fort_data$ants[[y]]$identifications[[1]]$targetAntID), time = fmTimeSinceEver())
  }
  ants <- fort_data$ants  # for each ant adjust the meta data if it is the queen or a treated worker
  for (x in ants) {
    if (x$identifications[[1]]$tagValue==0) {
      x$setValue("IsQueen", TRUE, time = fmTimeSinceEver())}
  }
  fort_data$save(paste0(directory,  'm_meta_keyed_c29.myrmidon')) ### adjust if more than one colony is done!
}

# move on to the ant orientation manuscript



