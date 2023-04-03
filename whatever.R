rm(list=ls())

# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()
library(Rcpp)         # contains sourceCpp (used for ant orientation)
library(circular)     # used for ant orientation
library(data.table)   # used to save files fwrite(list(myVector), file = "myFile.csv")
library(stringr)


directory <- "/home/gw20248/Documents/feeding_data/"
list.files(directory)

colonies <- c("c05", "c09", "c12", "c13", "c17", "c21")
basefile <- paste0(directory, "base_source.myrmidon")

for (colony in colonies) {
  tracking_data <- fmExperimentOpen(basefile)
  tracking_data$save(paste0(directory,"base_", colony, ".myrmidon"))
}
  
  