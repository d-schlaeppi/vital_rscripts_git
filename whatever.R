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



directory <- '/media/gw20248/gismo_hd5/trophy_data/' 
file <- "trophy_01_ants_oriented.myrmidon"
tracking_data <- fmExperimentOpen(paste(directory, file,sep=''))
tracking_data$ants[[2]]$identifications
tracking_data$ants[[2]]$identifications[[1]]$tagValue

number <- 255
hexadecimal <- sprintf("%x", number)
print(hexadecimal)

length(tracking_data$ants)
sprintf("0x%03x", tracking_data$ants[[2]]$identifications[[1]]$tagValue)
tracking_data$ants[[12]]$identifications[[1]]$tagValue

hex_values <- c()
for (i in 1:length(tracking_data$ants)) {
  # Check if the required elements exist
  if (!is.null(tracking_data$ants[[i]]$identifications) && 
      length(tracking_data$ants[[i]]$identifications) > 0 &&
      !is.null(tracking_data$ants[[i]]$identifications[[1]]$tagValue)) {
    # Extract the number and convert it to hexadecimal format
    number <- tracking_data$ants[[i]]$identifications[[1]]$tagValue
    hexadecimal <- sprintf("0x%03x", number)
    # Append the hexadecimal value to the vector
    hex_values <- c(hex_values, hexadecimal)
  } else {
    # Value is missing, skip this iteration
    next
  }
}
print(hex_values)

library(clipr)
write_clip(paste0(hex_values, collapse = "\n"))






  