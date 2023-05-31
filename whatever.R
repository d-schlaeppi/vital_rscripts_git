rm(list=ls())

# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()
library(Rcpp)         # contains sourceCpp (used for ant orientation)
library(circular)     # used for ant orientation
library(data.table)   # used to save files fwrite(list(myVector), file = "myFile.csv")
library(stringr)
library(utils)


#### just some base file generation ? ####
directory <- "/home/gw20248/Documents/feeding_data/"
list.files(directory)

colonies <- c("c05", "c09", "c12", "c13", "c17", "c21")
basefile <- paste0(directory, "base_source.myrmidon")
for (colony in colonies) {
  tracking_data <- fmExperimentOpen(basefile)
  tracking_data$save(paste0(directory,"base_", colony, ".myrmidon"))
}



#### messing around with the trophy data to get the right ant ID's ####
directory <- '/media/gw20248/gismo_hd5/trophy_data/' 
file <- "trophy_01_ants_oriented_old.myrmidon"
tracking_data <- fmExperimentOpen(paste(directory, file,sep=''))

#getting the right function to print the hexadecimal tagID
tracking_data$ants[[2]]$identifications
tracking_data$ants[[2]]$identifications[[1]]$tagValue
tracking_data$ants[[2]]$identifications[[1]]$targetAntID

number <- 255
hexadecimal <- sprintf("%x", number)
print(hexadecimal)
length(tracking_data$ants)
sprintf("0x%03x", tracking_data$ants[[2]]$identifications[[1]]$tagValue)
tracking_data$ants[[12]]$identifications[[1]]$tagValue

# loop to get all tag_ID's
{
tagID <- c()
ant_ID_old <- c()
for (i in 1:length(tracking_data$ants)) {
  if (!is.null(tracking_data$ants[[i]]$identifications) &&   # Check if the required elements exist 
      length(tracking_data$ants[[i]]$identifications) > 0 &&
      !is.null(tracking_data$ants[[i]]$identifications[[1]]$tagValue)) {
    # Extract the number and convert it to hexadecimal format
    hex_number <- tracking_data$ants[[i]]$identifications[[1]]$tagValue
    hexadecimal <- sprintf("0x%03x", hex_number)
    dec_number <- tracking_data$ants[[i]]$identifications[[1]]$targetAntID
    decimal <- sprintf("%03d", dec_number)
    # Append the hexadecimal value to the vector
    tagID <- c(tagID, hexadecimal)
    ant_ID_old <- c(ant_ID_old, decimal)
    
  } else {
    # Value is missing, skip this iteration
    next
  }
  ID_df <- data.frame(tagID, ant_ID_old)
}
print(ID_df)
}

#copy stuff to the clipboard to insert in excel
library(clipr)
write_clip(paste0(tagID, collapse = "\n"))
# copy stuff from clipboard into R
my_data <- read_clip() # read vector
#or
my_data <- read_clip_tbl() #read table
my_data



ID_vector <- read_clip() # copy vector of old antID sequences from clipboard to get corresponding tagIDs
# Create a new vector with the matching tagID's instead of antID's 
ID_vector_new <- character(length(ID_vector))  
for (i in 1:length(ID_vector)) {
  match_row <- which(ID_df$ant_ID_old == ID_vector[i])  # Find the row with matching ant_ID_old
  if (length(match_row) > 0) {
    ID_vector_new[i] <- ID_df$tagID[match_row]  # Copy the corresponding tagID value to the new vector
  }
}
write_clip(paste0(ID_vector_new, collapse = "\n"))


#### get ant ID s from new version and then check if things still fit ####

file_new <- "new_trophy_01_AntsCreated.myrmidon"
tracking_data <- fmExperimentOpen(paste(directory, file_new,sep=''))

{
  tagID <- c()
  ant_ID <- c()
  for (i in 1:length(tracking_data$ants)) {
    if (!is.null(tracking_data$ants[[i]]$identifications) &&   # Check if the required elements exist 
        length(tracking_data$ants[[i]]$identifications) > 0 &&
        !is.null(tracking_data$ants[[i]]$identifications[[1]]$tagValue)) {
      # Extract the number and convert it to hexadecimal format
      hex_number <- tracking_data$ants[[i]]$identifications[[1]]$tagValue
      hexadecimal <- sprintf("0x%03x", hex_number)
      dec_number <- tracking_data$ants[[i]]$identifications[[1]]$targetAntID
      decimal <- sprintf("%03d", dec_number)
      # Append the hexadecimal value to the vector
      tagID <- c(tagID, hexadecimal)
      ant_ID <- c(ant_ID, decimal)
      
    } else {
      # Value is missing, skip this iteration
      next
    }
    ID_df <- data.frame(tagID, ant_ID)
  }
  print(ID_df)
}

# now get the sequence of tag_id s from the annotated data and match the new ant ids, copy it into a vector and insert it into the original df
tagID_vector <- read_clip()
antID_vector <- character(length(tagID_vector))  # Initialize the new vector with empty strings
for (i in 1:length(tagID_vector)) {
  match_row <- which(ID_df$tagID == tagID_vector[i])  # Find the row with matching tagID
  if (length(match_row) > 0) {
    antID_vector[i] <- ID_df$ant_ID[match_row]  # Copy the corresponding ant_ID_old value to the new vector
  }
}
write_clip(paste0(antID_vector, collapse = "\n"))






  