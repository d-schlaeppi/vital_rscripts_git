# Directory path
dir_path <- "/media/cf19810/SeagateDesktopDrive/SICC_DATA/main_experiment/original_data/time_aggregation_info"
files <- list.files(path = dir_path, full.names = TRUE)


# Iterate through the files
for (file_path in files) {
  # Extract the base name of the file
  file_name <- basename(file_path)
  
  # Construct the new name by replacing "PreTreatment" with "PreTreatment_colony"
  #new_name <- gsub("Treatment", "Treatment_colony", file_name)
  input_vector <- strsplit(basename(file_name), "_")
  colony<-input_vector[[1]][1]
  desired_length=3
  num_zeros <- desired_length - nchar(colony)
  # Pad the values with leading zeros to the desired length
  colony <- str_pad(colony, width = desired_length, pad = "0")
  # Full path to the new file
  input_vector[[1]][1]<-colony
  new_name <- paste0(input_vector[[1]], collapse = "_")
  new_path <- file.path(dir_path, new_name)
  
  # Rename the file
  file.rename(file_path, new_path)
}

#for nested (random)
# List all files in the directory
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}
rename_files_in_dirs <- function(dirs) {
  
  for (dir_path in dirs) {
    # Get a list of all files in the directory
    files <- list.files(path = dir_path, full.names = TRUE)
    
    # Iterate through the files
    for (file_path in files) {
      # Extract the base name of the file
      file_name <- basename(file_path)
      
      # Extract the first part of the file name, pad it with zeros, and reconstruct the name
      input_vector <- strsplit(basename(file_name), "_")
      colony <- input_vector[[1]][1]
      desired_length = 3
      num_zeros <- desired_length - nchar(colony)
      # Pad the values with leading zeros to the desired length
      colony <- str_pad(colony, width = desired_length, pad = "0")
      # Full path to the new file
      input_vector[[1]][1] <- colony
      new_name <- paste0(input_vector[[1]], collapse = "_")
      new_path <- file.path(dir_path, new_name)
      
      # Rename the file
      file.rename(file_path, new_path)
    }
  }
}

dir_path <- "/media/cf19810/SeagateDesktopDrive/SICC_DATA/main_experiment/intermediary_analysis_steps/binned_interaction_lists/PreTreatment"

# list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(dir_path, n = 1)
rename_files_in_dirs(files_list)
