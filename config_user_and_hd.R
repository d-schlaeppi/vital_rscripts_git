
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Configure user and hard drive ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### Read Me ####

# Script written by DS
# Script containing a functions to save a bit of time when setting up the different directories across different computers and loading functions
# Can be sourced in other scripts and it will add the user name and the used hard drive to the environment which can then be used to quickly define generate
# directories that should be working across the different computers (Linux, mac & windows) and hard drives
# additionally it loads a couple of libraries and functions required in other scripts further down the analysis pipelines

#### Libraries ####
pacman::p_load(
  tcltk, igraph, ggraph, dplyr, tidyr, ggplot2, lme4, car, glmmTMB, multcomp, multcompView, crayon, plotrix, emmeans, 
  reshape2, circular, readr, data.table, measurements, viridis, plotwidgets , scales, parallel, epicontacts, gridExtra,
  Rcpp,   # contains sourceCpp()
  survival, #contains coxph() function
  moments, # contains skewness function
  DHARMa, blmeco, FortMyrmidon, progress #model testing comqqnorm
  )

# further libraries used by adriano but probably not required
# adehabitatHR, adehabitatLT, changepoint, e1071, gtools, entropy


#### get information on operation system, user and hard drive ####
getUserOptions <- function() {
  os <- Sys.info()["sysname"]
  usr <- ""
  hd <- ""
  if (os == "Linux") { 
    usr <- system("whoami", intern = TRUE)
    media_drives <- list.dirs(paste0("/media/", usr), full.names = TRUE, recursive = FALSE)
    if (length(media_drives) == 1) { 
      hd <- basename(media_drives)
    } else {
      hd <- basename(tk_choose.dir(default = paste0("/media/", usr), caption = "Select Hard Drive - double click HD and press ok"))  # select hd if more than one is loaded
    }
  } else if (os == "Darwin") {
    usr <- Sys.info()["user"]
    volumes <- list.dirs("/Volumes", full.names = TRUE, recursive = FALSE)
    volumes <- volumes[!grepl("/Volumes/Macintosh HD", volumes)]
    if (length(volumes) == 1) {
      hd <- basename(volumes)
    } else {
      hd <- basename(tk_choose.dir(default = "/Volumes", caption = "Select Hard Drive - double click HD and press ok"))
    }
    cat(red(paste("Warning: The code below has been coded for Linux and might not work as intended on this device")), "\n")
  } else if (os == "Windows") {
    usr <- Sys.getenv("USERNAME")
    drives <- system("wmic logicaldisk get name", intern = TRUE)
    drives <- drives[drives != ""]  # remove empty lines
    drives <- drives[!grepl("Name", drives)]  # remove the header
    drives <- gsub("\\s", "", drives)  # remove any whitespace
    if (length(drives) == 1) { 
      hd <- drives
    } else {
      hd <- tk_choose.dir(default = drives[1], caption = "Select Hard Drive - double click HD and press ok")
    }
  } else {
    cat(red(paste("Warning: Please define usr (User) and hd (hard drive) manually and adjust code to get data if necessary")), "\n")
  }
  # Output
  cat("Selected USER:", usr, "\n")
  cat("Selected HD:", hd, "\n")
  cat("Operating System:", os, "\n")
  assign("usr", usr, envir = .GlobalEnv)
  assign("hd", hd, envir = .GlobalEnv)
  assign("os", os, envir = .GlobalEnv)
}

getUserOptions()



#### define data paths based on user information ####
set_directories <- function(os, hd, usr) {
  if (os == "Linux") {
    DATADIR <- paste("/media", usr, hd, "vital/fc2", sep = "/")
    SCRIPTDIR <- paste("/media", usr, hd, "vital/vital_rscripts_git", sep = "/")
  } else if (os == "Darwin") {
    DATADIR <- paste("/Volumes", hd, "vital/fc2", sep = "/")
    SCRIPTDIR <- paste("/Volumes", hd, "vital/vital_rscripts_git", sep = "/")
  } else if (os == "Windows") {
    hd <- gsub("\\\\$", "", hd)  # Remove trailing backslash if present
    if (!grepl("^[A-Za-z]:", hd)) {
      stop("Invalid hard drive name. It should start with a drive letter, e.g., 'C:'.")
    }
    DATADIR <- paste(hd, "vital/fc2", sep = "/")
    SCRIPTDIR <- paste(hd, "vital/vital_rscripts_git", sep = "/")
  }
  # Set working directory and assign variables to global environment
  setwd(DATADIR)
  assign("DATADIR", DATADIR, envir = .GlobalEnv)
  assign("SCRIPTDIR", SCRIPTDIR, envir = .GlobalEnv)
  cat("DATADIR:", DATADIR, "\n")
  cat("SCRIPTDIR:", SCRIPTDIR, "\n")
  cat(blue(paste("Current working directory:", DATADIR)), "\n")
}

set_directories(os, hd, usr)


#### Additional Functions ####
plot_model_diagnostics <- function(model) {
  x <- compareqqnorm(model)   # Compare QQ plot
  cat("Press Enter to continue...")
  readline()
  cat(blue(paste("The plot showing the residuals of your model is:")), "\n")
  cat(red(paste(x), "\n"))
  cat("Press Enter to continue...")
  readline()
  cat("Press Enter to continue...")
  sim_res <- simulateResiduals(model)   # Simulate residuals using DHARMa
  plot(sim_res)
  readline()
  par(mfrow = c(1, 3))
  testUniformity(sim_res)  # Test for uniformity (KS test),
  testOutliers(sim_res)  # Test for outliers
  testDispersion(sim_res) # Test for homoscedasticity
  #text(1.5, 20, " Test for homoscedasticity against heteroscedacity", cex = 1)
  #text(1.5, 19,"(N.S.--> residuals have constant variance = homoscedasticity)", cex = 1)
  layout(1)
} # simple finction to display some model plots - not sure what models it works for... should be ok for lmer



choose_data_path <- function() {
  if (os == "Windows") {
    cat(red("Warning: Does not yet work on Windows\n"))
    stop("This is based on an older version of Adrianos scipts that does not support Windows but should work on Linux and Mac.")
  } else {
    if (os == "Darwin") { 
      base_path <- paste("/Volumes", hd, "vital/fc2/vital_experiment", sep = "/")
    }
    if (os == "Linux")  { 
      base_path <- paste("/media", usr, hd, "vital/fc2/vital_experiment", sep = "/")
    }
  }
  list(
    CLASSIC_INTERACTIONS = if (RUN_CLASSIC_INTERACTIONS) paste(base_path, "main_experiment", sep = "/") else NULL,
    GROOMING_INTERACTIONS = if (RUN_GROOMING_INTERACTIONS) paste(base_path, "main_experiment_grooming", sep = "/") else NULL,
    TROPHALLACTIC_INTERACTIONS = if (RUN_TROPHALLACTIC_INTERACTIONS) paste(base_path, "main_experiment_trophallaxis", sep = "/") else NULL
  )
}


# reading table with different delimiters
read_table_with_auto_delim <- function(file_path) {
  lines <- readLines(file_path, n = 5)   # Read first lines of  file
  potential_delims <- c("\t", " ") # Initialize potential delimiters (could be more than just two but I want it to be only space and tabs) c("\t", ",", ";", " ")
  delim_counts <- sapply(potential_delims, function(delim) { # Count delimiters
    sum(grepl(delim, lines))
  })
  chosen_delim <- potential_delims[which.max(delim_counts)] # select most likely delimiter based on count
  table <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, sep = chosen_delim) # read table
  return(table)
}

clean <- function(){
  rm(list=ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv)%in%to_keep], envir = .GlobalEnv)
  no_print <- gc(verbose=F) }

# read.tag <- function(tag_list){ #AW
#   tag <- paste0(tag_list,list.files(tag_list)[grep(colony,list.files(tag_list))])
#   tag <- read.table(tag,header=T,stringsAsFactors = F) #AW
#   return(tag)
# }

read.tag <- function(tag_list, colony){
  tag <- paste0(tag_list,list.files(tag_list)[grep(colony,list.files(tag_list))])
  tag <- read.table(tag,header=T,stringsAsFactors = F)
  return(tag)
} # used in subsequent scripts to read the tag list (generated earlier) for specific colonies


# test normal distribution of residuals in linear models:
test_norm <- function(mod_or_resids) { # function adjusted from Nathalie's test_norm # mod_or_resids <- model
  cat(blue("\nTesting normality of residuals for lm or lmer \n(assumptions likely different other types of models)\n"))
  
  # Check if the input is a model object or residuals
  if (inherits(mod_or_resids, "lm") || inherits(mod_or_resids, "lmerMod")) {
    resids <- residuals(mod_or_resids)  # Extract residuals from the model
    #cat("Model provided: Extracting residuals from the model.\n")
  } else if (is.numeric(mod_or_resids) || is.vector(mod_or_resids)) {
    resids <- mod_or_resids  # Directly use the provided residuals
    #cat("Residuals directly provided.\n")
  } else {
    stop("Invalid input: Provide either a model object (lm/lmer) or a numeric vector of residuals.")
  }
  
  resids <- residuals(object = mod_or_resids)
  
  if (length(resids) <= 300) {
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    shapiro_result <- shapiro.test(resids)
    print(shapiro_result)
    p_value <- shapiro_result$p.value
    if(p_value >= 0.05) {
      cat(green("Data is normally distributed\n --> Model is fine... go ahead \U1F44D \n"))
    } else {
      cat(red("Warning: p value below 0.05\n-->  Data significantly deviate from a normal distribution. \nAdjust or change model! \n"))
    }
  } else {
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero)")
    skew <- skewness(resids)
    print(paste0("Skewness = ", skew))
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    kurt <- kurtosis(resids)
    print(paste0("Kurtosis = ", kurt))
    if(skew >= -3 & skew <= 3 & kurt < 4){
      cat(green("Model is fine... go ahead \U1F44D \n"))
    } else {
      cat(red("Warning: Adjust or change model! \n"))
    }
  }
}

overdispersion_test <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson_chisq <- sum(rp^2)
  prat <- Pearson_chisq / rdf
  pval <- pchisq(Pearson_chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson_chisq, ratio = prat, rdf = rdf, p = pval)
}
