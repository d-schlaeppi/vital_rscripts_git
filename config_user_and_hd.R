
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
library(tcltk)
# library(adehabitatHR)
# library(adehabitatLT)
# library(changepoint)
# library(e1071)
library(igraph)
# library(gtools)
library(Rcpp)   # contains sourceCpp()
library(survival) #contains coxph() function
# library(entropy)
library(moments) # contains skewness function... 
library(crayon) # cloring messages in the terminal
library(tcltk)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(car)
library(glmmTMB)
library(multcomp)
library(plotrix)
library(blmeco) #model testing comqqnorm
library(DHARMa) # model testing

# getUserOptions <- function() {
#   # select user
#   options <- list(
#     "2A13_Office_Daniel" = "gw20248",
#     "Nath_office" = "bzniks",
#     "AEL-laptop" = "ael",
#     "mac_gismo" = "gismo",
#     "2A13_Office_Adriano" = "cf19810",
#     "other Linux computer" = "to_be_defined"
#   )
#   usr_index <- menu(names(options), title = "Select USER option:")
#   usr <- options[[names(options)[usr_index]]]
#   if (usr == "to_be_defined") {
#     usr <- basename(tk_choose.dir(default = paste0("~/"), caption = "Select User - just press ok"))
#   }
#   if (is.na(usr)) {
#     stop("User canceled the selection.")
#   }
#   # select hard drive
#   hd <- basename(tk_choose.dir(default = paste0("/media/", usr), caption = "Select Hard Drive - doubleclick HD and press ok")) # Opens a directory dialog for selecting the hard drive
#   if (is.na(hd) || hd == "") {
#     stop("User canceled the selection.")
#   }
#   
#   # output
#   cat("Selected USER:", usr, "\n")
#   cat("Selected HD:", hd, "\n")
#   assign("usr", usr, envir = .GlobalEnv)
#   assign("hd", hd, envir = .GlobalEnv)
# }


### old version opening user interface with buttons... but too much clicking when repeatedly rerunning the script ###
# getUserOptions <- function() {
#   ## old version using a little user interface to load user...
#   # win <- tktoplevel()
#   # tkwm.title(win, "Are you using Linux or is this Daniel messing around on his mac?")
#   # tkwm.geometry(win, "600x100")
#   # choice <- tclVar("") # Variable to store the user's choice
#   # on_button_click <- function(os) { # Function to handle button clicks
#   #   tclvalue(choice) <- os
#   #   tkdestroy(win)
#   # }
#   # # Create and position buttons
#   # button_linux <- tkbutton(win, text = "Linux", command = function() on_button_click("Linux"), font = tkfont.create(size = 15)) 
#   # button_mac <- tkbutton(win, text = "Mac", command = function() on_button_click("Mac"), font = tkfont.create(size = 15))
#   # button_else <- tkbutton(win, text = "else", command = function() on_button_click("else"), font = tkfont.create(size = 15))
#   # tkpack(button_linux, side = "left", padx = 10)
#   # tkpack(button_mac, side = "left", padx = 10)
#   # tkpack(button_else, side = "left", padx = 10)
#   # tkwait.window(win) # Wait for the user to click a button
#   os <- Sys.info()["sysname"]
#   usr <- ""
#   hd <- ""
#   mac <- FALSE
#   if (os == "Linux") { #use tclvalue(choice) isntead of os with the old version
#     usr <- system("whoami", intern = TRUE)
#     media_drives <- list.dirs(paste0("/media/", usr), full.names = TRUE, recursive = FALSE)
#     if (length(media_drives) == 1) { # if only one hd is plugged in it should be automatically detected
#       hd <- basename(media_drives)
#     } else {
#       hd <- basename(tk_choose.dir(default = paste0("/media/", usr), caption = "Select Hard Drive - double click HD and press ok"))  # select hd if more than one is loaded
#     }
#   } else if (os == "Darwin") {
#     usr <- Sys.info()["user"]
#     hd <- basename(tk_choose.dir(default = "/Volumes", caption = "Select Hard Drive - double click HD and press ok"))
#     cat(red(paste("Warning: The code below has been coded for Linux and might not work as intended on this device")), "\n")
#     mac <- TRUE
#   }  else if (os == "Windows") {
#     DATADIR <- "D:/r_data/"
#     SCRIPTDIR <- "D:/r_data/"
#   }else if (tclvalue(choice) == "else") {
#     cat(red(paste("Warning: Please define usr (User) and hd (hard drive) manually and adjust code to get data if necessary")), "\n")
#   }
#   # Output
#   cat("Selected USER:", usr, "\n")
#   cat("Selected HD:", hd, "\n")
#   assign("usr", usr, envir = .GlobalEnv);  assign("hd", hd, envir = .GlobalEnv);  assign("mac", mac, envir = .GlobalEnv)
# }
# getUserOptions()

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
    hd <- basename(tk_choose.dir(default = "/Volumes", caption = "Select Hard Drive - double click HD and press ok"))
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


#### functions #### 
choose_data_path <- function() {
  if (os == "Windows") {
    cat(red("Warning: Does not work on Windows\n"))
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
  no_print <- gc(verbose=F)
  Sys.sleep(1) }

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
