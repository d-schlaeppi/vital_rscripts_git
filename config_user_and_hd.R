
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Configure user and hard drive ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### Read Me ####

# Script written by DS
# This is just a little script containing a function to save a bit of time when setting all the different directories across different computers.
# It can be sourced in other scripts and it will add the user name and the used hard drive to the environment which can then be used to quickly define generate
# directories that should be working across the different computers (Linux and mac) and hard drives containing the data without much further issues.

# additionally it loads a couple of libraries and small functions already used by Adriano required for the main analysis and the scripts its sourcing

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

getUserOptions <- function() {
  ## old version using a little user interface to load user...
  # win <- tktoplevel()
  # tkwm.title(win, "Are you using Linux or is this Daniel messing around on his mac?")
  # tkwm.geometry(win, "600x100")
  # choice <- tclVar("") # Variable to store the user's choice
  # on_button_click <- function(os) { # Function to handle button clicks
  #   tclvalue(choice) <- os
  #   tkdestroy(win)
  # }
  # # Create and position buttons
  # button_linux <- tkbutton(win, text = "Linux", command = function() on_button_click("Linux"), font = tkfont.create(size = 15)) 
  # button_mac <- tkbutton(win, text = "Mac", command = function() on_button_click("Mac"), font = tkfont.create(size = 15))
  # button_else <- tkbutton(win, text = "else", command = function() on_button_click("else"), font = tkfont.create(size = 15))
  # tkpack(button_linux, side = "left", padx = 10)
  # tkpack(button_mac, side = "left", padx = 10)
  # tkpack(button_else, side = "left", padx = 10)
  # tkwait.window(win) # Wait for the user to click a button
  os <- Sys.info()["sysname"]
  usr <- ""
  hd <- ""
  mac <- FALSE
  if (os == "Linux") { #use tclvalue(choice) isntead of os with the old version
    usr <- system("whoami", intern = TRUE)
    media_drives <- list.dirs(paste0("/media/", usr), full.names = TRUE, recursive = FALSE)
    if (length(media_drives) == 1) { # if only one hd is plugged in it should be automatically detected
      hd <- basename(media_drives)
    } else {
      hd <- basename(tk_choose.dir(default = paste0("/media/", usr), caption = "Select Hard Drive - double click HD and press ok"))  # select hd if more than one is loaded
    }
  } else if (os == "Darwin") {
    usr <- Sys.info()["user"]
    hd <- basename(tk_choose.dir(default = "/Volumes", caption = "Select Hard Drive - double click HD and press ok"))
    cat(red(paste("Warning: The code below has been coded for Linux and might not work as intended on this device")), "\n")
    mac <- TRUE
  } else if (tclvalue(choice) == "else") {
    cat(red(paste("Warning: Please define usr (User) and hd (hard drive) manually and adjust code to get data if necessary")), "\n")
  }
  # Output
  cat("Selected USER:", usr, "\n")
  cat("Selected HD:", hd, "\n")
  assign("usr", usr, envir = .GlobalEnv);  assign("hd", hd, envir = .GlobalEnv);  assign("mac", mac, envir = .GlobalEnv)
}
getUserOptions()




#### functions #### 

choose_data_path <- function() {
  if (mac) { base_path <- paste("/Volumes", hd, "vital/fc2/vital_experiment", sep = "/")
  } else { base_path <- paste("/media", usr, hd, "vital/fc2/vital_experiment", sep = "/")
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
