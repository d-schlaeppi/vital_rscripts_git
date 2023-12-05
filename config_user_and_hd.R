
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Configure user and hard drive ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### Read Me ####

# Script written by DS
# This is just a little script containing a function to save a bit of time when setting all the different directories across different computers.
# It can be sourced in other scripts and it will add the user name and the used hard drive to the environment which can then be used to quickly define generate
# directories that should be working across the different computers (Linux and mac) and hard drives containing the data without much further issues.
# plus a couple of small functions already used by adriano.

library(tcltk)
getUserOptions <- function() {
  # select user
  options <- list(
    "2A13_Office_Daniel" = "gw20248",
    "Nath_office" = "bzniks",
    "AEL-laptop" = "ael",
    "mac_gismo" = "gismo",
    "2A13_Office_Adriano" = "cf19810"
  )
  usr_index <- menu(names(options), title = "Select USER option:")
  usr <- options[[names(options)[usr_index]]]
  if (is.na(usr)) {
    stop("User canceled the selection.")
  }

  # select hard drive
  hd <- basename(tk_choose.dir(default = "~/", caption = "Select Hard Drive")) # Opens a directory dialog for selecting the hard drive
  if (is.na(hd) || hd == "") {
    stop("User canceled the selection.")
  }
  
  # output
  cat("Selected USER:", usr, "\n")
  cat("Selected HD:", hd, "\n")
  assign("usr", usr, envir = .GlobalEnv)
  assign("hd", hd, envir = .GlobalEnv)
}



getUserOptions()

#### functions #### 
clean <- function(){
  rm(list=ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv)%in%to_keep], envir = .GlobalEnv)
  no_print <- gc(verbose=F)
  Sys.sleep(1) }


