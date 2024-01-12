
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
# library(survival)
# library(entropy)


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

read.tag <- function(tag_list){ #AW
  tag <- paste0(tag_list,list.files(tag_list)[grep(colony,list.files(tag_list))])
  tag <- read.table(tag,header=T,stringsAsFactors = F) #AW
  return(tag)
}
