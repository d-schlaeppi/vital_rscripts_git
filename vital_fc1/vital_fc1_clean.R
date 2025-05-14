# rm(list = ls())
rm(list = setdiff(ls(), "first_time_use_working_directory"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ###  VITAL - FC1 Analyses   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### READ ME ####
#' This script was written by Nathalie Stroeymeyt and Daniel Schläppi
#' Analysis of the binary food choice experiment conducted prior to the vital tracking experiment.
#' Feeding events of 16 colonies were colonies were could forage on a virus and/or a control food source are being analysed.
#' The data table contains annotations of ant feeding durations (Annotations till 1h until after the second food source has been discovered)

# 
#' INDEX
#' 1. Prerequisites
#' 



#________________________________________________________________________________________________________________________________________________
#### 1. Prerequisites ####


### directory

# NATHALIE DIRECTORY
# first_time_use_working_directory <- "~/Dropbox/SeniorLectureship_Bristol/Students_postdocs/Post-Docs/Daniel Schlaeppi/vital_rscripts_git-main/vital_fc1"

if (!exists("first_time_use_working_directory") || first_time_use_working_directory == "") {
  selected_dir <- tcltk::tk_choose.dir(default = "~/", caption = "Select Working Directory")
  if (is.null(selected_dir) || selected_dir == "") {
    cat("No directory selected. Exiting.\n")
    return()}
  setwd(selected_dir)
  first_time_use_working_directory <<- getwd()
  cat(crayon::blue(getwd()))
} else { setwd(first_time_use_working_directory)
  cat(crayon::blue(getwd())) }



### load data
dat_duration <- read.csv("vital_fc1_data_feedingdurations.csv", header = TRUE)
# remove 6 NA lines (artifact from annotations)
dat_duration <- dat_duration[!is.na(dat_duration$feeding_duration), ]



### libraries
# install.packages("pacman")
pacman::p_load(lubridate, plotrix, scales, car, lme4, Hmisc, 
               dplyr, blmeco, lmtest, lsmeans, lubridate,
               emmeans, multcompView, multcomp, viridis, crayon, 
               e1071, DHARMa, merTools, tidyr, pheatmap, grid,
               progress, tidyverse, ggplot2, glmmTMB)


# functions
sem <- function(x) {sd(x,na.rm=T)/sqrt(length(na.omit(x)))} # standard error of means
source("func_test_norm.R")  # adds test_norm() to the environment (it takes your model as an argument)

