#### vital  R script - all the steps ####

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###                                                                                                 ###
###   ### ### ### ###  ###          ###       ###  ### ### ### ###  ###       ###  ### ### ### ###  ###
###   ### ### ### ###  ###          ###       ###  ### ### ### ###  ###       ###  ### ### ### ###  ###
###   ###              ###          ###       ###  ###         ###  ###       ###  ###              ###
###   ###              ###          ###       ###  ###              ###       ###  ###              ###
###   ### ### ###      ###          ###       ###  ###              ###       ###  ### ### ### ###  ###
###   ### ### ###      ###          ###       ###  ###     ### ###  ###       ###  ### ### ### ###  ###
###   ###              ###          ###       ###  ###     ### ###  ###       ###              ###  ###
###   ###              ###          ###       ###  ###         ###  ###       ###              ###  ###
###   ###              ### ### ###    ###   ###    ### ### ### ###    ###   ###    ### ### ### ###  ###
###   ###              ### ### ###       ###       ### ### ### ###       ###       ### ### ### ###  ###
###                                                                                                 ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Script guiding through the firt steps of the Analysis of tracking data for the Flugus tracking experiment   ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Background information | Read me ####

# Useful information before starting the tracking experiment. 
# Tracking systems - don't change settings (camera height etc) throughout the experiment. Settings for one specific tracking system should stay the same to minimize manual orientation
# If you use separate systems for main tracking and treatment tracking: Don't mix them. Always the same tracking systems for main tracking and different systems for treatment tracking (Unless you can use exactly the same setup) 
# If you do short trackings (e.g. for treatments in addition to a main tracking increase the rate at which pictures are taken of each ant in the leto file)
# Use dedicated queen tags (0x000) if possible!

# Step by step processing of the tracking data for flugus:
# Step 1: For each tracking system setting the first colony is selected get the mean worker size (mm and pixcels used for data extrapolation and orientation of remaining colonies)
#   Step 1.1 Create base myrmidon files
#   Step 1.2 Automatically generate the ants for the selected tracking files using the "ant_generator"
#   Step 1.3 Manually orient files in fort myrmidon
#   Step 1.4 Get the mean worker size per tracking system using the "ant_ruler“

# Step 2: Run data extrapolation for all the data based on the mean body size in each tracking system (Nathalies script)

# Step 3: Use the extrapolated data for all the colonies do all the necessary post processing

# Step 3.1: Create all base myrmidon files
# Step 3.2: Generate all ants
# Step 3.3: Automatically create the metadata variables needed
# Step 3.4: Adjust replaced or re-glued tags and do required manual post processing
# Step 3.5: Automatically orient all ants in all extrapolated data using the ant_orient express (Includes capsule generation)
# Step 3.6: Post processing of queen meta data (Manual)

# Step 4 Data analyses
# Step 4.1 Define and apply different capsules to train the trophallaxis classifier
# Step 4.2 Train behaviour (Trophallaxis) classifier using Nathalies script ant computer
# Step 4.3 General data analyses i.e. Base Analysis and Next steps 


### ### ### ### ###
### Useful links ##
### ### ### ### ###

# for more information on fort-myrmidon and fort-studio see: 
# https://formicidae-tracker.github.io/myrmidon/latest/index.html

# Postprocessing tips
# https://uob.sharepoint.com/:w:/r/teams/grp-AntsEpidemiologyLab/_layouts/15/Doc.aspx?sourcedoc=%7B2562631B-A6E5-4289-907F-89502F6C27E6%7D&file=pre-processing_Adriano_June2022.docx&action=default&mobileredirect=true

# AEL Github repositories
# https://github.com/d-schlaeppi/vital_rscripts_git
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR
# https://github.com/Leckie-Bris/SICC
# https://github.com/EnricoGavagnin?tab=repositories
# Add Daniel 
# Add Linda

#### prerequisites ####
rm(list=ls())

# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()
library(Rcpp)         # contains sourceCpp (used for ant orientation)
library(circular)     # used for ant orientation
library(data.table)   # used to save files fwrite(list(myVector), file = "myFile.csv")
library(stringr)
library(reader)

### update the loading of the colonies etc nicer... 
# directories
directory_scripts <- "/home/gw20248/Documents/vital_rscripts_git/" # directory with the R scripts linked with github
directory_data    <- '/media/gw20248/gismo_hd2/vital/fc2/'  # directory with extrapolated data and the myrmidon files




















