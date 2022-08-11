rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1. ANT-POSE-CLONER  #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### Information ####
# If ants are used in two separate tracking systems it is usefull if the ant-pose, size and orientation are identical for both tracking files
# (e.g. tracking of a colony in one tracking system and then a few workers get samplet for treatment that is also tracked and then the workers get returned to their colony)
# especially if the treatment tracking is only short and there are no picutres saved yet that can be used for manual orientation
# (when doing short trackings make sure to adjudt the frequency at which pictures are taken - e.g every 10-15 minutes or so for a 2h tracking so you have enought to manually orient ants)
# This script will copy the ant orientation/pose and size from the main tracking and apply it the second tracking file with a scaling factor based on tag size (as not both tracking systems are set up identically)




