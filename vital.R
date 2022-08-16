#### vital  R script - all the steps ####

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###                                                                                             ###
###   ###               ###  ###  ### ### ### ### ###           ###           ###               ###
###    ###             ###   ###  ### ### ### ### ###          ## ##          ###               ###
###     ###           ###    ###          ###                 ### ###         ###               ###
###      ###         ###     ###          ###                ###   ###        ###               ###
###       ###       ###      ###          ###               ###     ###       ###               ###
###        ###     ###       ###          ###              ###  ###  ###      ###               ###
###         ###   ###        ###          ###             ### ### ### ###     ###               ###
###          ### ###         ###          ###            ###           ###    ###               ###
###           ## ##          ###          ###           ###             ###   ### ### ### ###   ###
###            ###           ###          ###          ###               ###  ### ### ### ###   ###
###                                                                                             ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Information before tracking
# Tracking systems - don t change settings for different colonies (only small focus adjustments)
# If you use separate systems for main tracking and treatment tracking: Don't mix them. Always the same tracking systems for main tracking and different systems for treatment tracking (Unless you can use exactly the same setup) 
# Use dedicated queen tags (0x000) if possible!

# this script contains:

# auto creation of 



# Step by step processing your tracking data:
# Step 1: For each tracking system used in tracking select one exemplary colony get mean worker size for each of them(you might need antID matcher and pose cloner for the treatment systems)
#   Step 1.1 Create base myrmidon files
#   Step 1.2 Automatically generate the ants for the selected tracking files using the "ant_generator"
#   Step 1.3 Manually orient files in fort myrmidon
#   Step 1.4 Get the mean worker size per tracking system using the "body_builder“

# Step 2: Run data extrapolation for all the data based on the mean body size in each tracking system (Nathalies script)

# Step 3: Use the extrapolated data to do all the necessary post processing

# Step 4: Data analysis 













