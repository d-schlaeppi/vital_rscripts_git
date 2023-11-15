library(igraph)
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)
# Metadata extraction file
rm(list=ls())
gc()
####################################################################################
#### THIS SCRIPT CONTAINS:
#### METADATA EXTRACTION FROM MYRMIDON FILES
####################################################################################

#### LIBRARIES
library(FortMyrmidon) ####R bindings
library(stringr)

#### FUNCTIONS
#list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}


#### ACCESS FILES
USER <- "Luke"

if (USER=="Luke") {
  WORKDIR <- "/media/ll16598/SeagateDesktopDrive/SICC_DATA" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR)#"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  SCRIPTDIR <- "/home/ll16598/Documents/SICC_ANALYSIS"
}
if (USER=="supercompAdriano") {
  WORKDIR <- "/media/cf19810/DISK4/ADRIANO" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  SCRIPTDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP1_analysis scripts"
}

###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"ANT_TASKS.R",sep="/"))


### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
#select REP folders
files_list <- files_list[grep('V',files_list)]
files_list <- files_list[!grepl('NW',files_list)]

# replicate folder
for (REP.n in 1:length(files_list)) {
  
  # REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = "withMetaData_QUEEN.myr")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
}




##extracting metadata - eventually in own loop

for (REP.FILES in REP.filefolder) {
  REP_TREAT <- stringr::str_extract(basename(REP.FILES), "[^_]*_[^_]*" )
  Metadata_Exp1 <- paste0(REP.folder,"/", "Metadata_Exp1_2022_",REP_TREAT, ".txt") 
  metadata <- read.table((Metadata_Exp1),header=T,stringsAsFactors = F, sep=",")
  exp<-fmExperimentOpen(REP.FILES)
  ###Using these, I will need to create a time dictionary to make a 3hr window and offset according to proper start time
  #From  <- fmQueryGetDataInformations(exp)$end + (TIME_HOURS - 24)*TimeWind - window_shift
  #To    <- fmQueryGetDataInformations(exp)$end + (TIME_HOURS - 21)*TimeWind - window_shift 
  From  <- fmQueryGetDataInformations(exp)$start#this is$end in adrianos
  To    <- fmQueryGetDataInformations(exp)$end
  start <- fmTimeCreate(offset=From) #end minus 48 hours plus incremental time
  end   <- fmTimeCreate(offset=To) #end minus 45 hours plus incremental time
  #gap<-fmTimeParse("10")$asPOSIXct()
  gap<-fmSecond(10)

  #compute_G <- function(exp, start, end, gap)
  G<-compute_G(exp, start, end, gap)#where does G go?
}
