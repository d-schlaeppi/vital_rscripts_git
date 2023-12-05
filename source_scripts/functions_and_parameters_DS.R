# 
# clean <- function(){
#   rm(list=ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv)%in%to_keep], envir = .GlobalEnv)
#   no_print <- gc(verbose=F)
#   Sys.sleep(1)
# }
# 
# closest_match <- function(x,y){
#   return(min(which(abs(x-y)==min(abs(x-y))),na.rm=T))
# }
# 
# read.tag <- function(tag_list){ #AW    ### probabaly not needed 
#   tag <- paste0(tag_list,list.files(tag_list)[grep(colony,list.files(tag_list))])
#   tag <- read.table(tag,header=T,stringsAsFactors = F) #AW
# return(tag)
# }
# 
# 
# 
# high_threshold <- 0.0411 * 0.5/0.3###Science paper: high threshold =0.0411 where 1 = load of treated
#                         ####In Adriano's experiment, spore concentration was the same but volume was 0.3 microliter instead of 0.5 microliter
# 
# info           <- read.table(paste(data_path,"/original_data/info.txt",sep=""),header=T,stringsAsFactors = F)
# treated        <- read.table(paste(data_path,"/original_data/treated_worker_list.txt",sep=""),header=T,stringsAsFactors = F)
# task_groups    <- read.table(paste(data_path,"original_data/task_groups.txt",sep="/"),header=T,stringsAsFactors = F)
# if (grepl("age",data_path)){
#   ages      <- read.table(paste(data_path,"original_data","ant_ages.txt",sep="/"),header=T,stringsAsFactors = F)
# }
# 
# if(file.exists(paste(data_path,"/original_data/info_dat.txt",sep=""))){
#   info_datfile   <- read.table(paste(data_path,"/original_data/info_dat.txt",sep=""),header=T,stringsAsFactors = F)
#   info_plume     <- read.table(paste(data_path,"/original_data/info_plume.txt",sep=""),header=T,stringsAsFactors = F)
# }
# 
# ### get time_aggregation_list
# if (!grepl("survival",data_path)){
#   input_aggregation_info        <- paste(data_path,"original_data/time_aggregation_info",sep="/")
#   setwd(input_aggregation_info)
#   split_list                    <- paste(input_aggregation_info,list.files(pattern="txt"),sep="/")
# }
# 

# ###get tag list ### probably not needed
# input_tag <- paste(data_path,"original_data/tag_files",sep="/")
# setwd(input_tag)
# tag_list <- paste(input_tag,list.files(pattern="tags"),sep="/")
# 
# 
# 
