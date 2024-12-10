
#### 3.3.6 Randomization II #### 
#' use randomization to see if difference between the two food sources is caused by first discovery.

# create same variables as used for first randomization

dat_duration_first <- dat_duration_first %>%
  group_by(colony_id) %>%
  mutate(
    first_feeding_virus     = min(feeding_start_seconds[food_source == "virus"], na.rm = TRUE),
    first_feeding_control   = min(feeding_start_seconds[food_source == "control"], na.rm = TRUE),
    first_source_fed_on = ifelse(
        first_feeding_virus < discovery_time_control, "virus",
        ifelse(discovery_time_control < first_feeding_virus, "control", NA_character_))) %>% 
  ungroup() %>% as.data.frame()

colony_metadata <- dat_duration_first %>%
  dplyr::select(colony_id, treatment, position_virus, block, 
         first_feeding_virus, first_feeding_control, first_source_fed_on) %>%
  distinct(colony_id, .keep_all = TRUE)

total_iterations <- 100 # number of iterations to run for the randomisation

if (RUN_ANALYSIS_AND_SIMULATIONS_DURATION) { # RUN_ANALYSIS_AND_SIMULATIONS_DURATION <- TRUE
  
  ### Get data calculate and plot cumulative exploitation and exploitation rate
  
  for (start_time in c( "discovery_first_food", "shifted_discovery", "start_experiment")) { # Loop over three different definitions of †0 # start_time = "discovery_first_food"
    time_origin <- start_time
    # get subsetted data
    if (start_time == "start_experiment") {
      subsetted_data <- dat_duration_first # complete data set (note: data has only been annotated until sixty minutes after the second food source has been discoverd or until 2h since start of the experiment have elapsed)
    }
    if (start_time == "shifted_discovery") { 
      summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id, food_source) %>%
        summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
      subsetted_data <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id", "food_source"))
      subsetted_data <- subsetted_data %>% # filter to only include feeding events starting within the 60 min (3600s) window since discovery (first actual feeding) of the food source of interest (different discovery times)
        mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
        mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>% 
        filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)
    }
    if (start_time == "discovery_first_food") {
      summary_time_first_feeding <- dat_duration_first %>% group_by(colony_id) %>%
        summarize(time_first_feeding = min(feeding_start_seconds, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
      subsetted_data <- dat_duration_first %>% left_join(summary_time_first_feeding, by = c("colony_id"))
      subsetted_data <- subsetted_data %>% # filter to only include first hour since discovery of first food. 
        mutate(feeding_start_seconds = feeding_start_seconds - time_first_feeding) %>%
        mutate(feeding_end_seconds = feeding_end_seconds - time_first_feeding) %>%
        filter(feeding_start_seconds >= 0 & feeding_start_seconds <= 3600)
    }
    
    
    # expand a full grid with second by second rows
    grid <- expand_grid(
      colony_id = unique(subsetted_data$colony_id), 
      food_source = unique(subsetted_data$food_source),
      time = seq(0, max(subsetted_data$feeding_end_seconds, na.rm = TRUE))
    ) %>% as.data.frame()
    
    
    # calculate cumulative exploitation at each second and exploitation rate. 
    cumulative_explotation_over_time <- NULL
    
    for (colony in unique(grid$colony_id)) { # loop over colony_id and food_source
      for (food in unique(grid$food_source)) {
        current_grid <- grid %>% filter(colony_id == colony, food_source == food) # filter grid and feeding data for the current colony and food source
        current_feedings <- subsetted_data %>% filter(colony_id == colony, food_source == food)
        
        
        # compute across all rows of current_grid using vectorized operations
        if (nrow(current_feedings) > 0) {
          current_feedings <- current_feedings %>% arrange(feeding_start_seconds) # Sort by start time
          current_grid <- current_grid %>%
            rowwise() %>%
            mutate(
              cumulated_exploitation_time = sum(
                pmin(current_feedings$feeding_duration_seconds, 
                     pmax(0, time - current_feedings$feeding_start_seconds))
              ),
              exploitation_rate = sum(
                time >= current_feedings$feeding_start_seconds &
                  time <= current_feedings$feeding_end_seconds
              )
            )
        } else {
          current_grid <- current_grid %>%
            mutate(cumulated_exploitation_time = 0,
                   exploitation_rate = 0)
        }
        
        # append to cumulative_exploitation_over_time
        cumulative_explotation_over_time <- bind_rows(cumulative_explotation_over_time, current_grid) 
      }
    }
    
    
    # aggregate data and calculate means of duration for plotting
    mean_exploitation_time <- cumulative_explotation_over_time %>%
      group_by(time, food_source) %>%
      summarize(mean_time = mean(cumulated_exploitation_time, na.rm = TRUE),
                sd = sd(cumulated_exploitation_time, na.rm = TRUE), 
                sem = sem(cumulated_exploitation_time))
    title_lab_1 <- paste0("Cumulative Food Exploitation (T0 = ", start_time, ")")
    print(
      ggplot(mean_exploitation_time, aes(x = time, y = mean_time, color = food_source)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = mean_time - sem, ymax = mean_time + sem, fill = food_source), alpha = 0.2) +
        labs(
          title = title_lab_1,
          x = "Time",
          y = "Mean Cumulated Exploitation Time (seconds)") +
        scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
        theme_bw() +
        theme(legend.title = element_blank())+
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    )
    
    
    # aggregate data and calculate means of exploitation rate for plotting
    mean_exploitation_rate <- cumulative_explotation_over_time %>%
      group_by(time, food_source) %>%
      summarize(mean_exploitation = mean(exploitation_rate, na.rm = TRUE),
                sd = sd(exploitation_rate, na.rm = TRUE))
    title_lab_2 <- paste0("Exploitation Rate (T0 = ", start_time, ")")
    print(
      ggplot(mean_exploitation_rate, aes(x = time, y = mean_exploitation, color = food_source)) +
        geom_line(size = 1) +
        labs(
          title = title_lab_2,
          x = "Time",
          y = "Mean Cumulated Exploitation Time (seconds)") +
        scale_color_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
        theme_bw() +
        theme(legend.title = element_blank())+
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    )
    
    
    ### Stats 
    
    # calculate means and SE
    dynamic_summary_dat_dur <- cumulative_explotation_over_time %>%
      group_by(food_source, time) %>%
      summarize(exploitation_rate_Mean   = mean(exploitation_rate, na.rm = TRUE),
                exploitation_rate_SE = sd(exploitation_rate, na.rm = TRUE) / sqrt(sum(!is.na(exploitation_rate))),
               .groups = 'drop') %>% 
      arrange(time) %>%
      mutate(across(contains("exploitation_rate_Mean"), as.numeric)) %>%
      mutate(across(contains("exploitation_rate_SE"), as.numeric)) %>% as.data.frame()
    
    # calculate one mean exploitation rate for each source, each colony and each feeding session
    summary_dat_duration <- aggregate(exploitation_rate   ~colony_id+food_source,FUN=mean,data=cumulative_explotation_over_time)
    summary_dat_duration <- summary_dat_duration %>% left_join(colony_metadata, by = "colony_id")
    
    # model
    print( paste("Statistics - T0 = ",start_time))
    model <- lmer(exploitation_rate ~ food_source + (1|colony_id) + (1|first_source_fed_on), data=summary_dat_duration )
    print(Anova(model))
    print(shapiro.test(residuals(model))) # (consider looking for alternative models because sometimes close. but this is not the test statistics to report so ignore for now
    
    # plot overall mean
    plot_data <- summary_dat_duration %>%
      group_by(food_source) %>%
      summarize(
        mean_exploitation_rate = mean(exploitation_rate, na.rm = TRUE),
        se_exploitation_rate = sd(exploitation_rate, na.rm = TRUE) / sqrt(n()), 
        .groups = 'drop') %>% as.data.frame()
    
    p <- ggplot(plot_data, aes(x = food_source, y = mean_exploitation_rate, fill = food_source)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = mean_exploitation_rate - se_exploitation_rate, 
                        ymax = mean_exploitation_rate + se_exploitation_rate),
                    width = 0.2) +
      labs(
        x = "Food Source",
        y = "Exploitation Rate (Mean ± SE)",
        title = "Mean Exploitation Rate by Food Source"
      ) +
      scale_fill_manual(values = c("virus" = "#FBB4AE", "control" = "#CCEBC5")) +
      theme_bw() +
      theme(legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(p)
    
    ### Plot exploitation over time data
    ymin <- min(dynamic_summary_dat_dur$exploitation_rate_Mean - dynamic_summary_dat_dur$exploitation_rate_SE, na.rm = TRUE)
    ymax <- max(dynamic_summary_dat_dur$exploitation_rate_Mean + dynamic_summary_dat_dur$exploitation_rate_SE, na.rm = TRUE)
    plot(exploitation_rate_Mean ~ time, type = "l", col = "red", lwd = 1,  # plot virus line
      xlab = "time [s]",
      data = dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), ],
      ylim = c(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)),
      bty = "l", 
      main = paste("Exploitation Rate -", time_origin, "t0"))
    polygon( # polygon virus
      x = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "time"],
        rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "time"])),
      y = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_Mean"] - dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_SE"],
        rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_Mean"] + dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "virus"), "exploitation_rate_SE"])),
      col = alpha("#FBB4AE", 0.5), border = NA)
    lines( # line control food 
      exploitation_rate_Mean ~ time, col = "green", lwd = 1,  
      data = dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), ])
    polygon( #  polygon - control food
      x = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "time"],
        rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "time"])),
      y = c(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_Mean"] - dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_SE"],
        rev(dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_Mean"] + dynamic_summary_dat_dur[which(dynamic_summary_dat_dur$food_source == "control"), "exploitation_rate_SE"])),
      col = alpha("#CCEBC5", 0.5), border = NA)
    legend("topright", legend = c("virus", "control"), col = c("red", "green"), lwd = 2, bty = "n")
 
    ### Randomisation test 
    # Aim: test whether the results are an artifact from virus food tending to be discovered earlier than healthy food
    # H0: there is no difference in ant behaviour towards either virus or healthy food, the only relevant point explaining exploitation rate is which food was discovered first
    
    # In the observed dataset the virus was discovered first a certain number of time
    # The randomisation will consist in randomly drawing the same number of colonies from the total and saying that in these colonies the first source had virus (irrespective of truth), and in all others the first source was healthy irrespective of truth)
    # If observed data is an artifact, we should observe similar exploitation rates in the number of ants between the two sources in randomised and real data
    
    ### Pre-calculations: nb of colonies that found the virus first in each session, and observed difference in number of ants across all colonies and all times

    
    nb_discoveries <- aggregate(colony_id~first_source_fed_on, FUN=length, data=colony_metadata)
    # mean difference in exploitation rate
    delta_exploitation_rate <- plot_data$mean_exploitation_rate[plot_data$food_source == "virus"] - plot_data$mean_exploitation_rate[plot_data$food_source == "control"]
    observed_Diff <- data.frame(
      status = "virus_minus_control",
      delta_exploitation_rate = delta_exploitation_rate)
    
    summary_dat_duration <- summary_dat_duration %>% left_join(colony_metadata, by = "colony_id")
    
    #initiate randomisation 
    random_data_dynamic_dur <- NULL
    mean_random_data_dynamic_dur <- NULL
    
    pb <- progress_bar$new(
      format = "Progress: :current/:total [:bar] :percent ETA: :eta",
      total = total_iterations,
      clear = FALSE,
      width = 60
    )
    
    for (i in 1 : total_iterations){ ### randomisation loop # i <- 1
      rand_dur   <- NULL
      nb_virus <- nb_discoveries[which(nb_discoveries$first_source_fed_on == "virus"), "colony_id"]
      subset_dur <- cumulative_explotation_over_time %>% left_join(colony_metadata, by = "colony_id") %>% as.data.frame() # which dataset is this  is this
      colony_list <- unique(subset_dur$colony_id)
      hypothetic_virus_first    <- sample(colony_list,size=nb_virus,replace = F)

      # now create first_source_discovered_RAND column according to hypothetical status 
      rand_subset_dur                                                                                            <- subset_dur
      rand_subset_dur$first_source_discovered_RAND                                                               <- "control"  # note: it is refering to first food source fed on and not first discoverd (slight difference in the actual data) 
      rand_subset_dur[which(rand_subset_dur$colony_id%in%hypothetic_virus_first),"first_source_discovered_RAND"] <- "virus" 

      # copy food_source column into a new column food_source_RAND, and switch the values for those colonies who have been assigned a different first source discovered than in the observed data
      rand_subset_dur$food_source_RAND  <- rand_subset_dur$food_source    
      rand_subset_dur[which(rand_subset_dur$first_source_discovered_RAND!=rand_subset_dur$first_source_fed_on &  rand_subset_dur$food_source=="virus")  , "food_source_RAND"] <- "control"
      rand_subset_dur[which(rand_subset_dur$first_source_discovered_RAND!=rand_subset_dur$first_source_fed_on &  rand_subset_dur$food_source=="healthy"), "food_source_RAND"]    <- "virus"
      
      # concatenate and store
      rand_dur <- rbind(rand_dur  , data.frame(RAND=i, rand_subset_dur))
      
      # get the mean number of ants at each time point in the hypothetical scenario, across all colonies
      mean_rand_dat_dur <- aggregate(exploitation_rate ~ food_source_RAND + time, FUN=mean, data=rand_dur)
      overall_mean_rand_dat_dur <- aggregate(exploitation_rate ~ food_source_RAND, FUN=mean, data=rand_dur)
    
      # concatenate
      random_data_dynamic_dur           <- rbind(random_data_dynamic_dur , data.frame(RAND=i, mean_rand_dat_dur))
      mean_random_data_dynamic_dur      <- rbind(mean_random_data_dynamic_dur , data.frame(RAND=i, overall_mean_rand_dat_dur))
      pb$tick()
    }
    
    mean_random_data_dynamic_dur_DIFF <- mean_random_data_dynamic_dur %>%
      pivot_wider(names_from = food_source_RAND, values_from = exploitation_rate) %>%
      mutate(delta_exploitation_rate = virus - control) %>% as.data.frame()
      
    
    ### Plot expected vs observed Exploitation rate
    xmin <- min (c(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate),observed_Diff[,"delta_exploitation_rate"])
    xmax <- max (c(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate),observed_Diff[,"delta_exploitation_rate"])
    hist(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate,col=alpha("grey",0.5),xlim=c(xmin,xmax),main=paste("Observed vs Expected,",capitalize(time_origin),"T0"), xlab=expression(paste(Delta, " exploitation rate")))
    arrows(x0=observed_Diff[,"delta_exploitation_rate"],
           y0=0,
           y1=-10,
           code=1,
           col="black",
           lwd=4,
           length=0.1)
    
    ### P value
    if (observed_Diff[,"delta_exploitation_rate"]>median(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate)){
      pval <- 2*length(which(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate>=observed_Diff[,"delta_exploitation_rate"]))/length(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate)
      mtext(paste("p=",pval),3, cex = 1.3)
    }else{
      pval <- 2*length(which(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate<=observed_Diff[,"delta_exploitation_rate"]))/length(mean_random_data_dynamic_dur_DIFF$delta_exploitation_rate)
      mtext(paste("p=",pval),3, cex = 1.3)
    }
  }
}
