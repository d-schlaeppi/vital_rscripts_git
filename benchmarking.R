### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Benchmarking #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#' Testing two different versions of a code to see which one is running faster / more efficient
#' The test code below is just setting the system to sleep (simulating some heavy calculation) multiple times once as regular code and once running parallel on multiple cores 

library(parallel) # running functions on multiple cores
library(microbenchmark) # testing efficiency of code


# Version 1
v1 <- expression({
  # Put code version 1 here (below is just some example code to replace with your code)
  n_iter <- 100; for (i in 1:n_iter) {Sys.sleep(0.01)}
})

# Version 2
v2 <- expression({
  # Put code version 2 here replacing the lines below
  n_iter <- 100; num_cores <- detectCores() - 1
  parallel_processing <- function(i) {Sys.sleep(0.01)}
  mclapply(1:n_iter,parallel_processing,mc.cores = num_cores)
})


compare_versions <- function(v1, v2, times = 10) {
  library(microbenchmark)
  
  # benchmarking the codes
  benchmark_results <- microbenchmark(
    Version_1 = eval(v1),
    Version_2 = eval(v2),
    times = times
  )
  
  # get mean execution time & convert time to milliseconds & identify which version is faster or slower
  benchmark_means <- aggregate(time ~ expr, data = benchmark_results, FUN = mean)
  benchmark_means$time_ms <- benchmark_means$time / 1e6
  faster_version <- benchmark_means$expr[which.min(benchmark_means$time)]
  slower_version <- benchmark_means$expr[which.max(benchmark_means$time)]
  speed_ratio <- max(benchmark_means$time) / min(benchmark_means$time)  # Calculate the speed ratio
  
  # Generate and print the message
  message <- sprintf("%s is %.1f times faster than %s.", faster_version, speed_ratio, slower_version)
  cat(blue(message, "\n"))
  
  # Return the benchmarking results
  return(benchmark_results)
}

compare_versions(v1,v2)



