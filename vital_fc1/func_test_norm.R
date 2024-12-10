#### Test Normal distribution of Residuals ####
#' Function to test the whether the distribution of residuals from linear models (lm and lmer) deviates from a normal distribution 
#' Glmer requires other testing as other assumptions are made... depending on the family argument
#' 
#' Literatur on the skewness and kurtosis approach: doi: 10.5395/rde.2013.38.1.52

test_norm <- function(mod) { # function from Nathalie # resids
  cat(blue("Testing normality of residuals for lm or lmer \n(assumptions likely different other types of models)\n"))
  resids <- residuals(object = mod)
  if (length(resids) <= 300) {
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    shapiro_result <- shapiro.test(resids)
    print(shapiro_result)
    p_value <- shapiro_result$p.value
    if(p_value >= 0.05) {
      cat(green("Data is normally distributed\n --> Model is fine... go ahead \U1F44D"))
    } else {
      cat(red("Warning: p value below 0.05\n-->  Data significantly deviate from a normal distribution. \nAdjust or change model!"))
    }
  } else {
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero)")
    skew <- skewness(resids)
    print(paste0("Skewness = ", skew))
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    kurt <- kurtosis(resids)
    print(paste0("Kurtosis = ", kurt))
    if(skew >= -3 & skew <= 3 & kurt < 4){
      cat(green("Model is fine... go ahead \U1F44D"))
    } else {
      cat(red("Warning: Adjust or change model!"))
    }
  }
}
