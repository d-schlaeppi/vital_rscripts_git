
#### create a vector with 9 colors as far apart as possible ####

library(viridis)
# Create a vector with 7 colors from viridis palette then add black and something close to white. 
colors <- viridis(7, begin = 0.1, end = 1, alpha = 1, option = "G")
colors <- c("#000000", colors, "#F7FDF9")
barplot(rep(1, 9), col = colors, space = 0, border = "black", lwd = 2, ylim = c(0, 1), yaxt = "n", names.arg = colors, cex.names = 0.8, las = 2)


#### find intermediate color ####
# Function to find intermediate color between two colors
intermediate_color <- function(color1, color2) {
  # Convert hexadecimal color codes to RGB values
  rgb1 <- col2rgb(color1)
  rgb2 <- col2rgb(color2)
  # Calculate the average RGB values between the two colors
  avg_rgb <- (rgb1 + rgb2) / 2
  # Convert RGB values to hexadecimal color code
  hex_color <- rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
  # Return the intermediate color as a hexadecimal color code
  return(hex_color)
}

# Find the intermediate color between two colors
color1 <- "#DEF5E5FF"
color2 <- "#FFFFFF"
intermediate <- intermediate_color(color1, color2)

# Print the intermediate color as a hexadecimal color code
intermediate
plot(1:10, type = "n", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", main = "Intermediate Color", xaxt='n', yaxt='n')
rect(xleft = 0, ybottom = 0, xright = 10, ytop = 10, col = intermediate, border = NA)



