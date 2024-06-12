

### ### ### ### ### ### ### ### ### ###
###     Colors in R with Daniel     ###
### ### ### ### ### ### ### ### ### ###



#### Pick a color from a pixcel: ####

### ### ### ### ### ### ### ### ### ###
### https://imagecolorpicker.com/   ###
### ### ### ### ### ### ### ### ### ###  






#### Function: Identify color | get color name ####
get_color <- function(input_color) {
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    install.packages("colorspace")}
  library(colorspace)
  hsb_to_rgb <- function(h, s, b) { # if necessary convert HSB to RGB
    hsb <- hex(HSV(h, s, b))
    t(col2rgb(hsb))}
  # List of named colors and their RGB values
  named_colors <- colors()
  named_rgb <- t(col2rgb(named_colors))
  # get the type of input used (function typically takes hexadecimal code as input but can also handle vectors with rgb or HSB values)
  if (length(input_color) == 1 && is.character(input_color) && grepl("^#", input_color)) {
    input_type <- "hex"
    input_rgb <- t(col2rgb(input_color))
  } else if (length(input_color) == 3 && is.numeric(input_color)) {
    input_type <- "rgb"
    input_rgb <- matrix(input_color, ncol = 3, byrow = TRUE)
  } else if (length(input_color) == 3 && is.numeric(input_color)) {
    input_type <- "hsb"
    input_rgb <- hsb_to_rgb(input_color[1], input_color[2], input_color[3])
  } else {
    stop("Input color should be a hex code, a vector of RGB values, or a vector of HSB values")
  }
  if (is.vector(input_rgb)) {  # Ensure input_rgb is a matrix
    input_rgb <- matrix(input_rgb, ncol = 3, byrow = TRUE)
  }
  
  # Find the nearest color based on eucledian distance
  dist <- apply(named_rgb, 1, function(x) sum((x - input_rgb)^2))
  nearest_index <- which.min(dist)
  nearest_named_color <- named_colors[nearest_index]
  
  # define function to the nearest main color hue
  find_main_color_hue <- function(rgb_color) {
    main_color_hues <- c("Red", "Orange", "Yellow", "Green", "Blue", "Violet", "Brown", "Black", "Grey", "White")
    main_color_ranges <- list(
      Red = list(min = c(160, 0, 0), max = c(255, 100, 100)),
      Orange = list(min = c(220, 45, 0), max = c(255, 150, 50)),
      Yellow = list(min = c(180, 80, 0), max = c(255, 255, 80)),
      Green = list(min = c(0, 90, 0), max = c(120, 255, 120)),
      Blue = list(min = c(0, 0, 150), max = c(100, 100, 255)),
      Violet = list(min = c(150, 0, 150), max = c(255, 80, 255)),
      Brown = list(min = c(100, 50, 0), max = c(200, 100, 50)),
      Black = list(min = c(0, 0, 0), max = c(100, 100, 100)),
      Grey = list(min = c(150, 150, 150), max = c(220, 220, 220)),
      White = list(min = c(220, 220, 220), max = c(255, 255, 255))
    )
    
    # Calculate distance to each main color hue
    distances <- sapply(main_color_ranges, function(range) {
      min_dist_to_range <- apply(cbind(t(rgb_color), range$min, range$max), 1, function(x) {
        min(max(x), max(min(x))) - min(min(x))
      })
      min(min_dist_to_range)
    })
    main_color_hues[which.min(distances)]
  }
  nearest_main_color_hue <- find_main_color_hue(input_rgb)
  
  # Output
  list(
    name = nearest_named_color,
    rgb = input_rgb,
    main_color_hue = nearest_main_color_hue
  )
}

# Example
get_color("#FF5733")







#### Other usefull color stuff... #### 

library(RColorBrewer)
display.brewer.all()
# Get colors from the palettes
pastel1_colors <- brewer.pal(9, "Pastel1")
bluecols <- brewer.pal(9, 'Blues')
print(pastel1_colors)
print(bluecols)

# plot
par(mfrow = c(1, 2))
barplot(1:9, col = pastel1_colors, main = "Pastel1 Colors")
barplot(1:9, col = bluecols, main = "Pastel2 Colors")

# What if you need more colors for some reason? --> use the function colorRampPalette() 
newcol <- colorRampPalette(bluecols)
ncols <- 1000
bluecols2 <- newcol(ncols) #apply the function to get 100 colours
par(mfrow = c(1, 2))
pie(rep(1, length(bluecols)), col = bluecols, labels = bluecols, main = "Blue Colors")
pie(rep(1, length(bluecols2)), col = bluecols2, border = bluecols2, labels = "", main = "Blue Colors")









#### viridis ####
library(viridis)
par(mfrow = c(2, 3))
palett_names <- c("viridis", "magma", "plasma", "inferno", "cividis", "mako") # other options I do not like: turbo and rocket
n_colors <- 100
for (palett_name in palett_names) { # palett_name <- "viridis"
  palett_func <- viridis_pal(option = palett_name)
  colors <- palett_func(n_colors)# Generate sequence of colors
  pie(rep(1, n_colors), col = colors, border = colors, labels = "", main = palett_name)
}

par(mfrow = c(1, 1))
