######################
### Random Numbers ###
######################

# Generate a random number between two values eg. 5.0 and 7.5 (xtremes 5 and 7.5 will not be generated) or ten values between these digits 
x1 <- runif(1, 5.0, 7.5)
x1

x2 <- runif(10, 5.0, 7.5)
x2

# Generate a or multiple random integer between 1 and 10 (replace function states if repeats are allowed)
x3 <- sample(1:10, 1)
x3

x4 <- sample(1:10, 5, replace=T)
x4

# Select 6 or X random numbers between two values eg. 1 and 40, without replacement (random number sequence)
x5 <- sample(1:40, 6, replace=F)
x5

# Select random items from a list eg 10 out of 50, the same function can be used to choos a random subsample of any given vector (you can do the same as above by first creating a vector)
sample(state.name, 10)

x6 <- c(1:60)
sample(x6, 60)

### create a vector and copy  it directly into the clipboard to past it eg into excel

anv1 <- sample(1:60, 60, replace=F)
writeClipboard(as.character(anv1))

anv2 <- sample(1:60, 60, replace=F)
writeClipboard(as.character(anv2))

anv3 <- sample(1:60, 60, replace=F)
writeClipboard(as.character(anv3))

anv4 <- sample(1:60, 60, replace=F)
writeClipboard(as.character(anv4))

anv5 <- sample(1:60, 60, replace=F)
writeClipboard(as.character(anv5))



