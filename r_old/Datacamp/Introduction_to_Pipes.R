#### Pipes in R - A tutorial for Beginners %>% %>% %>% ####

x <- c(0.109, 0.359, 0.63, 0.996, 0.515, 0.142, 0.017, 0.829, 0.907)
round(exp(diff(log(x))), 1) # nested code because R is a funcional language. Pipes can be used for easier reading. 

#alternatively with pipes: 
library(magrittr) # contains the pipe operator (but it is also contained in dplyr)

x %>% log() %>%
  diff() %>%
  exp() %>%
  round(1)

# why using pipes: Four reasons why you should be using pipes in R:
#You'll structure the sequence of your data operations from left to right, as apposed to from inside and out;
#You'll avoid nested function calls;
#You'll minimize the need for local variables and function definitions; And
#You'll make it easy to add steps anywhere in the sequence of operations.

### Additional Pipe operators ###

# 1. The compound assignment operator %<>%;
# Initialize `x` 
x <- rnorm(100)
# Update value of `x` and assign it to `x`
x %<>% abs %>% sort

# 2. The tee operator %T>%;
rnorm(200) %>%
  matrix(ncol = 2) %T>%
  plot %>% #short for { plot(.); . }
  colSums
# 3.The exposition pipe operator %$%.
data.frame(z = rnorm(100)) %$% 
  ts.plot(z)

### Basic Piping ###
# f(x) can be rewritten as x %>% f or function(argument) can be rewritten as follows: argument %>% function()
# Compute the logarithm of `x` 
log(x)
# Compute the logarithm of `x` 
x %>% log()

# What if a function needs more than one input variable?
#f(x, y) can be rewritten as x %>% f(y) or argument1 %>% function(argument2) 
# Round pi
round(pi, 6)
# Round pi 
pi %>% round(6)
# x %>% f %>% g %>% h can be rewritten as h(g(f(x))). Another example: 

# Import `babynames` data + Import `dplyr` library
library(babynames)
library(dplyr)
# Load the data
data(babynames)
# Count how many young boys with the name "Taylor" are born
sum(select(filter(babynames,sex=="M",name=="Taylor"),n))
# Do the same but now with `%>%`
babynames%>%filter(sex=="M",name=="Taylor")%>%
  select(n)%>%
  sum

# for some functions pipes do not work. E.g. Functions that use the current environment such as the 
# assignment function. 
# Assign `10` to `x`
assign("x", 10)
# Assign `100` to `x` then Return `x`
"x" %>% assign(100)
x

#The second call with the assign() function, in combination with the pipe, doesn't work.The value of x is not updated.
#The function assigns the new value 100 to a temporary environment used by %>%.
#So, if you want to use assign() with the pipe, you must be explicit about the environment:
# Define your environment
env <- environment()
# Add the environment to `assign()` then Return `x`
"x" %>% assign(100, envir = env)
x

# other functions that do not work: Functions with Lazy Evalution 
#e.g tryCatch(), try(), suppressMessages(), and suppressWarnings() in base R.


### Using Pipes as argument Placeholders ###

#f(x, y) can be rewritten as y %>% f(x, .)
"Ceci n'est pas une pipe" %>% gsub("une", "un", .)

pi %>% round(6) 
6 %>% round(pi, digits=.) #in this version we want the input to be the second argument of the function

### Re-using the placeholder for Attributes: ###
# f(x, y = nrow(x), z = ncol(x)) can be rewritten as x %>% f(y = nrow(.), z = ncol(.))

# Initialize a matrix `ma` 
ma <- matrix(1:12, 3, 4)
# Return the maximum of the values inputted
max(ma, nrow(ma), ncol(ma))
# Return the maximum of the values inputted
ma %>% max(nrow(ma), ncol(ma))


#f(y = nrow(x), z = ncol(x)) can be rewritten as x %>% {f(y = nrow(.), z = ncol(.))}
ma %>% {max(nrow(ma), ncol(ma))}



### Bilding unary functions ### (functions that only take one argument)
# Any pipeline that you might make that consists of a dot ., followed by functions and that is chained 
# together with %>% can be used later if you want to apply it to values. 
# Take a look at the following example of such a pipeline:

  . %>% cos %>% sin

# needs to be assigned to a to a variable to be used. 
# Unary function
f <- . %>% cos() %>% sin() 
f
# is equivalent to the basic ways of defining functions in R
f1 <- function(.) sin(cos(.))
f2 <- function(x) {sin(cos(x))}


### Compound Assignment Pipe Operations ###
#There are situations where you want to overwrite the value of the left-hand side, 
#just like in the example right below. Intuitively, you will use the assignment operator <- to do this.

# Load in the Iris data
iris <- read.csv(url("http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"), header = FALSE)
# Add column names to the Iris data
head(iris)
names(iris) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width", "Species")

# Compute the square root of `iris$Sepal.Length` and assign it to the variable
iris$Sepal.Length <- 
  iris$Sepal.Length %>%
  sqrt()
# However, there is a compound assignment pipe operator, which allows you to use a shorthand notation 
# to assign the result of your pipeline immediately to the left-hand side:
# Compute the square root of `iris$Sepal.Length` and assign it to the variable
iris$Sepal.Length %<>% sqrt
# Return `Sepal.Length`
iris$Sepal.Length

### Tee Operations with the Tee Operations ###
# The tee operator works exactly like %>%, but it returns the left-hand side value rather than 
# the potential result of the right-hand side operations.
# Is handy in situations where you have included functions that are used for their side effect, 
# such as plotting with plot() or printing to a file.
# In other words, functions like plot() typically don't return anything. 
#That means that, after calling plot(), for example, your pipeline would end. 
# However, in the following example, the tee operator %T>% allows you to continue your pipeline 
# even after you have used plot():

set.seed(123)
rnorm(200) %>%
  matrix(ncol = 2) %T>%
  plot %>% 
  colSums


### Exposing Data Variables with the Exposition Operator ###
# When you're working with R, you'll find that many functions take a data argument. 
#Consider, for example, the lm() function or the with() function. 
# These functions are useful in a pipeline where your data is first processed and then passed into the function.
# For functions that don't have a data argument, such as the cor() function, 
# it's still handy if you can expose the variables in the data. 
# That's where the %$% operator comes in. Consider the following example:

iris %>%
  subset(Sepal.Length > mean(Sepal.Length)) %$%
  cor(Sepal.Length, Sepal.Width)
0.336696922252551

# With the help of %$% you make sure that Sepal.Length and Sepal.Width are exposed to cor(). 
# Likewise, you see that the data in the data.frame() function is passed to the ts.plot() to plot several 
# time series on a common plot:

data.frame(z = rnorm(100)) %$%
  ts.plot(z)



### dplyr and magrittr ###
# this R package was built around five verbs, namely, "select", "filter", "arrange", "mutate" and "summarize". 
# which is most of what you want to do when you perform data manipulation 

#traditional example:
install.packages("hflights")
library(hflights)
grouped_flights <- group_by(hflights, Year, Month, DayofMonth)
flights_data <- select(grouped_flights, Year:DayofMonth, ArrDelay, DepDelay)
summarized_flights <- summarise(flights_data, 
                                arr = mean(ArrDelay, na.rm = TRUE), 
                                dep = mean(DepDelay, na.rm = TRUE))
final_result <- filter(summarized_flights, arr > 30 | dep > 30)
final_result

#vs the much clearer version of code
hflights %>% 
  group_by(Year, Month, DayofMonth) %>% 
  select(Year:DayofMonth, ArrDelay, DepDelay) %>% 
  summarise(arr = mean(ArrDelay, na.rm = TRUE), dep = mean(DepDelay, na.rm = TRUE)) %>% 
  filter(arr > 30 | dep > 30)

#In short, dplyr and magrittr are your dreamteam for manipulating data in R!
# you might want to insert keyboard shortcuts --> see https://rstudio.github.io/rstudioaddins/ for more informations


### When not to use Pipes or be careful ###

#Do not make them to long (eg. more than 10 steps)
# avoid them if you have multiple inputs and outputs
# if You are starting to think about a directed graph with a complex dependency structure.
# if You're doing internal package development