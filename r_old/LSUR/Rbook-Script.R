##########################
#### Learning R again ####
##########################

#### Introduction and first chapters ####

#first it is good to set the working directory
getwd()
setwd("/Users/gismo/Documents/DS-Bees/R/LSUR")
setwd("~/Documents/DS-Bees/R/LSUR")                    # ~ steht f�r standardpfad Users/gismo in this case
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/LSUR")
getwd()

# 1. Einem Objekt einen Wert zufügen

metallica <- c("Lars", "James", "Jason", "Kirk")

# Eckige Klammer auf Mac = alt + 5 und alt + 6 
# geschweifte Klammer auf Mac = alt + 8 und alt + 9 
# Mit einem Semikolon kann man zwei Befehle auf eine Zeile schreiben
# einen Wert aus einem Vektor löschen mit != 

metallica <- metallica[metallica != "Jason"]; metallica <- c(metallica, "Rob")     # in dieser Zeile wird in zwei Befehlen Jason durch Rob ersetzt
metallica


install.packages("DSUR")     # install a package (only once needed)
library(DSUR)                 # load a package (every timem needed aftere a restart if you want R to refer to it)

#if two packages have a different function saved under the same name you need to infer this to R. with  package::function()
car::recode()
Hmisc::recode()

#if you need help Google is your best friend other sources of information are: 
help(c)      # help(function)
?c           # ?function

# Quotes "" tell R that something is not numeric (string variables vs. numeric variables) eg.:
metallica.names <- c("Lars", "James", "Kirk", "Rob")
metallica.ages <- c(47,47,48,46)

metallica <- data.frame(name = metallica.names, age = metallica.ages)         #combining things in a dataframe: R creates a new object with called name with the values of metallica.names
                                                                              #and does the same for age. finally the two are combined in a dataframe called metallica

#refering to a variable with $ in the dataframe: dataframe$variablename
metallica$age

# creating a new variable straightforeward in a dataframe by creating a new object (vector)
metallica$child.ages <-  c(12,12,4,6)
names(metallica)         #gives us the the names of the variables in the dataframe

### 
#Data frames are very versatile and can combine different types of variables. But some functions need a list (list()) of objects and to combine
#numerical values the functions cbind() and rbind() can be used which combine rows or columns (cbind und rbind können nur einen variabel typ verwerten -->
#numerische daten werden zu Text, wenn mit variablen wie "namen" kombiniert)

metallicalist <- list(metallica.names, metallica.ages) 
metallicalist

metallicacolumn <- cbind(metallica.names, metallica.ages)
metallicacolumn

metallicarow <- rbind(metallica.names, metallica.ages)
metallicarow

###
# usefull operators in R
# + - / * normal calculations
# ^ exponential 
# < > >= <= 
# ==  bedeutet Exactly equal to!
# !=  not equal to
# !x  not x
# x|y x oder y                senkrechter strich auf mac ist alt+7 
# x&y e.g. name=="James" & age==47 means that variable age is equal to 47 and variable name is equal to James
# isTrue(x)     Test if x is true

metallica$fatherhood.age <- metallica$age - metallica$child.ages
metallica

# creating a strin Variable with teacher and student names
names <- c("Ben", "Martin", "Andy", "Paul", "Graham", "Carina", "Karina", "Dough", "Mark", "Zoe")

###Dates
# enter dates in R: enter date as string variable with a specific format (yyyy-mm-dd) and then tell R to recognize it as dates: 
# But when dates are entered the can be used for calculations:
husband <- c("1973-06-12", "1970-07-16","1948-11-11", "1969-05-24")
wife <- c("1984-11-12", "1973-08-02", "1948-11-11", "1983-07-23")
age.gap <- husband-wife     # ergibt eine fehlermeldung --> da nicht als datum erkannt
husband <- as.Date(husband)
wife <- as.Date(wife)
age.gap  <- husband-wife
age.gap                     #jetzt funktioniert die Rechnung und die Altersdifferenz wird in Jahren ausgespukt. 
# as.Date() wird ausserhalb der verktorfunktion c() geschrieben um einen stringvektor mit Daten zu erhalten 
birth.date <- as.Date(c("1977-07-03", "1969-05-24", "1973-06-21", "1970-07-16", "1949-10-10", "1983-11-05", "1987-10-08", "1989-09-16", "1973-05-20", "1984-11-12"))

###
# Creating coding variables (eg. Geschlecht --> 1 or 0, treatment vs control) to represent different groups
# Create a vector an tell R that it is a coding variable/factor: (eg. for the example of 5 teachers and 5 students)
job <- c(1,1,1,1,1,2,2,2,2,2)
# you can do this faster using the repetition function rep()
job <- c(rep(1,5), rep(2,5))
# this vector is now transformed the following general form: 
# factor(variable, levels = c(x,y, … z), labels = c("label1", "label2", … "label3"))“
job <- factor(job, levels = c(1:2), labels = c("Lecturer", "Student"))    #factor() will transform it into a nominal vector this can also be acheive with gl()
                                                                          # newFactor<-gl(number of levels, cases in each level, total cases, labels = c("label1", "label2"…))“
job <- gl(2,5, labels = c("Lecturer", "Student"))
levels(job)                                       #to see what levels a factor has use the function levels(variable)
levels(job) <- c("Medical lecturer", "Medical student")    # the levels can be modified 

friends <- c(5,2,0,4,1,10,12,15,12,17)
alcohol <- c(10,15,20,5,30,25,20,16,17,18)
income <- c(20000, 40000, 35000, 22000, 50000, 5000, 100, 3000, 10000, 10)
neurotic <- c(10, 17, 14, 13, 21, 7, 13, 9, 14, 13)

lecturer.data <- data.frame(name=names, birth.date=birth.date, job=job, nr.friends=friends,alcohol=alcohol, income=income, neuroticism=neurotic)
lecturer.data
print(lecturer.data)
neurotic <- c(10,17,NA,13,21,7,13,9,14,NA)    #missing values are  referred to as not available called NA in R

### Entering Data with r-Commander
install.packages("Rcmdr", dependencies=TRUE)
library(Rcmdr)

#load data into R if workingdirectory is set and saved as Textfile (Tabstobbsgetrennt)
# header = TRUE tells R that the first row are the variablenames
# use read.delim for .dat files or read.csv for csv.files
data <- read.table("filename.txt", header = TRUE)

#if you do not know the path or the location of the file use the function file.choose() which will open a dialogbox, the funciton can be entered in read table:
lecturer.data  <- read.table(file.choose(), header=TRUE)

#to save data as tab-delimited files use the command: write.table(), sep="\t" tells are to use tab as separator as needed!
write.table(metallica, "metallica data.txt", sep="\t", row.names = FALSE)



### manipulating the dataframes: 

#create new dataframe with only a few values using the basic command „newDataframe <- oldDataframe[rows, columns]“
lecturer.personality  <- lecturer.data[,c("nr.friends", "alcohol", "neuroticism")]
lecturer.personality

lecturer.only  <- lecturer.data[job=="Medical lecturer",]
lecturer.only

alcohol.personality  <- lecturer.data[alcohol > 10, c("nr.friends", "alcohol", "neuroticism")]
alcohol.personality

#selecting data with the subset() function: „newDataframe<-subset(oldDataframe, cases to retain, select = c(list of variables))“
lecturer.only2  <- subset(lecturer.data, job=="Medical lecturer" )
lecturer.only2

aclohol.personality2 <- subset(lecturer.data, alcohol > 10, select = c("nr.friends", "alcohol", "neuroticism"))
aclohol.personality2
all of the variables for those who drink 20 units or more or have a neuroticism score greater than 14.“

new.dataframe1  <- subset(lecturer.data, income >= 10000, select = c("name", "income", "job"))
new.dataframe2  <- subset(lecturer.data, alcohol <= 12, select = c("name", "job", "income", "nr.friends"))
new.datafrae

names(lecturer.data)

#some functions in R need matrixes instead of data frames, but matrixes can only contain numeric variables "newmatrix <- as.matrix(dataframe)"
alcohol.personality.matrix <- as.matrix(alcohol.personality)
#or directly from the original dataframe, by choosing the relevant numeric variables
alcohol.personality.matrix <- as.matrix(lecturer.data[alcohol > 10, c("nr.friends", "alcohol", "neuroticism")])


### Transforming between the long and the wide format:
# melt(), cast(), or for easy data sets stack() and unstack()

satisfaction <- read.delim("Honeymoon Period.dat", header=TRUE)
satisfaction

# stack function: „newDataFrame<-stack(oldDataFrame, select = c(variable list))“ creates a new Dataframe with the variables you want stacked
satisfaction.stacked <- stack(satisfaction, select = c("Satisfaction_Base","Satisfaction_6_Months","Satisfaction_12_Months","Satisfaction_18_Months"))
names(satisfaction.stacked) 

#unstack function: „newDataFrame<-unstack(oldDataFrame, scores ~ columns)“ with scores containing the variables (value) and column which specifies
#to which variable the value belonges to (ind) 

satisfaction.unstacked <- unstack(satisfaction.stacked, values ~ ind)
satisfaction.unstacked

# for more complex transformations use melt() and cast() which use the package reshape:
install.packages("reshape")
library(reshape)

# use melt() um die Daten aus dem weiten Format zu zusammenzuschmelzen in ein long format: „newDataFrame<-melt(oldDataFrame, id = c(constant variables), measured = c(variables that change across columns))“
# id: any variables that do not vary over time (e.g. gender or number identifiert...)
# measured: variables that do vary over time or are repeated measures (e.g scores).--> specify the names of variables that should be restructured in a new row!

restructured.data <- melt(satisfaction, id = c("Person","Gender"), measured = c("Satisfaction_Base","Satisfaction_6_Months","Satisfaction_12_Months","Satisfaction_18_Months"))
names(restructured.data)

# to cast data into an new dataframe from long to wide form use: „newData<-cast(moltenData, variables coded within a single column ~ variables coded across many columns, value = "outcome variable")“
# left: any variable that does not vary within an entity (previously used in id) like name or gender --> coding variables in the wide format
# right: anything that changes within the entity: variables that were set in measured in the functino melt (the variables that should be spread across multiple columns in the wide format) 
# value, enables you to specify a variable in the molten data that contains the actual scores 

satisfaction.wide.again <- cast(restructured.data, Person + Gender ~ variable, value="value" )
satisfaction.wide.again

###
#practice tasks
gender <- c(rep(1,10), rep(2,10))
gender  <- factor(gender, levels = c(1,2), labels = c("male", "female"))
treatment <- c(rep(1,5), rep(2,5), rep(1,5), rep(2,5))
treatment <- factor(treatment, levels = c(1,2), labels = c("shock", "nice"))
score <- c(15,14,20,13,13,10,9,8,8,7,6,7,5,4,8,12,10,7,8,13)
teaching  <- data.frame (gender = gender, treatment = treatment, score = score)
teaching

gender <- c(rep(1,24), rep(2,24))
gender  <- factor(gender, levels = c(1,2), labels = c("male", "female"))
treatment <- c(rep(1,12), rep(2,12), rep(1,12), rep(2,12))
treatment <- factor(treatment, levels = c(1,2), labels = c("partners face", "own face"))
number.of.bullets <- c(69,76,70,76,72,65,82,71,71,75,52,34,33,26,10,51,34,28,27,9,33,11,14,46,70,74,64,43,51,93,48,51,74,73,41,84,97,80,88,100,100,58,95,83,97,89,69,82)
infidelity  <-  data.frame(gender = gender, treatment=treatment, number.of.bullets = number.of.bullets)
write.csv(infidelity, "Infidelity.csv")



#### Chapter 4 - Discovering Data using Graphs ####


# aim: create simple plots with the best information / ink relationship

install.packages("ggplot2")
library(ggplot2)

#two ways of creating grpahs in ggplot 2: qplot() for quick easy plots and ggplot() for crazy graph shenaningans with great versatility
#ggplot will create plots using different layers with geometric objects (geoms) which can be aestetically modified (aes)

####possible geoms all creating a new layer:
#geom_bar() - layer with bars representic statistical properties
#geom_points() - layer showing data poits
#geom_lines() - layer connecting datapoints with a straight line
#geom_smooth() - layer with a smooth line summarizing the data
#geom_hist() - layer with a histogram
#geom_boxplot() - layer with a box whisker diagram
#geom_text() - layer with text
#geom_density() - layer with density plot
#geom_errorbar() - layer displaying errorbars
#geom_hline() / geom_vline() - layers with custom horizontal and vertical lines

# in the () of each geom the look is defined with some required and some optional astetics, eg ggplot needs to know which data to plot, but it is not necessary to define colours

### aestetics
# con be defined for the whole plot or just for specific layers
# possible aestetics: linetype, size, shape, color, alpha (making colours transparent)...

### anatomy of a plot
#1: create object that specifies your plot, set some basic aestetics -> which data to plot
# general version: myGraph <- ggplot(myData, aes(variable for x axis, variable for y axis))“
# at this toplevel additional aestetics like colour can already be set: myGraph <- ggplot(myData, aes(variable for x axis, variable for y axis, colour = gender))

# at this level also some options can be set for example the title: 
# + opts(title="Title)

#2: Now different layers (geoms) need to be added using the + symbol, after that the graph will contain graphic elements which will be plotted!
# mygraph + geom_bar() + geom points

#Withing the new layers aestetics can be defined which override the general settings
#myGraph <- ggplot(myData, aes(variable for x axis, „variable for y axis, colour = gender))
#myGraph + geom_bar() + geom_point(„shape = 17, colour = "Blue") + labels(x = "Text", y = "Text")   

# when you plot statistical graphs like boxplots ggplot needs to know the values for min, max, mean etz. which it mostly can get from the raw data automatically
# to avoid overplotting the function "position" can be used: position = "x" with x beeing one of the following five: 
# dodge - positions objects with no overlap at the side
# stack and fill. positions objects stacked with the largest in the back
# identity - no position adjustment
# jitter -  adds random offset to objects so they do not overlap

### other useful tools are ro create subplots using the commands: facet.grid() or facet.warp()
# general form: 
# + facet_warp ( ~ y, nrow = integer, ncol = integer)
# + facet_grip(x ~ y)
# x and y are the variables to plot, nrow and ncol are optional and give the numer of rows and columns where plots are putted

# saving graphs (not using the R-Studio help) but with a simple command: ggsave(filename)
# can save as tiff, jpeg, png, pdf and many more   eg. ggsave("Outlier Amazon.tiff")
# futher some additional seetings can be added: ggsave("Outlier Amazon.tiff", width = 2, height = 2) --> two inches high and two inches wide

###ggplot tutorial with practical stuff: 
getwd()
facebook.data <- read.delim("FacebookNarcissism.dat", header=TRUE)
facebook.data
names(facebook.data)

####
## 4.4 Plotting Graphic Realtionships

#scatterplots can be used to find outliers
#scatterplot of the example test-anxiety and test-result

exam.data <- read.delim("Exam Anxiety.dat", header=TRUE)
exam.data

#first level create plot calld scatter, and add aestetics with x und y achse 
scatter <- ggplot(exam.data, aes(Anxiety, Exam))
#now we ned to visualize it and add some dotts and nice labels
scatter 
scatter + geom_point()
scatter + geom_point() + labs(x="Exam Anxiety Score", y = "Exam Performance %")

#to add a regression line we add a "smoother" which smoothes out the raw data: -> gives a nonlinear curve with standard error shade
scatter + geom_point()  + geom_smooth() + labs(x="Exam Anxiety Score", y = "Exam Performance %")
#to add a straight line we change the method of the smoother to lm (linear model) or rlm (linear model which is less affected by outliers)
#add colour to line
scatter + geom_point()  + geom_smooth(method="lm", colour="red") + labs(x="Exam Anxiety Score", y = "Exam Performance %")
#standard error shade can be removed using se=F or given an other color using the aestetics "fill" and "alpha" (-->modifizieren Farbe und transparenz)
#total plot in one code
scatter <- ggplot(exam.data, aes(Anxiety, Exam))
scatter + geom_point()  + geom_smooth(method="lm", colour="red", se=F) + labs(x="Exam Anxiety Score", y = "Exam Performance %")
scatter + geom_point()  + geom_smooth(method="lm", colour="red", alpha=0.5, fill="yellow") + labs(x="Exam Anxiety Score", y = "Exam Performance %")

###
#If we want the scatterplot to differentiate betwee different groups we add gender as an asthetic
scatter <- ggplot(exam.data, aes(Anxiety, Exam, colour = Gender))
scatter + geom_point() + geom_smooth(method = "lm", aes(fill=Gender), alpha=0.2 ) + labs(x="Exam Anxiety Score", y = "Exam Performance %", colour=Gender)

###
# 4.6 Histogramms, a good way to spot problems
festi <- read.delim("DownloadFestival.dat", header=TRUE)
festi
names(festi)
summary(festi)

#first create object/plot, setting aes to the vriable we want to look at and specify the legend to be absent:
festi.histogram <- ggplot(festi, aes(day1), opts(legend.position="none"))
#now we add the histogramm layer (geom) to that plot:
festi.histogram + geom_histogram() 
#play aroung: change binwith, add labels
festi.histogram + geom_histogram(binwidth=0.4) + labs(y = "frequency", x = "hygiene (day 1)")
#problem: there is an outlier --> influences mean, blows up standarterror and makes graph ulgy! (regocnize ouliers: histogramm, boxplots or z-scores)

####
# 4.7 Boxplots (Box-Whiskers-Diagramme)
festi.boxplot <- ggplot(festi, aes(gender, day1))
festi.boxplot + geom_boxplot() + labs(x = "gender", y = "hygiene day 1")

#outlier is again the problem --> needs to be removed 
festi <- festi[order(festi$day1),]

festi.outlinefree <- festi[festi$day1 != 20.02,]
festi.outlinefree
festi.histogram.outlinefree <- ggplot(festi.outlinefree, aes(day1), opts(legend.position="none"))
festi.histogram.outlinefree + geom_histogram(binwidth=0.4) + labs(y = "frequency", x = "hygiene (day 1)")

festi.boxplot <- ggplot(festi.outlinefree, aes(gender, day1))
festi.boxplot + geom_boxplot() + labs(x = "gender", y = "hygiene day 1")


###
# 4.8 Density plots
#similar to histogram but with a smooth line rather than a bar: 

density <- ggplot(festi.outlinefree, aes(day1))
density + geom_density() + labs(x="Hygiene (Day 1 of Festival)", y="Density Estimates")


####
# 4.8 Graphic Means with errorbars

chickflick <- read.delim("ChickFlick.dat", header=TRUE)
bar <- ggplot(chickflick, aes(film, arousal))

#to get the mean we have to apply the mean function taking the following form: 
# stat_summary(function = x , geom = y)           (fill=farbe der Bars, colour,= Umrandung)
# fun.y =  mean für Mittelwerte
# errorbars: „stat_summary(fun.data = mean_cl_normal, geom = "pointrange")“   --> standard 95 % confidence intervall

bar + stat_summary(fun.y = mean, geom = "bar", fill = "White", colour = "Black") + stat_summary(fun.data = mean_cl_normal, geom = "pointrange") + labs(x="Film", y="Mean Arousal")

####
# 4.9 Bar charts for different independant variables
# separate the variables in the same graph         (other versions: set aestetigs with colour on gender or use faceting to create to different grpahs)

bar <- ggplot(chickflick, aes(film, arousal, fill=gender))    #fill=gender to have different colours for the two genders
#again add the means with stat_summary, but use position="dodge" to have the male and female bars be side by side instead of behind each other
# then as new layers add again error bars ans labels

bar + stat_summary(fun.y = mean, geom="bar", position = "dodge") + 
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", position=position_dodge(width=0.90), width = 0.2) + ### width =0.2 errorbar balken 20% of bar width
  labs(x="Film" , y="Mean Arousal", fill="Gender") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_text(color="yellow"))


#different plots for the variables males and females
bar <- ggplot(chickflick, aes(film, arousal, fill = film))
bar + stat_summary(fun.y = mean, geom = "bar") + 
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2) + facet_wrap(~gender) +
  labs(x="Film", y="Mean Arousal") +
  theme(axis.text=element_text(size=12),
      axis.title=element_text(size=16, face="bold"), 
      axis.title.x=element_text(color="yellow")) 

#to manually foll colours other than the default ones use scale_fill_manual
+ scale_fill_manual("Gender", c("Female" = "Blue", "Male" = "Green"))

###
# Line Graphs of single independant values
hiccups  <- read.delim("Hiccups.dat", header=TRUE)
names(hiccups)

hiccup <- stack(hiccups)
names(hiccup)  <-  c("Hiccups", "Intervention")
#for ggplot Intervention needs to be recognized as a factor
hiccup$intervention_factor  <- factor(hiccup$Intervention, levels = hiccup$Intervention)

line <- ggplot(hiccup, aes(intervention_factor, Hiccups))
line + stat_summary(fun.y = mean, geom = "point") + 
  stat_summary(fun.y = mean, geom = "line", aes(group = 1), colour="Blue", linetype = "dashed") +   # use geom line to plot a line and not points and define grout unter aesthetics to connect the means with a line
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2)   +                         #add errorbar using bootstrapping
  labs(x = "Intervention", y = "Mean Number of Hiccups") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

####
# line graph for several independant variables
text.data  <- read.delim("TextMessages.dat", header=TRUE)
text.data
names(text.data)

#transform the datafile to wideformat 
text.data$id = row(text.data[1])                                                                    # create a new variable to have a unique identifier for the reshape command
text.messages  <-  melt(text.data, id = c("id", "Group"), measured = c("Baseline", "Six_months"))   # then we create a new dataframe where we melt together the measured Points and keep the id ones
names(text.messages) <- c("id", "group", "time", "score")                                           # rename the variables 
text.messages$time  <- factor(text.messages$time, labels = c("Baseline", "6 Months"))               # transform time to a factor
print (text.messages)

#create Graph line with the dataframe text.messages, seeting aestethcs with time on x axis and score on y axis, and coulour set for groups giving different colours for the treatments
line  <- ggplot(text.messages, aes(time, score, colour = group)) 
line + stat_summary(fun.y = mean, geom = "point") + stat_summary(fun.y = mean, geom = "line", aes(group = group)) + # display means with sybmols and create lines with colours/grouping acording to variable group
 stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + labs(x="Time", y="Mean Grammar Score", colour="group") #add errorbars and labels

#the same again but with symbols and lines different for treatmentgroups
line + stat_summary(fun.y = mean, geom = "point", aes(shape=group), size=4) + 
  stat_summary(fun.y = mean, geom = "line", aes(group=group, linetype=group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + labs(x="Time", y="Mean Grammar Score", colour="group")
  
###
## Colors and Picture Options:

# the default theme is grey and there is a black white theme theme_bw() to set global options or single properties like colours fonts shapes of
# axes and titels grid lines and so on: 

# theme elements: 
line  all line elements (element_line)
rect	all rectangluar elements (element_rect)
text	all text elements (element_text)
title	all title elements: plot, axes, legends (element_text; inherits from text)
axis.title

label of axes (element_text; inherits from text)
axis.title.x	x axis label (element_text; inherits from axis.title)
axis.title.y	y axis label (element_text; inherits from axis.title)
axis.text	tick labels along axes (element_text; inherits from text)
axis.text.x	x axis tick labels (element_text; inherits from axis.text)
axis.text.y	y axis tick labels (element_text; inherits from axis.text)
axis.ticks	tick marks along axes (element_line; inherits from line)
axis.ticks.x	x axis tick marks (element_line; inherits from axis.ticks)
axis.ticks.y	y axis tick marks (element_line; inherits from axis.ticks)
axis.ticks.length	length of tick marks (unit)
axis.ticks.margin	space between tick mark and tick label (unit)
axis.line	lines along axes (element_line; inherits from line)
axis.line.x	line along x axis (element_line; inherits from axis.line)
axis.line.y	line along y axis (element_line; inherits from axis.line)
legend.background

background of legend (element_rect; inherits from rect)
legend.margin	extra space added around legend (unit)
legend.key	background underneath legend keys (element_rect; inherits from rect)
legend.key.size	size of legend keys (unit; inherits from legend.key.size)
legend.key.height	key background height (unit; inherits from legend.key.size)
legend.key.width	key background width (unit; inherits from legend.key.size)
legend.text	legend item labels (element_text; inherits from text)
legend.text.align	alignment of legend labels (number from 0 (left) to 1 (right))
legend.title	title of legend (element_text; inherits from title)
legend.title.align	alignment of legend title (number from 0 (left) to 1 (right))
legend.position	the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
legend.direction	layout of items in legends ("horizontal" or "vertical")
legend.justification	anchor point for positioning legend inside plot ("center" or two-element numeric vector)
legend.box	arrangement of multiple legends ("horizontal" or "vertical")
panel.background

background of plotting area (element_rect; inherits from rect)
panel.border	border around plotting area (element_rect; inherits from rect)
panel.margin	margin around facet panels (unit)
panel.grid	grid lines (element_line; inherits from line)
panel.grid.major	major grid lines (element_line; inherits from panel.grid)
panel.grid.minor	minor grid lines (element_line; inherits from panel.grid)
panel.grid.major.x	vertical major grid lines (element_line; inherits from panel.grid.major)
panel.grid.major.y	horizontal major grid lines (element_line; inherits from panel.grid.major)
panel.grid.minor.x	vertical minor grid lines (element_line; inherits from panel.grid.minor)
panel.grid.minor.y	horizontal minor grid lines (element_line; inherits from panel.grid.minor)
plot.background

background of the entire plot (element_rect; inherits from rect)
plot.title	plot title (text appearance) (element_text; inherits from title)
plot.margin	margin around entire plot (unit)
strip.background

background of facet labels (element_rect; inherits from rect)
strip.text	facet labels (element_text; inherits from text)
strip.text.x	facet labels along horizontal direction (element_text; inherits from strip.text)
strip.text.y	facet labels along vertical direction (element_text; inherits from strip.text)

### some playing around with the previous graph
line + theme_bw() + stat_summary(fun.y = mean, geom = "point", aes(shape=group), size=4) + 
  stat_summary(fun.y = mean, geom = "line", aes(group=group, linetype=group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + labs(x="Time", y="Mean Grammar Score", colour="group") +
  theme(panel.background = element_rect(colour = "yellow")) +
  labs(title = "whatever") +
  theme(plot.title = element_text(size = rel(2), colour = "blue"))

####
# Hometasks




####Chapter 5 - Exploring Assumptions   #####


# Drawing conclusions about data is based on some assumptions (Example of Doorbell that is allways assumed to be working). If Assumptions are violated then
# results are likely to become inacurate. If we use parametric tests they all are based on certain assumptions. 

###
# assumptions of parametric tests are: 
#
# 1) Normally distributed data - either the data or the errors in a model need to be normally distributed. 
# 2) Homogeneity of variance - This assumptions means that the variances should be tha same throughout the data. 
# 3) Interval data - Data should be measured at least at the interval level. 
# 4) Independance - Depends on the used test and means basically that measurements need to be independet from each other at a certain level.
#   in regression the erreors of the regressoin model need to be independant. 

#Packages for this chapter: 

install.packages("car")
install.packages("ggplot2")
install.packages("pastecs")
install.packages("psych")
install.packages("SparseM")
install.packages("Rcmdr")

library(car)
library(ggplot2)
library(pastecs)
library(psych)
library(Rcmdr)


# load the festival data and plot histograms with normaldistribution curve and look at qqplots to check for normality

dlf <- read.delim("DownloadFestival.dat", header=TRUE)
str(dlf)


dlf$day1 <- ifelse(dlf$day1 > 20, NA, dlf$day1)    #remove outlier

hist.day1 <- ggplot(dlf, aes(day1)) + theme(legend.position="none") + geom_histogram(aes(y=..density..), colour="black", fill="white") + labs(x="Hygiene score on day 1", y = "Density") 
hist.day1
hist.day1 + stat_function(fun = dnorm, args = list(mean(dlf$day1, na.rm = TRUE), sd = sd(dlf$day1, na.rm = TRUE)) , colour= "black", size = 1)

hist.day2 <- ggplot(dlf, aes(day2)) + theme(legend.position="none") + geom_histogram(aes(y=..density..), colour="black", fill="white") + labs(x="Hygiene score on day 1", y = "Density")
hist.day2
hist.day2 + stat_function(fun = dnorm, args = list(mean(dlf$day2, na.rm = TRUE), sd = sd(dlf$day2, na.rm = TRUE)) , colour= "black", size = 1)

hist.day3 <- ggplot(dlf, aes(day3)) + theme(legend.position="none") + geom_histogram(aes(y=..density..), colour="black", fill="white") + labs(x="Hygiene score on day 1", y = "Density")
hist.day3
hist.day3 + stat_function(fun = dnorm, args = list(mean(dlf$day3, na.rm = TRUE), sd = sd(dlf$day3, na.rm = TRUE)) , colour= "black", size = 1)

qqplot.day1 <- qplot(sample= dlf$day1, stat="qq") #oldschool way of doing it
qqplot.day1

pd1 <- ggplot(dlf, aes(sample = day1))              #new way of doing it!
pd1 + stat_qq()

pd2 <- ggplot(dlf, aes(sample = day2))            
pd2 + stat_qq()

pd3 <- ggplot(dlf, aes(sample = day3))             
pd3 + stat_qq()

### only the data from day 1 looks normally distributed the rest not really / but how much do they deviate from normallity?

###
# 5.5.2 Quantifying Normality with numbers
# we can use describe() from the package psych

describe(dlf$day1)
describe(dlf$day2)
describe(dlf$day3)

#or we use stat.desc(variable name, basic = TRUE, norm = FALSE) from the package pastec  --> basic stands for basic statistics, norm stands for normality statistics
# just choose between basic or norm

stat.desc(dlf$day1, basic = FALSE, norm = TRUE)
stat.desc(dlf$day2, basic = FALSE, norm = TRUE)
stat.desc(dlf$day3, basic = FALSE, norm = TRUE)

# use these functions for all variables at the same time

describe(cbind(dlf$day1, dlf$day2, dlf$day3))                               # cbind can be used to combine several variables by column
stat.desc(cbind(dlf$day1, dlf$day2, dlf$day3), basic = FALSE, norm = TRUE)

describe(dlf[, c("day1", "day2", "day3")])                                # we can also use [,] to choose [row, columns] --> wo chose all rows 
                                                                          #(nothing specified before the comma and selected the columns day 1-3
stat.desc(dlf[, c("day1", "day2", "day3")], basic = FALSE, norm = TRUE)

# under normal distribution the values for skew and kurtosis shuld be zero! Positive skew values indicate a pile of scores on the left of the distribution, 
# whereas negative ones indicate a pile to the right. Kurtosis: positive --> pointy and heavytailed curve, negative --> flat and light talied distribution
# the values are informative but need to be converted to a z-score (a score from a distribution with mean of 0 and SD = 1) then we can compare the 
# likelyhood a given skew or kurtosis is to come up.  

# an absolute  value of grater than 1.96 indicats significane at p < 0.05 , 2.58 --> p < 0.01, 3.29 --> 0.001
# In other words, if skew.2SE or kurt.2SE are greater than 1 (ignoring the plus or minus sign) then you have significant skew/kurtosis (at p < .05)
# in large samplesizes the visual check is better than testing, testing can be used in small samplesizes!

# we can use the functin round() to change the number of displayed decimal places and the scientific notation  to have better looking outputs
# round(object that we want to round, digits = x)

r <- stat.desc(dlf[, c("day1", "day2", "day3")], basic = FALSE, norm = TRUE)
round(r, digits = 3)






###5.5.3 Exploring groups of data
# often data can be split into different groups

rexam <-- read.delim("rexam.dat", header= TRUE)
head(rexam)  #look at the top 6 rows of the dataset
rexam$uni
rexam$uni <- factor(rexam$uni, levels = c(0,-1), labels = c("Duncetown University", "Sussex University")) #transform the variable uni from numbers to a factor representing the two universities
summary(rexam$uni)


#looking at the basic descriptive statistics of this new dataset
stat.desc(rexam, basic = FALSE, norm = TRUE)
par(mfrow=c(2,2))
a <- ggplot(rexam, aes(exam)) + theme(legend.position="none") + geom_histogram(aes(y=..density..), colour="black", fill="white")
a + stat_function(fun = dnorm, args = list(mean(rexam$exam, na.rm = TRUE), sd = sd(rexam$exam, na.rm= TRUE)), colour = "black", size = 1)
b <- ggplot(rexam, aes(computer)) + theme(legend.position="none") + geom_histogram(aes(y=..density..), colour="black", fill="white")
b + stat_function(fun = dnorm, args = list(mean(rexam$computer, na.rm = TRUE), sd = sd(rexam$computer, na.rm= TRUE)), colour = "black", size = 1)
c <- ggplot(rexam, aes(lectures)) + theme(legend.position="none") + geom_histogram(aes(y=..density..), colour="black", fill="white")
c + stat_function(fun = dnorm, args = list(mean(rexam$lectures, na.rm = TRUE), sd = sd(rexam$lectures, na.rm= TRUE)), colour = "black", size = 1)
d <- ggplot(rexam, aes(numeracy)) + theme(legend.position="none") + geom_histogram(aes(y=..density..), colour="black", fill="white")
d + stat_function(fun = dnorm, args = list(mean(rexam$numeracy, na.rm = TRUE), sd = sd(rexam$numeracy, na.rm= TRUE)), colour = "black", size = 1)

# 5.5.3.2 now we want to split the set into different groups: 

# the by function takes the general form: 
# by(data = dataFrame, INDICES = grouping variable, FUN = a function that you want to apply to the data)
# eg. 

by(data = rexam$exam, INDICES = rexam$uni, FUN = describe)
by(data = rexam$exam, INDICES = rexam$uni, FUN = stat.desc)

#if the order of the arguments are kept in order we do not have to write it all out and can just give the arguments without the description: 
by(rexam$exam, rexam$uni, describe)
by(rexam$exam, rexam$uni, stat.desc)

# add further informations eg. not to use basic statistics but have normality statistics. 
by(rexam$exam, rexam$uni, stat.desc, basic = FALSE, norm = TRUE)

#if we want the same for multiple variables we can use cbind() e.g.
by(cbind(data=rexam$exam, data=rexam$numeracy), rexam$uni, stat.desc, basic = FALSE, norm = TRUE)
by(cbind(data=rexam$exam,data=rexam$numeracy), rexam$uni, describe)
by(rexam[,c("exam", "numeracy")], rexam$uni, stat.desc, basic = FALSE, norm = TRUE)

#to plot histogramms and so on it is easyer to make subsets that only contain the data you want to have: 
dunce <- subset(rexam, rexam$uni == "Duncetown University")
sussex <- subset(rexam, rexam$uni == "Sussex University")

hist.dunce <- ggplot(dunce, aes(numeracy)) + theme(legend.position = "none") + geom_histogram(aes(y= ..density..), fill = "white", colour = "black", binwidth = 1) + labs(x = "Numeracy Score", y = "Density") + stat_function(fun=dnorm, args=list(mean = mean(dunce$numeracy, na.rm = TRUE), sd = sd(dunce$numeracy, na.rm = TRUE)), colour = "blue", size=1)
hist.dunce
hist.sussex <- ggplot(sussex, aes(numeracy)) + theme(legend.position = "none") + geom_histogram(aes(y= ..density..), fill = "white", colour = "black", binwidth = 1) + labs(x = "Numeracy Score", y = "Density") + stat_function(fun=dnorm, args=list(mean = mean(sussex$numeracy, na.rm = TRUE), sd = sd(sussex$numeracy, na.rm = TRUE)), colour = "blue", size=1)
hist.sussex
  
hist.dunce <- ggplot(dunce, aes(exam)) + theme(legend.position = "none") + geom_histogram(aes(y= ..density..), fill = "white", colour = "black", binwidth = 1) + labs(x = "Numeracy Score", y = "Density") + stat_function(fun=dnorm, args=list(mean = mean(dunce$exam, na.rm = TRUE), sd = sd(dunce$exam, na.rm = TRUE)), colour = "blue", size=1)
hist.dunce
hist.sussex <- ggplot(sussex, aes(exam)) + theme(legend.position = "none") + geom_histogram(aes(y= ..density..), fill = "white", colour = "black", binwidth = 1) + labs(x = "Numeracy Score", y = "Density") + stat_function(fun=dnorm, args=list(mean = mean(sussex$exam, na.rm = TRUE), sd = sd(sussex$exam, na.rm = TRUE)), colour = "blue", size=1)
hist.sussex  

#### 5.6 Testing wheter a distribution is normal ####


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### Chapter 14 Mixed models ####

install.packages("ez")
install.packages("ggplot2")
install.packages("nlme")
install.packages("pastecs")
install.packages("reshape")
install.packages ("WRS")

library(ez)
library(ggplot2)
library(nlme)
library(pastecs)
library(reshape)
library(WRS2)

dateData<-read.delim("LooksOrPersonality.dat", header = TRUE)
speedData<-melt(dateData, id = c("participant","gender"), 
                measured = c("att_high", "av_high", "ug_high", "att_some", "av_some", "ug_some", "att_none", "av_none", "ug_none"))
names(speedData)<-c("participant", "gender", "groups", "dateRating")

speedData$personality<-gl(3, 60, labels = c("Charismatic", "Average", "Dullard"))
speedData$looks<-gl(3,20, 180, labels = c("Attractive", "Average", "Ugly"))

#boxplots to have a first look at the data
dateBoxplot <- ggplot(speedData, aes(looks, dateRating, colour = personality))
dateBoxplot + geom_boxplot() + labs(x = "Attractiveness", y = "Mean Rating of Date", colour = "Charisma") + facet_wrap(~gender)
imageFile <- paste(imageDirectory,"14 Speed Date Boxplot.png",sep="/")
ggsave(file = imageFile)

by(speedData$dateRating, list(speedData$looks, speedData$personality, speedData$gender), stat.desc, basic = FALSE)

#setting contrasts
AttractivevsAv<-c(1, 0, 0)
UglyvsAv<-c(0, 0, 1)
contrasts(speedData$looks)<-cbind(AttractivevsAv, UglyvsAv)
HighvsAv<-c(1, 0, 0)
DullvsAv<-c(0, 0, 1)
contrasts(speedData$personality)<-cbind(HighvsAv, DullvsAv)

speedData$looks
speedData$personality

baseline <- lme(dateRating ~ 1, random = ~1|participant/looks/personality, data = speedData, method = "ML")
looksM<-lme(dateRating ~ looks, random = ~1|participant/looks/personality, data = speedData, method = "ML")

#or the shorter verision with the funktion "update()
looksM <- update(baseline, .~. +looks)
personalityM<-update(looksM, .~. + personality)
genderM<-update(personalityM, .~. + gender)

looks_gender<-update(genderM, .~. + looks:gender) 
personality_gender<-update(looks_gender, .~. + personality:gender) 
looks_personality<-update(personality_gender, .~. + looks:personality)
speedDateModel<-update(looks_personality, .~. + looks:personality:gender)

#compare the different models that were created by adding one factor or interaction at a time (always keeping the previus model)
anova(baseline, looksM, personalityM, genderM, looks_gender, personality_gender, looks_personality, speedDateModel)
summary(speedDateModel)

# use ggplot to draw a graph for means of the main effect gender
?
  





