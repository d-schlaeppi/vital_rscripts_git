# Linear Models and linear mixed models in R - a basic tutorial #
# http://www.bodowinter.com/tutorial/bw_LME_tutorial1.pdf

#question: is voice pitch determined/predicted by sex (male/female)
# pitch = dependent variable, sex = fixed effect (explanatory variable), other not controllable factors will be random factors

pitch = c(233,204,242,130,112,142)
sex = c(rep("female",3),rep("male",3))

my.df = data.frame(sex, pitch)

xmdl = lm(pitch ~ sex, my.df)
summary(xmdl)

#Model output

# Multiple R-squared = the	 statistic	 R2 which	 is	 a	 measure	 of "variance	 explained"	(ranging from 0 to 1)
# showing	that	92.1%	of	the	stuff	that's	happening	in	our	dataset	is	"explained" by	our	model.
# in this case the gixed effec (only one) explains 92% of the variance in the dataset

#F-statistic: 46.61 on 1 and 4 DF,  p-value: 0.002407
#p value gives probability that the obtained data is collected if the null hypothesis (sex has no effect on pitch) is true... 
# in this case it is very unlikely (low p value). If it is lower than a treshhold (most often 0.05 --> 95 % probability) 
# thus we assume that the alternative hypothesis "sex affects pitch is more likely und thus this result is statistically significant
# to report: "We	 constructed	a	linear	model	of	 pitch	as	a	 function	 of	 sex.	This	model	was	significant	(F(1,4)=46.61,	p<0.01).	(.)"

#Coefficient table
# as we only had one fixe effect p value for sexmale is the same as the overall effect
# Estimate intercept: = Mean of the variable of the reference group
# Estimate sexmale: in this case negative showing that the mean for men lies -98.33 below the mean of the reference group = estimate of the difference between males and females
# because its a linear model the model looks at the difference as a slope (line form famles to males) und thus needs to go down -98.33
# Females are on x axis = zero and males are on x coordinate = 1
# the p-values to the right correspond to tests wheter the coefficents under estimate are "non-zero"
# Obviously,	226.33	Hz	is	different	from	zero,	so	the	intercept	is	"significant"	
# with	a	very	low	p-value.	The	slope	-98.33	is	also	different	 from	zero (but	in	 the negative	direction),	and	so	this	is	significant	as	well.

# the thinking with slopes is practical as it works also if the expolaining variable is not categorical but continuous like zb age, that could also explain the pitch of a voice

age = c(14,23,35,48,52,67)
pitch = c(252,244,240,233,212,204)
my.df = data.frame(age,pitch)
xmdl = lm(pitch ~ age, my.df)
summary(xmdl)

# Look at coefficient again: 
# The intercept now is 267.08 which is the predicted pitch value for people with age 0
# age: the estimate readS: for every increase of age by 1 the voice pitch is decreased by 0.9099
plot(my.df$pitch ~ my.df$age)
abline(xmdl)
# the line represents the coefficients form the model: it goes at 267.08 through zero and decreases by 0.9 by each year age increases

# meaning full and meaningless intercepts
# create a new variable that centers the data 
my.df$age.c = my.df$age - mean(my.df$age)
xmdl = lm(pitch ~ age.c, my.df)
summary(xmdl)
mean(my.df$age.c)
# the slopw ist stille the same but the intercept is now the mean voice pitch 

# However sex and age could both affect the voice pitch and maybe even dialect as well
# --> p value at the bottom will be for the overall model, and the ones in the coefficient table for the fixed effects



# Assumptions: 
# A model has some assumptions and thus , there are conditions that need to be satisfied in order for the liner model to be meaningful
# 1) Linearity. 
# to check for that look at the residual plot (residuals =  the amount each samples varries from the expeced value of the model)
plot(fitted(xmdl),residuals(xmdl))
abline(a=0, b=0)
# if residuals deviat in a specific pattern the data violates the assumtption of linearity 

# 2) absence of collinearity
# if two fixed effects are correlated they are said to be collinear
# thus the model becomes unstable --> beforhand think of which explanatory variable is most meaningfull

# 3) Homoskedasticity (or absence of herteroskedasticity)
# 	It	says	that	that the variance should be more or less equal across the range of the predicted values
# residuals need to roughly have a similar amount of deviation from a mean --> a good plot looks like a blob 

plot(rnorm(100),rnorm(100))
# This	creates	two	sets	of	100	normally	distributed	random	numbers	with	a	mean	
# of	 0	 and	 a	 standard	 deviation	 of	 1.	 If	 you	 type	 this	 in	 multiple	 times	 to	 create	
# multiple	plots,	you	can	get	a	feel	of	how	a	"normal"	residual	plot	should	look	like.

# if any of the assumptions are violated consider transforming the data (eg. logtransforming of response or squaring explanatory variable)
# (4)	Normality	of	residuals (Least important) --> qqplot
hist(residuals(xmdl))
qqnorm(residuals(xmdl))
# the histogramm would need to be bellshaped and in the QQ plot the dots should fall on a straight line

#(5)	Absence	of	influential	data	points
# influental datapoints can drasticly change the interpretation of data and thus should be aressed
# use dfbeta()
dfbeta(xmdl)
#for each coeffecient the function returns the dfbeta values (leave one out diagnostics) if a specific datapoint is left out
# The	first	row	means	that	the	coefficient	for	age	(which was	 -0.9099)	 has	 to	 be	 adjusted	 by	 0.06437573 if	 data	 point	 1	 is	
# excluded.	 That	 means	 that	 the	 coefficient	 of	 the	 model	 without	 the	 data	 point
# would	 be	 -0.9742451	 (which	 is	 -0.9099	 minus	 0.0643757)
# any value that changes the slope drasticly (or even the sign of the slop) then consider it to be influetial
# if you have such a datapoint it is not legit to just remove it (unless it is an obvious error), but you can run the analyses with and without the point and report both results

# 6) INDEPENDANCE!!
# in our case the subjects were all different und thus the measured variables were independet from each other

### if data is not independant you will need mixed effect models!!! ### 



#### A very basic tutorial for performing linear mixed effect analyses (TUTORIAL 2) #### 
# http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf

# if variable are not independant (eg. multiple measurements form the same subject) a random effect needs to be added to a model, 
# this will allow it tho assume a different baseline value of the response variable for each subject
# in the case of the pitch study a baselin pitch per subject that might then vary depending on politeness and sex
# eg. compare the mean pitch per subject
# now we no longer generalize the error term across the board with a general error term "??" but add a structur to this error term 
# by adding a random effect for subject that characterizes the idiosymcratic variation due to individual differences
# pitch ~ politeness + sex + (1|subject) + ??
# (1|subject): 1 = interrcept, will get mutliple responses per subject and each subject will have different baseline
# we can add another one (1|item) which will be that different response items each might account for variation of pitch (eg. excusing for coming too late might always resutl in higher pitch)
# thus the same response items might not be independant
# We now "resolved" those non-independencies (ourmodel knows that there are multiple responses per subject and per item), and we accounted for by-subject and by-item variation in overall pitch levels
# pitch ~ politeness + sex + (1|subject) + (1|item) + ??

install.packages("lme4") # contains the function lmer() (mixed model equvalent for lm)
library(lme4)
#read data from server
politeness=
  read.csv("http://www.bodowinter.com/tutorial/politeness_data.csv")
head(politeness)
which(is.na(politeness$frequency))
which(!complete.cases(politeness))

boxplot(frequency ~ attitude*gender,
        col=c("white","lightgray"),politeness)
#start with first model
politeness.model = lmer(frequency ~ attitude +
                          (1|subject) + (1|scenario), data=politeness)

summary(politeness.model)

# what does the output mean: 
# random effects
# standard deviation = how much variability in the dependent is there due to scenario and subject (more for subject) and residual (which stands for the variability thats due to neither of the random factors)
# fixxed effect = estimate -mean form inpolite, and -19.694 to go to the mean of the polite
# that mean is the average of the data between males and females... --> we need to include gender

politeness.model = lmer(frequency ~ attitude +
                          gender + (1|subject) +
                          (1|scenario), data=politeness)

summary(politeness.model)

#random effects
# the variability due to subject dropped considerably because before it was confounded with gender the variance is now more explained by the fixed effects and a bit less by random terms

#fixed effects: 
# Intercept is now for the mean for femals infromal, we need to decrease by 19 to go to the formal state and decrease by 109 to switch to males!

#but are these differences statistically relevant? p values are not straight foreward in mixed models

#one possibility is to focus on the Likelyhood Ratio Test as means to attain p-values
# thus we compare our modle with a model that does not contain the factor we are interested in. 

politeness.null = lmer(frequency ~ gender +
                         (1|subject) + (1|scenario), data=politeness,
                       REML=FALSE)
# REML = FALSE is necessary to do when models are compared using likelyhood ratio tests
politeness.model = lmer(frequency ~ attitude +
                          gender + (1|subject) + (1|scenario),
                        data=politeness, REML=FALSE)

anova(politeness.null, politeness.model)

# this output can then be reported as: politeness affected pitch (??2(1)=11.62, p=0.00065), lowering it by about 19.7 Hz ± 5.6 (standard errors) 

# what happens if you have an interaction?
# this can be compared by comparing
# full model: frequency ~ attitude*gender
# reduced model: frequency ~ attitude + gender
# if the two models differ significantly then you know that attitude ande gender are significantly interdependent


# super crucial: random slopes and random interceps
coef(politeness.model)
# so far the slope is the same for all subjects (rondom intercept model), but that is not always valid  
# each subject could have different slopes --> random slope model 
politeness.model = lmer(frequency ~ attitude +
                          gender + (1+attitude|subject) +
                          (1+attitude|scenario),
                        data=politeness,
                        REML=FALSE)

# The notation "(1+attitude|subject)" means that you tell the model to expect differing baseline-levels of frequency (the intercept,
# represented by 1) as well as differing responses to the main factor in question, which is "attitude" in this case. 
coef(politeness.model)
#the coefiicient for politeness differs now between subject but the variability is similar among participants indicating that the voice goes down wehn speaking politely

#get p-value
politeness.null = lmer(frequency ~ gender +
                         (1+attitude|subject) + (1+attitude|scenario),
                       data=politeness, REML=FALSE)

anova(politeness.null,politeness.model)

# If possible always include random lopes that are justified (for all fixed effects that are importatnt for the overall interpretation) 
# and  keep the model maximal 

#Assumptions: 
# linearity, collinearity, homoscedasticity, normality of residual (all the same as for lm)


