

#### Create Barplot for Biology talk ####

getwd()
setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/ANV")
dat <- read.table("biology18_bartable.txt", header = TRUE)
dat <- subset(dat, sample != "bee_pupae") #remove positive controls
dat <- subset(dat, sample != "single_ants")

dat
head(dat)

DWV <- subset(dat, virus=="DWV")
ABPV <- subset(dat, virus=="ABPV")

head(ABPV)
#colors
#control: darkolivegreen3
#DWV: greenblue - deepskyblue3
#ABPVred - 

col1 <- c("deepskyblue3","darkolivegreen3","deepskyblue3","darkolivegreen3")
col2 <- c("red","darkolivegreen3","red","darkolivegreen3")
black <- c("black")
names <- c("colonies treatment (N=18)", "colonies control (N=2)", "queens treatment (N=16)","queens control (N=1)")
names.neg <- c("colonies treatment (N=18)", "colonies control (N=2)", "queens treatment (N=16)","queens control (N=1)")

#DWV

par(las=1)
par(mar=c(5,18,4,2)) # adjust as needed
barplot(rev(DWV$pos), main = "DWV analyses", col = rev(col1), xlab= "% of samples tested positive for DWV", ylab= "", 
        names.arg = rev(names), xlim = c(0,100),
        cex.main=2, cex.lab=1.5, horiz=TRUE, cex.names = 1.5, cex.axis = 1.5)
mtext(side = 2, text = "tested samples", cex=1.5, at = 5.5, line = 1.5, font = 2)
box()
grid(nx=NULL, ny=NA)

barplot(rev(DWV$neg_strand), main = "DWV negative strand analyses", col = rev(col1), xlab= "% of samples tested positive", ylab= "", 
        names.arg = rev(names.neg), xlim = c(0,100),
        cex.main=2, cex.lab=1.5, horiz=TRUE, cex.names = 1.5, cex.axis = 1.5)
mtext(side = 2, text = "tested samples", cex=1.5, at = 5.5, line = 1.5, font = 2)
box()
grid(nx=NULL, ny=NA)

#ABPV
barplot(rev(ABPV$pos), main = "ABPV analyses", col = rev(col2), xlab= "% of samples tested positive for ABPV", ylab= "", 
        names.arg = rev(names), xlim = c(0,100),
        cex.main=2, cex.lab=1.5, horiz=TRUE, cex.names = 1.5, cex.axis = 1.5)
mtext(side = 2, text = "tested samples", cex=1.5, at = 5.5, line = 1.5, font = 2)
box()
grid(nx=NULL, ny=NA)

barplot(rev(ABPV$neg_strand), main = "ABPV negative strand analyses", col = rev(col2), xlab= "% of samples tested positive", ylab= "", 
        names.arg = rev(names.neg), xlim = c(0,100),
        cex.main=2, cex.lab=1.5, horiz=TRUE, cex.names = 1.5, cex.axis = 1.5)
mtext(side = 2, text = "tested samples", cex=1.5, at = 5.5, line = 1.5, font = 2)
grid(nx=NULL, ny=NA)
box()


# graph eurobee
par(mfrow = c(2,2)) # wieviele graphen aufs mal geplottet 
par(cex=1.7)
# par(oma = c(0.1, 2, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
# par(mar = c(0, 0, 0, 0)) # make the plots be closer together

par(las=1)
par(mar=c(3,12,4.5,2)) # adjust as needed
par(xpd=FALSE)
barplot(rev(DWV$pos), main = "DWV", col = rev(black), xlab= "", ylab= "", 
        names.arg = rev(names), xlim = c(0,100), horiz=TRUE) #cex.lab=1.8, cex.names = 1.8, cex.axis = 1.8, cex.main = XX
mtext(side = 2, text = "tested samples", cex= 1.7, at = 6, line = 2.5, font = 2)
mtext(side = 1, text = "% of samples tested positive", cex = 1.7, line = 2)
grid(nx=NULL, ny=NA)
box()

barplot(rev(ABPV$pos), main = "ABPV", col = rev(black), xlab= "% of samples tested positive", ylab= "", 
        names.arg = rev(names), xlim = c(0,100), horiz=TRUE)
mtext(side = 2, text = "tested samples", cex= 1.7, at = 6, line = 2.5, font = 2)
mtext(side = 1, text = "% of samples tested positive", cex = 1.7, line = 2)
grid(nx=NULL, ny=NA)
box()

barplot(rev(DWV$neg_strand), main = "", col = rev(black), xlab= "% of samples tested positive", ylab= "", 
        names.arg = rev(names.neg), xlim = c(0,100),
        cex.main=2.5, horiz=TRUE)
mtext(side = 1, text = "% of samples tested positive", cex = 1.7, line = 2)
grid(nx=NULL, ny=NA)
box()

barplot(rev(ABPV$neg_strand), main = "", col = rev(black), xlab= "% of samples tested positive", ylab= "", 
        names.arg = rev(names.neg), xlim = c(0,100),
        cex.main=2.5, horiz=TRUE)
mtext(side = 1, text = "% of samples tested positive", cex = 1.7, line = 2)
grid(nx=NULL, ny=NA)
box()
par(xpd=TRUE)




par(op) # At end of plotting, reset to previous settings










#### Mortality in the 11 week feeding regime ####

dat <- read.table("ANVlong_full.txt", header = TRUE)
date1 <- subset(dat, date == "16.07.2017")
date2 <- subset(dat, date == "01.10.2017")
date1 <- subset(date1, treatment_1=="control") #(without neonic treatments)
date2 <- subset(date2, treatment_1=="control")

date1$diff_adults <- date2$adults-date1$adults
date1$diff_pupae <- date2$pupae-date1$pupae

date1$diff_adults

boxplot(date1$diff_adults ~ date1$treatment_2, xlab= "treatments", ylab = "nr. of adult ants", main = "Increase of colonysize over 11 weeks")
boxplot(date1$diff_pupae ~ date1$treatment_2)
boxplot(date1$pupae ~ date1$treatment_2)
boxplot(date1$adults ~ date1$treatment_2)
boxplot(date2$adults ~ date2$treatment_2)
boxplot(date2$pupae ~ date2$treatment_2)



# t test for adults
shapiro.test(date1$diff_adults) #non significant --> Ho: Normaldistributed data can be assumed
hist(date1$diff_adults)
t.test(date1$diff_adults ~ date1$treatment_2)

# t test for pupae
shapiro.test(date1$diff_pupae) #non significant --> Ho: Normaldistributed data can be assumed
hist(date1$diff_pupae)
t.test(date1$diff_pupae ~ date1$treatment_2)

# significant difference for the difference in numbers of adults in the two groups 





        