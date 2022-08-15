############################
### Kaspar Diplomarbeit  ###
############################

getwd()
setwd("H:/DS-Bees/R/kaspar")

crit <- read.table("critidia.txt", header = TRUE)
nosema <- read.table("nosema.txt", header = TRUE)
nosemaold <- read.table("nosema_old.txt", header = TRUE)

crit$logcps <- log10(crit$cps)
nosema$logcps <- log10(nosema$cps)
nosemaold$logcps <- log10(nosemaold$cps)



crit$treat <- as.factor(crit$treat) 
crit$pathogen <- as.factor(crit$pathogen)

nosema$treat <- as.factor(nosema$treat) 
nosema$pathogen <- as.factor(nosema$pathogen)

nosemaold$treat <- as.factor(nosemaold$treat) 
nosemaold$pathogen <- as.factor(nosemaold$pathogen)

####
#deskriptive statistik 
names(nosema)
controlnosema <- subset(nosema, treat=="control")
mean(controlnosema$logcps)
sd(controlnosema$logcps)

nosematreat <- subset(nosema, treat!="control")
mean(nosematreat$logcps)
sd(nosematreat$logcps)

#nosema alt 
cnosalt <- subset(nosemaold, treat=="control")
mean(cnosalt$logcps)
sd(cnosalt$logcps)

wnosalt <- subset(nosemaold, treat=="high" & product=="wax")
mean(wnosalt$logcps)
sd(wnosalt$logcps)

pnosalt <- subset(nosemaold, treat!="control" & product=="pollen")
mean(pnosalt$logcps)
sd(pnosalt$logcps)

#crithidia
critc <- subset(crit, treat=="control")
mean(crit$logcps)
sd(crit$logcps)

critx <- subset(crit, treat!="control" & product=="pollen" | treat!="control" & product=="wax")
critx
mean(critx$logcps)
sd(critx$logcps)

critt <-  subset(crit, treat!="control" & product=="honey")
mean(critt$logcps)
sd(critt$logcps)

#Bienenprodukte
nfresht <- subset(nosema, treat!="control")
mean(nfresht$logcps)
sd(nfresht$logcps)

noldh <- subset(nosemaold, treat!="control" & product == "honey")
mean(noldh$logcps)
sd(noldh$logcps)

noldp <- subset(nosemaold, treat!="control" & product == "pollen")
mean(noldp$logcps)
sd(noldp$logcps)

crith <- subset(crit, treat!="control" & product == "honey")
mean(crith$logcps)
sd(crith$logcps)

critr <- subset(crit, treat!="control" & product == "pollen" | treat!="control" & product=="wax")
mean(critr$logcps)
sd(critr$logcps)




####
# Graphen mit deutscher beschriftung

### Critidia

#create subsets with products
crit$treat <- factor(crit$treat, levels = c("control", "low", "medium", "high"))
honeyc <- subset(crit, product=="honey")
pollenc <- subset(crit, product=="pollen")
waxc <- subset(crit, product=="wax")

par(mfrow=c(2,2))

boxplot(honeyc$logcps~honeyc$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Honig", ylim=c(3.5,9.5))
text(1,9, c("a"), cex=1)
text(2,9, c("b"), cex=1)
text(3,9, c("b"), cex=1)
text(4,9, c("b"), cex=1)

boxplot(pollenc$logcps~pollenc$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Pollen", ylim=c(3.5,9.5))
text(1,9, c("a"), cex=1)
text(2,9, c("a"), cex=1)
text(3,9, c("a"), cex=1)
text(4,9, c("a"), cex=1)

boxplot(waxc$logcps~waxc$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Wachs", ylim=c(3.5,9.5))
text(1,9, c("a"), cex=1)
text(2,9, c("a"), cex=1)
text(3,9, c("a"), cex=1)
text(4,9, c("a"), cex=1)

### Nosema

par(mfrow=c(2,2))

nosema$treat <- factor(nosema$treat, levels = c("control", "low", "medium", "high"))
honeyn <- subset(nosema, product=="honey")
pollenn <- subset(nosema, product=="pollen")
waxn <- subset(nosema, product=="wax")

boxplot(honeyn$logcps~honeyn$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Honig", ylim=c(4.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("b"), cex=1)
text(3,9.8, c("c"), cex=1)
text(4,9.8, c("d"), cex=1)

boxplot(pollenn$logcps~pollenn$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Pollen", ylim=c(4.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("b"), cex=1)
text(3,9.8, c("b"), cex=1)
text(4,9.8, c("b"), cex=1)

boxplot(waxn$logcps~waxn$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Wachs", ylim=c(4.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("b"), cex=1)
text(3,9.8, c("c"), cex=1)
text(4,9.8, c("c"), cex=1)


### Nosema old

par(mfrow=c(2,2))

nosemaold$treat <- factor(nosemaold$treat, levels = c("control", "low", "medium", "high"))
honeyno <- subset(nosemaold, product=="honey")
pollenno <- subset(nosemaold, product=="pollen")
waxno <- subset(nosemaold, product=="wax")

boxplot(honeyno$logcps~honeyno$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Honig", ylim=c(4.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("a"), cex=1)
text(3,9.8, c("a"), cex=1)
text(4,9.8, c("b"), cex=1)

boxplot(pollenno$logcps~pollenno$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Pollen", ylim=c(4.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("b"), cex=1)
text(3,9.8, c("b"), cex=1)
text(4,9.8, c("b"), cex=1)

boxplot(waxno$logcps~waxno$treat, col = "brown1", outline=TRUE, range=7, xlab="Behandlungsgruppe", ylab="Anzahl Parasiten pro Biene [log]", names=c("Kontrolle", "Tief", "Mittel", "Hoch"), main="Wachs", ylim=c(4.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("a"), cex=1)
text(3,9.8, c("a"), cex=1)
text(4,9.8, c("b"), cex=1)




############################
### Infectivity of Bee products
############################

# Nosema fresh
#subset ohne Kontrolle
nosematreat <- subset(nosema, treat!="control")

par(mfrow=c(1,3))
par(cex=1)

boxplot(nosematreat$logcps~nosematreat$product, col = "brown1", outline=TRUE, range=1.5, xlab="Produkt", ylab="Anzahl Parasiten pro Biene [log]", names=c("Honig", "Pollen", "Wachs"), main="Nosema frisch", ylim=c(3.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("a"), cex=1)
text(3,9.8, c("a"), cex=1)


# Nosema old
nosemaoldtreat <- subset(nosemaold, treat!="control")
boxplot(nosemaoldtreat$logcps~nosemaoldtreat$product, col = "brown1", outline=TRUE, range=1.5, xlab="Produkt", ylab="Anzahl Parasiten pro Biene [log]", names=c("Honig", "Pollen", "Wachs"), main="Nosema alt", ylim=c(3.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("b"), cex=1)
text(3,9.8, c("c"), cex=1)


# Critidia
crittreat <- subset(crit, treat!="control")
boxplot(crittreat$logcps~crittreat$product, col = "brown1", outline=TRUE, range=1.5, xlab="Produkt", ylab="Anzahl Parasiten pro Biene [log]", names=c("Honig", "Pollen", "Wachs"), main="Crithidia", ylim=c(3.5, 10))
text(1,9.8, c("a"), cex=1)
text(2,9.8, c("b"), cex=1)
text(3,9.8, c("b"), cex=1)






citation()
