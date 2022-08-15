### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Possible effects of Networkstructure on disease transmission    ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

setwd("G:/BIENEN/Personnel/Daniel_Schlaeppi/DS-Bees_new/R/Nest")

library(dplyr)
library(igraph)

mono <- read.table("monodomous_interactions.txt", header = TRUE)
poly <- read.table("polydomous_interactions.txt", header = TRUE)

# create new variable to determine interaction-duration as weight for the network analyses 
mono$duration <- (mono$Stoptime-mono$Starttime)+0.25
poly$duration <- (poly$Stoptime-poly$Starttime)+0.25
# cut superlong interactions that might arise from ants resting close to each other to a max. of 120 seconds
mono$duration[mono$duration>120] <- 120
poly$duration[poly$duration>120] <- 120
# delete some no longer needed information
mono[3:6] = NULL
poly[3:6] = NULL

# sum up the interaction durations for each pair of ants
mono_weighted <- mono %>% 
  group_by(Ant1, Ant2) %>% 
  summarise(duration = sum(duration))

poly_weighted <-poly %>% 
  group_by(Ant1, Ant2) %>% 
  summarise(duration = sum(duration))

# create weighted graph opjects for the package i-graph
# define duration as the weight?
g1 <- graph_from_data_frame(d = mono_weighted, directed = FALSE) 
g2 <- graph_from_data_frame(d = poly_weighted, directed = FALSE)

#plot the graphs to see how the networks look like 
V(g1)$vertex_degree <- degree(g1)
plot(g1, 
     vertex.label=NA,
     edge.width=log10((E(g1)$duration)+1)*0.1,
     edge.color = 'black',
     vertex.size = log10((V(g1)$vertex_degree)+1),
     layout = layout_with_fr(g1))

V(g2)$vertex_degree <- degree(g2)
plot(g2, 
     vertex.label.cex = 0.8,
     edge.width=log10((E(g2)$duration)+1)*0.1,
     edge.color = 'black',
     vertex.size = log10((V(g2)$vertex_degree)+1),
     layout = layout_with_fr(g1))

### Calculation of network characteristics and derive assumptions on disease transmission ### 
# based on known effects of global network properties on disease transmission 
# (comp. Tab1 in Stroeymeyt et al. (2018) Science, 362(6417), 941-945)
# Network caracteristics that enhance or decrease disease transmission
# ENHANCiNG: Density, Network efficiency, Degree centrality
# DECREASING: Diameter, Modularity, Clustering coefficient

# Density - proportion of realized connections among all possible connections 
ed1 <- edge_density(g1)
ed2 <- edge_density(g2)

# Network efficiency - average connection efficiency of all pairs of nodes (inverse average shortest distance between nodes)
ef1 <- 1/mean_distance(g1)
ef2 <- 1/mean_distance(g2)

# Degree centrality - Average number of edges connecting a node to other network nodes normalized by dividing with the max. possible, i.e. (N-1) to account for network size
dc1 <- mean(centr_degree(g1)$`res`/(gorder(g1)-1))
dc2 <- mean(centr_degree(g2)$`res`/(gorder(g2)-1))

# Diameter - Maximum length of the shortest paths between all pairs of nodes
dia1 <- diameter(g1)
dia2 <- diameter(g2)

# Modularity
ebc1 <- edge.betweenness.community(g1, weights = E(g1)$duration)
sizes(ebc1)
mod1 <- modularity(g1,membership(ebc1))

ebc2 <- edge.betweenness.community(g2, weights = E(g2)$duration)
sizes(ebc2)
mod2 <- modularity(g2,membership(ebc2))

# Transitivity measures the probability that the adjacent vertices of a vertex are connected 
# somethimes refered to as clustering coefficient
# Average tendency of  neighboring nodes of any given node to also be connected
  
cc1 <- mean(transitivity(g1, type = "weighted", weights = E(g1)$duration))
cc2 <- mean(transitivity(g2, type = "weighted", weights = E(g2)$duration))  


# Based on the direction of effects the variables are entered for enhancing effects as positive values and for decreasing effects as negative values in a new vector
# Diameter is 
monodomous <- c(ed1, ef1, dc1, -(dia1), -(mod1), -(cc1))
polydomous <- c(ed2, ef2, dc2, -(dia2), -(mod2), -(cc2))
#remove any variable in the vectors if they are the same, as no effect of nest architecture on transmission dynamicson can be deduced from them
y <- monodomous != polydomous

# compare the remaining variables of the two vectors
# If the output is TRUE, then variable favours disease transmission in the monodomous nest 

output <- monodomous[y] > polydomous[y]
names <- c("Density", "Network efficiency", "Degree centrality", "Modularity", "Clustering coefficient")
z <- cbind(output, names)
z

### ### Conclusion ### ### 
# There might be an effect of nest architecture on disease dynamics. 
# In the network of the polydomous nest most network characteristics appear to disfavor 
# the transmission of diseases compared to the network in the monodomous nest.
# Thus, one might hypothesize that polydomy is a favorable trait when disease transmission is considered.  










