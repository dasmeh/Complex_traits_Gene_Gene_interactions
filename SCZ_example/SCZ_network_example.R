# Load libraries
library("readr")
library("STRINGdb")
library("igraph")
library("rbioapi")
library("fitdistrplus")


# We will store the AIC values of three different degree distributions in this vector
fit_aic=rep(NA, 3)


## Getting MAGMA p-values 
SCZ=read.table("magma.genes_SCZ.out", header = TRUE)


## Selecting the genes whose p-values are below a threshold (~1e-6)
proteins=SCZ$SYMBOL[which(SCZ$P < 1e-6)]

## Now we map our protein IDs
proteins_mapped <- rba_string_map_ids(ids = proteins,
                                      species = 9606)

## What we need and will use for the rest of this vignette is the `stringId` column
proteins_ids=proteins_mapped$stringId

int_net <- rba_string_interactions_network(ids = proteins_ids,
                                           species = 9606,
                                           required_score = 900)

## Extract the names of the first and second genes, along with their scores    
matrix=int_net[,c(3,4,6)]


# Create igraph graph object
graph <- graph_from_data_frame(matrix, directed=FALSE)

# Set edge weights to connectivity values
E(graph)$weight <- matrix$score


# Remove duplicates    
graph <- simplify(graph)

# Plot the graph

# Set node and edge colors
node_color <- "#FFCCCC"  # Light red
edge_color <- "#CCCCCC"  # Gray

plot(graph, 
     vertex.size = 8, 
     vertex.label.cex = 0.5 , main="SCZ", vertex.label.color="black", 
     vertex.color = node_color, edge.color = edge_color, edge.width=5)

# Calculate the degree distribution
degree_distribution_SCZ <- degree(graph)


# Plot the degree distribution
hist(degree_distribution_SCZ)



# Fit three types of distributions to the degree distribution of the SCZ network and determine the
# distribution with the lowest value of Akaike Information Criteria (AIC)
fit_normal <- fitdist(degree_distribution_SCZ,  "norm")
fit_lognormal <- fitdist(degree_distribution_SCZ,  "lnorm")
fit_exp <- fitdist(degree_distribution_SCZ,  "exp")


fit_aic[1]= fit_normal$aic
fit_aic[2]= fit_lognormal$aic
fit_aic[3]= fit_exp$aic

names(fit_aic)=c("Normal Dist", "LogNorm Dist", "Exponential")


