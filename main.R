library(tidyverse)
library(readxl)
library(MASS)
library(ape)
library(CCA)
library(CCP)

inequality <- read_excel("data/main.xlsx")

ineq.matrix <- as.matrix(inequality[2:7])
row.names(ineq.matrix) <- inequality$country

# just for comparison we'll do PCA on the dataset
pca <- prcomp(ineq.matrix, scale=T, center=T)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)

set.seed(123)
pca.km <- kmeans(pca$x[,1:2], centers = 4)

plot(pca$x[,1], pca$x[,2],
  xlab = paste0("PC1 - ", pca.var.per[1], "%"),
  ylab = paste0("PC2 - ", pca.var.per[2], "%"),
  xlim = range(pca$x[,1])*1.2,
  type = "n",
  main = "PCA Graph"
)
text(pca$x[,1], pca$x[,2], labels = row.names(pca$x), cex=0.6, col=pca.km$cluster)

# now it is time for MDS
distance.matrix <- dist(scale(ineq.matrix, center=T, scale=T), method = "euclidean")

ineq.mds <- cmdscale(distance.matrix, eig=T)
mds.var.per <- round(ineq.mds$eig/sum(ineq.mds$eig)*100,1)

set.seed(123)
mds.km <- kmeans(ineq.mds$points, centers = 4)

plot(ineq.mds$points[,1],ineq.mds$points[,2],
     xlab = paste0("Coordinate 1 - ", mds.var.per[1], "%"),
     ylab = paste0("Coordinate 2 - ", mds.var.per[2], "%"),
     xlim = range(ineq.mds$points[,1])*1.2,
     type = "n",
     main = "MDS plot using Euclidean distance"
)
text(ineq.mds$points[,1], ineq.mds$points[,2], labels = row.names(ineq.mds$points), cex=0.6,
     col = mds.km$cluster)

# Minimum spanning tree
st <- mst(distance.matrix)
plot(ineq.mds$points[,1], ineq.mds$points[,2],
     xlab = "Coordinate 1", ylab = "Coordinate 2",
     xlim = range(ineq.mds$points[,1])*1.2, type = "n",
     main = "Minimum spanning tree for the inequality dataset")
for (i in 1:nrow(inequality)) {
  w1 <- which(st[i, ] == 1)
  segments(ineq.mds$points[,1][i], ineq.mds$points[,2][i], ineq.mds$points[,1][w1], ineq.mds$points[,2][w1])
  }
text(ineq.mds$points[,1], ineq.mds$points[,2], labels = inequality$country, cex = 0.6)

# Non-metric MDS
ineq.nmds <- isoMDS(d=distance.matrix, k=2)

set.seed(123)
nmds.km <- kmeans(ineq.nmds$points, centers=4)

plot(ineq.nmds$points[,1], ineq.nmds$points[,2],
     xlab = "Coordinate 1",
     ylab = "Coordinate 2",
     xlim = range(ineq.nmds$points[,1])*1.2,
     type = "n",
     main = "Non-metric MDS"
)
text(ineq.nmds$points[,1], ineq.mds$points[,2], labels = row.names(ineq.nmds$points), 
     cex=0.6, col=nmds.km$cluster)

ineq.sh <- Shepard(distance.matrix, ineq.nmds$points)
plot(ineq.sh, pch = ".", xlab = "Dissimilarity",  ylab = "Distance",
     main = "Shepard diagram")
lines(ineq.sh$x, ineq.sh$yf, type = "S")

# Canonical correlation analysis
# load first dataset
dietary <- read.csv2("data/dietary.csv")
names(dietary)[1] <- "country"
for(i in c(2:ncol(dietary))){
  dietary[,i] <- as.numeric(dietary[,i])
}

X <- as.matrix(dietary[,2:10])
rownames(X) <- dietary$country

# load second dataset
sectors <- read.csv2("data/sectors.csv")
names(sectors)[1] <- "country"
for(i in c(2:ncol(sectors))){
  sectors[,i] <- as.numeric(sectors[,i])
}

Y <- as.matrix(sectors[,2:10])
rownames(Y) <- sectors$country

# canonical correlation
res.cc <- cc(X,Y)
barplot(res.cc$cor, xlab = "Component", ylab = "Canonical correlation values", names.arg = 1:9)
plt.cc(res.cc, var.label = TRUE)

# canonical correlations
res.cc$cor

# raw canonical coefficients (interpreted in a way that is analogous to interpreting regression coefficients)
res.cc[3:4]

# compute canonical loadings
cc2 <- comput(X, Y, res.cc)
# display canonical loadings
cc2[3:6]

# tests of canonical dimensions
rho <- res.cc$cor
## Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(X)[1]
p <- ncol(X)
q <- ncol(Y)

## Calculate p-values using the F-approximations of different test statistics:
# the first test of the canonical dimensions tests whether all 9 dimensions are significant (and they are significant)
# the next test tests whether dimensions 2 to 9 combined are significant (they are not)
# and so on until the 9th dimension (all not significant)
p.asym(rho, n, p, q, tstat = "Wilks")

