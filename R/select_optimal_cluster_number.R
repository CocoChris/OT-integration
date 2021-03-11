library(NbClust)
library(ggplot2)
library(RColorBrewer)

data.seqwell <- read.table("~/Workspace/code/OT-integration/data/pbmc/cca/seqwell.txt")
data.seqwell <- as.matrix(data.seqwell)
data.10x <- read.table("~/Workspace/code/OT-integration/data/pbmc/cca/10x.txt")
data.10x <- as.matrix(data.10x)

data <- rbind(data.seqwell, data.10x)

res <- NbClust(data.10x[,1:20], diss=NULL, distance = "euclidean", min.nc=2, max.nc=10, 
               method = "kmeans", index = "all") 
res$All.index
res$Best.nc
res$All.CriticalValues
res$Best.partition
  


