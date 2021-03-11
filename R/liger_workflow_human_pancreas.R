library(liger)
library(cowplot)
library(pracma)

data.human1 <- read.table('~/Workspace/code/OT-integration/data/human_pancreas/human1.txt',
                          sep="\t",stringsAsFactors=F,header=T,row.names = 1)
data.human2 <- read.table('~/Workspace/code/OT-integration/data/human_pancreas/human2.txt',
                          sep="\t",stringsAsFactors=F,header=T,row.names = 1)
data.human3 <- read.table('~/Workspace/code/OT-integration/data/human_pancreas/human3.txt',
                          sep="\t",stringsAsFactors=F,header=T,row.names = 1)
data.human4 <- read.table('~/Workspace/code/OT-integration/data/human_pancreas/human4.txt',
                          sep="\t",stringsAsFactors=F,header=T,row.names = 1)

data.human1 <- as(as.matrix(data.human1), "dgCMatrix")
data.human2 <- as(as.matrix(data.human2), "dgCMatrix")
data.human3 <- as(as.matrix(data.human3), "dgCMatrix")
data.human4 <- as(as.matrix(data.human4), "dgCMatrix")

data.human1@Dimnames[[2]] <- paste('human1', data.human1@Dimnames[[2]], sep = '_')
data.human2@Dimnames[[2]] <- paste('human2', data.human2@Dimnames[[2]], sep = '_')
data.human3@Dimnames[[2]] <- paste('human3', data.human3@Dimnames[[2]], sep = '_')
data.human4@Dimnames[[2]] <- paste('human4', data.human4@Dimnames[[2]], sep = '_')

data = list(human1 = data.human1, human2 = data.human2, human3 = data.human3, human4 = data.human4)

# Create liger object
human.liger <- createLiger(data)

human.liger <- normalize(human.liger)
# Can pass different var.thresh values to each dataset if one seems to be contributing significantly
# more genes than the other
human.liger <- selectGenes(human.liger, var.thresh = c(0.3, 0.9), do.plot = F)

human.liger <- scaleNotCenter(human.liger)

# running suggestK on multiple cores can greatly decrease the runtime
# k.suggest <- suggestK(a.pbmc, num.cores = 5, gen.new = T, return.results = T, plot.log2 = F,
#                       nrep = 5)

# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
human.liger <- optimizeALS(human.liger, k=20, thresh = 5e-5, nrep = 3)

human.liger <- quantile_norm(human.liger)
human.liger <- louvainCluster(human.liger, resolution = 0.5)

# save clustering result
clusters <- human.liger@clusters
clusters <- as.numeric(clusters)
write.table(clusters, file = strcat(c('~/Workspace/code/OT-integration/data/human_pancreas/clusters/', 'liger_human_clusters.txt')), sep = '\t', quote = FALSE)

# save tsne result
human.liger <- runTSNE(human.liger)
tsne.coords <- human.liger@tsne.coords
write.table(tsne.coords, file = strcat(c('~/Workspace/code/OT-integration/data/human_pancreas/tsne/', 'liger_human.txt')), sep = '\t', quote = FALSE)

