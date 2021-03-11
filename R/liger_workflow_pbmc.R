library(liger)
library(cowplot)

pbmc.10x <- read.table('~/Workspace/code/OT-integration/data/pbmc/pbmc_10x.txt',
                       sep="\t",stringsAsFactors=F,header=T,row.names = 1)
pbmc.seqwell <- read.table('~/Workspace/code/OT-integration/data/pbmc/pbmc_seqwell.txt',
                           sep="\t",stringsAsFactors=F,header=T,row.names = 1)
pbmc.data = list(tenx=pbmc.10x, seqwell=pbmc.seqwell)

# Create liger object
a.pbmc <- createLiger(pbmc.data)

a.pbmc <- normalize(a.pbmc)
# Can pass different var.thresh values to each dataset if one seems to be contributing significantly
# more genes than the other
a.pbmc <- selectGenes(a.pbmc, var.thresh = c(0.3, 0.875), do.plot = F)

# In this case, we have a precomputed set of variable genes
s.var.genes <- readRDS('~/Downloads/data/pbmc_alignment/var_genes.RDS')
a.pbmc@var.genes <- s.var.genes
a.pbmc <- scaleNotCenter(a.pbmc)

# running suggestK on multiple cores can greatly decrease the runtime
# k.suggest <- suggestK(a.pbmc, num.cores = 5, gen.new = T, return.results = T, plot.log2 = F,
#                       nrep = 5)

# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
a.pbmc <- optimizeALS(a.pbmc, k=22, thresh = 5e-5, nrep = 3)

# a.pbmc <- quantileAlignSNF(a.pbmc, resolution = 0.4, small.clust.thresh = 20)
# Could also forgo small cluster extraction if desired
# a.pbmc <- quantileAlignSNF(a.pbmc, resolution = 0.4)

a.pbmc <- quantile_norm(a.pbmc)
a.pbmc <- louvainCluster(a.pbmc, resolution = 0.16)

# save clustering result
clusters <- a.pbmc@clusters
clusters <- as.numeric(clusters)
write.table(clusters, file = strcat(c('~/Workspace/Git/OT-integration/data/techs/pbmc/clusters/', 'liger_clusters.txt')), sep = '\t', quote = FALSE)

