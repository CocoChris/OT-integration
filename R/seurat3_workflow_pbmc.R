library(Seurat)
library(Seurat2)
library(ggplot2)
library(cowplot)
library(pracma)

prefix <- '~/Workspace/code/OT-integration/data/pbmc/'
data.seqwell <- read.delim(strcat(c(prefix, 'pbmc_seqwell.txt')), header = TRUE, row.names = 1, check.names = FALSE)
data.seqwell <- as(as.matrix(data.seqwell), 'dgTMatrix')
pbmc.seqwell <- CreateSeuratObject(data.seqwell)
pbmc.seqwell@meta.data$tech <- 'seqwell'
data.10x <- read.delim(strcat(c(prefix, 'pbmc_10x.txt')), header = TRUE, row.names = 1, check.names = FALSE)
data.10x <- as(as.matrix(data.10x), 'dgTMatrix')
pbmc.10x <- CreateSeuratObject(data.10x)
pbmc.10x@meta.data$tech <- '10x'

pbmc.list <- list()
pbmc.list[[1]] <- pbmc.seqwell
pbmc.list[[2]] <- pbmc.10x

for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] <- Seurat::NormalizeData(pbmc.list[[i]], verbose = FALSE)
  pbmc.list[[i]] <- Seurat::FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", 
                                             nfeatures = 3000, verbose = FALSE)
}

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:30)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)

DefaultAssay(pbmc.integrated) <- "integrated"

pbmc.integrated <- Seurat::ScaleData(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- Seurat::RunPCA(pbmc.integrated, npcs = 30, verbose = FALSE)
pbmc.integrated <- Seurat::RunTSNE(pbmc.integrated, reduction = "pca", dims = 1:30)

pbmc.integrated <- Seurat::FindNeighbors(pbmc.integrated)
pbmc.integrated <- Seurat::FindClusters(pbmc.integrated, resolution = 0.4)

# save clustering result
clusters <- pbmc.integrated@meta.data$seurat_clusters
clusters <- as.numeric(clusters)
write.table(clusters, file = strcat(c(prefix, 'clusters/seurat3_clusters.txt')), sep = '\t', quote = FALSE)

# save tsne result
tsne_coords <- pbmc.integrated@reductions$tsne@cell.embeddings
write.table(tsne_coords, file = strcat(c(prefix, 'tsne/', 'pbmc_seurat3.txt')), sep = '\t', quote = FALSE)