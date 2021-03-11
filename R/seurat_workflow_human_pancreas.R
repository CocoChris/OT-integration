################################################################################
### Alignment workflow for human and mouse pancreatic islet datasets
################################################################################

library(Seurat)
library(Seurat2)
library(Matrix)
library(RColorBrewer)
library(pracma)

# Read in human expression matrices
human1.data <- read.delim("~/Workspace/code/OT-integration/data/human_pancreas/human1.txt", header = TRUE, row.names = 1, check.names = FALSE)
human2.data <- read.delim("~/Workspace/code/OT-integration/data/human_pancreas/human2.txt", header = TRUE, row.names = 1, check.names = FALSE)
human3.data <- read.delim("~/Workspace/code/OT-integration/data/human_pancreas/human3.txt", header = TRUE, row.names = 1, check.names = FALSE)
human4.data <- read.delim("~/Workspace/code/OT-integration/data/human_pancreas/human4.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Convert to sparse matrices for efficiency
human1.data <- as(as.matrix(human1.data), "dgCMatrix")
human2.data <- as(as.matrix(human2.data), "dgCMatrix")
human3.data <- as(as.matrix(human3.data), "dgCMatrix")
human4.data <- as(as.matrix(human4.data), "dgCMatrix")

# Create and setup Seurat objects for each dataset
human1 <- CreateSeurat2Object(raw.data = human1.data)
human1 <- NormalizeData(human1)
human1 <- FindVariableGenes(human1, do.plot = F, display.progress = F)
human1 <- ScaleData(human1)
human1@meta.data$human <- "human1"
colnames(human1@data) <- paste('human1', make.names(colnames(human1@data), unique = T), sep = '_')
colnames(human1@raw.data) <- paste('human1', make.names(colnames(human1@raw.data), unique = T), sep = '_')
colnames(human1@scale.data) <- paste('human1', make.names(colnames(human1@scale.data), unique = T), sep = '_')
rownames(human1@meta.data) <- paste('human1', make.names(rownames(human1@meta.data), unique = T), sep = '_')
human1@cell.names <- paste('human1', make.names(human1@cell.names, unique = T), sep = '_')

human2 <- CreateSeurat2Object(raw.data = human2.data)
human2 <- NormalizeData(human2)
human2 <- FindVariableGenes(human2, do.plot = F, display.progress = F)
human2 <- ScaleData(human2)
human2@meta.data$human <- "human2"
colnames(human2@data) <- paste('human2', make.names(colnames(human2@data), unique = T), sep = '_')
colnames(human2@raw.data) <- paste('human2', make.names(colnames(human2@raw.data), unique = T), sep = '_')
colnames(human2@scale.data) <- paste('human2', make.names(colnames(human2@scale.data), unique = T), sep = '_')
rownames(human2@meta.data) <- paste('human2', make.names(rownames(human2@meta.data), unique = T), sep = '_')
human2@cell.names <- paste('human2', make.names(human2@cell.names, unique = T), sep = '_')

human3 <- CreateSeurat2Object(raw.data = human3.data)
human3 <- NormalizeData(human3)
human3 <- FindVariableGenes(human3, do.plot = F, display.progress = F)
human3 <- ScaleData(human3)
human3@meta.data$human <- "human3"
colnames(human3@data) <- paste('human3', make.names(colnames(human3@data), unique = T), sep = '_')
colnames(human3@raw.data) <- paste('human3', make.names(colnames(human3@raw.data), unique = T), sep = '_')
colnames(human3@scale.data) <- paste('human3', make.names(colnames(human3@scale.data), unique = T), sep = '_')
rownames(human3@meta.data) <- paste('human3', make.names(rownames(human3@meta.data), unique = T), sep = '_')
human3@cell.names <- paste('human3', make.names(human3@cell.names, unique = T), sep = '_')

human4 <- CreateSeurat2Object(raw.data = human4.data)
human4 <- NormalizeData(human4)
human4 <- FindVariableGenes(human4, do.plot = F, display.progress = F)
human4 <- ScaleData(human4)
human4@meta.data$human <- "human4"
colnames(human4@data) <- paste('human4', make.names(colnames(human4@data), unique = T), sep = '_')
colnames(human4@raw.data) <- paste('human4', make.names(colnames(human4@raw.data), unique = T), sep = '_')
colnames(human4@scale.data) <- paste('human4', make.names(colnames(human4@scale.data), unique = T), sep = '_')
rownames(human4@meta.data) <- paste('human4', make.names(rownames(human4@meta.data), unique = T), sep = '_')
human4@cell.names <- paste('human4', make.names(human4@cell.names, unique = T), sep = '_')

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(human1, human2, human3, human4)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 2))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

# Run multi-set CCA
pancreas.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 30)

# calc alignment metric
# cca.nums = 0:20
# alignment.metric = c(0)
# 
# pancreas.cca <- AlignSubspace(pancreas.integrated, reduction.type = "cca",
#                               grouping.var = "human", dims.align = 1:30)
# 
# for (num in 1:20) {
#   val <- CalcAlignmentMetric(pancreas.cca, reduction.use = "cca.aligned",
#                              dims.use = 1:num, grouping.var =  "human")
#   alignment.metric <- append(alignment.metric, val)
# }
# 
# df <- data.frame(nums = cca.nums, vals = alignment.metric)
# ggplot(data = df, mapping = aes(x = nums, y = vals)) + geom_line(color = brewer.pal(7, "Blues")[3], size = 0.8) + geom_point(color = brewer.pal(7, "Blues")[4], size = 1.2) + xlab("Numbers of vectors n") + ylab("Alignment") + ggtitle("Alignment metric n = 20") + theme(plot.title = element_text(hjust = 0.5))


# CC Selection
# MetageneBicorPlot(pancreas.integrated, grouping.var = "human", dims.eval = 1:30)

# for (id in 1:4) {
#   human.id <- strcat(c('human', id))
#   human.cells <- pancreas.integrated@cell.names[pancreas.integrated@meta.data$human == human.id]
#   human.cell.embeddings <- pancreas.integrated@dr$cca@cell.embeddings[human.cells, ]
#   write.table(human.cell.embeddings, file = strcat(c('~/Workspace/Git/OT-integration/data/species/cca/', human.id, '.txt')), sep = '\t', quote = FALSE)
# }

# Run rare non-overlapping filtering
# pancreas.integrated <- CalcVarExpRatio(object = pancreas.integrated, reduction.type = "pca",
#                                        grouping.var = "human", dims.use = 1:10)
# pancreas.integrated <- SubsetData(pancreas.integrated, subset.name = "var.ratio.pca",
#                                   accept.low = 0.5)

# Alignment
pancreas.integrated <- AlignSubspace(pancreas.integrated,
                                     reduction.type = "cca",
                                     grouping.var = "human",
                                     dims.align = 1:13)

# t-SNE and Clustering
pancreas.integrated <- FindClusters(pancreas.integrated, reduction.type = "cca.aligned",
                                    dims.use = 1:13, save.SNN = T, resolution = 0.75)
pancreas.integrated <- RunTSNE(pancreas.integrated,
                               reduction.use = "cca.aligned",
                               dims.use = 1:13,
                               check_duplicates = FALSE)

# Visualization
# TSNEPlot(pancreas.integrated, do.label = T)

# save tsne result
tsne_coords <- pancreas.integrated@dr$tsne@cell.embeddings
write.table(tsne_coords, file = '~/Workspace/code/OT-integration/data/human_pancreas/tsne/seurat_human.txt', sep = '\t', quote = FALSE)

# save clustering result
clusters <- pancreas.integrated@ident
clusters <- as.numeric(clusters)
write.table(clusters, file = strcat(c('~/Workspace/code/OT-integration/data/human_pancreas/clusters/', 'seurat_human_clusters.txt')), sep = '\t', quote = FALSE)

# save cca result
human1.cells <- pancreas.integrated@cell.names[pancreas.integrated@meta.data$human == 'human1']
human1.cca <- pancreas.integrated@dr$cca@cell.embeddings[human1.cells, ]
write.table(human1.cca, file = '~/Workspace/code/OT-integration/data/human_pancreas/cca/human1.txt', sep = '\t', quote = FALSE)

human2.cells <- pancreas.integrated@cell.names[pancreas.integrated@meta.data$human == 'human2']
human2.cca <- pancreas.integrated@dr$cca@cell.embeddings[human2.cells, ]
write.table(human2.cca, file = '~/Workspace/code/OT-integration/data/human_pancreas/cca/human2.txt', sep = '\t', quote = FALSE)

human3.cells <- pancreas.integrated@cell.names[pancreas.integrated@meta.data$human == 'human3']
human3.cca <- pancreas.integrated@dr$cca@cell.embeddings[human3.cells, ]
write.table(human3.cca, file = '~/Workspace/code/OT-integration/data/human_pancreas/cca/human3.txt', sep = '\t', quote = FALSE)

human4.cells <- pancreas.integrated@cell.names[pancreas.integrated@meta.data$human == 'human4']
human4.cca <- pancreas.integrated@dr$cca@cell.embeddings[human4.cells, ]
write.table(human4.cca, file = '~/Workspace/code/OT-integration/data/human_pancreas/cca/human4.txt', sep = '\t', quote = FALSE)

# save pca result
human1 <- RunPCA(human1)
human1.pc <- human1@dr$pca@cell.embeddings
write.table(human1.pc, file = '~/Workspace/code/OT-integration/data/human_pancreas/pca/human1.txt', sep = '\t', quote = FALSE)

human2 <- RunPCA(human2)
human2.pc <- human2@dr$pca@cell.embeddings
write.table(human2.pc, file = '~/Workspace/code/OT-integration/data/human_pancreas/pca/human2.txt', sep = '\t', quote = FALSE)

human3 <- RunPCA(human3)
human3.pc <- human3@dr$pca@cell.embeddings
write.table(human3.pc, file = '~/Workspace/code/OT-integration/data/human_pancreas/pca/human3.txt', sep = '\t', quote = FALSE)

human4 <- RunPCA(human4)
human4.pc <- human4@dr$pca@cell.embeddings
write.table(human4.pc, file = '~/Workspace/code/OT-integration/data/human_pancreas/pca/human4.txt', sep = '\t', quote = FALSE)


