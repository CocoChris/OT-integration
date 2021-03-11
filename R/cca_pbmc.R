################################################################################
### cca workflow for pbmc datasets
################################################################################

library(Seurat)
library(Seurat2)
library(Matrix)
library(RColorBrewer)
library(pracma)

# Read in human expression matrices
pbmc.seqwell.data <- read.delim("~/Workspace/code/OT-integration/data/pbmc/pbmc_seqwell.txt", header = TRUE, row.names = 1, check.names = FALSE)
pbmc.10x.data <- read.delim("~/Workspace/code/OT-integration/data/pbmc/pbmc_10x.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Convert to sparse matrices for efficiency
pbmc.seqwell.data <- as(as.matrix(pbmc.seqwell.data), "dgCMatrix")
pbmc.10x.data <- as(as.matrix(pbmc.10x.data), "dgCMatrix")

# Create and setup Seurat objects for each dataset
pbmc.seqwell <- CreateSeurat2Object(raw.data = pbmc.seqwell.data)
pbmc.seqwell <- NormalizeData(pbmc.seqwell)
pbmc.seqwell <- FindVariableGenes(pbmc.seqwell, do.plot = F, display.progress = F)
pbmc.seqwell <- ScaleData(pbmc.seqwell)
pbmc.seqwell@meta.data$tech <- "seqwell"

pbmc.10x <- CreateSeurat2Object(raw.data = pbmc.10x.data)
pbmc.10x <- NormalizeData(pbmc.10x)
pbmc.10x <- FindVariableGenes(pbmc.10x, do.plot = F, display.progress = F)
pbmc.10x <- ScaleData(pbmc.10x)
pbmc.10x@meta.data$tech <- "10x"

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(pbmc.seqwell, pbmc.10x)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 3000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

# Run multi-set CCA
pbmc.integrated <- RunCCA(object = pbmc.seqwell, object2 = pbmc.10x, genes.use = genes.use, num.cc = 30)

# calc alignment metric
cca.nums = 0:30
alignment.metric = c(0)

pbmc.cca <- AlignSubspace(pbmc.integrated, reduction.type = "cca",
                          grouping.var = "tech", dims.align = 1:30)

for (num in 1:30) {
  val <- CalcAlignmentMetric(pbmc.cca, reduction.use = "cca.aligned",
                             dims.use = 1:num, grouping.var =  "tech")
  alignment.metric <- append(alignment.metric, val)
}

df <- data.frame(nums = cca.nums, vals = alignment.metric)
ggplot(data = df, mapping = aes(x = nums, y = vals)) + geom_line(color = brewer.pal(7, "Blues")[3], size = 0.8) + geom_point(color = brewer.pal(7, "Blues")[4], size = 1.2) + xlab("Numbers of vectors n") + ylab("Alignment") + ggtitle("Alignment metric n = 20") + theme(plot.title = element_text(hjust = 0.5))

# save cca results
pbmc.seqwell.cells <- pbmc.integrated@cell.names[pbmc.integrated@meta.data$tech == 'seqwell']
pbmc.seqwell.cca <- pbmc.integrated@dr$cca@cell.embeddings[pbmc.seqwell.cells, ]
write.table(pbmc.seqwell.cca, file = '~/Workspace/code/OT-integration/data/pbmc/cca/seqwell.txt', sep = '\t', quote = FALSE)

pbmc.10x.cells <- pbmc.integrated@cell.names[pbmc.integrated@meta.data$tech == '10x']
pbmc.10x.cca <- pbmc.integrated@dr$cca@cell.embeddings[pbmc.10x.cells, ]
write.table(pbmc.10x.cca, file = '~/Workspace/code/OT-integration/data/pbmc/cca/10x.txt', sep = '\t', quote = FALSE)
