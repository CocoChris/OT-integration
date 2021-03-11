library(riverplot)

prop = readRDS('gamma_human_pancreas.Rda')
prop = as.matrix(prop)

prop = prop[,15:42]

sum = apply(prop, 1, sum)
sum = rep(sum, 28)
dim(sum) = c(14,28)

prop = prop / sum

max = apply(prop, 2, max)
max = rep(max, each = 14)
dim(max) = c(14,28)

prop[prop < 0.8*max] = 0

p = as.numeric(prop) * 20

edges = data.frame(N1 = c(paste0('human2_', rep(c('schwann', 'endothelial', 'ductal', 'quiescent_stellate', 'beta', 'activated_stellate', 'delta', 'alpha', 'T_cell', 'gamma', 'macrophage', 'epsilon', 'acinar', 'mast'), each = 14)),
                          paste0('human3_', rep(c('schwann', 'endothelial', 'ductal', 'quiescent_stellate', 'beta', 'activated_stellate', 'delta', 'alpha', 'T_cell', 'gamma', 'macrophage', 'epsilon', 'acinar', 'mast'), each = 14))),
                   N2 = c(paste0('C', rep(1:14, 28))),
                   Value = p,
                   stringsAsFactors = F)

nodes = data.frame(ID = unique(c(edges$N1, edges$N2)), stringsAsFactors = FALSE)

nodes$x = c(rep(65, 14), rep(67, 14), rep(66, 14))
nodes$y = c(65:78, 65:78, 65:78)
nodes$y = -2*nodes$y

rownames(nodes) = nodes$ID

# 添加颜色
library(RColorBrewer)
palette = c(paste0(brewer.pal(8, "Set2"), "60"), paste0(brewer.pal(8, "Pastel2"), "60"))

# 对每个节点生成相应的格式
styles = lapply(nodes$y, function(n) {
  list(col = palette[-(n/2)-64], lty = 0, textcol = "black", srt = "0")
})
names(styles) = nodes$ID


# 以list的结构保存一遍调用
rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
plot(rp, plot_area = 0.95, yscale=0.06)
