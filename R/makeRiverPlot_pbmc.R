library(riverplot)

prop = read.table('~/Workspace/Git/OT-integration/data/techs/pbmc/ot/X.txt')
prop = as.matrix(prop)

sum = apply(prop, 1, sum)
sum = rep(sum, 15)
dim(sum) = c(6,15)

prop = prop / sum

max = apply(prop, 2, max)
max = rep(max, each = 6)
dim(max) = c(6,15)

prop[prop < 0.5*max] = 0
p = as.numeric(prop) * 10
# 构造连接节点的数据框
edges = data.frame(N1 = c(rep(c('Bcell', 'CD4', 'CD8', 'DC', 'NK', 'Myeloid'), each = 6), paste0(rep('10X_', each = 54), rep(1:9, each = 6))),
                   N2 = rep(1:6, 15),
                   Value = p,
                   stringsAsFactors = F)

nodes = data.frame(ID = unique(c(edges$N1, edges$N2)), stringsAsFactors = FALSE)
nodes$x = c(rep(65, 6), rep(67, 9), rep(66, 6))
nodes$y = c(65:70, 65:73, 65:70)
nodes$y = -nodes$y
#
rownames(nodes) = nodes$ID

# 添加颜色
library(RColorBrewer)
# 后面加调淡颜色
palette = paste0(brewer.pal(9, "Set1"), "60")

# 对每个节点生成相应的格式
styles = lapply(nodes$y, function(n) {
  list(col = palette[-n-64], lty = 0, textcol = "black", srt = "0")
})
names(styles) = nodes$ID

# 以list的结构保存一遍调用
rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")

plot(rp, plot_area = 0.95, yscale=0.06)

