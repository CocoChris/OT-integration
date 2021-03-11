library(mclust)
library(pracma)
library(ggplot2)
library(RColorBrewer)

# data.human <- read.table("~/Workspace/Git/OT-integration/data/species/cca/human.txt")
# data.human <- as.matrix(data.human)
# data.mouse <- read.table("~/Workspace/Git/OT-integration/data/species/cca/mouse.txt")
# data.mouse <- as.matrix(data.mouse)
# 
# data <- rbind(data.human, data.mouse)

# data.human1 <- read.table("~/Workspace/Git/OT-integration/data/species/cca/human1.txt")
# data.human1 <- as.matrix(data.human1)
# data.human2 <- read.table("~/Workspace/Git/OT-integration/data/species/cca/human2.txt")
# data.human2 <- as.matrix(data.human2)
# data.human3 <- read.table("~/Workspace/Git/OT-integration/data/species/cca/human3.txt")
# data.human3 <- as.matrix(data.human3)
# data.human4 <- read.table("~/Workspace/Git/OT-integration/data/species/cca/human4.txt")
# data.human4 <- as.matrix(data.human4)
# 
# data <- rbind(data.human1, data.human2, data.human3, data.human4)

data.seqwell <- read.table("~/Workspace/Git/OT-integration/data/techs/pbmc/cca/seqwell.txt")
data.seqwell <- as.matrix(data.seqwell)
data.10x <- read.table("~/Workspace/Git/OT-integration/data/techs/pbmc/cca/10x.txt")
data.10x <- as.matrix(data.10x)

data <- rbind(data.seqwell, data.10x)

mod = Mclust(data[,1:20], G = 2:15, , modelNames = "VVV")
plot(mod$BIC)

clusters.nums = 2:15
bic = as.vector(mod$BIC)

df <- data.frame(nums = clusters.nums, vals = bic)

ggplot(data = df, mapping = aes(x = nums, y = vals)) + geom_line(color = brewer.pal(7, "Blues")[3], size = 0.8) + geom_point(color = brewer.pal(7, "Blues")[4], size = 1.2) + xlab("Numbers of clusters C") + ylab("BIC") + ggtitle("BIC metric   C = 9") + theme(plot.title = element_text(hjust = 0.5))
