library(igraph)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = "data/TFs_list.RData",
                            help = "Path to TFs_list RData file [default %default]", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
TFs_list_RData <- opt$input

load(TFs_list_RData)

## 1. 取出边列表 -------------------------------------------------
el <- TFs_list$FOSF_RegMTF_Cor_EnTFs[ , c("TF", "Target", "Weight", "Regulation")]
colnames(el) <- c("from", "to", "weight", "regulation")

## 2. 建图 -------------------------------------------------------
g <- graph_from_data_frame(el, directed = TRUE)

## 3. 给边加颜色：Positive = 绿色，Negative = 红色 --------------
E(g)$color <- ifelse(E(g)$regulation == "Positive",
                     "forestgreen", "firebrick")

## 4. 给节点加类型属性（TF vs Target） --------------------------
V(g)$type <- ifelse(V(g)$name %in% unique(el$from), "TF", "Target")

## 5. 简单布局 & 绘图 -------------------------------------------
set.seed(123)
# layout <- layout_randomly(g) # Places the vertices completely randomly
# plot(g,
#      layout          = layout,
#      vertex.size     = 5,
#      vertex.label.cex= .7,
#      vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
#      edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
#      edge.arrow.size = .3,
#      vertex.label    = V(g)$name,
#      vertex.label.dist = .5,
#      main = "TF-Target Regulatory Network(randomly)")

layout <- layout_in_circle(g) # Deterministic layout that places the vertices on a circle
plot(g,
     layout          = layout,
     vertex.size     = 5,
     vertex.label.cex= .7,
     vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
     edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
     edge.arrow.size = .3,
     vertex.label    = V(g)$name,
     vertex.label.dist = .5,
     main = "TF-Target Regulatory Network(circle)")

# layout <- layout_on_sphere(g) # Deterministic layout that places the vertices evenly on the surface of a sphere
# plot(g,
#      layout          = layout,
#      vertex.size     = 5,
#      vertex.label.cex= .7,
#      vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
#      edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
#      edge.arrow.size = .3,
#      vertex.label    = V(g)$name,
#      vertex.label.dist = .5,
#      main = "TF-Target Regulatory Network(sphere)")

# layout <- layout_with_drl(g) # The Drl (Distributed Recursive Layout) algorithm for large graphs
# plot(g,
#      layout          = layout,
#      vertex.size     = 5,
#      vertex.label.cex= .7,
#      vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
#      edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
#      edge.arrow.size = .3,
#      vertex.label    = V(g)$name,
#      vertex.label.dist = .5,
#      main = "TF-Target Regulatory Network(drl)")

# layout <- layout_with_fr(g) # Fruchterman-Reingold force-directed algorithm
# plot(g,
#      layout          = layout,
#      vertex.size     = 5,
#      vertex.label.cex= .7,
#      vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
#      edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
#      edge.arrow.size = .3,
#      vertex.label    = V(g)$name,
#      vertex.label.dist = .5,
#      main = "TF-Target Regulatory Network(fr)")

# layout <- layout_with_kk(g) # Kamada-Kawai force-directed algorithm
# plot(g,
#      layout          = layout,
#      vertex.size     = 5,
#      vertex.label.cex= .7,
#      vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
#      edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
#      edge.arrow.size = .3,
#      vertex.label    = V(g)$name,
#      vertex.label.dist = .5,
#      main = "TF-Target Regulatory Network(kk)")

# layout <- layout_with_lgl(g) # The LGL (Large Graph Layout) algorithm for large graphs
# plot(g,
#      layout          = layout,
#      vertex.size     = 5,
#      vertex.label.cex= .7,
#      vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
#      edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
#      edge.arrow.size = .3,
#      vertex.label    = V(g)$name,
#      vertex.label.dist = .5,
#      main = "TF-Target Regulatory Network(lgl)")

layout <- layout_as_tree(g) # Reingold-Tilford tree layout, useful for (almost) tree-like graphs
plot(g,
     layout          = layout,
     vertex.size     = 5,
     vertex.label.cex= .7,
     vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
     edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
     edge.arrow.size = .3,
     vertex.label    = V(g)$name,
     vertex.label.dist = .5,
     main = "TF-Target Regulatory Network(tree)")

layout <- layout_nicely(g)
plot(g,
     layout          = layout,
     vertex.size     = 5,
     vertex.label.cex= .7,
     vertex.color    = ifelse(V(g)$type == "TF", "gold", "skyblue"),
     edge.width      = abs(E(g)$weight) * 50,   # 权重映射到线宽
     edge.arrow.size = .3,
     vertex.label    = V(g)$name,
     vertex.label.dist = .5,
     main = "TF-Target Regulatory Network(nicely)")

dev.off()