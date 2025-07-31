library(tidyverse)
library(ggsignif)
library(Seurat)

setwd('/data/work/Single-Cell-Pipeline/output/GRN-gene2role/test')
set.seed(123)

my_colors <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
		"#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
		"#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
		"#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

glb_cluster_distance <- read.table('/data/work/Single-Cell-Pipeline/output/GRN-gene2role/test/glb_distance_0.csv', sep = '\t', header = TRUE)

plotEmbedding <- function(embedding){
	p1 <- ggplot(data = embedding, mapping = aes(x = nan_ratio, y = mean, color = as.factor(community), size = size)) + 
			 				 geom_point() + 
			 				 scale_color_manual(values = my_colors) + 
			 				 theme_bw() +
			 				 theme(panel.grid.major = element_blank(),
			 				 	   panel.grid.minor = element_blank(),
			 				 	   legend.position = "right"
			 				 	   ) +
			 				 labs(x = 'NA%', y = 'Mean distance',size = "No. of genes", color = "Gene module")
	return(p1)
}

p1 <- plotEmbedding(glb_cluster_distance)

ggsave('glb_gene_modules.pdf',p1,dpi = 300, width = 6, height = 5)


# library(clusterProfiler)
# library('org.Hs.eg.db')
# glb_index <- read.table('/data/work/GRN-gene2role/test3/result/splitMatrix/index_tracker.tsv', sep = '\t', header = TRUE)
# rownames(glb_index) <- glb_index$X0h

# glb_lovain <- read.table('/data/work/GRN-gene2role/plot/rep_data/rep_data/result4/0h_edgelist_lovain.csv', sep = '\t', header = TRUE)
# glb_lovain_0_5 <- glb_lovain[glb_lovain$community %in% c(0,5),]
# glb_lovain_0_5$gene_name <- glb_index[glb_lovain_0_5$gene,]$X

# glb_lovain_0_5_list <- split(glb_lovain_0_5$gene_name, glb_lovain_0_5$community)

# glb_lovain_0_5_BP <- compareCluster(glb_lovain_0_5_list,
#                              fun=enrichGO,
#                              OrgDb = 'org.Hs.eg.db',
#                              keyType = 'SYMBOL',
#                              ont = 'BP', 
#                              pAdjustMethod = 'BH',
#                              qvalueCutoff = 0.05)

# p1 <- dotplot(glb_lovain_0_5_BP, showCategory=5) + 
#     scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + 
#     theme(axis.title.x = element_blank())

# ggsave('result3/glb_0_5_BP.pdf', p1, dpi = 200, units = 'cm', width = 15, height = 10)