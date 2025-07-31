# 绘制diagram，小提琴图，直方图等，缺少绘制韦恩图

library(ggVennDiagram)
library(tidyverse)
library(ggsignif)
library(Seurat)
library(ggrepel)
library(optparse)
option_list <- list(
  make_option(c("-d", "--distance_csv"),
              type = "character",
              default = "glb_eeisp_distance.csv",
              help = "path to distance CSV file",
              metavar = "path"),
  make_option(c("-t", "--dtg_threshold"),
              type = "numeric",
              default = 0.90,
              help = "DTG threshold [default %default]"),
  make_option(c("-i", "--input_rds"),
              type = "character",
              default = "../checked.rds",
              help = "input Seurat RDS file",
              metavar = "path"),
  make_option(c("-e", "--deg_threshold"),
              type = "numeric",
              default = 0.01,
              help = "DEG adjusted p-value threshold [default %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
distance_csv  <- opt$distance_csv
dtg_threshold <- opt$dtg_threshold
input_rds     <- opt$input_rds
deg_threshold <- opt$deg_threshold


set.seed(123)
################################################################################################

glb_distance <- read.table(distance_csv, sep = '\t', header = TRUE)
seu <- readRDS(input_rds)
Idents(seu) <- seu$celltype
plotEmbedding <- function(embedding){
    p1 <- ggplot(data = embedding, mapping = aes(x = distance_avg, y = distance_sd, color = color)) + 
                             geom_point(alpha = 0.9,size =3) + 
                             theme_bw() +
                             theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   legend.position = "none"
                                   ) +
                             labs(x = 'Average distance', y = 'Standard deviation of distance')
    return(p1)
}

if ("distance" %in% colnames(glb_distance)){    
    print("----------------------- plot: pair -----------------------------")
    threshold <- quantile(glb_distance$distance, dtg_threshold, na.rm = TRUE) # 前10%作为DTG
    glb_distance <- glb_distance[!is.na(glb_distance$distance),]
    glb_distance$color <- ifelse(glb_distance$distance > threshold, "DTGs", "non-DTGs")

    p1 <- ggplot(glb_distance,aes(x=distance, fill = color)) + 
                geom_histogram(binwidth=0.05, alpha=0.9)+
                scale_fill_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2")) +
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.position = c(0.9, 0.9),
                      legend.background = element_rect(fill = NA)) +
                labs(x = 'Distance', y = 'Count') +
                guides(fill = guide_legend(title = NULL))

    ggsave('hist.pdf',p1,dpi = 300, width = 6, height = 5)

    glb_dtgs <- glb_distance[glb_distance$color == 'DTGs',]
    colnames(glb_dtgs)[1] <- 'symbol'

    non_dtgs <- glb_distance[glb_distance$color == 'non-DTGs',]

    # DEGs
    # glb_degs <- read.table('../data/singleCell/glioblastoma/merge_seurat.markers.txt', 
    # 											 sep = ' ', 
    # 											 header = TRUE)
    # symbol_ensembl <- read.table('../../metadata/formatted_genes.txt')
    # colnames(symbol_ensembl) <- c('ensembl', 'symbol')
    # rownames(symbol_ensembl) <- symbol_ensembl$ensembl
    # glb_degs$symbol <- symbol_ensembl[glb_degs$gene,]$symbol
    cluster_value <- unique(seu$celltype)
    markers <- FindAllMarkers(seu, only.pos = TRUE); head(markers)
    length(markers$gene)
    markers <- markers[markers$p_val_adj < deg_threshold,]
    length(markers$gene)
    glb_degs <- markers
    glb_degs$symbol <- glb_degs$gene

    gene_lists <- list(
      DTGs = unique(glb_dtgs$symbol),
      DEGs = unique(na.omit(glb_degs$symbol))
    )
    save(gene_lists, file="gene_lists.Rdata")
    p <- ggVennDiagram(gene_lists,
                       label = "count",     # 显示交集数量
                       label_alpha = 1)     # 填充标签背景
    ggsave("DTGs_vs_DEGs_venn.pdf", plot = p, dpi = 300)


    ##### Comparing DTGs with DEGs

    # Seurat Object that have 
    merge_seurat <- seu
    merge_metadata <- merge_seurat@meta.data
    head(merge_metadata)
    # merge_metadata$CD164 <- merge_seurat@assays$RNA$data['ENSG00000135535',]
    # merge_metadata$DKK3 <- merge_seurat@assays$RNA$data['ENSG00000050165',]
    # merge_metadata$PBK <- merge_seurat@assays$RNA$data['ENSG00000168078',]
    # merge_metadata$MCM4 <- merge_seurat@assays$RNA$data['ENSG00000104738',]

    top_dtgs <- head(glb_distance$Unnamed..0, n = 5)
    for (one_dtg in top_dtgs) {
        merge_metadata[[one_dtg]] <- merge_seurat@assays$RNA$data[one_dtg, ]
        p1 <- ggplot(merge_metadata,
                                 aes(x = celltype,
                                         y = .data[[one_dtg]])) +
            geom_violin(aes(fill = celltype)) +
            geom_jitter(width = 0.15, size = 0.1) +
            geom_signif(comparisons = list(cluster_value),
                                    test = "wilcox.test") +
            labs(x = "",
                     y = "Normalized expression level") +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) +
            guides(fill = "none")
        ggsave(paste0(one_dtg, ".pdf"), p1, dpi = 300, width = 6, height = 3)
    }
} else {
    print("----------------------- plot: multiple -----------------------------")
    pbmc_distance <- glb_distance
    pbmc <- seu
    # pbmc_distance <- read.table('result2/pbmc_all_distance.txt', sep = '\t', header = TRUE)
    colnames(pbmc_distance)[1] <- 'name'
    pbmc_distance$color <- ifelse((pbmc_distance$distance_avg > 
                                   quantile(pbmc_distance$distance_avg, dtg_threshold, na.rm = TRUE)  | pbmc_distance$distance_sd > 
                                   quantile(pbmc_distance$distance_sd, dtg_threshold, na.rm = TRUE)), 'DTGs', 'non-DTGs')
    
    pbmc_DTGs <- pbmc_distance[pbmc_distance$color == 'DTGs',]$name
    pbmc_distance$label <- ifelse(pbmc_distance$name  == head(pbmc_DTGs, n=5), pbmc_distance$name, NA)
    pbmc_distance$label <- ifelse(pbmc_distance$name %in% pbmc_DTGs, pbmc_distance$name, "")
    
    p1 <- plotEmbedding(pbmc_distance) +
                                 geom_hline(yintercept = 0.4, linetype = "dashed", color = "red", size = 0.5) +
                                 scale_color_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2"))  +
                                 geom_vline(xintercept=3.35,linetype = "dashed", color = "red", size = 0.5) +
                                 geom_text_repel(aes(label = label), 
                                                 size = 2,
                                                 max.overlaps = Inf,
                                                 nudge_x = 0.5,
                                                 nudge_y = 0.5,
                                                 box.padding = 0.5, 
                                                 segment.size = 0.1) + 
                                theme(legend.position = "right")
    
    ggsave('mean_sd.pdf',p1,dpi = 300, width = 4.8, height = 4)
    
    pbmc.marker <- FindAllMarkers(pbmc, only.pos = TRUE)
    pbmc.marker <- pbmc.marker[pbmc.marker$p_val_adj < deg_threshold,]
    gene_lists <- list(DTGs = unique(pbmc_DTGs),DEGs = unique(pbmc.marker$gene))
    save(gene_lists, file="gene_lists.Rdata")
    p <- ggVennDiagram(gene_lists,label = "count",label_alpha = 1)
    ggsave("DTGs_vs_DEGs_venn.pdf", plot = p, dpi = 300)
    # write.table(pbmc.marker,'markers.txt', quote = FALSE, row.names = FALSE)
}