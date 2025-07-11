# Reference: https://github.com/swaruplabUCI/Relapse-to-cocaine-seeking-is-regulated-by-medial-habenula-Nr4a2/blob/main/hdWGCNA/hdWGCNA.Rmd

# conda activate voyager
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(ggrastr)
library(patchwork)
theme_set(theme_cowplot())
set.seed(12345)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(harmony)

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# utility script with plotting code
source('/dfs7/swaruplab/smorabit/analysis/scWGCNA/bin/spatial_functions.R')

# set up directories
setwd('/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/Nurr2c_vs_GFP/revision/hdWGCNA')
fig_dir <- 'figures/'
data_dir <- 'data/'

# load seurat object
seurat_obj <- readRDS(file='/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/Nurr2c_vs_GFP/revision/data/harmony_annotated_integration.rds')

# re-load hdWGCNA seurat object (after analysis is complete)
seurat_obj <- readRDS(file = paste0(data_dir, 'harmony_annotated_hdWGCNA.rds'))




#---------------------------------------------------------#
# Part 1: Set up seurat object & create metacells
#---------------------------------------------------------#

# subset by medial habenula neurons only:
seurat_obj <- subset(seurat_obj, cell_type %in% c('MHb-Neuron') & Group %in% c('Nurr2c', 'GFP'))
seurat_obj$annotation <- droplevels(seurat_obj$annotation)
seurat_obj$Sample <- droplevels(seurat_obj$Sample)

# set up the hdWGCNA experiment
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  #group.by = 'cell_type',
  wgcna_name = "MHb_metacell"
)
length(GetWGCNAGenes(seurat_obj))

# compute metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c('annotation', "cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'annotation',
  reduction='harmony',
  target_metacells=100,
  min_cells = 25
)

# normalize metacells
seurat_obj <- NormalizeMetacells(seurat_obj)

#---------------------------------------------------------#
# Part 2: Construct co-expression network
#
# Figure S7B (dendrogram)
#---------------------------------------------------------#

# setup the gene expression matrix for coex network analysis
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "MHb-Neuron",
  group.by = "cell_type"
)

# test the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)
plot_list <- PlotSoftPowers(seurat_obj)

#assemble with patchwork
pdf(paste0(fig_dir, 'softpower.pdf'), width=12, height=8)
wrap_plots(plot_list, ncol=2)
dev.off()

# build the coex net
seurat_obj <- ConstructNetwork(
    seurat_obj, 
    tom_name='MHb_metacell', 
    overwrite_tom=TRUE
)

# plot the dendrogram 
pdf(paste0(fig_dir, "dendro.pdf"),height=3, width=6)
PlotDendrogram(seurat_obj, main='MHb hdWGCNA Dendrogram')
dev.off()

# get the modules table
modules <- GetModules(seurat_obj)
table(modules$module)

# checking some individual genes
subset(modules, gene_name == 'Nr4a2')
subset(modules, gene_name == 'Gabra2')
subset(modules, gene_name == 'Sox5')

#---------------------------------------------------------#
# Part 2.1: Rename modules and set up color scheme
#---------------------------------------------------------#

# run RenameModules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "M"
)

library(MetBrewer)

modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mod_colors <- dplyr::select(modules, c(module, color)) %>%
  distinct %>% arrange(module) %>% .$color
n_colors <- length(mod_colors) -1

new_colors <- paste0(met.brewer("Cross", n=n_colors, type='discrete'))
seurat_obj <- ResetModuleColors(seurat_obj, new_colors)


#---------------------------------------------------------#
# Part 3: Calculate Module Eigengenes
#
# Module Feature plots (Figure 6B)
# ME violin plots (Figure 6D)
#---------------------------------------------------------#

# compute module eigengenes
seurat_obj$Assignment <- droplevels(seurat_obj$Assignment)
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="Assignment" # snRNAseq batch
)

# compute module connectivity:
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'MHb-Neuron'
)

# plot the module feature plots
plot_list <- ModuleFeaturePlot(
  seurat_obj, 
  order=TRUE, 
  raster=TRUE, 
  restrict_range=FALSE, 
  raster_dpi=400, 
  raster_scale=0.5, 
  point_size=1,
  reduction = 'neuron_umap'
  )

plot_list <- lapply(plot_list, function(x){
  x + NoLegend() + theme(plot.margin=margin(0,0,0,0))
})

pdf(paste0(fig_dir, "MHb_featureplot_MEs.pdf"),height=6, width=12)
wrap_plots(plot_list, ncol=5)
dev.off()

# module violin plots
MEs <- GetMEs(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- custom_vln(
    seurat_obj,
    features = mods,
    group.by = 'annotation',
    split.by = 'Group',
    groups = c('MHb-1', 'MHb-2', 'MHb-3', 'MHb-4', 'MHb-5'),
    add_boxplot=FALSE,
    split_colors=c('darkgoldenrod3', 'hotpink3'),
    add_colorbar=TRUE,
    plot_ymin = NA
  )

pdf(paste0(fig_dir, 'MHb_hME_vln_stack.pdf'), width=3, height=5)
p
dev.off()


#---------------------------------------------------------#
# Part 4: GO term enrichment 
#
# Enrichment dotplot (Figure 6C)
#---------------------------------------------------------#

library(enrichR)

dbs <-c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021','WikiPathway_2021_Human', 'KEGG_2021_Human')

# compute GO terms:
seurat_obj <- RunEnrichr(seurat_obj, dbs=dbs, max_genes=Inf)

# inspect the enrichr table:
enrichr_df <- GetEnrichrTable(seurat_obj) %>% subset(P.value < 0.05)

write.table(enrichr_df, quote=FALSE, row.names=FALSE, sep='\t', file=paste0(data_dir, 'MHb_hdWGCNA_enrichr.tsv'))

# save the significant results
subset(enrichr_df, P.value < 0.05) %>% 
write.table(quote=FALSE, row.names=FALSE, sep='\t', file=paste0(data_dir, 'MHb_hdWGCNA_enrichr_signif.tsv'))

# get module info and colors
modules <- GetModules(seurat_obj)
color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  rename(c(group=module, colour=color))

# plot selected GO terms
combined_output <- GetEnrichrTable(seurat_obj)
combined_output $ngenes <- unlist(lapply(strsplit(combined_output $Genes, ';'), function(x){length(x)}))
combined_output <- subset(combined_output , ngenes >= 3)

selected_terms <- read.delim('data/MHb_hdWGCNA_enrichr_selected.txt', sep='\t', header=1)

# subset selected terms
combined_output$ngenes <- unlist(lapply(strsplit(combined_output$Genes, ';'), function(x){length(x)}))
selected_terms <- subset(combined_output, Term %in% selected_terms$Term & P.value < 0.05 & ngenes > 2)
selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mods
)

# set max pval for plotting
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove GO Term ID to shorten the term
selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")
selected_terms <- selected_terms %>%
  arrange(group)

# order the terms
selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

# GO Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = Term, color =logp, size=log(Combined.Score))) +
  geom_point() +
  scale_color_stepsn(colors=rev(magma(256))) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )+ labs(
    color = bquote("-log"[10]~"(P)"),
    size= bquote("log"[10]~"(Enrich)")
  )


# # make the colorbar as its own heatmap
color_df$var <- 1
cp <- color_df$colour; names(cp) <- color_df$group
colorbar <- color_df %>%
  subset(group %in% unique(selected_terms$group)) %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_fill_manual(values=cp) +
  coord_equal() +
  NoLegend() + RotatedAxis() +
  theme(
    plot.title=element_blank(),
    axis.line=element_blank(),
    axis.ticks.y =element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.margin=margin(0,0,0,0),
  )

pdf(paste0(fig_dir, 'selected_go_terms_bubbleplot.pdf'), width=7.3, height=7)
p / colorbar
dev.off()


#---------------------------------------------------------#
# Part 5: Network plotting 
#
# Figure 6A and S7A
#---------------------------------------------------------#

library(igraph)
library(reshape2)

# individual module networks
ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  #label_center=TRUE,
  outdir = paste0(fig_dir, 'MHb_hubNetworks/')
)

# co-expression UMAP
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 3,
  n_neighbors=10,
  min_dist=0.3,
  spread=3,
  supervised=TRUE,
  target_weight=0.5
)


# additional genes to label aside from hub genes
label_genes <- c('Gabra2', 'Nr4a2')

pdf(paste0(fig_dir, 'MHb_hubgene_umap_igraph.pdf'), width=10, height=10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, 
  label_hubs=3,
  label_genes = label_genes
)
dev.off()


#---------------------------------------------------------#
# Part 6: Save the results
#---------------------------------------------------------#

saveRDS(seurat_obj, file = paste0(data_dir, 'harmony_annotated_hdWGCNA.rds'))

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
write.csv(modules, file=paste0(data_dir, 'modules.csv'), quote=FALSE, row.names=FALSE)





#---------------------------------------------------------#
# Perform DME analysis to compare NURR2C & GFP
#---------------------------------------------------------#

# loop over each MHb cluster
groups <- as.character(unique(seurat_obj$annotation))
DMEs <- data.frame()
for(cur_group in groups){

  # get barcodes for each conditiohn
  group1 <- seurat_obj@meta.data %>% subset(annotation == cur_group & Group == "Nurr2c") %>% rownames
  group2 <- seurat_obj@meta.data %>% subset(annotation == cur_group & Group == "GFP") %>% rownames

  # perform DME analysis
  cur_DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    pseudocount.use=0.01
  )
  cur_DMEs$group <- cur_group 
  DMEs <- rbind(DMEs, cur_DMEs)

}

#---------------------------------------------------------#
# Plot DME heatmap (Figure 6E)
#---------------------------------------------------------#

# plot the result as a heatmap
maxval <- 0.5; minval <- -0.5
plot_df <- DMEs
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)

plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'black')
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)

# set factor levels:
plot_df$module <- factor(as.character(plot_df$module), levels=mods)

p <- plot_df %>% 
  ggplot(aes(x=group, y=fct_rev(module), fill=avg_log2FC)) +
  geom_tile() +
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) +
  scale_fill_gradient2(low='hotpink3', mid='grey95', high='goldenrod3') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') 

# Plot the result
pdf(paste0(fig_dir, 'DME_heatmap.pdf'),height=2.5, width=4)
p 
dev.off()




#---------------------------------------------------------------------------#
# Make a table of primary and secondary Nr4a2 target genes 
#---------------------------------------------------------------------------#

cur_tf <- 'Nr4a2'

tf_nets <- dir('../tf_net/data/tf_nets/')
tf_nets <- tf_nets[grepl('importance', tf_nets)]

# parameters for regulons
n_tfs <- 5
importance_thresh <- 0.001
combined_output <- data.frame()

primary_genes <- list()
secondary_genes <- list()

for(cur_net_file in tf_nets){

  print(cur_net_file)

  # load the tf-gene table
  importance_df <- read.csv(paste0('../tf_net/data/tf_nets/', cur_net_file ))

  tmp <- strsplit(cur_net_file, '_')[[1]]
  cur_celltype <- tmp[2]
  cur_group <- tmp[3]

  #---------------------------------------------------------------------------#
  # Define the TF regulons
  #---------------------------------------------------------------------------#

  regulons <- importance_df %>% 
    subset(Gain > importance_thresh) %>% 
    group_by(gene) %>%
    slice_max(order_by=Gain, n=n_tfs) %>% 
    ungroup()

  # compute the degree for each TF:
  tf_degrees <- table(regulons$tf)

  #---------------------------------------------------------------------------#
  # Get the primary & secondary targets of Nr4a2
  #---------------------------------------------------------------------------#

  # primary target genes 
  cur_primary<- regulons %>% 
    subset(tf == cur_tf) 

  # which of these pimary target genes are tfs?
  cur_primary_tfs <- cur_primary %>% 
    subset(gene %in% unique(regulons$tf)) %>% .$gene

  cur_tfs <- unique(c(cur_tf, cur_primary_tfs))

  # get the regulons for these TFs:
  cur_secondary <- subset(regulons, tf %in% cur_primary_tfs)
  cur_secondary_tfs <- cur_primary %>% 
    subset(gene %in% unique(regulons$tf)) %>% .$gene

  primary_genes[[paste0(cur_celltype, '_', cur_group)]] <- cur_primary$gene 
  secondary_genes[[paste0(cur_celltype, '_', cur_group)]] <- cur_secondary$gene

}

#---------------------------------------------------------------------------#
# Test overlaps between Nr4a2 targets and modules
#---------------------------------------------------------------------------#

# get nr4a2 targets
primary <- primary_genes[['MHb-Neuron_Nurr2c']]
secondary <- secondary_genes[['MHb-Neuron_Nurr2c']]
secondary <- setdiff(secondary, primary)

# get modules
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# genome size for the overlap test
genome.size <- nrow(seurat_obj)

# run overlap for each module
overlap_df <- data.frame()
plot_df <- data.frame()
for(cur_mod in mods){

  cur_genes <- subset(modules, module == cur_mod) %>% .$gene_name

  cur_primary <- intersect(cur_genes, primary)
  cur_secondary <- intersect(cur_genes, secondary)
  cur_other <- setdiff(cur_genes, unique(c(cur_primary, cur_secondary)))

  cur_df <- data.frame(
    module = cur_mod,
    target_type = c('primary', 'secondary', 'other'),
    n = c(
      length(cur_primary),
      length(cur_secondary),
      length(cur_other)
    )
  )
  cur_df$percentage <- cur_df$n / sum(cur_df$n)
  plot_df <- rbind(plot_df, cur_df)

  cur_overlap_primary <- testGeneOverlap(newGeneOverlap(
      primary,
      cur_genes,
      genome.size=genome.size
  ))
  cur_overlap_secondary <- testGeneOverlap(newGeneOverlap(
      secondary,
      cur_genes,
      genome.size=genome.size
  ))

  cur_overlap <- data.frame(
    'odds.ratio' = c(cur_overlap_primary@odds.ratio, cur_overlap_secondary@odds.ratio),
    'pval' = c(cur_overlap_primary@pval, cur_overlap_secondary@pval),
    'Jaccard' = c(cur_overlap_primary@Jaccard, cur_overlap_secondary@Jaccard),
    'size_intersection' = c(length(cur_overlap_primary@intersection), length(cur_overlap_secondary@intersection)),
    'target_type' = c('primary', 'secondary'),
    'module' = cur_mod
  )
  overlap_df <- rbind(overlap_df, cur_overlap)

}
overlap_df$fdr <- p.adjust(overlap_df$pval, method='fdr')
overlap_df$shape <- ifelse(overlap_df$fdr < 0.05, 21, 4)

#---------------------------------------------------------------------------#
# Make a stacked bar chart to show the overlaps (Fig 6F)
#---------------------------------------------------------------------------#

# plot the stacked bar chart 
plot_df$target_type <- factor(as.character(plot_df$target_type), levels=c('other', 'secondary', 'primary'))

p1 <- plot_df %>%
  ggplot(aes(x = percentage, y=fct_rev(module), fill=target_type)) + 
  geom_bar(position='stack', stat='identity') + 
  # scale_fill_manual(values=amylo_cp) + 
  geom_text(aes(label=abs(n)), position = position_stack(vjust=0.5)) +
  #xlim(-plot_max, plot_max) +
#  scale_x_continuous(expand = c(0, 0), limits = c(-plot_max, plot_max)) + 
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_blank()
  ) + 
  xlab("Proportion") + NoLegend() + 
   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) 
  #facet_wrap(~direction, ncol=2) + NoLegend()
#      ylab(bquote("-log"[10]~"(Adj. P-value)")) +

pdf(paste0(fig_dir, 'Nr4a2_target_modules_barplot.pdf'), width=3, height=4.5, useDingbats=FALSE)
p1
dev.off()

library(ggbreak)

p2 <- overlap_df %>% 
  ggplot(aes(x=fct_rev(module), y = odds.ratio, size=-log10(fdr), color=target_type, fill=target_type)) + 
  geom_hline(yintercept=1, color='lightgrey', linetype='dashed') +
     geom_linerange(aes(x=fct_rev(module), xmax=fct_rev(module), ymin=0, ymax=odds.ratio), size=0.5, position=position_dodge(width=1)) +
  geom_point(position=position_dodge(width=1), shape=overlap_df$shape, color='black') +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + ylab('') + xlab('')

pdf(paste0(fig_dir, 'Nr4a2_target_modules_lollipop_legend.pdf'), width=4, height=4.5, useDingbats=FALSE)
p2
dev.off()





modules <- GetModules(seurat_obj, wgcna_name)
mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct()
mod_cp <- mod_colors$color; names(mod_cp) <- as.character(mod_colors$module)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

grey_genes <- subset(modules, module == 'grey') %>% .$gene_name

# get a list of all TFs:
importance_df <- read.csv(paste0('../tf_net/data/tf_nets/TFnet_MHb-Neuron_Nurr2c_importance.csv'))
tf_list <- unique(importance_df$tf)

grey_tfs <- tf_list[tf_list %in% grey_genes]

# get just the kME values from the modules table
kMEs <- modules[,4:ncol(modules)]
kMEs <- kMEs[,colnames(kMEs) != "kME_grey"]


reassigned <- sapply(grey_tfs, function(cur_gene){
  cur_kMEs <- kMEs[cur_gene,]
  max_kME <- max(cur_kMEs)
  if(max_kME < 0){
    return('kME_grey')
  }
  colnames(kMEs)[which(cur_kMEs == max_kME)]
})

# add the reassigned modules to the modules table
reassigned <- do.call(rbind, strsplit(reassigned, 'kME_'))[,2]
reassigned <- factor(as.character(reassigned), levels=levels(modules$module))

# new colors:
reassigned_colors <- as.character(mod_cp[as.character(reassigned)])

# reassign modules and colors
modules[grey_tfs,'module'] <- reassigned
modules[grey_tfs,'color'] <- reassigned_colors

subset(modules, gene_name %in% grey_tfs) %>% .$module %>% table

seurat_obj <- SetModules(seurat_obj, modules)





modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mod_colors <- dplyr::select(modules, c(module, color)) %>%
  distinct %>% arrange(module) %>% .$color
cp <- mod_colors; names(cp) <- mods
mods <- mods[mods != 'grey']

# get a list of all TFs:
importance_df <- read.csv(paste0('../tf_net/data/tf_nets/TFnet_MHb-Neuron_Nurr2c_importance.csv'))
tf_list <- unique(importance_df$tf)

cur_tf <- 'Nr4a2'

tf_nets <- dir('../tf_net/data/tf_nets/')
tf_nets <- tf_nets[grepl('importance', tf_nets)]
tf_nets <- c('TFnet_MHb-Neuron_Nurr2c_importance.csv')
n_tfs <- 10
importance_thresh <- 0.001
network_df <- data.frame()

for(cur_net_file in tf_nets){

  print(cur_net_file)

  tmp <- strsplit(cur_net_file, '_')[[1]]
  cur_celltype <- tmp[2]
  cur_group <- tmp[3]

  importance_df <- read.csv(paste0('../tf_net/data/tf_nets/', cur_net_file ))

  #---------------------------------------------------------------------------#
  # Define the TF regulons
  #---------------------------------------------------------------------------#

  regulons <- importance_df %>% 
    subset(Gain > importance_thresh) %>% 
    group_by(gene) %>%
    slice_max(order_by=Gain, n=n_tfs) %>% 
    ungroup()

  # compute the degree for each TF:
  tf_degrees <- table(regulons$tf)

  #---------------------------------------------------------------------------#
  # Get the primary & secondary targets of Nr4a2
  #---------------------------------------------------------------------------#

  # primary target genes 
  cur_primary<- regulons %>% 
    subset(tf == cur_tf) 

  # which of these pimary target genes are tfs?
  cur_primary_tfs <- cur_primary %>% 
    subset(gene %in% unique(regulons$tf)) %>% .$gene

  cur_tfs <- unique(c(cur_tf, cur_primary_tfs))

  # get the regulons for these TFs:
  cur_secondary <- subset(regulons, tf %in% cur_primary_tfs)
  cur_secondary_tfs <- cur_primary %>% 
    subset(gene %in% unique(regulons$tf)) %>% .$gene

  # combine the primary and secondary into one table 
  cur_network <- rbind(cur_primary, cur_secondary)
  cur_network$Gain <- cur_network$Gain * sign(cur_network$Cor)

  cur_genes <- unique(cur_network$gene)
  length(cur_genes)

  # make an igraph network from the nr4a2 regulon:
  cur_network <- cur_network %>%
    dplyr::rename(c(source=tf, target=gene)) %>%
    mutate(Score = sign(Cor) * Gain)

  primary_genes <- unique(cur_primary$gene)
  secondary_genes <- unique(cur_network$target[! cur_network$target %in% primary_genes])

  cur_network$target_type <- ifelse(cur_network$target %in% primary_genes, 'primary', 'secondary')
  cur_network$cell_group <- cur_celltype 
  cur_network$group <- cur_group

  network_df <- rbind(network_df, cur_network)
}


# get a list of co-expression modules that 

modules <- GetModules(seurat_obj)
module_tfs <- subset(modules, gene_name %in% tf_list)

length(unique(network_df$source))

umap_df <- GetModuleUMAP(seurat_obj)

cur_network <- subset(network_df, group == 'Nurr2c' & cell_group == 'MHb-Neuron')
cur_network <- subset(cur_network, source %in% module_tfs$gene_name & target %in% module_tfs$gene_name)
cur_tfs <- unique(c(cur_network$source, cur_network$target))
cur_tfs <- cur_tfs[cur_tfs %in% modules$gene_name]

#---------------------------------------------------------------------------#
# Plot with ggraph
#
# Figure 6G network plot
#---------------------------------------------------------------------------#

library(ggraph)
library(tidygraph)

graph <- as_tbl_graph(cur_network) %>% 
  activate(nodes) %>% 
  mutate(degree  = centrality_degree())  

tmp <- tf_degrees[names(V(graph))]; tmp <- tmp[!is.na(tmp)]
V(graph)[names(tmp)]$degree <- as.numeric(tmp)

V(graph)$gene_type <- ifelse(names(V(graph)) %in% unique(regulons$tf), 'TF', 'Gene')
V(graph)$gene_type <- ifelse(names(V(graph)) == cur_tf, cur_tf, V(graph)$gene_type)

# make the layout table using the umap coords:
umap_layout <- umap_df[names(V(graph)),] %>% dplyr::rename(c(x=UMAP1, y = UMAP2, name=gene))
rownames(umap_layout) <- 1:nrow(umap_layout)
lay <- create_layout(graph, umap_layout)

# add extra info
lay$tf_name <- ifelse(lay$name %in% c(cur_tf, cur_primary_tfs), lay$name, NA)
lay$size <- ifelse(lay$name %in% unique(regulons$tf), 5, 2)
lay$type <- ifelse(lay$name %in% primary_genes, 'Primary', 'Secondary')
lay$type <- ifelse(lay$name == cur_tf, cur_tf, lay$type)
lay$type <- factor(lay$type, levels = c(cur_tf, 'Primary', 'Secondary'))

# shape layout:
cur_shapes <- c(18, 17, 16); names(cur_shapes) <- c(cur_tf, 'TF', 'Gene')
cur_shapes <- c(23, 24, 25); names(cur_shapes) <- c(cur_tf, 'Primary', 'Secondary')

p <- ggraph(lay) + 
  ggrastr::rasterise(
    geom_point(inherit.aes=FALSE, data=umap_df, aes(x=UMAP1, y=UMAP2), color=umap_df$color, alpha=0.1, size=3),
    dpi=500
   ) +
  geom_edge_fan(
    aes(color=Cor, alpha=abs(Cor)),
    arrow = arrow(length = unit(2, 'mm'), type='closed'), 
    end_cap = circle(3, 'mm')
  ) + 
  geom_node_point(data=subset(lay, (! name %in% regulons$tf) | name == cur_tf ), aes(fill=module, shape=type, size=degree), color='black') +
  geom_node_point(data=subset(lay, name %in% regulons$tf & name != cur_tf), aes(fill=module, size=degree, shape=type), color='black') +
  geom_node_label(aes(label=tf_name, color=module), repel=TRUE, max.overlaps=Inf, fontface='italic') +
  scale_edge_colour_gradient2(high='orange2', mid='white', low='dodgerblue')  + 
  scale_colour_manual(values=cp) + 
  scale_fill_manual(values=cp) + 
  scale_shape_manual(values=cur_shapes) 

pdf(paste0(fig_dir,'nr4a2_tf_network_umap.pdf'), width=8, height=7)
print(p)
dev.off()

#---------------------------------------------------------------------------#
# Quantify the number of TF links between modules and plot as a heamap
#
# Figure 6I heatmaps
#---------------------------------------------------------------------------#

tf_regulons <- subset(regulons, tf %in% module_tfs$gene_name & gene %in% module_tfs$gene_name)
ix <- match(tf_regulons$tf, modules$gene_name)
tf_regulons$source_module <- modules$module[ix]
ix <- match(tf_regulons$gene, modules$gene_name)
tf_regulons$target_module <- modules$module[ix]
tf_regulons %<>% dplyr::rename(source=tf, target=gene) %>% mutate(Gain = Gain * sign(Cor))
cur_network <- tf_regulons

# add module info to the network:
ix <- match(cur_network$source, modules$gene_name)
cur_network$source_module <- modules$module[ix]
ix <- match(cur_network$target, modules$gene_name)
cur_network$target_module <- modules$module[ix]

# make empty matrices to store the number of links between mods
pos_mat <- matrix(0, length(mods), length(mods))
rownames(pos_mat) <- mods; colnames(pos_mat) <- mods
neg_mat <- matrix(0, length(mods), length(mods))
rownames(neg_mat) <- mods; colnames(neg_mat) <- mods

# loop through combos of mods and identify the number of +/- links between mods
combos <- expand.grid(mods, mods)
for(i in 1:nrow(combos)){
  m1 <- as.character(combos[i,'Var1'])
  m2 <- as.character(combos[i, 'Var2'])

  # how many positive connections from m1 to m2:
  cur_pos <- subset(cur_network, target_module == m1 & source_module == m2 & Gain >= 0)
  cur_neg <- subset(cur_network, target_module == m1 & source_module == m2 & Gain < 0)

  if(m1 == m2){
    pos_mat[m1,m2] <- 0
    neg_mat[m1,m2] <- 0
  } else{
    pos_mat[m1,m2] <- nrow(cur_pos)
    neg_mat[m1,m2] <- nrow(cur_neg)
  }

}

max_val <- 1
plot_df <- reshape2::melt(neg_mat)
plot_df$count <-plot_df$value
plot_df$label <- ifelse(plot_df$value > 5, plot_df$value, "")
tmp <-table(module_tfs$module)
ix <- match(plot_df$Var1, names(tmp))
plot_df$value <- plot_df$value / as.numeric(tmp)[ix]
plot_df1 <- plot_df

plot_df$value <- ifelse(plot_df$value > max_val, max_val, plot_df$value)
p1 <- plot_df %>% 
  ggplot(aes(x=Var1, y=fct_rev(Var2), fill=value)) + 
  geom_tile() + 
  geom_text(aes(label=label)) +
  scale_fill_gradient(low='grey95', high='dodgerblue') + 
  coord_fixed() + RotatedAxis() + 
  xlab('Target') + ylab('Source') + 
  ggtitle('Repressive interactions') + 
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(linewidth=1, fill=NA, color='black')
  )

plot_df <- reshape2::melt(pos_mat)
plot_df$count <-plot_df$value
plot_df$label <- ifelse(plot_df$value > 5, plot_df$value, "")

tmp <-table(module_tfs$module)
ix <- match(plot_df$Var1, names(tmp))
plot_df$value <- plot_df$value / as.numeric(tmp)[ix]
plot_df2 <- plot_df

plot_df$value <- ifelse(plot_df$value > max_val, max_val, plot_df$value)
p2 <- plot_df %>% 
  ggplot(aes(x=Var1, y=fct_rev(Var2), fill=value)) + 
  geom_tile() + 
  geom_text(aes(label=label)) +
  scale_fill_gradient(low='grey95', high='orange2') + 
  coord_fixed() + RotatedAxis()  + 
  xlab('Target') + ylab('Source') + 
  ggtitle('Activating interactions') + 
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(linewidth=1, fill=NA, color='black')
  )


pdf(paste0(fig_dir,'test_network_interactions_heatmap_norm.pdf'), width=7, height=4)
p1 + p2 + plot_layout(guides='collect')
dev.off()

#---------------------------------------------------------------------------#
# Delta of Activating / Repressive interactions
#---------------------------------------------------------------------------#

plot_df1$delta <-  plot_df2$value - plot_df1$value

p3 <- plot_df1 %>% 
  ggplot(aes(x=Var1, y=fct_rev(Var2), fill=delta)) + 
  geom_tile() + 
  #geom_text(aes(label=label)) +
  scale_fill_gradient2(low='dodgerblue', mid='grey95', high='orange2') + 
  coord_fixed() + RotatedAxis()  + 
  xlab('Target') + ylab('Source') + 
  ggtitle('Activating - repressing') + 
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(linewidth=1, fill=NA, color='black')
  )


pdf(paste0(fig_dir,'test_network_interactions_heatmap_norm_delta.pdf'), width=4, height=10)
p1 / p2 / p3 
dev.off()

#---------------------------------------------------------------------------#
# Dot Plot to show expression of these TFs
#---------------------------------------------------------------------------#

plot_genes <- subset(cur_network, target_type == 'primary') %>% .$target %>% unique
plot_genes <- c('Nr4a2',plot_genes)
plot_genes <- modules %>% arrange(module) %>% subset(gene_name %in% plot_genes)

# make dotplot
p <- DotPlot(
  subset(seurat_obj, Group == 'Nurr2c'),
  group.by='annotation',
  features = rev(plot_genes$gene_name)
#  scale=FALSE, col.max=5
) + coord_flip() + RotatedAxis() +
  scale_color_gradient(high='darkorchid3', low='lightgrey') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text.y = element_text(color=rev(plot_genes$color), face='italic')
  )

pdf(paste0(fig_dir, 'tf_expression_dotplot.pdf'), width=4.25, height=6)
p
dev.off()

#---------------------------------------------------------------------------#
# Bar plot showing the top and bottom TF targets of nr4a2 
#
# Figure 6H barplot
#---------------------------------------------------------------------------#

plot_genes <- subset(cur_network, target_type == 'primary') %>% .$target %>% unique
plot_genes <- c('Nr4a2',plot_genes)
plot_genes <- modules %>% arrange(module) %>% subset(gene_name %in% plot_genes)


plot_df <- subset(cur_network, source == 'Nr4a2' & target %in% plot_genes$gene_name)
plot_df %<>% mutate(rank = dense_rank(dplyr::desc(Gain))) %>% 
  arrange(Gain, desc=TRUE)


p <- plot_df %>%
  ggplot(aes(y=rev(factor(rank)), x = Gain, fill=Gain)) + 
  geom_bar(stat='identity', width=1) +
  geom_vline(xintercept=0, color='black') + 
  geom_text(aes(label=target), color='black', size=3.5, hjust='center') +
  scale_fill_gradient2(high="orange2", mid='grey95', low='dodgerblue') +
  NoLegend() + xlab('') + ylab('') +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

pdf(paste0(fig_dir, 'Nr4a2_target_ranks.pdf'), width=3, height=4)
p
dev.off()





# load the hashikawa dataset 
seurat_hashikawa <- readRDS('/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/analysis/test_harmony/data/hashikawa_seurat.rds')
seurat_hashikawa$cell_type <- ifelse(grepl("MHb", seurat_hashikawa$celltype_neuron), 'MHb', seurat_hashikawa$celltype_neuron)

seurat_hashikawa <- ProjectModules(
  seurat_obj = seurat_hashikawa,
  seurat_ref = seurat_obj,
  wgcna_name_proj="MHb_projected",
  wgcna_name = "MHb_metacell"
)

# plot modules in hashikawa data:
plot_list <- ModuleFeaturePlot(
  seurat_hashikawa, order="shuffle", raster=TRUE,
  restrict_range=FALSE, point_size=0.5, raster_scale=0.5
)

# remove legend, title, then plot
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + NoLegend() + ggtitle('') + theme(plot.margin = margin(0,0,0,0))
}

pdf("figures/MHb_featureplot_hMEs_hashikawa.pdf",height=6, width=12)
wrap_plots(plot_list, ncol=5)
dev.off()

################################################################################
# Module preservation
################################################################################

seurat_hashikawa <- SetDatExpr(
  seurat_hashikawa,
  group_name = "MHb",
  group.by = "cell_type",
  use_metacells = FALSE
)

# run module preservation function
seurat_hashikawa <- ModulePreservation(
  seurat_hashikawa,
  seurat_ref = seurat_obj,
  name="MHb_projected",
  verbose=3,
  n_permutations=200
)

plot_list <- PlotModulePreservation(
  seurat_hashikawa,
  name="MHb_projected",
  statistics = "summary"
)

pdf(paste0(fig_dir, 'hashikawa_module_preservation_summary.pdf'), width=10, height=5)
wrap_plots(plot_list, ncol=2)
dev.off()

saveRDS(seurat_hashikawa, file='data/hashikawa_hdWGCNA.rds')
seurat_hashikawa <- readRDS('data/hashikawa_hdWGCNA.rds')

################################################################################
# Preservation in the naive dataset:
################################################################################

seurat_full <- readRDS(file='/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/Nurr2c_vs_GFP/revision/data/harmony_annotated_integration.rds')
seurat_naive <- subset(seurat_full, cell_type == 'MHb-Neuron' & Group %in% c('NN', 'NGFP'))

seurat_naive <- ProjectModules(
  seurat_obj = seurat_naive,
  seurat_ref = seurat_obj,
  wgcna_name_proj="MHb_projected",
  wgcna_name = "MHb_metacell"
)

################################################################################
# Module preservation
################################################################################

seurat_naive <- SetDatExpr(
  seurat_naive,
  group_name = "MHb-Neuron",
  group.by = "cell_type",
  use_metacells = FALSE
)

# run module preservation function
seurat_naive <- ModulePreservation(
  seurat_naive,
  seurat_ref = seurat_obj,
  name="MHb_projected",
  verbose=3,
  n_permutations=200
)

plot_list <- PlotModulePreservation(
  seurat_naive,
  name="MHb_projected",
  statistics = "summary"
)

pdf(paste0(fig_dir, 'naive_module_preservation_summary.pdf'), width=10, height=5)
wrap_plots(plot_list, ncol=2)
dev.off()


saveRDS(seurat_naive, file='data/naive_hdWGCNA.rds')
seurat_naive <- readRDS(file='data/naive_hdWGCNA.rds')

