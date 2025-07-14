# Part 3: Regulatory network analysis and visualization
# filtered_regulatory_relationships <- filtered_regulatory # added by yd

pdf("network.pdf")
plot_tf_network(TFs_list, layout = "grid")
plot_tf_network(TFs_list, layout = "sphere") 
plot_tf_network(TFs_list, layout = "circle")
plot_tf_network(TFs_list, layout = "random")
dev.off()



enrich_module <- function(Kmeans_result, org.db, enrich.db ,fun_num = 5,
                          pvalueCutoff = 0.05, use_internal_data = TRUE, organism = NULL) {
  all_gene <- Kmeans_result
  all_gene<-all_gene[order(all_gene$KmeansGroup),]
  le<-levels(as.factor(all_gene$KmeansGroup))
  for (i in le) {
    acc11<-rownames(all_gene[all_gene$KmeansGroup == i,])
    # gene1 <- clusterProfiler::bitr(acc11, fromType = "ENSEMBL",
    #               toType = c("SYMBOL", "ENTREZID"),
    #               OrgDb = org.db)
    gene1 <- data.frame(ENTREZID = acc11, row.names = acc11)
    if (enrich.db =='KEGG') {
      k1 <- clusterProfiler::enrichKEGG(gene = gene1$ENTREZID,
                                        pvalueCutoff = pvalueCutoff
                                        ,organism = organism
                                        ,use_internal_data = use_internal_data)
    }else if(enrich.db =='GO'){
      k1 = clusterProfiler::enrichGO(gene = gene1$ENTREZID,
                    OrgDb = org.db,
                    keyType = "GID",
                    ont = "BP",
                    pvalueCutoff = pvalueCutoff)
    }
    acc2 <- k1@result
    acc2$'-log10(q-value)' <- -log10(acc2$qvalue)
    acc2 <- acc2[order(-acc2$`-log10(q-value)`),]
    if (i=='1' | i == 1) {
      acc21 <- acc2[1:fun_num,]
      acc21$module <- rep(i,fun_num)
    }else{
      acc22 <- acc2[1:fun_num,]
      acc22$module<-rep(i,fun_num)
      acc21 <- rbind(acc21,acc22)
    }
  }
  acc21 <- acc21[,c(1,2,11,10,3:9)]
  return(acc21)
}
enrichment_GO <- enrich_module(Kmeans_clustering_ENS, org.db='org.Ga.eg.db', enrich.db = 'GO', fun_num = 5, pvalueCutoff = 0.05, use_internal_data = TRUE, organism = 'Ga')

### select functions that you want to present in the figure
enrichment_GO_select <- enrichment_GO[c(1,11,21),]
### plotting
plot_intramodular_network(TFs_list,enrichment_GO_select,layout = 'circle')
dev.off()