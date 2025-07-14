# bamfilepath1 <- 'E1_N2_filter.bam'
# bamfilepath2 <- 'E1_N1_filter.bam'
# bamfilepath3 <- 'E1_0_filter.bam'
# bamfilepath4 <- 'E1_P1_filter.bam'
# ### set parameter 'workers' to make this function run in parallel
# cuts1 <- cal_footprint_cuts(bamfilepath = bamfilepath1,bedfile = list2[[1]],workers = 32,index_bam = T)
# cuts2 <- cal_footprint_cuts(bamfilepath = bamfilepath2,bedfile = list2[[1]],workers = 32,index_bam = T)
# cuts3 <- cal_footprint_cuts(bamfilepath = bamfilepath3,bedfile = list2[[1]],workers = 32,index_bam = T)
# cuts4 <- cal_footprint_cuts(bamfilepath = bamfilepath4,bedfile = list2[[1]],workers = 32,index_bam = T)
# cut_list <- list(cuts1,cuts2,cuts3,cuts4)

### get related genes of footprints with high FOS
potential_regulation <- Footprints_FOS(cut_list,list2[[2]], FOS_threshold = 0.1); length(potential_regulation$TF)
### Use information of footprints with high FOS to refine regulatory relationships
# filtered_regulatory <- filter_ATAC(potential_regulation,regulatory_relationships)


filter_ATAC <- function (FOS, regulary_relationships) 
{
    # validInput(regulary_relationships, "regulary_relationships", "df")
    # if (grepl("ENS", FOS[1, 1])) {
    #     TfIndex <- 1
    # } else {
    #     TfIndex <- 2
    # }
    # if (grepl("ENS", FOS[1, 2])) {
    #     TargetIndex <- 4
    # } else {
    #     TargetIndex <- 5
    # }
    TfIndex <- 1; TargetIndex <- 4 # added by yd
    pair1 <- paste(FOS[, 1], FOS[, 2])
    pair2 <- paste(regulary_relationships[, TfIndex], regulary_relationships[, 
        TargetIndex])
    filtered <- regulary_relationships[pair2 %in% pair1, ]
    return(filtered)
}
head(potential_regulation, n=2)
head(regulatory_relationships, n=2)
filtered_regulatory <- filter_ATAC(potential_regulation,regulatory_relationships); head(filtered_regulatory)
length(filtered_regulatory$TF)

TFs_list <- network_analysis(filtered_regulatory,Kmeans_clustering_ENS,TFFDR1 = 2,TFFDR2 = 2)
# TFs_list <- network_analysis(filtered_regulatory,Kmeans_clustering_ENS,TFFDR1 = 1,TFFDR2 = 1)
str(TFs_list)

save(TFs_list, filtered_regulatory, potential_regulation, Kmeans_clustering_ENS, motif1, file = paste0(basename(genie3_rdata),"_regulatory.RData"))
message(paste0("Saved results(TFs_list, filtered_regulatory, potential_regulation, Kmeans_clustering_ENS, motif1) in: ", basename(genie3_rdata), "_regulatory.RData"))