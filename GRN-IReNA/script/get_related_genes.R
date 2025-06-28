# get_related_genes
function (footprints, motif, Species, txdb, tssRegion = c(-3000, 
    3000)) 
{
    validInput(footprints, "footprints", "df")
    validInput(motif, "motif", "df")
    validInput(Species, "Species", "character")
    validInput(txdb, "txdb", "txdb")
    validInput(tssRegion, "tssRegion", "vector")
    footprintslist <- merge_extent_footprints(footprints, motif)
    merged_footprints <- footprintslist[[2]]
    if (Species == "Hs") {
        annodb <- "org.Hs.eg.db"
    }
    else if (Species == "Mm") {
        annodb <- "org.Mm.eg.db"
    }
    else if (Species == "Zf") {
        annodb <- "org.Dr.eg.db"
    }
    else if (Species == "Ch") {
        annodb <- "org.Gg.eg.db"
    }
    reference_GRange <- GenomicRanges::GRanges(seqnames = merged_footprints[, 
        1], IRanges::IRanges(start = as.numeric(merged_footprints[, 
        2]), end = as.numeric(merged_footprints[, 3])), strand = merged_footprints[, 
        4])
    peakAnno <- ChIPseeker::annotatePeak(reference_GRange, tssRegion = tssRegion, 
        TxDb = txdb, annoDb = annodb)
    region <- peakAnno@anno@elementMetadata$annotation
    gene <- peakAnno@anno@elementMetadata$ENSEMBL
    start1 <- peakAnno@anno@ranges@start
    merged_footprints2 <- merged_footprints[merged_footprints$V2 %in% 
        start1, ]
    exon1 <- grep("exon", region)
    Intron1 <- grep("Intron", region)
    Intergenic1 <- grep("Intergenic", region)
    Downstream1 <- grep("Downstream", region)
    Promoter1 <- grep("Promoter", region)
    UTR3 <- grep("3' UTR", region)
    UTR5 <- grep("5' UTR", region)
    region2 <- rep(NA, length(region))
    region2[exon1] = "Exon"
    region2[Intron1] = "Intron"
    region2[Downstream1] = "Downstream"
    region2[Promoter1] = "Promoter"
    region2[UTR3] = "3' UTR"
    region2[UTR5] = "5' UTR"
    region2[Intergenic1] = "Intergenic"
    table(region2)
    peak_region1 <- paste(as.character(peakAnno@anno@seqnames), 
        as.character(peakAnno@anno@ranges), sep = ":")
    peak_region2 <- paste0(merged_footprints[, 1], ":", merged_footprints[, 
        2], "-", merged_footprints[, 3])
    merged_footprints2 <- merged_footprints[peak_region2 %in% 
        peak_region1, ]
    merged_footprints2$gene <- gene
    merged_footprints2$region <- region2
    merged_footprints2 <- merged_footprints2[, c(9, 8, 1:7)]
    colnames(merged_footprints2) <- c(paste0("V", 1:9))
    footprintslist[[2]] <- merged_footprints2
    footprintslist[[1]] <- footprintslist[[1]][peak_region2 %in% 
        peak_region1, ]
    return(footprintslist)
}