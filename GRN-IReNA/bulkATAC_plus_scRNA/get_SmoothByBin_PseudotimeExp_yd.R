# Ref:
get_SmoothByBin_PseudotimeExp_yd <- function(seurat_object, FC = TRUE, Bin_key = "seurat_clusters",
                                          method = 'ydgenomics',FcType = "Q95"){
    hvg_object <- seurat_object
    # hvg_object[["RNA"]]<-as(object=hvg_object[["RNA"]],Class="Assay")
    n_cell <- length(unique(hvg_object@meta.data[[Bin_key]]))
    list_cell <- unique(hvg_object@meta.data[[Bin_key]])

    Exp1 <- as.matrix(hvg_object@assays$RNA@data)
    Exp2 <- array(0, dim = c(nrow(Exp1), n_cell))

    for (i in 1:n_cell) {
      Cells1 <- subset(hvg_object, seurat_clusters==list_cell[i])
      # print(colnames(Cells1))
      if (nrow(Cells1) > 1) {
        Exp2[, i] <- rowMeans(Exp1[, match(colnames(Cells1), colnames(Exp1))])
      } else if (nrow(Cells1) == 1) {
        Exp2[, i] <- Exp1[, match(colnames(Cells1), colnames(Exp1))]
      }
    }
    rownames(Exp2) <- rownames(Exp1)
    Exp1 <- Exp2
    FcType1 <- FcType
    if (FC == TRUE) {
      if (FcType1 == "Q90") {
        Prob1 <- c(0.1, 0.9)
      } else {
        Prob1 <- c(0.05, 0.95)
      }
      Fc1 <- apply(Exp1, 1, function(x1) {
        x12 <- quantile(x1, probs = Prob1)
        x2 <- x12[2] - x12[1]
        return(x2)
      })
      Exp2 <- cbind(Fc1, Exp1)
      # colnames(Exp2)[1:ncol(Exp2)] <- c(paste0("FoldChange", FcType1),
      #                                 paste0("SmExp_", 1:ncol(Exp1)))
      colnames(Exp2)[1:ncol(Exp2)] <- c(paste0("FoldChange", FcType1),
                                      paste0("SmExp_", list_cell))
      var1 <- as.data.frame(Exp2)
      return(var1)
    } else {
      colnames(Exp1)[1:ncol(Exp1)] <- c(paste0("SmExp", 1:ncol(Exp1)))
      var1 <- as.data.frame(Exp1)
      return(var1)
    }
}