# Date: 250710
# Ref: https://rdrr.io/github/jiang-junyao/IReNA/src/R/smooth.R

smoothByBin <- function(Exp1, Pseudotime1, PseudotimeRange1 = NULL,
                        SmoothLength1 = 100, ByBin1 = c("Equal.Pseudotime", "Equal.Cells")) {
  if (ncol(Exp1) != nrow(Pseudotime1)) {
    stop("The length of pseudotime is not equal to the number of cells")
  } else if (!is.element("Pseudotime", colnames(Pseudotime1))) {
    stop("No pseudotime inforamtion in variable Pseudotime1")
  }

  if (ByBin1[1] == "Equal.Pseudotime") {
    if (is.null(PseudotimeRange1)) {
      PseudotimeBin1 <- seq(min(Pseudotime1$Pseudotime), max(Pseudotime1$Pseudotime),
                            length.out = SmoothLength1 + 1)
    } else {
      PseudotimeBin1 <- seq(PseudotimeRange1[1], PseudotimeRange1[2],
                            length.out = SmoothLength1 + 1)
    }
  } else {
    Pseudotime1 <- Pseudotime1[order(Pseudotime1$Pseudotime), ]
    Bin1 <- ceiling(nrow(Pseudotime1) / SmoothLength1)
    PseudotimeBin1 <- seq(1, nrow(Pseudotime1), by = Bin1)
    PseudotimeBin1[length(PseudotimeBin1) + 1] <- nrow(Pseudotime1)
  }

  Exp2 <- array(0, dim = c(nrow(Exp1), SmoothLength1))
  for (i in 1:(length(PseudotimeBin1) - 1)) {
    if (ByBin1[1] == "Equal.Pseudotime") {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] &
                                Pseudotime1$Pseudotime <= PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] &
                                Pseudotime1$Pseudotime < PseudotimeBin1[i + 1], ]
      }
    } else {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:(PseudotimeBin1[i + 1] - 1), ]
      }
    }

    if (nrow(Cells1) > 1) {
      Exp2[, i] <- rowMeans(Exp1[, match(rownames(Cells1), colnames(Exp1))])
    } else if (nrow(Cells1) == 1) {
      Exp2[, i] <- Exp1[, match(rownames(Cells1), colnames(Exp1))]
    }
  }
  rownames(Exp2) <- rownames(Exp1)

  return(Exp2)
}




smoothByState <- function(seurat_obj, Pseudotime1, each_state_bin1){
  ### order cells
  Pseudotime1 <- Pseudotime1[order(Pseudotime1$Pseudotime),]
  Pseudotime1$State <- as.character(Pseudotime1$State)
  seurat_exp <- as.matrix(seurat_obj@assays$RNA@data)
  seurat_exp <- seurat_exp[,rownames(Pseudotime1)]
  ### order states
  order_state <- as.character(Pseudotime1$State[!duplicated(Pseudotime1$State)])
  ### divide cells in each state to bins
  each_state_bin <- each_state_bin1
  for (j in order_state) {

    Pseudotime2 <- Pseudotime1[Pseudotime1$State==j,]
    seurat_exp2_state <- seurat_exp[,rownames(Pseudotime2)]
    PseudotimeBin1 <- round(seq(1,nrow(Pseudotime2),length.out=each_state_bin+1))
    Exp2 <- array(0, dim = c(nrow(seurat_exp), each_state_bin))
    for (i in 1:(length(PseudotimeBin1) - 1)) {
        Exp2[, i] <- rowMeans(seurat_exp2_state[, PseudotimeBin1[i]:PseudotimeBin1[i+1]])
    }
    if (j==1) {
      Exp3 = Exp2
    }else{Exp3 <-cbind(Exp3,Exp2)}
  }
  rownames(Exp3) <- rownames(seurat_exp)
  return(Exp3)
}



#' Smooth cells into bins based on pseudotime or State
#' @description Divide cells into bins across pseudotime and return expression
#' profile
#' @param seurat_object seurat object, where meta.data should contain
#' Pseudotime
#' @param FC logic, indicating whether to add FoldChangeQ95 to the first column,
#' default is TRUE
#' @param Bin numeric, indicating the numbers of bins which divide the pseudotime,
#' default is 50
#' @param FcType 'Q95' or 'Q90', FoldChange threshold, default is 'Q95'
#' @param method 'Pseudotime' or 'State'. Pseudotime method firstly orders all
#' cells according to pseudotime, and then smooth these cells into ‘Bin' number
#' of bins. State method firstly orders all cells according to pseudotime, and
#' then smooth cells in each state into ('Bin'/number of state) number of bins.
#' @importFrom stats quantile
#' @return return a expression profile
#' @export
#'
#' @examples load(system.file("extdata", "test_seurat.rda", package = "IReNA"))
#' monocle_object = get_pseudotime(test_seurat)
#' seurat_object_pseudotime=add_pseudotime(seurat_object = test_seurat,monocle_object = monocle_object)
#' get_SmoothByBin_PseudotimeExp(seurat_object_pseudotime, Bin = 50, FcType = "Q95")
get_SmoothByBin_PseudotimeExp <- function(seurat_object, FC = TRUE, Bin = 50,
                                          method = 'Pseudotime',FcType = "Q95") {
  validInput(seurat_object,'seurat_object','seurat')
  validInput(FC,'FC','logical')
  validInput(Bin,'Bin','numeric')
  if (!method %in% c('Pseudotime','State')) {
    stop('method should be Pseudotime or State')
  }
  FcType1 <- FcType
  ByBin1 <- c("Equal.Pseudotime")
  pbmc1 <- seurat_object
  if (method == 'Pseudotime') {
    TotalBin1 <- Bin
    Exp1 <- smoothByBin(as.matrix(pbmc1@assays$RNA@data), pbmc1@meta.data,
                      SmoothLength1 = TotalBin1, ByBin1 = ByBin1[1])
  }
  if (method == 'State') {
    TotalBin1 <- round(Bin/length(levels(as.factor(pbmc1@meta.data$State))))
    Exp1 <- smoothByState(pbmc1, pbmc1@meta.data,
                          TotalBin1)
  }
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
    colnames(Exp2)[1:ncol(Exp2)] <- c(paste0("FoldChange", FcType1),
                                      paste0("SmExp", 1:ncol(Exp1)))
    var1 <- as.data.frame(Exp2)
    return(var1)
  } else {
    colnames(Exp1)[1:ncol(Exp1)] <- c(paste0("SmExp", 1:ncol(Exp1)))
    var1 <- as.data.frame(Exp1)
    return(var1)
  }
}

#' filter expression profile
#' @description remove empty columns and genes with low foldchange in expression profile
#' @param expression_profile expression profile, generated by
#' \code{\link{get_SmoothByBin_PseudotimeExp}}
#' @param filterfc logic, indicating whether you want to filter genes
#' based on logFC
#' @param FC numeric, indicating logFC threshold to filter genes
#'
#' @return return filtered expression profile
#' @export
#'
#' @examples load(system.file("extdata", "test_clustering.rda", package = "IReNA"))
#' expression_profile = test_clustering[,-1]
#' filter_expression_profile(expression_profile)
filter_expression_profile <- function(expression_profile, filterfc = TRUE,
                                       FC = 0.1) {
  validInput(expression_profile,'expression_profile','df')
  validInput(filterfc,'filterfc','logical')
  validInput(FC,'FC','numeric')
  pro <- expression_profile
  acc1 <- colSums(pro)
  acc2 <- c()
  for (i in 1:length(acc1)) {
    if (acc1[i] != 0) {
      acc2 <- c(acc2, i)
    }
  }
  var1 <- pro[, acc2]
  if (filterfc == TRUE) {
    var1 <- var1[var1[, 1] > FC, ]
    return(var1)
  }
  else {
    return(var1)
  }
}

#' clustering DEGs by K-means
#' @description  This function is used to cluster the DEGs in the expression
#' profile by K-means algorithmn
#' @param RNA1 expression profile
#' @param K1 numbers of bin
#' @param ColumnGroup1 vector, if it is not NULL, this function wil sepate columns
#' according to varible ColumnGroup1, default is NULL
#' @param Scale1 character, indicating if the values should be centered and scaled
#' in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param Range1 numeric, indicating the range of expression values, defalut is c(-Inf, Inf)
#' @param Reorder1 logic, indicating whether to order gene through hclustering,
#' defalut is TRUE
#' @param RevOrder1 logic, indicating whether to reverser the order
#' @param NAcolum1 vector, indicating which columns need to be removed due to NA values
#' @importFrom stats kmeans
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @return return a data.frame, and the first column is K-means group
#' @export
#'
#' @examples load(system.file("extdata", "test_clustering.rda", package = "IReNA"))
#' expression_profile = test_clustering[,-1];expression_profile<-na.omit(expression_profile)
#' clustering_Kmeans(expression_profile, K1=4, Range1=c(-3,3))
clustering_Kmeans <- function(RNA1, K1 = 1, ColumnGroup1 = NULL, Scale1 = "row",
                              Range1 = c(-Inf, Inf), Reorder1 = TRUE,
                              RevOrder1 = -1, NAcolum1 = NULL) {
  validInput(K1,'K1','numeric')
  validInput(Range1,'Range1','numeric')
  RowGroup1 = NULL;NumColumnBlank1 = NULL
  if (is.numeric(K1)) {
    K2 <- K1
  } else{K2 <- nrow(K1)}
  if (Scale1 == "row") {
    RNA1 <- t(scale(t(RNA1)))
  }
  RNA1[is.na(RNA1)] = 0
  RNA1[is.nan(RNA1)] = 0
  if (!is.null(ColumnGroup1)) {
    print("sepate columns according to varible ColumnGroup1")
    ColumnBlank1 <- array(0, dim = c(nrow(RNA1), NumColumnBlank1))
    uColumnGroup1 <- unique(ColumnGroup1)
    if (length(uColumnGroup1) == 1) {
    } else {
      for (i in 1:length(uColumnGroup1)) {
        RNA12 <- RNA1[, ColumnGroup1 == uColumnGroup1[i]]
        if (i == 1) {
          RNA21 <- cbind(RNA12, ColumnBlank1)
        } else if (i < length(uColumnGroup1)) {
          RNA21 <- cbind(RNA21, RNA12, ColumnBlank1)
        } else {
          RNA21 <- cbind(RNA21, RNA12)
        }
      }
    }
  } else {
    RNA21 <- RNA1
  }

  print("Perform k-means")
  if (is.null(RowGroup1)) {
    if (!is.null(NAcolum1)) {
      cRNA1 <- kmeans(RNA1[, -NAcolum1], K1)
      RNA02 <- cbind(cRNA1$cluster, RNA1[, -NAcolum1])
      Cluster1 <- cRNA1$cluster
    } else {
      if (is.numeric(K1)) {
        if (K1 == 1) {
          Cluster1 <- rep(1, nrow(RNA1))
        } else {
          cRNA1 <- kmeans(RNA1, K1)
          Cluster1 <- cRNA1$cluster
        }
      }else{
        cRNA1 <- kmeans(RNA1, K1)
        Cluster1 <- cRNA1$cluster
      }
      RNA02 <- cbind(Cluster1, RNA1)
    }
    RNA2 <- cbind(Cluster1, RNA21)
    NameInd1 <- 1:K2
    colnames(RNA2)[1] <- c("KmeansGroup")
    print(table(RNA2[, "KmeansGroup"]))
  } else {
    if (is.factor(RowGroup1)) {
      Name1 <- levels(RowGroup1)
      Name2 <- sort(levels(RowGroup1))
      NameInd1 <- match(Name1, Name2)
      RowGroup1 <- as.numeric(RowGroup1)
      uRowGroup1 <- NameInd1
    } else {
      uRowGroup1 <- sort(unique(RowGroup1))
    }

    if (!is.null(NAcolum1)) {
      RNA02 <- cbind(RowGroup1, RNA1[, -NAcolum1])
    } else {
      RNA02 <- cbind(RowGroup1, RNA1)
    }
    RNA2 <- cbind(RowGroup1, RNA21)
    K1 <- length(uRowGroup1)
    colnames(RNA2)[1] <- c("KmeansGroup")
    tRNA2 <- table(RNA2[, "KmeansGroup"])
    if (is.factor(RowGroup1)) {
      names(tRNA2) <- Name2
      tRNA2 <- tRNA2[NameInd1]
    }
    print(tRNA2)
  }

  print("Sort genes")
  RNA22 <- c()
  for (i in 1:K2) {
    if (is.null(RowGroup1)) {
      RNA20 <- RNA2[RNA2[, "KmeansGroup"] == i, ]
      RNA03 <- RNA02[RNA02[, 1] == i, ]
    } else {
      RNA20 <- RNA2[RNA2[, "KmeansGroup"] == uRowGroup1[i], ]
      RNA03 <- RNA02[RNA02[, 1] == uRowGroup1[i], ]
    }
    if (Reorder1 == TRUE) {
      if (nrow(RNA03) > 1) {
        Hier1 <- hclust(as.dist((1 - cor(t(RNA03[, 2:ncol(RNA03)]))) / 2))
        if (RevOrder1[1] != -1) {
          RevOrder2 <- FALSE
          for (j in 1:length(RevOrder1)) {
            if (RevOrder1[j] == i) {
              RevOrder2 <- TRUE
              break
            }
          }
          if (RevOrder2 == TRUE) {
            Ind1 <- rev(Hier1$order)
          } else {
            Ind1 <- Hier1$order
          }
        } else {
          Ind1 <- Hier1$order
        }
      } else {
        Ind1 <- 1:nrow(RNA20)
      }
    } else {
      Ind1 <- 1:nrow(RNA20)
    }
    RNA22 <- rbind(RNA22, RNA20[Ind1, ])

  }

  print("Revise outlier")
  RNA23 <- RNA22[, 2:ncol(RNA22)]
  print(paste("Number of outlier:", c(length(RNA23[RNA23 < Range1[1]]),
                                      length(RNA23[RNA23 > Range1[2]]))))
  RNA23[RNA23 < Range1[1]] <- Range1[1]
  RNA23[RNA23 > Range1[2]] <- Range1[2]
  RNA3 <- cbind(RNA22[, 1], RNA23)
  colnames(RNA3)[1] <- "Module"
  RNA4 <- RNA22[RNA22[, 1] != 0, ]
  RNA5 <- cbind(RNA4[match(rownames(RNA1), rownames(RNA4)), "KmeansGroup"], RNA1)
  colnames(RNA5)[1] <- "KmeansGroup"
  RNA6 <- RNA5[order(RNA5[, "KmeansGroup"]), ]
  RNA6<-as.data.frame(RNA6)
  return(RNA6)
}

#' add ENSID
#' @description  Add ENSEMBLE ID to the first column of your data.frame
#' @param Kmeans_result Kmeans result data.frame, row names should be Symbol ID
#' @param GeneInf1 data.frame, correspondence file of gene ID, row names should
#' be ENSEMBLE ID, first column should be Symbol ID, if GeneInf1 = NULL this
#' function will be built-in corresponding gene ID file.
#' @param Spec1 If you don’t have a gene ID corresponding file, you can also use
#' our built-in corresponding gene ID file, 'Mm' for mus musculus
#' @return add ENSEMBL ID to the first column of data frame
#' @export
#'
#' @examples
#' load(system.file("extdata", "test_clustering.rda", package = "IReNA"))
#' add_ENSID(test_clustering, Spec1 = "Hs")
add_ENSID <- function(Kmeans_result, GeneInf1 = NULL, Spec1 = "") {
  validInput(Kmeans_result,'Kmeans_result','df')
  if (is.null(GeneInf1)) {
    if (Spec1 == "") {
      stop('Please input a gene names corresponding table or species name')
    }
  }
  con2 <- Kmeans_result
  con1 <- Converse_GeneIDSymbol(rownames(Kmeans_result), GeneInf1, Spec1 = Spec1)
  con2 <- con2[!is.na(match(rownames(con2), con1[, 2])), ]
  con2$ENS <- con1[, 1]
  con2$Symbol <- rownames(con2)
  rownames(con2) <- con2$ENS
  con2 <- con2[, c(ncol(con2), 1:(ncol(con2) - 2))]
  return(con2)
}

#' plot kmeans pheatmap
#' @description  Display the result of clustering based on pheatmap, the input
#' data.frame should contain KmeansGroup column.
#' @param Kmeans_result expression profile which should contain column "KmeansGroup"
#' @param start_column numeric, indicating the start column of expression value,
#' defalut is 3
#' @param Gene1 custom labels for rows that are used instead of rownames
#' @param NumRowBlank1 the blank between each group
#' @param show_colnames logic, indicating whether show column names, default is FALSE
#' @param clustering_distance_rows distance measure used in clustering rows.
#' Possible values are "correlation" for Pearson correlation and all the
#' distances supported by dist, such as "euclidean", etc. If the value
#' is none of the above it is assumed that a distance matrix is provided
#' @param ModuleColor1 each color for each module
#' @param ModuleScale1 change the relative proportion of module bar and expression heatmap
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object.
#' @param fontsize font size
#' @param legend1 logic, indicating whether show the legend, default is TRUE
#' @param border_color color of cell borders on heatmap, use NA if no border
#' should be drawn.
#' @param ByPanel logic, if TRUE, reorder columns of heatmap according to the
#' parameter ColumnGroup1, default is FALSE
#' @param Color1 colors used in heat map
#' @param Show.Module logic, indicating whether show the module, defalut is TURE
#' @param ColumnGroup1 vector, if it is not NULL, this function wil sepate columns
#' according to varible ColumnGroup1, default is NULL
#' @param Range1 numeric, indicating the range of expression values which will
#' affect the range of heat map, defalut is c(-3, 3)
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange
#' @export
#' @return Kmeans clustering figure
#' @examples
#' col1 <- c('#67C1E3','#EF9951','#00BFC4','#AEC7E8','#C067A9','#E56145','#2F4F4F')
#' load(system.file("extdata", "test_clustering.rda", package = "IReNA"))
#' plot_kmeans_pheatmap(test_clustering, ModuleColor1 = col1,Range1=c(-3,3),NumRowBlank1=1,ModuleScale1 = 20)
#'
plot_kmeans_pheatmap <- function(Kmeans_result, start_column = 3, Gene1 = NULL,
                                 NumRowBlank1 = 20, show_colnames = FALSE,
                                 clustering_distance_rows = "correlation",
                                 ModuleColor1 = NULL, ModuleScale1 = 20,
                                 cluster_cols = FALSE, fontsize = 10,
                                 legend1 = TRUE, border_color = "white",
                                 ByPanel = FALSE, Color1 = NULL, Show.Module = TRUE,
                                 Range1 = c(-3, 3), ColumnGroup1 = NULL) {
  validInput(Kmeans_result,'Kmeans_result','df')
  validInput(start_column,'start_column','numeric')
  validInput(NumRowBlank1,'NumRowBlank1','numeric')
  validInput(ModuleScale1,'ModuleScale1','numeric')
  validInput(fontsize,'fontsize','numeric')
  validInput(show_colnames,'show_colnames','logical')
  validInput(ByPanel,'ByPanel','logical')
  validInput(Show.Module,'Show.Module','logical')
  validInput(legend1,'legend1','logical')

  RNA1 <- Kmeans_result
  g1<-unique(RNA1$KmeansGroup)
  for (i in g1) {
    if (i==g1[1]) {
      c1<-RNA1[RNA1$KmeansGroup==i,]
      RowBlank1 <- array(0, dim=c(NumRowBlank1, ncol(RNA1)))
      colnames(RowBlank1)<-colnames(RNA1)
      c2<-rbind(c1,RowBlank1)
    }else {
      c1<-RNA1[RNA1$KmeansGroup==i,]
      RowBlank1 <- array(0, dim=c(NumRowBlank1, ncol(RNA1)))
      colnames(RowBlank1)<-colnames(RNA1)
      c2<-rbind(c2,c1,RowBlank1)
    }
  }
  RNA1<-c2
  RowGroup1 <- as.factor(RNA1$KmeansGroup)
  K1 <- length(levels(as.factor(RowGroup1)))
  NameInd1 <- 1:K1
  RNA3 <- RNA1[start_column:ncol(RNA1)]
  RNA3 <- t(scale(t(RNA3)))
  RNA3[RNA3 < Range1[1]] <- Range1[1]
  RNA3[RNA3 > Range1[2]] <- Range1[2]
  if (is.null(Gene1)) {
    show_rownames <- FALSE
    Gene2 <- ""
  } else if (Gene1[1] == "all") {
    show_rownames <- TRUE
    Gene2 <- rownames(RNA3)
  } else {
    show_rownames <- TRUE
    Gene2 <- rownames(RNA3)
    Gene2[!is.element(Gene2, Gene1)] <- ""
    Gene3 <- rownames(RNA3)[is.element(rownames(RNA3), Gene1)]
    print(sort(Gene3))
  }
  RNA3[is.nan(RNA3)] = 0
  if (is.null(Color1)) {
    Color1 <- c(rgb(72 / 255, 85 / 255, 167 / 255), rgb(255 / 255, 255 / 255, 255 / 255),
                rgb(239 / 255, 58 / 255, 37 / 255))
  }

  if (Show.Module == TRUE) {
    if (is.null(ModuleColor1)) {
      ModuleColor1 <- c(rgb(255 / 255, 255 / 255, 255 / 255), rgb(152 / 255, 152 / 255, 152 / 255),
                        rgb(72 / 255, 85 / 255, 167 / 255))
    } else {
      if (!is.null(RowGroup1) & is.factor(RowGroup1)) {
        ModuleColor1<-c('white',ModuleColor1)
      }
    }

    if (K1 == 1) {
      ph1 <- pheatmap(RNA1$KmeansGroup, clustering_distance_rows = clustering_distance_rows,
                      show_rownames = FALSE, show_colnames = show_colnames, border_color = "white",
                      cluster_rows = FALSE, cluster_cols = cluster_cols, color = ModuleColor1,
                      legend = FALSE, silent = TRUE, breaks = seq(1, nrow(RNA3)))
    } else {
      ph1 <- pheatmap(RNA1$KmeansGroup, clustering_distance_rows = clustering_distance_rows,
                      show_rownames = FALSE, show_colnames = show_colnames, border_color = "white",
                      cluster_rows = FALSE, cluster_cols = cluster_cols, color = ModuleColor1,
                      legend = FALSE, silent = TRUE)
    }
    if (ByPanel == TRUE) {
      NumColumnBlank1 <- 10
      uColumnGroup1 <- unique(ColumnGroup1)
      RNA31 <- RNA3[, 2:(length(ColumnGroup1[ColumnGroup1 == uColumnGroup1[1]]) + 1)]
      ph12 <- pheatmap(RNA3[, 2:(length(ColumnGroup1[ColumnGroup1 == uColumnGroup1[1]]) + 1)],
                       clustering_distance_rows = clustering_distance_rows,
                       clustering_method = clustering_method, show_rownames = show_rownames,
                       show_colnames = show_colnames, labels_row = Gene2, border_color = border_color,
                       cluster_rows = FALSE, cluster_cols = cluster_cols,
                       color = colorRampPalette(Color1)(100),
                       fontsize = fontsize, legend = legend1, silent = TRUE)
      ph2 <- pheatmap(RNA31, clustering_distance_rows = clustering_distance_rows,
                      clustering_method = clustering_method, show_rownames = show_rownames,
                      show_colnames = show_colnames, labels_row = Gene2, border_color = border_color,
                      cluster_rows = FALSE, cluster_cols = cluster_cols,
                      color = colorRampPalette(Color1)(100),
                      fontsize = fontsize, legend = legend1, silent = TRUE)
      plot_list <- list()
      plot_list[[1]] <- ph1[[4]]
      plot_list[[2]] <- ph2[[4]]
      ModuleScale2 <- rep(2, ceiling(ncol(RNA31) / ncol(RNA3) * ModuleScale1))
      if (length(uColumnGroup1) > 1) {
        for (i in 2:length(uColumnGroup1)) {
          GroupLeng1 <- length(ColumnGroup1[ColumnGroup1 < uColumnGroup1[i]])
          GroupLeng2 <- length(ColumnGroup1[ColumnGroup1 == uColumnGroup1[i]])
          RNA32 <- RNA3[, (GroupLeng1 + NumColumnBlank1 * (i - 1) + 2):(GroupLeng1 +
                                                                          NumColumnBlank1 * (i - 1) + GroupLeng2 + 1)]
          ph3 <- pheatmap(RNA32, clustering_distance_rows = clustering_distance_rows,
                          clustering_method = clustering_method, show_rownames = show_rownames,
                          show_colnames = show_colnames, labels_row = Gene2, border_color = border_color,
                          cluster_rows = FALSE, cluster_cols = cluster_cols,
                          color = colorRampPalette(Color1)(100),
                          fontsize = fontsize, legend = legend1, silent = TRUE)
          plot_list[[i + 1]] <- ph3[[4]]
          ModuleScale2 <- c(ModuleScale2, rep(i + 1, ceiling(ncol(RNA32) / ncol(RNA3) * ModuleScale1)))
        }
      }
      layout_ncol <- length(uColumnGroup1) + 1
      layout_matrix <- matrix(c(1, ModuleScale2), nrow = 1)
    } else {
      ph2 <- pheatmap(RNA3[, 1:ncol(RNA3)], clustering_distance_rows = clustering_distance_rows,
                      clustering_method = clustering_method, show_rownames = show_rownames,
                      show_colnames = show_colnames, labels_row = Gene2, border_color = border_color,
                      cluster_rows = FALSE, cluster_cols = cluster_cols,
                      color = colorRampPalette(Color1)(100),
                      fontsize = fontsize, legend = legend1, silent = TRUE)
      plot_list <- list(ph1[[4]], ph2[[4]])
      layout_ncol <- 2
      layout_matrix <- matrix(c(1, rep(2, ModuleScale1)), nrow = 1)
    }

    g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = plot_list, ncol = layout_ncol,
                                                        layout_matrix = layout_matrix))
  } else {
    pheatmap(RNA3, clustering_distance_rows = clustering_distance_rows,
             clustering_method = clustering_method,
             show_rownames = show_rownames, show_colnames = show_colnames, labels_row = Gene2,
             border_color = border_color, cluster_rows = FALSE, cluster_cols = cluster_cols,
             color = colorRampPalette(Color1)(100), fontsize = fontsize, legend = legend1)
  }
}