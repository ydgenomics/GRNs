# Sort TF data of different species

[PlantTFDb](https://planttfdb.gao-lab.org/download.php)

IReNA需要的motif1对象格式(5列：Accession, ID, Name, TFs, EnsembleID)
```shell
head /data/work/test0715/input/Ath_TF_binding_motifs_information2.txt
# Accession       ID      Name    TFs     EnsemblID
# MP00119 MYB_related     MYB_related     AT1G01060       AT1G01060
# MP00120 ERF     ERF     AT1G01250       AT1G01250
# MP00100 bHLH    bHLH    AT1G01260       AT1G01260
# MP00121 NAC     NAC     AT1G01720       AT1G01720
# MP00090 SBP     SBP     AT1G02065       AT1G02065
# MP00122 NAC     NAC     AT1G02230       AT1G02230
# MP00123 NAC     NAC     AT1G02250       AT1G02250
# MP00124 ERF     ERF     AT1G03800       AT1G03800
# MP00125 C2H2    C2H2    AT1G03840       AT1G03840
```

```R
# 基于PlantTFDB下载的TF_blinding_motif内容构建motif1
# wget https://planttfdb.gao-lab.org/download/motif/Ath_TF_binding_motifs_information.txt
data <- read.table('/data/work/test0715/input/Ath_TF_binding_motifs_information.txt', sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(data)

data$ID <- data$Family
data$EnsemblID <- data$Gene_id
data <- data[c("Matrix_id", "ID", "Family", "Gene_id", "EnsemblID")]
colnames(data) <- c("Accession","ID","Name","TFs","EnsemblID")
length(data$TFs)
length(unique(data$TFs))
write.table(data, file = "/data/work/test0715/input/Ath_TF_binding_motifs_information2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# motif1 <- read.table(motif_txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE) # added by yd
# head(motif1) # added by yd

# If your data is un-modle species, you should do blast to change value of EnsemblID column
```