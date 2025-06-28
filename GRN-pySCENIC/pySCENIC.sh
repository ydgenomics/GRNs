# Author: lili, ydgenomics
# Date: 20250525
# Prepore files: input.rds; TF_gene_PlantTFDB.txt; ITAG4.1.regions_vs_motifs.rankings.feather; ITAG4.1_MOTIF_PlantTFDB.tbl
# Output files: scenic.data.csv; scenic.loom; grn.tsv; ctx.csv; aucell.loom
# What are .feather and .tbl files?

scenic_loom="/data/work/0.peanut/GRN/peanut/H2014/scenic.loom"
tf_list="/data/work/0.peanut/GRN/peanut/TF_gene_maped.txt"
feather_file="/data/work/0.peanut/GRN/peanut/peanut.regions_vs_motifs.rankings.feather"
tbl_file="/data/work/0.peanut/GRN/peanut/TF_motifs_many2many.tbl"
n_cpus=10
rank_threshold=20000
auc_threshold=0.05

# step1 grn
pyscenic grn --num_workers $n_cpus --output grn.tsv --method grnboost2 $scenic_loom $tf_list

# step2 ctx Image: pySCENIC
pyscenic ctx grn.tsv $feather_file --annotations_fname $tbl_file --expression_mtx_fname $scenic_loom \
--mode "dask_multiprocessing" --output ctx.csv --num_workers $n_cpus --mask_dropouts \
--rank_threshold $rank_threshold --auc_threshold $auc_threshold
# You could use `pyscenic ctx -h` to check its parameters.
#rank_threshold：控制进入排名的基因或区域数量，默认为5000。较高的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。
#auc_threshold：控制 AUC 值的阈值，默认为0.05。较低的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。

# step3 AUCell
pyscenic aucell $scenic_loom ctx.csv --output aucell.loom --num_workers $n_cpus
pyscenic aucell $scenic_loom ctx.csv --output aucell.csv --num_workers $n_cpus

# Rscript 03.plot.R \
# --aucell_loom /data/work/output/pySCENIC/EFM/aucell.loom \
# --input_rds /data/users/yuantingting/yuantingting_fd18e4ec5adc46dd89e2e5c607e9e7e3/online/output/PlantPhone/Seurat/Merge/EFM/RNA_T_0.5/Annotated_EFM_RNA_T_0.5.rds \
# --cluster_key assign.ident # Image: GRN-SCENIC-R--03 /opt/conda/bin/R