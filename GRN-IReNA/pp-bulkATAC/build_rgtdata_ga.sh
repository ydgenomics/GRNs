sudo mkdir ~/rgtdata/ga
cd ~/rgtdata/ga
sudo cp /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genome_ga.fa .
sudo cp /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/chrom.sizes.ga .
sudo cp /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genes_ga.bed .
sudo cp /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/ga.gtf .
sudo cp /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/alias_ga.txt .
sudo cp /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genome_ga.fa.fai .

cd /data/work/test
rgt-hint footprinting --atac-seq --paired-end --organism=ga \
/data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/merge/four_samples.bam \
/data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/8.peak/peaks.bed