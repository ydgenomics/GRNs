algorithm_distance_py="/GRNs/GRN-gene2role/wdl/tools/SignedS2V/src/algorithm_distance.py"
algorithm_py="/GRNs/GRN-gene2role/wdl/tools/SignedS2V/src/algorithm.py"
graph_py="/GRNs/GRN-gene2role/wdl/tools/SignedS2V/src/graph.py"
main_py="/GRNs/GRN-gene2role/wdl/tools/SignedS2V/src/main.py"
signeds2v_py="/GRNs/GRN-gene2role/wdl/tools/SignedS2V/src/signeds2v.py"
utils_py="/GRNs/GRN-gene2role/wdl/tools/SignedS2V/src/utils.py"

eeisp_py="/GRNs/GRN-gene2role/wdl/codes/eeisp.py"
merge_edgelist_py="/GRNs/GRN-gene2role/wdl/codes/merge_edgelist.py"
spearman_py="/GRNs/GRN-gene2role/wdl/codes/spearman.py"
split_cell_py="/GRNs/GRN-gene2role/wdl/codes/split_cells.py"

pipeline_py="/GRNs/GRN-gene2role/wdl/pipeline.py"

mkdir ./tools
mkdir ./tools/SignedS2V
mkdir ./tools/SignedS2V/emb
mkdir ./tools/SignedS2V/graph
mkdir ./tools/SignedS2V/pickles
mkdir ./tools/SignedS2V/src

cp $algorithm_distance_py ./tools/SignedS2V/src
cp $algorithm_py ./tools/SignedS2V/src
cp $graph_py ./tools/SignedS2V/src
cp $main_py ./tools/SignedS2V/src
cp $signeds2v_py ./tools/SignedS2V/src
cp $utils_py ./tools/SignedS2V/src

mkdir ./codes
cp $eeisp_py ./codes
cp $merge_edgelist_py ./codes
cp $spearman_py ./codes
cp $split_cell_py ./codes

cp $pipeline_py .

n_workers=$1

export PATH="/opt/software/miniconda3/envs/gene2role/bin:$PATH"
if [ -f "sample_merged.edgelist" ]; then
    echo "sample_merged.edgelist existed"
    /opt/software/miniconda3/envs/gene2role/bin/python pipeline.py 1 2 2 sample_merged.edgelist \
    --dimensions 128 --workers $n_workers --cell_metadata metadata.csv
else
    echo "sample_merged.edgelist unexisted"
    /opt/software/miniconda3/envs/gene2role/bin/python pipeline.py 3 2 2 count.csv --workers $n_workers --cell_metadata metadata.csv
fi
