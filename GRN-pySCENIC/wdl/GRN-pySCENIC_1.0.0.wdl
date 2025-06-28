version 1.0
#You need to declaration version information(version 1.0)
workflow grn_pyscenic{
  input{
    File rds_h5ad
    String layers="RNA"
    File tf_txt
    File tbl_file
    File feather_file
    Int rank_threshold=5000
    String cluster_key
    Int pyscenic_cpu=10
    Int pyscenic_mem=50
    Int mem_scdatacg=10
    Int mem_get_loom=30
    Int mem_plot=20
    String url_scdatacg="stereonote_hpc/yangdong_6c3a7cd28b5d4861ad87065a5644f7ca_private:latest"
    String url_allscenic="stereonote_hpc/yangdong_4a8497ef4b8d46eea0a6c572c8389923_private:latest"
    String url_pyscenic="stereonote_hpc_external/yangdong_168e88e158bc4b66bda4b279ce72d6b3_private:latest"
  }
  call scdatacg{
    input:
    rds_h5ad=rds_h5ad,
    layers=layers,
    cpu=1,
    mem=mem_scdatacg,
    url=url_scdatacg,
  }
  call get_loom{
    input:
    rds=scdatacg.rds,
    cpu=2,
    mem=mem_get_loom,
    url=url_pyscenic,
  }
  call pyscenic{
    input:
    scenic_loom=get_loom.loom,
    tf_list=tf_txt,
    feather_file=feather_file,
    tbl_file=tbl_file,
    rank_threshold=rank_threshold,
    auc_threshold=0.05,
    cpu=pyscenic_cpu,
    mem=pyscenic_mem,
    url=url_pyscenic,
  }
  call plot{
    input:
    rds=scdatacg.rds,
    h5ad=scdatacg.h5ad,
    aucell_loom=pyscenic.aucell_loom,
    aucell_csv=pyscenic.aucell_csv,
    cluster_key=cluster_key,
    cpu=2,
    mem=mem_plot,
    url=url_allscenic,
  }
  output{
    File pyscenic_plot=plot.pyscenic_plot
    File loom=pyscenic.aucell_loom
    File csv=pyscenic.aucell_csv
    File ctx=pyscenic.ctx
  }
}

task scdatacg{
  input {
    File rds_h5ad
    String layers
    Int cpu
    Int mem
    String url
  }
  command <<<
    #!/bin/bash
    input_file=~{rds_h5ad}
    ext="${input_file##*.}"
    echo "input file extension is: $ext"
    if [ "$ext" == "rds" ]; then
        echo "Converting rds to h5ad..."
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers ~{layers}
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        echo "Copying rds file..."
        cp "$input_file" ./
    elif [ "$ext" == "h5ad" ]; then
        echo "Converting h5ad to rds..."
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers ~{layers}
        echo "Copying h5ad file..."
        cp "$input_file" ./
    else
        echo "Error: Unsupported file extension '$ext'. Only 'rds' and 'h5ad' are supported."
    fi
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi" 
  }
  output {
    File h5ad = glob("*.h5ad")[0]
    File rds = glob("*.rds")[0]
  }
}

task get_loom{
  input {
    File rds
    Int cpu
    Int mem
    String url
  }
  command {
    #get_csv
    Rscript /script/grn_pyscenic/01.csv.R \
    --input_rds ~{rds} --output_csv scenic.data.csv # Image: GRN-SCENIC-R--03 /opt/conda/bin/R
    #csv2loom
    python /script/grn_pyscenic/02.csv2loom.py \
    --input_csv scenic.data.csv --output_loom scenic.loom # Image:GRN-SCENIC-R--03 /opt/conda/bin/python
  }
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File loom = "scenic.loom"
  }
}

task pyscenic{
  input {
    File scenic_loom
    File tf_list
    File feather_file
    File tbl_file
    Int rank_threshold
    Float auc_threshold
    Int cpu
    Int mem
    String url
  }
  Int n_cpus = cpu
  command {
    # step1 grn
    pyscenic grn --num_workers ~{n_cpus} --output grn.tsv --method grnboost2 ~{scenic_loom} ~{tf_list}
    
    # step2 ctx Image: pySCENIC
    pyscenic ctx grn.tsv ~{feather_file} --annotations_fname ~{tbl_file} --expression_mtx_fname ~{scenic_loom} \
    --mode "dask_multiprocessing" --output ctx.csv --num_workers ~{n_cpus} --mask_dropouts \
    --rank_threshold ~{rank_threshold} --auc_threshold ~{auc_threshold}
    # You could use `pyscenic ctx -h` to check its parameters.
    #rank_threshold：控制进入排名的基因或区域数量，默认为5000。较高的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。
    #auc_threshold：控制 AUC 值的阈值，默认为0.05。较低的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。
    
    # step3 AUCell
    pyscenic aucell ~{scenic_loom} ctx.csv --output aucell.loom --num_workers ~{n_cpus}
    pyscenic aucell ~{scenic_loom} ctx.csv --output aucell.csv --num_workers ~{n_cpus}
  }
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File aucell_loom = "aucell.loom"
    File aucell_csv = "aucell.csv"
    File ctx = "ctx.csv"
  }
}

task plot{
  input {
    File rds
    File h5ad
    File aucell_loom
    File aucell_csv
    String cluster_key
    Int cpu
    Int mem
    String url
  }
  command {
    mkdir pyscenic_plot
    cd pyscenic_plot
    Rscript /script/grn_pyscenic/03.plot.R --aucell_loom ~{aucell_loom} \
    --input_rds ~{rds} --cluster_key ~{cluster_key} # Image: GRN-SCENIC-R--03 /opt/conda/bin/R
    
    python /script/grn_pyscenic/Regulon_specificity_Z−scores.py \
    --input_h5ad ~{h5ad} --aucell_csv ~{aucell_csv} --cluster_key ~{cluster_key}
  }
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File pyscenic_plot = "pyscenic_plot"
  }
}