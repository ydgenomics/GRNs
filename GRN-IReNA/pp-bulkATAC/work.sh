for c in ~{sep=" " bams}; do
    echo $c >> bams.txt
done
for n in ~{sep=" " names}; do
    echo $n >> names.txt
done

bams_txt="bams.txt"
names_txt="names.txt"
genome_size=1444625381
sh pp_bulkATAC.sh $bams_txt $names_txt $genome_size
sh footprints.sh 