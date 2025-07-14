echo -e "Query_ID\tRefer_ID\tIdentity(%)\tAlignment_Length\tMismatches\tGap_Openings\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score" > header.tsv
n=0
for i in $(ls */*.txt)
do 
  cat header.tsv $i > "$n".txt
  awk -F '\t' '$3 >= 70' "$n".txt > "$n"_filter.txt
  awk '!seen[$1]++' "$n"_filter.txt > "$n"_unique.tsv
  let n++
done

awk 'NR==FNR{a[$2"_"$1]=$1}NR!=FNR{if(a[$1"_"$2])print $1"\t"a[$1"_"$2]}' 0_unique.tsv 1_unique.tsv > reciprocal_best.txt