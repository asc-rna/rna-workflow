ANS_DIR=$1
#ANS_DIR=/mnt/treasure/asc25/jiazhaopeng/rna/results/THU_rlen12000_p16_docker

CASE_ID0=SRR23538290
CASE_ID1=SRR23538291
CASE_ID2=SRR23538292

ANS_FILE0=${ANS_DIR}/${CASE_ID0}.filtered.tsv
FILTERED_ANS_FILE0=./${CASE_ID0}.ans.filtered.tsv
awk -F'\t' -v OFS='\t' '$8=="true" && $7+0<=1e-6' ${ANS_FILE0} > ${FILTERED_ANS_FILE0}
ANS_FILE1=${ANS_DIR}/${CASE_ID1}.filtered.tsv
FILTERED_ANS_FILE1=./${CASE_ID1}.ans.filtered.tsv
awk -F'\t' -v OFS='\t' '$8=="true" && $7+0<=1e-6' ${ANS_FILE1} > ${FILTERED_ANS_FILE1}
ANS_FILE2=${ANS_DIR}/${CASE_ID2}.filtered.tsv
FILTERED_ANS_FILE2=./${CASE_ID2}.ans.filtered.tsv
awk -F'\t' -v OFS='\t' '$8=="true" && $7+0<=1e-6' ${ANS_FILE2} > ${FILTERED_ANS_FILE2}

awk -F'\t' -v OFS='\t' 'NR==FNR {a[$1,$2,$3]=1; b[$1,$2,$3]=$4; c[$1,$2,$3]=$5; next} ($1,$2,$3) in a {print $1,$2,$3,$4+b[$1,$2,$3],$5+c[$1,$2,$3]}' "${FILTERED_ANS_FILE0}" "${FILTERED_ANS_FILE1}" > temp.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR {a[$1,$2,$3]=1; b[$1,$2,$3]=$4; c[$1,$2,$3]=$5; next} ($1,$2,$3) in a {print $1,$2,$3,$4+b[$1,$2,$3],$5+c[$1,$2,$3]}' temp.tsv "${FILTERED_ANS_FILE2}" > $2

rm ${FILTERED_ANS_FILE0} ${FILTERED_ANS_FILE1} ${FILTERED_ANS_FILE2} temp.tsv

echo "the m5C site count in $1 is:"
wc -l < $2

