# ANS_DIR=/mnt/treasure/asc25/jiazhaopeng/rna/results/THU_rlen100000_original
ANS_DIR=/mnt/ramdisk/rna/jzp/output
STD_FILE=/mnt/treasure/asc25/jiazhaopeng/rna/results/THU_rlen12000_p16_test


bash ./intersect.sh ${ANS_DIR} ans.tsv
bash ./intersect.sh ${STD_FILE} std.tsv



precision=$(awk 'NR==FNR {a[$1,$2,$3]=1; next} ($1,$2,$3) in a' std.tsv ans.tsv | wc -l | awk -v total=$(wc -l < ans.tsv) '{printf "%.2f", ($1/total)*100}')

echo "ANS_DIR=${ANS_DIR}"
echo "STD_FILE=${STD_FILE}"
echo "presicion:"
echo ${precision}

echo "Correlation:"
awk -F'\t' -v OFS='\t' 'NR==FNR {a[$1,$2,$3]=$4/$5; next} ($1,$2,$3) in a {printf "%.5f\t%.5f\n",a[$1,$2,$3],$4/$5}' std.tsv ans.tsv > cmp_ur.tsv
awk '{
    sumXY += $1 *$2;
    sumX += $1;
    sumY += $2;
    sumX2 += $1 *$1;
    sumY2 += $2 *$2;
    n++;
} END {
    print(n * sumXY - sumX * sumY) / (sqrt(n * sumX2 - sumX^2) * sqrt(n * sumY2 - sumY^2))
}' cmp_ur.tsv

rm cmp_ur.tsv ans.tsv std.tsv

