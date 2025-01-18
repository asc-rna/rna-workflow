LOG_TAG=$(date "+%m-%d_%H-%M-%S")

snakemake --cores 64 --config case_id=SRR23538291 2>&1 | tee snake-$LOG_TAG.log
snakemake --cores 64 --config case_id=SRR23538292 2>&1 | tee snake-$LOG_TAG.log