LOG_TAG=$(date "+%m-%d_%H-%M-%S")
snakemake --snakefile Snakefile-stage2-disk all --cores 64  | tee stage-2-$LOG_TAG.log