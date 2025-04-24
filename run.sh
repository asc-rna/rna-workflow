export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH && snakemake --snakefile Snakefile-stage1 all --cores 64 --config case_id=SRR23538290
