#!/bin/bash

# 函数用于在远程机器上执行命令
executeremotecommand() {
  local host=$1
  local command=$2
  ssh $host "$command"
}

VENV_ACTIVATE_PATH=~/venv-rna/bin/activate
NUM_CORES=64
CASE_ID1=SRR23538290
CASE_ID2=SRR23538291
CASE_ID3=SRR23538292

executeremotecommand "i1" "tmux new-session -d -s rna 'source $(VENV_ACTIVATE_PATH) && cd $(pwd) && snakemake --snakefile Snakefile-stage1 all --cores $NUM_CORES --config case_id=$CASE_ID1'"
executeremotecommand "i2" "tmux new-session -d -s rna 'source $(VENV_ACTIVATE_PATH) && cd $(pwd) && snakemake --snakefile Snakefile-stage1 all --cores $NUM_CORES --config case_id=$CASE_ID2 && snakemake --snakefile Snakefile-stage2 all --cores $NUM_CORES'"
executeremotecommand "i3" "tmux new-session -d -s rna 'source $(VENV_ACTIVATE_PATH) && cd $(pwd) && snakemake --snakefile Snakefile-stage1 all --cores $NUM_CORES --config case_id=$CASE_ID3'"

# snakemake --snakefile Snakefile-stage1 all --cores 64 --config case_id=SRR23538290
# snakemake --snakefile Snakefile-stage2 all --cores 64

