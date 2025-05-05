#!/bin/bash

echo "note: You should run this script on i2(corresponding to the largest input)"

# execute command in remote host
executeremotecommand() {
  local host=$1
  local command=$2
  echo "Executing on $host:$command"
  ssh $host "$command"
}

RESULT_PATH=/mnt/treasure/asc25/jiazhaopeng/rna/output-p4_p12
VENV_ACTIVATE_PATH=$HOME/venv-rna/bin/activate
NUM_CORES=64
CASE_ID1=SRR23538290
CASE_ID2=SRR23538291
CASE_ID3=SRR23538292
PWD=$(pwd)

echo pwd is $PWD

# run with tmux to help monitor their behavior"
executeremotecommand "i1" "echo running on i1 && tmux new-session -d -s rna 'source $VENV_ACTIVATE_PATH && cd $PWD && time snakemake --snakefile Snakefile-stage1 all --cores $NUM_CORES --config case_id=$CASE_ID1'"
executeremotecommand "i3" "echo running on i3 && tmux new-session -d -s rna 'source $VENV_ACTIVATE_PATH && cd $PWD && time snakemake --snakefile Snakefile-stage1 all --cores $NUM_CORES --config case_id=$CASE_ID3'"


# run the largest case
echo running on i2 && source $VENV_ACTIVATE_PATH && time snakemake --snakefile Snakefile-stage1 all --cores $NUM_CORES --config case_id=$CASE_ID2

# check if all cases are done
while [ ! -f "$RESULT_PATH/${CASE_ID1}_done" ]; do
  echo "waiting..."
  sleep 0.5
done
while [ ! -f "$RESULT_PATH/${CASE_ID2}_done" ]; do
  echo "waiting..."
  sleep 0.5
done
while [ ! -f "$RESULT_PATH/${CASE_ID3}_done" ]; do
  echo "waiting..."
  sleep 0.5
done

# run stage 2
time snakemake --snakefile Snakefile-stage2 all --cores $NUM_CORES

# snakemake --snakefile Snakefile-stage1 all --cores 64 --config case_id=SRR23538290
# snakemake --snakefile Snakefile-stage2 all --cores 64

