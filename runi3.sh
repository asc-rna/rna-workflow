LOG_TAG=$(date "+%m-%d_%H-%M-%S")

docker run --cap-add SYS_NICE \
-v /mnt/ramdisk/rna/input:/data/input:ro \
-v /mnt/ramdisk/docker-output/scy24:/data/output:rw \
-v /mnt/treasure/asc25/scy24/result:/data/result:rw \
--rm -i asc-rna:0219 bash -c "snakemake --snakefile Snakefile-stage1 all --cores 64 --config case_id=SRR23538292" 2>&1 \
| tee i3-docker-$LOG_TAG.log