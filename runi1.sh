
LOG_TAG=$(date "+%m-%d_%H-%M-%S")

docker run --cap-add SYS_NICE -v /mnt/treasure/asc25/zhangyang/rna:/data/input:ro \
-v /mnt/ramdisk/docker-output/scy:/data/output:rw \
-v /mnt/treasure/asc25/scy24/result:/data/result:rw \
--rm -it asc-rna bash -c "snakemake all --cores 64 --config case_id=SRR23538290" | tee i1-docker-$LOG_TAG.log
