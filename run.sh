rm -f /mnt/ramdisk/rna/yyyyxh/output/*
rm -f /mnt/ramdisk/rna/yyyyxh/result/*
rm -f /mnt/ramdisk/rna/yyyyxh/tmp/*

snakemake --snakefile Snakefile-tmp all --cores 64 --config case_id=SRR23538290

snakemake --snakefile Snakefile-tmp all --cores 64 --config case_id=SRR23538291

snakemake --snakefile Snakefile-tmp all --cores 64 --config case_id=SRR23538292

snakemake --snakefile Snakefile-stage2 all --cores 64