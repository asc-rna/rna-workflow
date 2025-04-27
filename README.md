# RNA m5C

## File structure

```raw
src
├── build.sh            // A script used in Dockerfile
├── Dockerfile
├── eval                // Evaluation scripts
├── hisat-3n
├── htslib
├── m5C-UBSseq
├── README.md
├── runi1.sh            // "docker run" script for the first case
├── runi2.sh            // "docker run" script for the second case
├── runi3.sh            // "docker run" script for the third case
├── samtools
├── Snakefile-stage1    // Snakefile for stage-1
├── Snakefile-stage2    // Snakefile for stage-2
└── Umicollapse
```


## Preparation and running without using docker

1. Clone the Repository:

```sh
git clone git@github.com:asc-rna/rna-workflow.git
cd rna-workflow && git submodule update --init --recursive
```

2. Python env creating

python venv(recommended)

```sh
python -m venv venv-rna
source ~/venv-rna/bin/activate
pip install --upgrade pip
pip install cutseq snakemake polars scipy
```

[miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install)

```sh
conda create --name rna
conda info --envs
conda activate rna

pip install --upgrade pip
pip install cutseq snakemake polars scipy
```

3. Compiling the tools:

```sh
./build.sh
```

You may need to handle the errors manually when running the `build.sh`, e.g. installing `libbz2-dev` and `liblzma-dev`.

4. (Don't do this! Set INPUT DIR to `/mnt/treasure/asc25/rna/input` or `/mnt/ramdisk/rna/input` in the running phase to skip this step) Build indexes for the reference genome and non-coding RNA (ncRNA) if files like `Homo_sapiens.GRCh38.dna.primary_assembly.fa.3n.CT.1.ht2` are missing.

```sh
export INPUT_DIR=/mnt/treasure/asc25/rna/input
export OUTPUT_DIR=/mnt/treasure/asc25/jiazhaopeng/rna/output
make
```

5. When all things are prepared, run the workflow as follows:

Modify the paths in Snakefile-stage1 and Snakefile-stage2 to make sure you have the authority to read/write them. It's suggested to put the input files and intermediate files into ramdisk (`/mnt/ramdisk` in i1-i3) to accelerate IO.

An example of Snakefile-stage1 using `treasure`(disk):

```makefile
OUTPUT_DIR = "/mnt/treasure/asc25/jiazhaopeng/rna/output"
TEMP_DIR = config.get('temp_dir') or "/mnt/treasure/asc25/jiazhaopeng/rna/tmp"
CASE_ID = config.get('case_id') or "SRR23538290"

INPUT_DIR = "/mnt/treasure/asc25/rna/input"
REF_DIR = INPUT_DIR
FASTQ_DIR = "/mnt/treasure/asc25/rna/fastq"

RESULT_DIR = "/mnt/treasure/asc25/jiazhaopeng/rna/result"
```

Then run `snakemake` for every case:

```sh
snakemake --snakefile Snakefile-stage1 all --cores 64 --config case_id=SRR23538290
snakemake --snakefile Snakefile-stage1 all --cores 64 --config case_id=SRR23538291
snakemake --snakefile Snakefile-stage1 all --cores 64 --config case_id=SRR23538292
snakemake --snakefile Snakefile-stage2 all --cores 64
```

You can find the three `.filtered.tsv` file in `RESULT_DIR`.

6. Evaluate. Change the `ANS_DIR` in `eval/eval.sh` to your `RESULT_DIR` and run the script.

```sh
bash eval.sh
```

The output should be like:

```raw
presicion:
97.67
Correlation:
0.99796
```

It is required that $precision >= 95\%$ and $correlation >= 90\%$

## Build and run docker

1. Build docker on three different nodes, for example, i1, i2, i3.

```sh
docker build -t asc-rna:0219 . --network=host
```

2. Modify `runi{1-3}.sh`, change mount directory of <input_dir>, <tmp_dir>, <result_dir>.

3. Prepare input fastq and index file in <input_dir>, including:

```
Homo_sapiens.GRCh38.dna.primary_assembly.fa*
Homo_sapiens.GRCh38.ncrna.fa*
SRR2353829{0,1,2}.fastq
```

For our experiment setup, we create 100G RAM disk named "ramdisk" each for i{1-3}, storing index file and corresponding SRR2353829{0,1,2}.fastq file on each RAM disk.

For our experiment setup, "treasure" is a disk shared by all three nodes, where the results from the three cases are merged and further processed.

For our experiment setup, case2 takes the longest time. so we run stage2 on scripts for i2(case2), immediately after case2-stage1 is finished. It should be reconsiderated if the cases or their processing times change.

4. Run docker simultaneously on three different nodes.

```.sh
# seperately on i{1-3}
./runi{1,2,3}.sh
```


## 静态链接 htslib

（暂时没搞定）

方便起见，我们使用静态链接方法。

```sh
cd htslib
autoreconf -i
export CFLAGS=-flto -O3
export LDFLAGS=-flto
./configure --enable-static --disable-shared
make -j 16
```

经过（漫长的）等待后，得到一个 `libhts.a` 文件。直接把这个文件扔到仓库里面。

但运行的时候需要各种连接，暂时没搞定这个问题

---

## Older README versions

### workflow summary

Updated on 2025.1.16

<!-- （以下均为对Problem_Description，下称文档，的提取）

分为5个阶段：构建索引、清洗数据、基因组比对、排序去重、检测和过滤，
每个阶段用不同软件去处理

初始输入文件在 文档中P8给了链接，建议用SRA Toolkit下载和提取

P8-P12详细介绍了每个环节的输入输出，用到的软件和常用的命令。
P14-P17 供参考的软件使用的脚本代码
-->

克隆仓库后，参考下面内容安装除了 sra-tools 以外的依赖软件，输入数据和产生的中间文件放在 /mnt/treasure/ 里，参考 Makefile 运行。

#### 依赖软件

+ sra-tools（只用于下载数据）: <https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz>

+ hisat3n: <http://daehwankimlab.github.io/hisat2/hisat-3n/>
```
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout -b hisat-3n origin/hisat-3n
make
```

+ samtools: <https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2>
```
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
# make install
```

+ m5C-UBSseq
```
git clone https://github.com/y9c/m5C-UBSseq
```

+ cutseq

```
python3 -m venv cutseq
cd cutseq
source ./bin/activate.sh
pip3 install cutseq
```

+ umicollapse: <https://raw.githubusercontent.com/Daniel-Liu-c0deb0t/UMICollapse/refs/heads/master/umicollapse.jar>

dependency:

```
mkdir lib
cd lib
curl -O -L https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar
curl -O -L https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
cd ..
```

#### 下载

+ 下载 GRCh38

目前在网上找到的较新版本：

dna: <http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz>

ncrna: <http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz>

+ 下载输入文件

```
prefetch <ID>
fasterq-dump --fasta <ID>
```

输入文件：ascii 字符串，都是 ACGTU 状物，fastq 格式 120+G，这题共有 3 个。

小样例：SRR23538290，文件约 20G，是题目 workflow example 里面的例子，应该保证能跑通流程。

#### 复现步骤

Homo\_sapien 开头的数据一般不用动，如果要改可以通过 make 两个 index\_phony 生成。

单个 case 的 make 写在了单独的 Makefile\_case 里，通过外面指定 CASE\_ID 来决定处理哪个 case，应该也只需要优化这一部分。

直接在这个目录 make all 就可以得到所有答案，make -j3 可以让三个 CASE 并行。

### RNA m5C Modification Site Detection and Performance Optimization Challenge

Updated on 2025.1.8

<!-- （以下均为对Problem_Description，下称文档，的提取）

分为5个阶段：构建索引、清洗数据、基因组比对、排序去重、检测和过滤，
每个阶段用不同软件去处理


初始输入文件在 文档中P8给了链接，建议用SRA Toolkit下载和提取
TODO1:下载一下输入文件，用自己的话描述一下这是什么。（是不是文本文件，能不能直接读？）

P8-P12详细介绍了每个环节的输入输出，用到的软件和常用的命令。
P14-P17 供参考的软件使用的脚本代码
疑似除了最后一步的软件在文档中给出链接的m5C-UBSseq仓库里有外，别的都软件都只给出了名字，得自己去找。

TODO2:下载一下工具，跑一下，用自己的话讲一个每个阶段在做什么，输入输出是啥，记录脚本 -->

克隆仓库后，参考下面内容安装除了 sra-tools 以外的依赖软件，输入数据和产生的中间文件放在 /mnt/treasure/ 里，参考 Makefile 运行。

### 依赖软件

+ sra-tools（只用于下载数据）: <https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz>

+ hisat3n: <http://daehwankimlab.github.io/hisat2/hisat-3n/>
```
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout -b hisat-3n origin/hisat-3n
make
```

+ samtools: <https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2>
```
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
# make install
```

+ m5C-UBSseq
```
git clone https://github.com/y9c/m5C-UBSseq
```

+ cutseq

```
python3 -m venv cutseq
cd cutseq
source ./bin/activate.sh
pip3 install cutseq
```

+ umicollapse: <https://raw.githubusercontent.com/Daniel-Liu-c0deb0t/UMICollapse/refs/heads/master/umicollapse.jar>

dependency:

```
mkdir lib
cd lib
curl -O -L https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar
curl -O -L https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
cd ..
```

### 下载

+ 下载 GRCh38

目前在网上找到的较新版本：

dna: <http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz>

ncrna: <http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz>

+ 下载输入文件

```
prefetch <ID>
fasterq-dump --fasta <ID>
```

输入文件：ascii 字符串，都是 ACGTU 状物，fastq 格式 120+G，这题共有 3 个。

小样例：SRR23538290，文件约 20G，是题目 workflow example 里面的例子，应该保证能跑通流程。

### 复现步骤

Homo\_sapien 开头的数据一般不用动，如果要改可以通过 make 两个 index\_phony 生成。

单个 case 的 make 写在了单独的 Makefile\_case 里，通过外面指定 CASE\_ID 来决定处理哪个 case，应该也只需要优化这一部分。

直接在这个目录 make all 就可以得到所有答案，make -j3 可以让三个 CASE 并行。
