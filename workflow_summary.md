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

## 依赖软件

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

## 下载

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

## 复现步骤

Homo\_sapien 开头的数据一般不用动，如果要改可以通过 make 两个 index\_phony 生成。

单个 case 的 make 写在了单独的 Makefile\_case 里，通过外面指定 CASE\_ID 来决定处理哪个 case，应该也只需要优化这一部分。

直接在这个目录 make all 就可以得到所有答案，make -j3 可以让三个 CASE 并行。

