<!-- （以下均为对Problem_Description，下称文档，的提取）

分为5个阶段：构建索引、清洗数据、基因组比对、排序去重、检测和过滤，
每个阶段用不同软件去处理


初始输入文件在 文档中P8给了链接，建议用SRA Toolkit下载和提取
TODO1:下载一下输入文件，用自己的话描述一下这是什么。（是不是文本文件，能不能直接读？）

P8-P12详细介绍了每个环节的输入输出，用到的软件和常用的命令。
P14-P17 供参考的软件使用的脚本代码
疑似除了最后一步的软件在文档中给出链接的m5C-UBSseq仓库里有外，别的都软件都只给出了名字，得自己去找。

TODO2:下载一下工具，跑一下，用自己的话讲一个每个阶段在做什么，输入输出是啥，记录脚本 -->

## 下载

+ 下载 sra tools: <https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz>

+ 下载输入文件

```
prefetch <ID>
fasterq-dump --fasta <ID>
```

输入文件：ascii 字符串，都是 ACGTU 状物，解压缩后一个数据点 fasta 格式 60+G，fastq 格式 120+G，这题共有 3 个。

小样例：SRR000001，但不知道是否能得到预期的结果。

## 各步骤

+ 建立 reference index

hisat3n-build: <http://daehwankimlab.github.io/hisat2/hisat-3n/>
```
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout -b hisat-3n origin/hisat-3n
make
```

samtools: <https://github.com/samtools/samtools>
```
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
# make install
```

需要 .fasta 格式的输入文件。

运行：hisat3n-build -p <jobs> --base-change C,T <input1>,<input2> <output-basename>

Reading reference sizes 阶段观察到只有一个核在工作，过程中（反反复复）产生多个 .ht2l 文件，IO 非常高，

中途观测到近 300G 的内存占用（还没开始并行）。

然后是（并行）的 sorting 过程。

最后 Returning from GFM constructor, 奇慢无比，cpu 跑不满 100% 也没观察到 IO 瓶颈。

