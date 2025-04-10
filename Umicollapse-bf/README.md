# README

编译方法

cd htslib
autoreconf -i
./configure --prefix=$HOME/local
make -j 8
make install


就可以编译、安装好 htslib 啦。这样会安装到 `~/local` 文件夹里面，后面的 makefile 都是按照这个来的。

使用的话就直接 make 就好。目前使用方法是传入一个输入 BAM 文件，一个输出 BAM 文件名，然后把输入的 BAM 做一点基本的解析。还没有真正实现 UMICollapse 的逻辑。

例子：（记得先加载动态链接库）

```bash
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
./umicollapse-bf /mnt/treasure/asc25/jiazhaopeng/rna/output/SRR23538290.mRNA.genome.mapped.sorted.dedup.unsorted.bam foo | head -n 200
```

