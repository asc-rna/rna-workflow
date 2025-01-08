# RNA m5C Modification Site Detection and Performance Optimization Challenge

## What is
这道题主要围绕 RNA 中的 m5C（5-甲基胞嘧啶）修饰位点 的检测与性能优化展开。它要求在高通量测序（HTS）数据基础上，通过一条多工具集成的生物信息学流程（workflow），准确地识别出 RNA 上的 m5C 位点，并对该流程的运行效率进行优化。

- 背景
    - m5C 是 RNA 中重要的化学修饰类型，和基因表达调控、转录后修饰以及蛋白质翻译密切相关。
	- 现有检测方法（RNA-Bis-seq、UBS-seq 等）多基于 C→T 转换信号，但容易带来较高的假阳性率。
	- 本题的目标是综合多种软件（如 hisat3n、samtools、UMICollapse 等），从原始测序数据开始，一步步过滤、比对、统计，最终得到精准且假阳性低的 m5C 位点列表。

## What to do


（1）实现一条多步骤的分析流程

题目明确给出了从原始 FASTQ 到最终 m5C 位点列表的各个主要步骤，包括：
	1.	参考基因组和参考 rRNA/tRNA 的索引构建 : 对参考序列（基因组 / rRNA / tRNA）做 C→T 的“base-change”处理，然后用 hisat3n-build 建索引。
	2.	数据预处理（cutseq）: 去接头、去低质量碱基、去 poly(A) 等操作，得到干净的 FASTQ 数据。
	3.	rRNA / tRNA 过滤 : 按照先后顺序对清洗后的 reads 先比对到 rRNA，再对剩余 reads 比对到 tRNA，以去除这些污染。
	4.	基因组比对（C→T 索引）: 将 rRNA 和 tRNA 都未比对上的 reads，再比对到 C→T 处理的基因组索引上。
	5.	BAM 文件整理、排序、去重 : 利用 samtools sort、UMICollapse 等对比对结果进行排序与重复去除，以便后续统计。
	6.	m5C 位点检测和过滤
	    -	在去重且经过严格过滤的 BAM 文件上，利用 hisat3n-table 或自带的脚本检测可能的 C→T 转换位点。
	    -	根据 unfiltered / filtered、uniq / multi 等不同条件分类统计。
	7.	合并统计，计算重要指标，筛选候选位点
	    -	join_pileup.py: 合并同一样本在不同条件下的碱基统计信息
	    -	group_pileup.py: 合并多个样本（同组）的信息，计算覆盖度、转化率等关键指标
	    -	select_sites.py: 对 group 级别的数据进行初步阈值筛选，得到候选位点
	    -	filter_sites.py: 进一步进行背景校正、统计检验（如二项分布检验），输出最终的 m5C 位点
	8.	跨重复的整合
	    -	最后只保留在多个生物学重复（默认为三次重复）中都显著（p 值 < 10^-6）通过的位点。

（2）评估指标和结果提交
	•	精确度（Precision）$\text{Precision} = \frac{TP}{TP + FP}$
        给定标准答案 (true.tsv) 和你的检测结果 (detected.tsv)，要至少达到 95%。
	•	相关性（Correlation）
	•	计算检测到的位点和标准答案中 “unconverted ratio (ur)” 的 Pearson 相关系数（要求 ≥ 90%）。
	•	提交内容:
    
	1.	Workflow 描述文件（包括中间文件命名、QC 和比对统计结果、每一步的运行时间）
	2.	m5C sites 文件（每个测序数据集对应一个 .filtered.tsv）
	3.	软件或容器打包（如 conda 环境等，保证可复现）
	4.	执行流程的时间截图（从 cutseq 开始到最后结束的总运行时间）

（3）性能优化
	•	当准确度达标后（Precision ≥ 95%、Correlation ≥ 90%），需要考虑如何提升流程的运行效率：
	•	优化多线程设置（-p, -@ 等参数）；
	•	调整内存用量 (samtools sort -m)；
	•	合理安排 I/O 操作；
	•	可能的脚本改进或并行化策略等。
