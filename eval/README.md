# Evaluation scripts

To evaluate the accuracy according to `GSE225614_HeLa-WT_sites.tsv`: Change your results' path in `eval.sh` and then run it.

```sh
bash eval.sh
```

To evaluate the accuracy according to another three .filtered.tsv files: Change the results' path in `compare.sh` and then run it.

```sh
bash compare.sh
```

## method

From email correspondence with ASC25 Committee:

For each of the three files after the workflow is completed, perform the following operations respectively: Filter out the data items where "passed = true" and "pval <= 1e-6"; then take the intersection、calculate the precision according the standard set reference. 
The supplementary file of this paper `GSE225614_Hela - WT_sites.tsv.gz` is very close to our reference set for correctness, but not exactly the same. This is because this file does not provide clear screening criteria. Before submit the proposal you can send us your three result files, we will helping to check the correctness of the results for you, but we won’t distribute the standard answer. 

Step by step:

1. Take the intersection of the three filtered datasets and then calculate the precision based on the standard set reference.
2. You can calculate the ur value by using `u_total` / `d_total`.
3. The order of the detected sites is not important. It is recommended to simply sort them in the way of "ref + pos".

-----------------

We may adjust the thresholds for correctness evaluation metrics. Additionally, the m5C site counts must not show significant deviation from our reference results. The specific numerical thresholds will be officially disclosed on-site during the final competition.
We will appropriately adjust the correctness thresholds for final competition task to avoid requiring excessive time and effort in studying the minor impacts of application runtime parameters on result accuracy. In general, you can directly reuse the parameters from the preliminary round.

