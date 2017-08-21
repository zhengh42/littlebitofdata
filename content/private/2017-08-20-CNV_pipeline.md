---
title: "Detection of Copy Number Variation in Targeted Sequencing Samples"
date: '2017-08-20'
slug: cnv_pipeline
---

## Introduction
## Methods
## Results

### Ground-truth 1000 Genome CNV data

The processing can be found in /home/zhengh/projects/CNV/105genes/1000g/.

```{r}
########
download
########
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz

########
preprocess
########
pre=estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline
zcat $pre.vcf.gz | egrep -v '^#' | sed 's/;/\t/g' | grep -v CIPOS | awk 'OFS="\t"{print $1,$2,$13,$12,$10,$9}' | sed 's/END=//;s/SAMPLE=//;s/SVTYPE=//;s/CALLID=//' > $pre.preciseCNV.txt
zcat $pre.vcf.gz | egrep -v '^#' | sed 's/;/\t/g' | grep CIPOS | awk 'OFS="\t"{print $1,$2,$15,$12,$10,$9}' | sed 's/END=//;s/SAMPLE=//;s/SVTYPE=//;s/CALLID=//' > $pre.impreciseCNV.txt
cat $pre.preciseCNV.txt $pre.impreciseCNV.txt > $pre.CNV.txt
```


### XHMM
The demo is based on targeted sequencing data of 105 genes in 118 samples.

__1) Calculate depth of coverage from BAM files__

__Script__: /home/zhengh/projects/CNV/105genes/run_DepthOfCoverage.sh

__Input__:  

*	A list of BAM files, generated properly from the established pipeline for variant calling.
*	The BED file of targeted regions in proper format. The chromosome names should match with the BAM files. The chromosome names and starting coordinates should be sorted. Redundant or overlapping regions should be merged (bedtools merge). The choice of regions in the BED file is fairly important, as it is the basic unit for coverage calculation and CNV calling.  

```{r}
work_dir=/home/zhengh/projects/CNV/105genes
target=105genespure
interval_file=/home/database/targetBed/$target.bed

/home/pipelines/DepthOfCoverage.sh \
       $work_dir/bam.list $interval_file $work_dir/QC $target \
       1> $work_dir/logs/$target.DepthOfCoverage.log 2>&1
```

__Output__：

/home/zhengh/projects/CNV/105genes/QC/105genespure.bam.DP.sample\_interval\_summary

```
Target  total_coverage  average_coverage        HG00707_S2_total_cvg    HG00707_S2_mean_cvg     HG00707_S2_granular_Q1 HG00707_S2_granular_median      HG00707_S2_granular_Q3  HG00707_S2_%_above_15
1:17345334-17345466     662773  4983.26 12917   97.12   94      99      104     100.0
1:17349079-17349268     878128  4621.73 16233   85.44   78      87      97      100.0
1:17350434-17350584     499734  3309.50 9242    61.21   58      64      69      100.0
1:17354219-17354360     669584  4715.38 13648   96.11   90      96      105     100.0
1:17355074-17355252     1200155 6704.78 22278   124.46  109     129     146     100.0
1:17359524-17359658     849562  6293.05 15855   117.44  108     123     134     100.0
```

__2) XHMM calling__  

__Script__: /home/zhengh/projects/CNV/105genes/xhmm/run_XHMM.sh

__2.1__) Change data format, remove samples and target regions that do not fulfill coverage requirments, and mean-center the coveage.

a) The"chr" are added to chromosome names to facilitate plots.  

b) Here the default thresholds for coverage are used.   

target mean coverge: between 10 and 500    
sample mean coverage: bewteen 25 and 200    

In this case, 17 samples and 15 targets are removed due to low coverage. They are found in files ending with "filtered\_centered.RD.txt.filtered\_targets.txt" and "filtered\_centered.RD.txt.filtered\_samples.txt".  

Importantly, if the experimental setup is expected to result in mean depths that deviate from these generic values, then these values need to be changed to reflect the dataset and only remove true outlier targets and samples. The default values work well for the current dataset.

c) Center each column of the matrix (RD.txt) so that the mean target value is 0 in preparation for PCA based normalization in the next step.


```{r}
########
echo "change data format. Good luck."
########
less $covdata | sed '2,$s/^/chr/' > $work_dir/xhmm/$pre.DP.sample_interval_summary
$xhmm --mergeGATKdepths -o $work_dir/xhmm/$pre.RD.txt --GATKdepths $work_dir/xhmm/$pre.DP.sample_interval_summary

#######
echo "mean-centers the targets. Good luck."
#######
$xhmm --matrix -r $work_dir/xhmm/$pre.RD.txt --centerData --centerType target \
-o $work_dir/xhmm/$pre.filtered_centered.RD.txt \
--outputExcludedTargets $work_dir/xhmm/$pre.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples $work_dir/xhmm/$pre.filtered_centered.RD.txt.filtered_samples.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150
```

Raw coverage (RD.txt) and the output mean-centered coverage after filtering (filtered\_centered.RD.txt) are shown in __Figure 1A and 1B__. 

__2.2__) Run principal component analysis (PCA) on mean-centered data, and then normalizes the data using PCA information.

PCA analysis is used to determine the strongest independent ways (principal components) in which the data varies and remove the strongest signals that are driven by factors that do not reflect true CNV. 

In this case, the first 2 PCs are removed from the datasets. 

```{r}
########
echo "Runs PCA on mean-centered data. Good luck."
#######
$xhmm --PCA -r $work_dir/xhmm/$pre.filtered_centered.RD.txt --PCAfiles $work_dir/xhmm/$pre.RD_PCA

#######
echo "Normalizes mean-centered data using PCA information. Good luck."
#######
$xhmm --normalize -r $work_dir/xhmm/$pre.filtered_centered.RD.txt --PCAfiles $work_dir/xhmm/$pre.RD_PCA \
--normalizeOutput $work_dir/xhmm/$pre.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
```

In the above codes, the first step outputs 3 files to be used in the next step:

DATA.RD\_PCA.PC.txt: the data projected into the principal components  
DATA.RD\_PCA.PC\_LOADINGS.txt: the loadings of the samples on the principal components  
DATA.RD\_PCA.PC\_SD.txt: the variance of the input read depth data in each of the principal components  

In the second step, the number of principal components to be removed is determined by the `--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7` argument, which indicates that any principal component in which the data variance is greater than 0.7 times the mean variance over all components will be removed. This argument can be changed to retain additional components or remove more components for normalization.

The coverage after PCA normalization is shown in __Figure 1C__.


__2.3__)   Z-score normalization and further filtering

The PCA-normalized coverage data might still have some targets that have very high variance that may be reflective of failed normalization, thus, they are removed in this step. Note, however, that even after this filtering step, exon targets with very high variance in the original read depths still remain with relatively high (though attenuated) variance in depth after normalization, since large differences between samples and targets are by necessity maintained in order to preserve existing CNV signal as well. Then, for each sample, a z-score is calculated by centering relative to all target depths in that sample.

Next, the same targets and samples that are removed during the normalization process are also removed in the original dataset. This matrix will be used for annotation purposes in the subsequent CNV discovery steps.

```{r}
#######
echo "Filters and z-score centers (by sample) the PCA-normalized data. Good luck."
#######
$xhmm --matrix -r $work_dir/xhmm/$pre.PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o $work_dir/xhmm/$pre.PCA_normalized.filtered.sample_zscores.RD.txt \
--outputExcludedTargets $work_dir/xhmm/$pre.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples $work_dir/xhmm/$pre.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
--maxSdTargetRD 30

#######
echo "Filters original read-depth data to be the same as filtered, normalized data. Good luck."
#######
$xhmm --matrix -r $work_dir/xhmm/$pre.RD.txt \
--excludeTargets $work_dir/xhmm/$pre.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets $work_dir/xhmm/$pre.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples $work_dir/xhmm/$pre.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples $work_dir/xhmm/$pre.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o $work_dir/xhmm/$pre.same_filtered.RD.txt
```

The coverage after Z score normalization is shown in __Figure 1D__.

__2.4__)  Genotyping CNVs

A hidden Markov model (HMM) with certain properties that were derived from [a large-scale trio data set](https://www.ncbi.nlm.nih.gov/pubmed/23040492) is applied in this step. The parameters, which determine the rate and length of the CNV called, are specified in params.txt (1e-8    6       70      -3      1.00    0       1.00    3       1.00)

The 9 values in this file correspond to quantities that respectively define the:

- Exome-wide CNV rate 
- Mean number of targets in a CNV call
- Mean distance between targets within a CNV (in KB)
- Mean of DELETION z-score distribution
- Standard deviation of DELETION z-score distribution
- Mean of DIPLOID z-score distribution
- Standard deviation of DIPLOID z-score distribution
- Mean of DUPLICATION z-score distribution
- Standard deviation of DUPLICATION z-score distribution

These parameters are used in the HMM for CNV discovery and genotyping. They may need modification to suit the needs of particular project. For example, to increase the number of CNV calls, increase the first parameter from its default value of 1e-8, e.g., to 1e-7. 

The parameters used for the current datasets are (1e-8    6       70      -3      1.1     0       0.9     3       1.1)

```{r}
#######
echo "Discovers CNVs in normalized data. Good luck."
#######
$xhmm --discover -p $param_adjusted \
-r $work_dir/xhmm/$pre.PCA_normalized.filtered.sample_zscores.RD.txt -R $work_dir/xhmm/$pre.same_filtered.RD.txt \
-c $work_dir/xhmm/$pre.xcnv -a $work_dir/xhmm/$pre.aux_xcnv -s $work_dir/xhmm/$pre
```

Output file (.xcnv):
```
SAMPLE  CNV     INTERVAL        KB      CHR     MID_BP  TARGETS NUM_TARG        Q_EXACT Q_SOME  Q_NON_DIPLOID   Q_START Q_STOP  MEAN_RD MEAN_ORIG_RD
NA12889_S61     DUP     chr7:6017194-6018339    1.15    chr7    6017766 412..413        2       67      88      88      20      25      6.33    135.48
NA18873_S38     DUP     chr7:6017194-6018339    1.15    chr7    6017766 412..413        2       50      50      50      6       29      5.46    59.16
HG00451_S34     DUP     chr7:6012999-6048725    35.73   chr7    6030862 411..425        15      10      99      99      1       53      4.79    157.75
NA18609_S49     DUP     chr7:6017194-6018339    1.15    chr7    6017766 412..413        2       25      36      36      22      7       4.96    107.23
HG00692_S10     DEL     chr7:6012999-6022639    9.64    chr7    6017819 411..414        4       7       99      99      23      5       -4.84   16.56
HG01872_S48     DUP     chr22:29115357-29121386 6.03    chr22   29118371        1698..1700      3       14      70      70      14      8       4.52    136.17
HG00684_S13     DUP     chr7:6012999-6018339    5.34    chr7    6015669 411..413        3       5       99      99      5       21      5.10    79.65
NA18622_S55     DUP     chr2:233198528-233208282        9.76    chr2    233203405       212..218        7       30      99      99      14      50      3.72    76.54
NA18558_S39     DEL     chr15:80454562-80460700 6.14    chr15   80457631        1162..1164      3       35      39      39      4       27      -4.04   29.29
HG00629_S27     DEL     chr7:6012999-6018339    5.34    chr7    6015669 411..413        3       21      99      99      21      5       -5.13   12.96
NA18624_S59     DEL     chr7:6012999-6018339    5.34    chr7    6015669 411..413        3       3       99      99      3       7       -4.88   21.31
HG02522_S57     DUP     chr7:6012999-6018339    5.34    chr7    6015669 411..413        3       20      41      41      21      13      4.11    67.96
NA11995_S51     DEL     chr7:6012999-6018339    5.34    chr7    6015669 411..413        3       9       99      99      9       27      -5.21   0.15
NA07357_S45     DUP     chr19:45855446-45873862 18.42   chr19   45864654        1621..1643      23      5       99      99      10      5       3.03    51.15
NA18628_S59     DUP     chr7:6012999-6018339    5.34    chr7    6015669 411..413        3       9       37      37      9       8       3.96    38.32
HG00634_S26     DUP     chr16:23614761-23615012 0.25    chr16   23614886        1247..1247      1       30      30      30      30      23      7.65    155.58
NA18572_S44     DUP     chr9:271606-289600      18.00   chr9    280603  560..562        3       48      80      80      48      15      4.80    105.55
```
Each line is one CNV called in an individual. Column meaning: 

- SAMPLE	sample name in which CNV was called
- CNV	type of copy number variation (DEL or DUP)
- INTERVAL	genomic range of the called CNV
- KB	length in kilobases of called CNV
- CHR	chromosome name on which CNV falls
- MID_BP	the midpoint of the CNV
- TARGETS	the range of the target indices over which the CNV is called
- NUM_TARG	# of exome targets of the CNV
- Q_EXACT	Phred-scaled quality of the exact CNV event along the entire interval
- Q_SOME	Phred-scaled quality of some CNV event in the interval
- Q\_NON_DIPLOID	Phred-scaled quality of not being diploid, i.e., DEL or DUP event in the interval
- Q_START	Phred-scaled quality of “left” breakpoint of CNV
- Q_STOP	Phred-scaled quality of ‘right” breakpoint of CNV
- MEAN_RD	Mean normalized read depth (z-score) over interval
- MEAN\_ORIG_RD	Mean read depth (# of reads) over interval


__3)__  Compare with ground truth



### CODEX

## Conclusion 

