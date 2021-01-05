---
title: "Strandness in RNASeq"
date: '2017-08-17'
slug: Strandness_in_RNASeq
categories: ["data at fingertips"]
tags: ["RNASeq","Bioinformatics"]
---

How to tell whether the paired-end sequencing reads in an RNASeq library are strand-specific or not? According to how read 1 and read 2 align to DNA and RNA sequences, there are three types of RNASeq libraries:

<img src="https://github.com/zhengh42/myfiles/blob/master/research/RNASeq_strandness.png?raw=true" width="480" height="280" />

- If sequences of read 1 align to the RNA strand, the library is "stranded".
- If sequences of read 2 align to the RNA strand, the library is "reversed stranded".
- Sometimes sequences of read 1 align to the RNA strand; the other times sequences of read 2 align to the RNA strand. The library is "unstranded".

Different tools have different names for stranded setting. To name a few:


|	            |Case A	| Case B	| Case C |
| ---         | ---         | ---         | ---         |
|METHODS<br> or KITS	|Ligation<br>Standard SOLiD         |	dUTP<br> IlluminaTruSeq Stranded |	Standard Illumina             |
|TopHat	      |--library-type<br>fr-secondstrand	  | --library-type<br>fr-firststrand	| --library-type<br>fr-unstranded  |
|HTSeq	      |stranded=yes                     |	stranded=reverse              |	stranded=no                   |
|FeatureCounts|-s 1                             |	-s 2                          |	-s 0                          |
|RSEM	        |--forward-prob 1                 |	--forward-prob 0              |	--forward-prob 0.5            |
|Kallisto     |	--fr-stranded	                  | --rf-stranded                 | 	                            |
|Salmon       |	-l ISF                          |	-l ISR                        |	-l IU                         |
|collectRnaSeq<br>Metrics |	FIRST\_READ\_<br>TRANSCRIPTION\_<br>STRAND |	SECOND\_READ\_<br>TRANSCRIPTION\_<br>STRAND	|              |
|Trinity	    |   --SS\_lib\_type FR            |	--SS\_lib\_type RF	          |                               |



More information regarding RNA-Seq analysis pipelines and the parameter settings for the tools can be found:

[GitHub: RNASeq pipeline](https://github.com/zhengh42/RNASeq_pipeline/tree/develop)

[GigaScience: Benchmark of RNA-Seq pipelines](https://academic.oup.com/gigascience/article/8/12/giz145/5663671?guestAccessKey=c350886b-32ec-416a-b941-a5cf68840cb8)

---

What is the most convenient way to tell the strandness? Usually I take about 50 thousands reads from both read 1 and read 2, naming it as test\_1.fg.gz and test\_2.fq.gz, and run one of the following tools:


- <a href="https://github.com/alexdobin/STAR" target="_blank">STAR</a>

```
STAR --genomeDir <reference index directory> --runThreadN <assigned threads> \
      --readFilesIn test_1.fq.gz test_2.fq.gz --readFilesCommand zcat \
      --outFileNamePrefix <out dir and prefix> --outSAMtype BAM  SortedByCoordinate \
      --twopassMode Basic --sjdbOverhang <readlength - 1 > \
      --quantMode TranscriptomeSAM GeneCounts \
      --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD
```

The `--quantMode GeneCounts` option will output a file with suffix "ReadsPerGene.out.tab", which counts the number of reads mapping to each gene. For example:

```
ENSG00000198804.2       555285  290702  281079
ENSG00000198938.2       471151  238612  232541
ENSG00000198886.2       466966  236203  230767
ENSG00000210082.2       404889  203289  201602
ENSG00000198712.1       359278  175022  184256
ENSG00000198727.2       297601  150393  147574
ENSG00000198763.3       288383  149779  138612
ENSG00000156508.17      189858  96202   93663
```

2rd column: counts for unstranded RNASeq (condition C)  
3rd column: counts for the stranded RNASeq (condition A)  
4rd column: counts for the reverse stranded RNASeq (condition B)  

In the case above, the library is unstranded. If 4rd column counts are almost zero, the library is stranded, and vice versa.

- <a href="https://github.com/pachterlab/kallisto" target="_blank">Kallisto</a>

I use the following script to get the strand information. The reads were processed with Kallisto under three library type settings. The output of the three settings were compared and the the library type can be found in the "test.libtype" file.

```
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test.un ${SEQ_DIR}/test_1.fg.gz ${SEQ_DIR}/test_2.fq.gz
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test.rf ${SEQ_DIR}/test_1.fg.gz ${SEQ_DIR}/test_2.fq.gz --rf-stranded
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test.fr ${SEQ_DIR}/test_1.fg.gz ${SEQ_DIR}/test_2.fq.gz --fr-stranded

paste ${OUT_DIR}/test.fr/abundance.tsv ${OUT_DIR}/test.rf/abundance.tsv ${OUT_DIR}/test.un/abundance.tsv | cut -f1,4,9,14  | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' > ${OUT_DIR}/test.libtypetesting
less ${OUT_DIR}/test.libtypetesting | awk '{print $2/$1,$3/$1,$3/$2}' | awk '{if($1<0.3 && $3>3)print "stranded";else if($1>3 && $2>3)print "reverse";else print "unstranded"}' > ${OUT_DIR}/test.libtype
```


- <a href="https://combine-lab.github.io/salmon/getting_started/" target="_blank">Salmon</a>

This tool will automatically detect the strandness of the library.

```
salmon quant -i <transcriptome index> --libType A \
             -o <out dir and prefix> -1 test_1.fq.gz -2 test_2.fq.gz
             -p <assigned threads>
```

The `--libType A` option will allow Salmon to automatically infer the library type. Check the running log for the strand information.
