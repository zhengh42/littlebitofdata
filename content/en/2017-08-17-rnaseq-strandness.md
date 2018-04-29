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


|	            |Condition A	| Condition B	| Condition C |
| ---         | ---         | ---         | ---         |
|METHODS/KITS	|Ligation, Standard SOLiD         |	dUTP, IlluminaTruSeq Stranded |	Standard Illumina             |
|TopHat	      |--library-type fr-secondstrand	  | --library-type fr-firststrand	| --library-type fr-unstranded  |
|HTSeq	      |stranded=yes                     |	stranded=reverse              |	stranded=no                   |
|FeatureCounts|-s 1                             |	-s 2                          |	-s 0                          |
|RSEM	        |--forward-prob 1                 |	--forward-prob 0              |	--forward-prob 0.5            |
|Kallisto     |	--fr-stranded	                  | --rf-stranded                 | 	                            |
|Salmon       |	-l ISF                          |	-l ISR                        |	-l IU                         |
|collectRnaSeqMetrics |	FIRST\_READ\_TRANSCRIPTION\_STRAND |	SECOND\_READ\_TRANSCRIPTION\_STRAND	|              |
|Trinity	    |   --SS\_lib\_type FR            |	--SS\_lib\_type RF	          |                               |


What is the most convenient way to tell the strandness? Usually I take several thousands reads from both read 1 and read 2, naming it as test\_1.fg.gz and test\_2.fq.gz, and run one of the following tools:


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



- <a href="https://combine-lab.github.io/salmon/getting_started/" target="_blank">Salmon</a>

This tool will automatically detect the strandness of the library.

```
salmon quant -i <transcriptome index> --libType A \
             -o <out dir and prefix> -1 test_1.fq.gz -2 test_2.fq.gz 
             -p <assigned threads>
```

The `--libType A` option will allow Salmon to automatically infer the library type. Check the running log for the strand information. 
