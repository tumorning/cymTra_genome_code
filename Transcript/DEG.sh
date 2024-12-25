#!/bin/bash

#reads mapping to genome
#index
hisat2-build cymTra.fa ./index/cymTra
#mapping
for i in *; do hisat2 -p4 --dta -t --phred33 -x ./index/cymTra -1 ${i}.1.fq -2 ${i}.2.fq -S ${i}.sam; done
#sam to bam
for i in *; do samtools view -bS ${i}.sam > ${i}.bam; done
for i in *; do samtools sort ${i}.bam -o ${i}.sorted.bam; done

#count
featureCounts T 4 -F GTF -t exon -g gene_id -s 0 -Q 10 -C -B -p -a cymTra.gtf -o counts.txt *.bam


#DEG
Rscript DEG.R
