#!/bin/bash


#blast
/home/liunyw/biosoft/ncbi-blast-2.14.0+/bin/makeblastdb -in cymTra.pep.fa -dbtype prot -parse_seqids -out cymTra
/home/liunyw/biosoft/ncbi-blast-2.14.0+/bin/blastp -query cymTra.pep.fa.fa -db cymTra -out cym2cym.blast.out -evalue 1e-10 -num_threads 16 -outfmt 6 -num_alignments 5

#gff
awk -F '[\t;]' '{if($3=="mRNA")print "ct"$1"\t"$9"\t"$4"\t"$5}' cymTra.gff | grep Chr | sed 's/Chr//g' | sed 's/ID=//g' > ct.gff

#run MCScanX
/home/liunyw/biosoft/MCScanX/MCScanX ct


