#!/bin/bash

awk -v OFS="\t" '{$NF=null;if($1=="Orthogroup") print "Desc",$0; else print "(null)",$0}' ~/OrthoFinder/Results_Jan02/Orthogroups/Orthogroups.GeneCount.tsv > orchid.genefamily.txt
cp ~/OrthoFinder/Results_Jan02/Orthogroups/Orthogroups.tsv .


#filter
python3 CAFE5/tutorial/clade_and_size_filter.py -i orchid.genefamily.txt -o orchid.genefamily.filter.cafe5.txt -s

#run cafe4
cafe5 -t orchid.tree -i orchid.genefamily.filter.cafe5.txt
mv results/ base_results

#gene for expansion
grep "$(awk '{if($15>0)print $1}' ../base_results/Base_change.tab | grep OG)" ../Orthogroups.tsv | awk -v FS="\t" '{print $9}' | awk '{for(i = 1; i <= NF; i++) printf("%s%s", $i,"\n")}' | sed 's/,//g' > Increase.geneid.txt

grep "$(grep "$(grep -v 'FamilyID' ../base_results/Base_family_results.txt | awk '{if($3=="y")print $1}')" ../base_results/Base_change.tab | awk '{if($15>0)print $1}')" ../Orthogroups.tsv | awk -v FS="\t" '{print $9}' | awk '{for(i = 1; i <= NF; i++) printf("%s%s", $i,"\n")}' | sed 's/,//g' > Increase.geneid.pvalue.txt

#contracted gene
grep "$(awk '{if($15<0)print $1}' ../base_results/Base_change.tab | grep OG)" ../Orthogroups.tsv | awk -v FS="\t" '{print $9}' | awk '{for(i = 1; i <= NF; i++) printf("%s%s", $i,"\n")}' | sed 's/,//g' > Decrease.geneid.txt
grep "$(grep "$(grep -v 'FamilyID' ../base_results/Base_family_results.txt | awk '{if($3=="y")print $1}')" ../base_results/Base_change.tab | awk '{if($15<0)print $1}')" ../Orthogroups.tsv | awk -v FS="\t" '{print $9}' | awk '{for(i = 1; i <= NF; i++) printf("%s%s", $i,"\n")}' | sed 's/,//g' > Decrease.geneid.pvalue.txt



#enrichment
Rscript enrichment_expansion.R
Rscript enrichment_contracted.R