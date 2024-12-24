#!/bin/bash

#run hmmer
#C terminal
~/hmmer-3.1b2/bin/hmmsearch Terpene_synth_C.hmm cymTra.pep.fa > cymTra_Terpene_synth_C.out

#filter E-value < 0.00001
head -n 237 cymTra_Terpene_synth_C.out | tail -n +17 | awk '{if($1<=0.00001)print $9}' > cymTra_Terpene_synth_C.filter.id.txt

#N terminal
#/hmmer-3.1b2/bin/hmmsearch Terpene_synth.hmm step05.Cym.pep > cymTra_Terpene_synth.out
#filter E-value < 0.00001
head -n 208 cymTra_Terpene_synth.out | tail -n +17 | awk '{if($1<=0.00001)print $9}' > cymTra_Terpene_synth.filter.id.txt

#intersection
grep -f cymTra_Terpene_synth_C.filter.id.txt cymTra_Terpene_synth.filter.id.txt > cymTra_tps.id.txt




for i in apoShe araTha cymTra denNob orySat phaEqu popTri selMoe solLyc; 
	do 
		~/hmmer-3.1b2/bin/hmmsearch Terpene_synth.hmm ~/${i}.pep.fa > ${i}_Terpene_synth.out 
		~/hmmer-3.1b2/bin/hmmsearch /home/zhangshibao/tumengling/tps_family_2024/Terpene_synth_C.hmm /home/zhangshibao/tumengling/tps_family_2024/clean_data/${i}.pep.fa > /home/zhangshibao/tumengling/tps_family_2024/20240107_result/${i}_Terpene_synth_C.out
	done


for i in apoShe araTha cymTra denNob orySat phaEqu popTri selMoe solLyc; 
	do 
		grep -f ${i}_Terpene_synth_C.filter.id.txt ${i}_Terpene_synth.filter.id.txt > ${i}_tps.id.txt
	done


#multiple sequence alignment
~/mafft/bin/mafft --maxiterate 1000 --localpair TPS_phylogenetic_tree_sequences > TPS.mafft.txt 

#trimming
~/trimal-1.4.1/source/trimal -in TPS.mafft.txt -out TPS.mafft.trimal.txt -automated1

#building phylogeny tree
~/iqtree-2.2.2.7-Linux/bin/iqtree2 -s TPS.mafft.trimal.txt -m MFP -pre TPS.tree -nt AUTO -B 1000 -wbt --seqtype AA

