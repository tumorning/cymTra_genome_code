#!/bin/bash

#align
for i in OG*; do /home/liunyw/biosoft/mafft/bin/mafft --maxiterate 1000 --localpair $i > ../align/$i 2>> ../3.align.log; done

#gblock
mkdir gblock && cd align

for i in *fa; do /home/liunyw/biosoft/Gblocks_0.91b/Gblocks $i -t=d &>> ../4.gblock.log; sed -i 's/ //g' ${i}-gb; mv *gb ../gblock; mv *-gb.htm ../gblock/; done

cd ..

python3 ../script/filter_gblock.py OG_cds/ gblock/ OG_cds_gblock_filter.txt

mkdir OG_cds_gblock_filter
for i in $(cat OG_cds_gblock_filter.txt); do cp gblock/${i}-gb ./OG_cds_gblock_filter/${i}; done


mkdir iqtree && cd OG_cds_gblock_filter

#concatenate alignment 
/home/liunyw/biosoft/catfasta2phyml/catfasta2phyml.pl * -f > ../iqtree/all.fa 2> ../iqtree/partitions.txt

cd ../iqtree
#iqtree building phylogeny
/home/liunyw/biosoft/iqtree-2.2.2.7-Linux/bin/iqtree2 -s all.fa -m MFP -pre orchid -nt AUTO -B 1000 -wbt


#calculate divergence time
~/newick_utils/bin/nw_reroot ../iqtree/orchid.treefile ambTri > orchid.rooted.tree
#fasta convert to phylip
perl ../../script/Fasta2Phylip.pl ../iqtree/all.fa all.phy

#delet branch length and bootstrap
#(ambTri,(nymCol,(lirTul,((((araTha,rosChi),vitVin),solLyc),(acoGra,((((anaCom,orySat),elaGui),((apoShe,((((cymMan,cymTra),phaEqu),denNob),gasEla)),aspOff)),zosMar))))));
#add the Fossil correction time into the tree
#A. thaliana and R. chinensis (102.0-112.5 Mya), A. comosus and O. sativa (94.1-117.0 Mya), C. mannii and C. tracyanum (9.3-45.0 Mya), C. tracyanum and D. nobile (12.3-51.0 Mya), C. tracyanum and A. officinalis (92.5-118.5 Mya), C. tracyanum and A. shenzhenica (72.2-114.1 Mya), C. tracyanum and A. trichopoda (179.9-205.0 Mya), C. tracyanum and N. colorata (168.4 -191.6 Mya).
#and save as species.tree

#run mcmctree
~/paml-4.10.7/bin/mcmctree mcmctree.ctl
#get divergence time
#(ambTri: 1.819414, (nymCol: 1.779875, (lirTul: 1.735647, ((((araTha: 1.085434, rosChi: 1.085434): 0.005674, vitVin: 1.091108): 0.031307, solLyc: 1.122415): 0.581851, (acoGra: 1.546898, ((((anaCom: 1.076287, orySat: 1.076287): 0.276442, elaGui: 1.352729): 0.019176, ((apoShe: 1.006690, ((((cymMan: 0.039787, cymTra: 0.039787): 0.072416, phaEqu: 0.112203): 0.028617, denNob: 0.140820): 0.838808, gasEla: 0.979628): 0.027062): 0.003754, aspOff: 1.010444): 0.361462): 0.016736, zosMar: 1.388642): 0.158256): 0.157367): 0.031381): 0.044228): 0.039539);

