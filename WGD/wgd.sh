#!/bin/sh
set -vex
date

hostname

## make blast db
ncbi-blast-2.15.0/bin/makeblastdb  -in C.tracyanum.pep -dbtype prot -title C.tracyanum.pep
## blastp 
ncbi-blast-2.15.0/bin/blastp -query C.tracyanum.pep -db C.tracyanum.pep -evalue 1e-05 -out C.tracyanum.blas

##cal Ks
export PERL5LIB=''
__conda_setup="$('/export/personal/software/software/miniconda/v4.8.2/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate /export/personal/software/software/evolution/

    export PATH=/export/personal/software/software/NCBI_blast/bin:$PATH
    export PATH=/export/personal/software/Pipeline/evolution/v1.0_test2/./script/WGD/:$PATH
    
###gff
less C.tracyanum.longest.gff3 |sed 's/;/\t/g;s/ID=//g'| awk '$3=="mRNA"' | awk '{print $1"\t"$9"\t"$4"\t"$5}'  > C.tracyanum.gff 

###run MCScanX
    /export/personal/software/Pipeline/evolution/v1.0_test2/./script/WGD//MCScanX/MCScanX C.tracyanum
###KaKs_Calculator calculate ks
    perl  /export/personal/software/Pipeline/evolution/v1.0_test2/./script/WGD//cal_4DTV_ks.pl C.tracyanum.cds C.tracyanum.cds C.tracyanum.collinearity C.tracyanum.out


less C.tracyanum.out.kaks |sed -n '2,$p' |cut -f4 |grep -v "NA" |awk -v spe=${i} '{print spe"\t"$1}' > C.tracyanum.ks

###plot


export R_LIBS_USER=/export/personal/software/software/R_LIB/:$R_LIBS_USERS
/export/pipeline/RNASeq/Software/R/R-3.6.1/bin/Rscript /sugon/personal1/xuxiaoman/projects/hongzhang/06_evolution_v2/04_WGD/KS_plot.R --ks all.ks

date
