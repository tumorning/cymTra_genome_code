#!/bin/bash

virtualenv -p python3 mcscan_venv
source mcscan_venv/bin/activate
python3 -m pip install jcvi

#Convert the GFF to BED file and rename them.
python -m jcvi.formats.gff bed --type=mRNA --key=ID cymTra.gff -o ct.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID cymMan.gff -o cm.bed
python3 make_mcscan.py denNob.gff > dn.bed


#Pairwise synteny search with last alignment
cd ../ct_cm
python -m jcvi.compara.catalog ortholog cm ct --no_strip_names

cd ../ct_dn
python -m jcvi.compara.catalog ortholog dn ct --no_strip_names

cd ../ct_ct
python -m jcvi.compara.catalog ortholog ct ct --no_strip_names
