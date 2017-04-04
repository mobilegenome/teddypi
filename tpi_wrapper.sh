#!/bin/bash

# Start per sample filtering and intersection across TE/SV callers
./teddypi.py -c test/teddypi.yaml &&

# Find orthologous variants across samples for
# Ref+ insertions
./tpi_ortho.py  -i test/output/testA.TE.cls.bed \
								-i test/output/testB.TE.cls.bed \
								-i test/output/testC.TE.cls.bed \
								-i test/output/testD.TE.cls.bed \
								-n A \
								-n B \
								-n C \
								-n D \
								--nexus \
								--inverse \ 
								--suffix REF \
								-c test/teddypi.yaml

# and Ref- insertions
./tpi_ortho.py  -i test/output/testA_retroseq.04_REFTEs.vcf \
								-i test/output/testB_retroseq.04_REFTEs.vcf \
								-i test/output/testC_retroseq.04_REFTEs.vcf \
								-i test/output/testD_retroseq.04_REFTEs.vcf \
								-n A \
								-n B \
								-n C \
								-n D \
								--nexus \
								--suffix NONREF \
								-c test/teddypi.yaml

# Unite both datasets
./tpi_unite.py  -i test/output/TPI-DEMO.NONREF.tsv \
								-i test/output/TPI-DEMO.REF.tsv \
								--names "A,B,C,D,outgroup" \
								--nexus \
								-c test/teddypi.yaml
