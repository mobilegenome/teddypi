#!/usr/bin/bash

# tpi_maskmissing.sh
# This script screens a TeddyPi BED matrix for intersection with inappropriatly mapped
# regions in using a mask (*.coverage_exclude.bed.midpass).
# Identified loci are then altered in the the BED matrix, so that samples with 'bad'
# coverage are moved to an additional column.
#
# run as ./tpi_maskmissing.sh INPUT OUTPUT

# Declarations
taxon_list="UamMi Mur14 Uma01 Uma02 UarFi Hma13 Uth16 Tor17 Tor18"

TE_MATRIX=$1

# Step 1.
for taxon in $taxon_list;
do
  outfname="$taxon.is_missing.bed"
  outfname_list="$outfname_list $outfname"
  # change the path for -b to your coverage masks
  bedtools intersect -u -a $TE_MATRIX -b coverage/$taxon.coverage_exclude.bed.midpass > $taxon.is_missing
  sed "s/$taxon//g; s/\,\,/,/g" $taxon.is_missing | awk -v t="$taxon" -F "\t" '{OFS="\t"; print $0,t}' > $taxon.is_missing.bed
done


cat $outfname_list | bedtools sort -i - > all_missing.sorted.bed
bedtools merge -o distinct,distinct,distinct -c 4,5,6 -i  all_missing.sorted.bed > all.is_missing.bed
bedtools intersect -v -a $TE_MATRIX -b all.is_missing.bed | awk -F "\t" '{OFS="\t"; print $0,"NA"}' > pruned_matrix.bed
cat all.is_missing.bed pruned_matrix.bed | bedtools sort > all_with_missing.bed

# Step 2.
./tpi_maskmissing.py $2

# cleanup
rm pruned_matrix.bed all_missing.sorted.bed all.is_missing.bed

# EOF
