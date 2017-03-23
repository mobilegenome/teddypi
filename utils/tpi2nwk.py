#!/usr/bin/env python2.7

'''
tpi2nwk.py

This script converts BED files with information on present and missing TEs from TeddyPi to bipartions tree in newick format.
Such list of genetrees can be used as input for Multispecies-Coalescent methods. 

Usage:
./tpi2nwk.py INPUT OUTPUT
'''
from sys import argv

loci = []  # list of loci strings (rows from tsv file)
trees = []  # list of trees

# define list of taxa here
taxa = "UmaRef", "Mur14", "Uma02", "Uma01", "Uth16", "UamMi", "UarFi", "Hma13", "Tor17", "Tor18","OUT"
# define nwk template string

tree_template = "(({pres}),({abs}));"
tree_template_1p = "({pres},({abs}));"  # template for loci with only 1 present-taxon because Astral seems to require single taxon denoted without parentheses
# tree_template_0 = "({group});"
print "start tpi2nwk.py ...",
# open input and store lines as list
with open(argv[1]) as fin:
    for line in fin.readlines():
        locus = [e.strip() for e in line.split("\t")]
        loci.append(locus)

# iterate over loci, extract present and missing taxa and create absent taxa
for l in loci:
    if len(l) > 4:
        present = set([e.strip() for e in l[3].split(",") if e !=  ""])
        if len(l) == 6:  # has missing
            missing = set([e.strip() for e in l[5].split(",") if e not in ("NA","")])
            absent = set(taxa) - present - missing
            present -= missing
        else:
            absent = set(taxa) - present

        if len(absent) == 0: # skip loci without "absent" taxa
            if len(present) >= 3:
                pass
                #trees.append(tree_template_0.format(group=",".join(present)))
                #print tree_template_0.format(group=",".join(present))
            else:
                continue
        elif len(present) == 0: # skip loci without "present" taxa
            if len(absent) >= 3:
                pass
               # trees.append(tree_template_0.format(group=",".join(absent)))
               # print tree_template_0.format(group=",".join(absent))
            else:
                continue

        elif (len(absent) + len(present)) < 3: # skip loci with <3 taxa, such gene trees are not informative for Astral
            continue
        elif len(present) == 1:
            trees.append(tree_template_1p.format(pres=",".join(present), abs=",".join(absent)))
        elif len(present) > 1:
            trees.append(tree_template.format(pres=",".join(present), abs=",".join(absent)))  # write outfile
        else:
            continue

with open(argv[2], "w") as fout:
    fout.writelines([e + "\n" for e in trees])

print "finished"
# EOF
