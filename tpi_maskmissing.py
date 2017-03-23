#!/usr/bin/env python2.7


from tpi_helpers import create_nexus
from pybedtools import BedTool
from sys import argv
in_table = []
out_table = []

# Read new generated output table from Step 1
with open(argv[1]) as fin:
    for l in fin.readlines():
        in_table.append([e.strip() for e in l.split("\t")])

# tidy up present field and generate new table
for row in in_table:
    present = set([e for e in row[3].split(",") if e != ""])
    missing = set(row[5].split(","))
    present = present - missing
    if present:
        new_row = row[:3] + [",".join(present)] + [row[4]] + [",".join(missing)]
        out_table.append("\t".join(new_row ))
    elif not present:
        print "skipping locus, no TE present"
        continue

# save table in BED format
b = BedTool("\n".join(out_table), from_string=True).saveas(argv[2])


# EOF
