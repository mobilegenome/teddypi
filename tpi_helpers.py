import random
import string
import os
import sys
from pybedtools import BedTool
import vcf
import tempfile

modulename = ""
config = ""

def create_out_path(dirname):
    if not os.path.exists(os.path.dirname(dirname)):
        try:
            os.makedirs(os.path.dirname(dirname))
        except:
            print "Error creating output path"
            sys.exit()
    return dirname


def random_string(length=32):
    rand_string = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(length)])
    return rand_string

def main_config(fname):
    import yaml
    with open(fname) as fin:
        return yaml.load(fin)

def count(iter): # count elements in generator
    return sum(1 for _ in iter)


def make_BED_fromVCF(fname, breakpoint=False):
    """

    :param fname: filename of VCF file
    :param breakpoint: indicates whether the POS field is the breakpoint or ?
    :return: BedTool object with BED formatted variants
    """
    fin = open(fname)
    variants = vcf.Reader(fin)
    variant_list = []
    reversed_count = 0
    for record in variants:
        if record.INFO.has_key("END") and record.INFO.has_key("SVLEN"):
            chrom, start, end, size, ori = record.CHROM, record.POS, int(record.INFO["END"][0]), int(record.INFO["SVLEN"][0]), "+"
        else:
            chrom, start, end, size, ori = record.CHROM, record.POS-1, record.POS, "NA", "+"
        if breakpoint:
            end = start  # 1nt length
        if start > end:
            print start, end,
            reversed_count += 1
            ori = "-"
            start, end = end, start
            print start, end
        variant_list.append((chrom, start, end, size))
    if reversed_count > 0:
        print "Warning: reversed %i entries" % reversed_count
    temp_file = tempfile.NamedTemporaryFile()
    bt_list = BedTool(variant_list).sort().saveas(fname.replace(".vcf", ".bed"))

#        bt_list = BedTool(variant_list).sort().saveas(self.fname.replace(".vcf", ".filtered.bed"))
    fin.close()
    return bt_list

def create_nexus(bt_matrix, sample_list, absence_code = 0, has_missing = False):
    '''
    create NexusWriter object from BedTool matrix
    BedTool object apparently can contain empty lines, a simple length-check skips these

    :param has_missing indicate whether input BED file should be parsed for uncertainly called loci, that will be coded as "?"
    requires pybedtools, nexus
    '''
    from nexus import NexusWriter
    n = NexusWriter()
    # if isinstance(bt_matrix, BedTool):
    #     matrix_lines = str(bt_matrix).split("\n")  # create lines from BedTool object
    # else:
    #     matrix_lines = str(bt_matrix).split("\n")
    matrix_lines = str(bt_matrix).split("\n")
    current_chrom = ""
    if has_missing:
        print u"[ {} ] adding missing characters".format(modulename)

    for line in matrix_lines:
        line = line.split("\t")

        if len(line) < 4:  # check for incomplete line
            print u"[ {} ]".format(modulename),
            print u"skipping incomplete or empty line:  %s" % ",".join(line)
            continue

        if not has_missing:

            chrom, start, end, samples_present = line[0], line[1], line[2], line[3].split(
                ",")  # create locus and present_sample objects
            locus = "%s_%i_%i" % (chrom, int(start), int(end))  # define locus

            if current_chrom != chrom:
                print u"[ {} ] Operating on".format(modulename),
                print chrom
                current_chrom = chrom
            # add locus to nexus object
            for taxon in sample_list:  # iterate over taxons
                char = 1 if taxon in samples_present else absence_code
                n.add(taxon, locus, char)

        elif has_missing:
            chrom, start, end, samples_present, type, samples_missing = line[0], line[1], line[2], line[3].split(
                ","), line[4], line[5].split(",")  # create locus and present_sample objects
            locus = "%s_%i_%i" % (chrom, int(start), int(end))  # define locus

            if current_chrom != chrom:
                print u"[ {} ] Operating on".format(modulename),
                print chrom
                current_chrom = chrom

            # add locus to nexus object
            for taxon in sample_list:  # iterate over taxons
                if taxon in samples_present:
                    char = 1
                elif taxon in samples_missing:
                    char = "?"
                else:
                    char = absence_code
                n.add(taxon, locus, char)

    return n