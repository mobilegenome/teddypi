#!/usr/bin/env python

__copyright__ = """
based on breakdancer2vcf.py
https://github.com/ALLBio/allbiotc2/blob/master/breakdancer/breakdancer2vcf.py
modified by  Fritjof Lammers (2015-2018)

"""

__desc__ = """Convert mobster output to pseudo .vcf file format."""

import argparse
import csv
import os.path
import sys
from Bio import SeqIO
from collections import OrderedDict
import datetime

now = datetime.datetime.now()


def main(tsvfile, vcffile):
    '''
    :param tsvfile: filename of input file.tsv
    :type tsvfile: string
    :param vcffile: filename of output file.vcf
    :type vcffile: string
    '''
    with open(tsvfile) as reader:
        # Parse file
        dictreader = _parse_tsvfile(reader)
        print dictreader.fieldnames

        # Write out file
        _format_vcffile(dictreader, vcffile)

    # Quick output
    with open(vcffile) as reader:
        print reader.read(1000)


def _create_out_path(dirname):
    if not os.path.exists(os.path.dirname(dirname)) and os.path.dirname(dirname):
        try:
            os.makedirs(os.path.dirname(dirname))
        except:
            print "Error creating output path"
            sys.exit()
    return dirname


def _parse_tsvfile(readable):
    '''
    Read readable using csv.Sniffer and csv.DictReader
    :param readable: open file.tsv handle to read with csv.DictReader
    :type readable: file
    '''
    prev, curr = 0, 0
    while True:
        line = readable.readline()
        if not line.startswith('#'):
            # lets start from prev # line, without the hash sign
            _tsv_fields = line.split("\t")
            # readable.readline()
            break
        else:
            prev = curr
            curr = readable.tell()

    # Determine dialect
    curr = readable.tell()
    # dialect = csv.Sniffer().sniff(readable.read(3000))
    dialect = 'excel-tab'
    readable.seek(curr)

    # Read file
    dictreader = csv.DictReader(readable, dialect=dialect, fieldnames=_tsv_fields)
    return dictreader


# 'Chr1': '1',
# 'Pos1': '269907',
# 'Orientation1': '39+39-',
# 'Chr2': '1',
# 'Pos2': '270919',
# 'Orientation2': '39+39-',
# 'Type': 'DEL',
# 'Size': '99',
# 'Score': '99',
# 'num_Reads': '38',
# 'num_Reads_lib': '/home/allbio/ERR031544.sort.bam|38',
# 'ERR031544.sort.bam': 'NA'


# def _get_sample_info(fieldnames , field):
# sample_names = fieldnames[11:]
# sample_info_dic = {}
# item_dic = {}
# field = field.split(":")
# items_per_field = [{item.split("|")[0]:item.split("|")[1] for item in field]
# for item in field:
# item_dic[item.split("|")[0]] = item.split("|")[1]
#
#     for sample_name in sample_names:
#         sample_info_dic
#         if sample_name in item_dic.keys():
#             sample_dic = {"sample":sample_name,
#             "read-depth":item_dic.get(item_dic),
#             "genotype":"./."}
#         else:
#             sample_dic = {"sample":sample_name,
#             "read-depth":,
#             "genotype":"./."}
#         formats.append(sample_name)
#
#     format_string = ["%s:%s"%(str(gt) , str(rd)) for gt,rd in zip(genotypes,read_depths)]
#     return "\t".join(format_string)

def _vcf_fields(tsv_fieldnames):
    # sample_names = tsv_fieldnames[11:]
    fieldnames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', "FORMAT",
                  args.sample_name]  # + sample_names #  'FORMAT'
    return fieldnames


def write_vcf(vcf_reader, vcf_list, out_fname):
    with open(out_fname, mode='w') as writer:
        # Write header
        _metadata = OrderedDict([("fileformat", "VCFv4.2"), ("fileDate", now.strftime("%Y-%m-%d %H:%M")),
                                 ("source", "mobster"), ("reference", "Uma")])
        _INFO = [OrderedDict(
            [("ID", "SVTYPE"), ("Number", "1"), ("Type", "String"), ("Description", "\"Type of structural variant\"")]),
                 OrderedDict([("ID", "SVLEN"), ("Number", "1"), ("Type", "Integer"),
                              ("Description", "\"Length of structural variant\"")]),
                 OrderedDict([("ID", "PROGRAM"), ("Number", "1"), ("Type", "String"),
                              ("Description", "\"Program used to call variant\"")]),
                 OrderedDict([("ID", "END"), ("Number", "1"), ("Type", "Integer"),
                              ("Description", "\"End position of the variant described in this record\"")]),
                 OrderedDict([("ID", "DP"), ("Number", "1"), ("Type", "Integer"),
                              ("Description", "\"Total read depth at the locus\"")]),
                 OrderedDict([("ID", "SCORE"), ("Number", "1"), ("Type", "Integer"),
                              ("Description", "\"Quality score\"")]),
                 OrderedDict([("ID", "MEINFO"), ("Number", "4"), ("Type", "String"),
                              ("Description",
                               "\"String,Description=Mobile element info of the form NAME,START,END,POLARITY\"")])]
        _FORMAT = [	OrderedDict([("ID", "GT"), ("Number", "1"), ("Type", "String"), ("Description", "\"Genotype\"") ]),
					OrderedDict([("ID", "DP"), ("Number", "1"), ("Type", "Integer"), ("Description", "\"Depth\"") ])]
        for key, value in _metadata.iteritems():
            writer.write("##{}={}\n".format(key, value))
        for elem in _INFO:
            _INFO_line = ["%s=%s" % (key, value) for key, value in elem.iteritems()]
            writer.write("##INFO=<{}>\n".format(",".join([e for e in _INFO_line])))
        for elem in _FORMAT:
            _FORMAT_line = ["%s=%s" % (key, value) for key, value in elem.iteritems()]
            writer.write("##FORMAT=<{}>\n".format(",".join([e for e in _FORMAT_line])))
        _tsv_fields = vcf_reader.fieldnames
        _no_samples = len(_tsv_fields[11:])

        writer.write('#{}\n'.format('\t'.join(_vcf_fields(_tsv_fields))))
        # Write record
        writer.write(vcf_list)

    return True


def _format_vcffile(dictreader, vcffile):
    '''
    Create a pseudo .vcf file based on values read from DictReader instance.
    :param dictreader: DictReader instance to read data from
    :type dictreader: csv.DictRedaer
    :param vcffile: output file.vcf filename
    :type vcffile: string
    EXAMPLE VCF HEADER
    ##fileformat=VCFv4.0
    ##fileDate=20100101
    ##source=pindel
    ##reference=Uma
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
    ##INFO=<ID=PF,Number=1,Type=Integer,Description="The number of samples carry the variant">
    ##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
    ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=TSD,Number=1,Type=Integer,Description="TSD presence inferred by Mobster (0/1)">
    ##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
    ##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference depth, how many reads support the reference">
    ##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth, how many reads support this allele">

    '''

    _create_out_path(vcffile)

    ref = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    output_vcf = []
    discard_vcf = []
    discard_vcffile = vcffile.replace(".vcf", ".discard.vcf")

    for line in dictreader:
        if "scaffold" in line['Chr']:
            CHROM1 = line['Chr'].replace("chr", "").strip()
        else:
            CHROM1 = line['Chr'].strip()
        # TODO Figure out whether we have zero or one based positioning
        POS = int(line['Insert Point'])

        FILTER = ""
        if POS >= (len(ref.get(CHROM1)) - 500) or POS <= 500:
            print "Prediction at very beginning or end of scaffold, skipping. %i/%i" % (POS, len(ref.get(CHROM1)))
            FILTER = ",".join([FILTER, "START/END"]) if FILTER else  "START/END"
        if "NA" in (line['cluster5 hits'], line['cluster3 hits']):
            DP = [d for d in line['cluster5 hits'], line['cluster3 hits'] if type(d) == int]
            if not DP: DP = "."
            FILTER = ",".join([FILTER, "NO_HIT_SIDE"]) if FILTER else  "NO_HIT_SIDE"
        else:
            DP = sum([int(i) for i in (line['cluster5 hits'], line['cluster3 hits'])])
            if DP < args.qual_threshold:
                FILTER = ",".join([FILTER, "LOW_COV"]) if FILTER else  "LOW_COV"

        ID = "."
        if not FILTER:
            FILTER = "PASS"

        if "START/END" in FILTER:
            REF = "."
        else:
            REF = ref.get(CHROM1)[POS]
        ALT = "<INS>"
        IMPRECISE_START = int(line['border5'])
        IMPRECISE_END = int(line['border3'])
        SVLEN = "."  # int(line['border3']) - int(line['border5'])
        SVTYPE = "INS"  # line['Mobile Element'].strip()
        TSD = "1" if line['target site duplication\n'] == "duplication" else 0
        MEITYPE = line['Mobile Element'].strip()
        SCORE = "."
        QUAL = "."  #int(line['.'])

        INFO ='PROGRAM=mobster;SVTYPE={};SVLEN={};SCORE={};TSD={};MEINFO={},{},{},{};DP={}'.format(SVTYPE,SVLEN,SCORE,TSD,MEITYPE,IMPRECISE_START,IMPRECISE_END,".",DP)
        FORMAT = "GT:DP\t{}:{}".format("1/1", DP)


        # Create record
        if FILTER != "PASS":
            discard_vcf.append([CHROM1, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
        elif FILTER == "PASS":
            output_vcf.append([CHROM1, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
        else:
            print "Error with FILTER flag!"

    # Sort all results
    output_vcf.sort()
    output = "\n".join(["\t".join(map(str, vcf_row)) for vcf_row in output_vcf])
    write_vcf(dictreader, output, vcffile)
    discard_vcf.sort()
    discard = "\n".join(["\t".join(map(str, vcf_row)) for vcf_row in discard_vcf])
    write_vcf(dictreader, discard, discard_vcffile)

    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--mobstertsv', type=str,
                        help='Mobster TSV outputfile', required=True)
    parser.add_argument('-o', '--outputvcf', dest='outputvcf', type=str,
                        help='Output vcf to', required=True)
    parser.add_argument('-q', '--qual_threshold', dest='qual_threshold', type=int,
                        help='Minimum number of reads', default=5, required=False)
    parser.add_argument('-s', '--sample_name', dest='sample_name',
                        help='Minimum number of reads', required=True)
    parser.add_argument('-r', '--reference', dest='reference',
                        help='Reference fasta', required=True)

    args = parser.parse_args()
    _tsv_fields =  ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', "FORMAT",
                  args.sample_name]

    main(args.mobstertsv, args.outputvcf)
