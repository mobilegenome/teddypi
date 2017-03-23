#!/usr/bin/env python2.7

'''
1. merge tab-seperated matrices
2. create NEXUS object and save
'''

import os
import sys
import argparse
import tpi_helpers
from pybedtools import BedTool


def get_samplelist(input_data, arg_snames=False):
    """
    Create list of samples from input data
     or take given command line argument
    :param: Inputdata from TSV file
    :return: List of samples
    """
    sample_list = set()
    if arg_snames:
        sample_list = [e.strip() for e in arg_snames.split(",")]
    else:
        for line in input_data.split("\n"):
            if len(line.split("\t")) > 3:
                sample_list |= set(line.split("\t")[3].split(","))
    print u"[ {} ] Operating on sample:".format(tpi_helpers.modulename),
    print u",".join(list(sample_list))

    return sample_list


def parse_args(args):
    parser = argparse.ArgumentParser(description='Arguments passed to the program')
    parser.add_argument('-i', '--input', action="append", required=True, help='insertion_matrix')
    #  parser.add_argument('-o', '--output', required=True,help='Output path')
    parser.add_argument('-c', '--config', required=True, help='TeddyPi configuration file')
#    parser.add_argument('-r', '--refname', required=False, default='OUTGROUP', help='Name of reference taxon')
    parser.add_argument('--nexus', help='Create NEXUS output file?', default=False, action="store_true")
    parser.add_argument('--missing', required=False, default=False, action="store_true", help='Are missing states declared?')

    parser.add_argument('-n', '--names', help='List of taxa-names, as comma-separated list. If not given, tries to get list from data.', required=False, default=False)
    return parser.parse_args()


def main(args):
    options = parse_args(args)
    fcontents = ""
    tpi_helpers.modulename = "Unify"

    config = tpi_helpers.main_config(options.config)
    tpi_helpers.config = config
    # Create output path if not exists
    tpi_helpers.create_out_path(config['out_dir'])

    # Load and concatenate input files
        with open(fname) as fin:
            fcontents = fcontents + fin.read().strip() + "\n"
    print u"[ {} ]".format(tpi_helpers.modulename),
    print u"loaded {} entries".format(len(fcontents.split("\n")))

    print u"[ {} ]".format(tpi_helpers.modulename),
    out_tsv_unite = "{p}.unite.tsv".format(p=config['project_name'])
    print u"write concatenated datasets to {} ...".format(os.path.join(config['out_dir'], out_tsv_unite)),
    with open(os.path.join(config['out_dir'], out_tsv_unite), "w") as fout:
        for l in sorted(fcontents.split("\n")):
            fout.write(l + "\n")
    print u"done"

    # If argument --nexus given, do NEXUS conversion
    if options.nexus:
        sample_list = get_samplelist(fcontents, options.names)
        print u"[ {} ]".format(tpi_helpers.modulename),
        print u"create NEXUS object"
        nex_object = tpi_helpers.create_nexus(fcontents, sample_list, has_missing=options.missing)
        out_nex_unite = "{p}.unite.nex".format(p=config['project_name'])
        print u"[ {} ]".format(tpi_helpers.modulename),
        print u"saving to file %s" % os.path.join(config['out_dir'], out_nex_unite)
        nex_object.write_to_file(filename=os.path.join(config['out_dir'], out_nex_unite), interleave=True,
                                 charblock=True)
        print u"[ {} ] Finished. \n".format(tpi_helpers.modulename),

    return 1


if __name__ == '__main__':
    main(sys.argv[1:])
