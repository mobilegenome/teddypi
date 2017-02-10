#!/usr/bin/env python2.7

import os
import argparse
import sys
# from tpi_helpers import make_BED_fromVCF, main_config, create_nexus
import tpi_helpers
from pybedtools import BedTool


def load_variants(flist, nlist):
    """
    This function opens supplied VCF files and saves metadata and variants in the dictionary spec_data.
    """

    datadic = {}  # create empty list that will be populated which dictionaries of VCF file information

    for fname, sname in zip(flist, nlist):
        # DO VCF to BED conversion

        if fname.endswith(".vcf"):
            tpi_helpers.make_BED_fromVCF(fname)
            fname = fname.replace(".vcf", ".bed")  # ?
        # else:
        # copyfile(os.path.join(path, fname), os.path.join(options.out_path, fname))

        bed_data = BedTool(fname)
        print u"[ {} ]".format(tpi_helpers.modulename),
        print u"loading %s" % fname,
        print u"\t %i entries" % (len(bed_data)),
        bed_data = bed_data.sort(header=True).merge(header=True)
        print u"\t %i after sorting and merging overlapping calls" % (len(bed_data))
        if len(bed_data) < 1:  # filter out empty BedTool objects
            print u"[ {} ] ".format(tpi_helpers.modulename),
            print u"Warning: There seems to be no data for {}. It will be skipped".format(sname)
            continue
        datadic[sname] = bed_data

        # print u"\t performing analyses: "
        # SORT AND MERGE
    return datadic


def get_intersection_matrix(data):
    """
    Create intervals from variant coordinates.
    Assume only on TE type is selected and present in the data.
    """
    # data = get_file_dic(options.input, options)
    matrix = BedTool().multi_intersect(i=[bed_data.fn for sample_name, bed_data in data.iteritems()],
                                       names=[sample_name for sample_name, bed_data in data.iteritems()],
                                       header=False)
    return matrix


def parse_args(args):
    parser = argparse.ArgumentParser(description='Arguments passed to the program')
    parser.add_argument('-i', '--input', required=True, action='append', help='Input files, each added by -i FILE')
    parser.add_argument('-c', '--config', required=True, help='Main TeddyPi config file')

    parser.add_argument('-n', '--names', required=True, action='append',
                        help='Sample names for the input files, each added by -i NAME')
    parser.add_argument('-r', '--refname', required=False, default='OUTGROUP', help='Name of reference taxon')
    parser.add_argument('-d', '--distance', default=50, required=False, help='Distance for merge (in bp)', type=int)
    parser.add_argument('--nexus', action='store_true', default=False, required=False,
                        help='Create nexus output?')
    parser.add_argument('-v', '--inverse', action='store_true', default=False, required=False,
                        help='Inverse calls (when using Ref+ calls)')
    # parser.add_argument('-o', '--output', required=True, help='Output path')
    parser.add_argument('--suffix', required=True, help='Suffix for output files. Prefix will be the project name. ')

    return parser.parse_args(args)


def inverse_samples(spresent, slist):
    inverted_samples = slist - set(spresent)
    return ",".join(inverted_samples)


def main(args):
    tpi_helpers.modulename = "ISECT"
    print u"[ {} ] Initiating.".format(tpi_helpers.modulename)

    options = parse_args(args)
    inverse_flag = ""

    config = tpi_helpers.main_config(options.config)
    projectname = config['project_name']

    # Create output path if not exists
    tpi_helpers.create_out_path(config['out_dir'])
    # 1. Load files
    variants = load_variants(options.input, options.names)
    # 2. Do optional VCF > BED conversion
    out_tsv_fname = "{p}.{s}{v}.tsv".format(p=projectname, s=options.suffix, v=inverse_flag)

    print u"[ {} ] Create Multi-intersection matrix for presence / absence patterns".format("ISECT")
    isect_matrix = get_intersection_matrix(variants)
    isect_matrix = isect_matrix.sort().merge(d=options.distance, c=5, o="collapse").saveas(
        os.path.join(config['out_dir'], out_tsv_fname))
    print u"[ {} ] The list of sample names is {}".format("ISECT", ",".join(options.names))
    print u"[ {} ] The outgroup is {}".format("ISECT", options.refname)

    sample_list = set(options.names) | {options.refname}
    # 3. If DEL, do 0<>1 conversion

    if options.inverse:
        inverse_flag = ".inv"
        print u"[ {} ] Performing 0<>1 conversion for Ref+ calls".format("ISECT")

        inverse_matrix = []
        for line in str(isect_matrix).split("\n"):
            line = line.split("\t")
            if len(line) < 4: continue  # filter empty/ incomplete lines lines
            inverse_matrix.append(tuple(line[0:3] + [inverse_samples(line[3].split(","), sample_list)]))  #

        out_tsv_fname = "{p}.{s}{v}.tsv".format(p=projectname, s=options.suffix, v=inverse_flag)
        print u"[ {} ] Save presence absence data in {}".format("ISECT", os.path.join(config['out_dir'], out_tsv_fname))
        isect_matrix = BedTool(inverse_matrix).saveas(os.path.join(config['out_dir'], out_tsv_fname))
    else:
        pass

    if options.nexus:
        out_nex_fname = out_tsv_fname = "{p}.{s}{v}.nex".format(p=projectname, s=options.suffix, v=inverse_flag)
        print u"[ {} ] Save NEXUS files {}".format("ISECT", os.path.join(config['out_dir'], out_nex_fname))
        nex = tpi_helpers.create_nexus(isect_matrix, sample_list)
        nex.write_to_file(filename=os.path.join(config['out_dir'], out_nex_fname), interleave=True,
                          charblock=True)
    print u"[ {} ] Finished.\n".format("ISECT")

    return 1


if __name__ == '__main__':
    main(sys.argv[1:])
