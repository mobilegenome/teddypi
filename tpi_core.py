#!/usr/bin/env python2.7

'''
(c) Fritjof Lammers

fritjoflammers@gmail.com

RetroSeq: assume filename of type {TAXON}.PE.filtered_gaps.{MEITYPE}.{FILTER_STEPS}.vcf

'''

import os, argparse, sys
import numpy as np
from nexus import NexusWriter
from pybedtools import BedTool
from helpers import create_out_path
import subprocess32
from shutil import copyfile


species_dict_ref = {"sample": "UmaBG",
                    "fname": False,
                    "ftype": "fasta",
                    "meitype": None,
                    "f_bedtools": None,
                    "fname_sm": None}

#
# DEFINTIONS
#


def bt_sort(vcf):
    return vcf.sort(header=True)


def bt_merge(vcf, dist):
    return vcf.merge(d=dist, header=True)


def get_file_dic(path, options):
    '''
    This function opens supplied VCF files and saves metadata and variants in the dictionary spec_data.
    '''
    global sample_list

    spec_data = []  # create empty list that will be populated which dictionaries of VCF file information
    file_list = [f for f in os.listdir(path) if
                 options.file_pattern in f]  # get file list from input path and fname pattern

    for fname in file_list:

        if not options.inverse_deletions:  # retroseq results

            # DO VCF to BED conversion
            if fname.endswith(".vcf"):
                subprocess32.call(["vcf2bed", "-d"], stdin=open(os.path.join(path, fname)), stdout=open(os.path.join(options.out_path, fname.replace(".vcf", ".bed")),"w"))
                fname = fname.replace(".vcf", ".bed")
            else:
                copyfile(os.path.join(path, fname), os.path.join(options.out_path,fname))

            fname_splitted = fname.split(".")

            species_dict = {"sample": fname_splitted[0],
                            "fname": fname,
                            "ftype": fname_splitted[-1],
                            "meitype": fname_splitted[3],
                            "f_bedtools": BedTool(os.path.join(options.out_path, fname)),
                            "fname_sm": fname.replace(".gq.bed", ".sorted.merged.gq.bed")}

        elif options.inverse_deletions:

            fname_splitted = fname.split(".")
            species_dict = {"sample": fname_splitted[0],
                            "fname": fname,
                            "ftype": fname_splitted[-1],
                            "meitype": fname_splitted[-2],  #fname_splitted[2],
                            "f_bedtools": BedTool(os.path.join(path, fname)).saveas(os.path.join(options.out_path,fname)),
                            "fname_sm": fname.replace(".bed", ".sorted.merged.bed")}

        print "\n loading %s" % fname,
        print "\t BedTool object length: %i" % (len(species_dict.get("f_bedtools"))),

        if len(species_dict.get("f_bedtools")) < 3 or species_dict.get(
                "meitype") == "DNA":  # filter out empty BedTool object and DNA insetions
            continue
        print "\t performing analyses: ",
        for analyses in prep_analyses:  # perform initial analyses
            print "\t %s" % analyses.__name__,
            species_dict["f_bedtools"] = analyses(species_dict.get("f_bedtools")).saveas(os.path.join(options.out_path, species_dict.get("fname_sm"))) #.saveas(os.path.join(options.out_path, species_dict.get("fname_sm"))) # save again to dictionary

        # species_dict.get("f_bedtools").saveas(
        #     os.path.join(options.out_path, species_dict.get("fname_sm")))  # save to file
        spec_data.append(species_dict)  # append to list
    sample_list = set([l.get("sample") for l in spec_data])

    return spec_data


def estimate_intervals(options):
    '''
    Create intervals from variant coordinates.
    '''
    data = get_file_dic(options.input, options)
    matrices = {}
    for mei_type in set([d.get("meitype") for d in data]):
        data_mei = [d for d in data if d.get("meitype") == mei_type]
        bt_new = BedTool()
        matrices[mei_type] = bt_new.multi_intersect(i=[l.get("f_bedtools").fn for l in data_mei],
                                                    names=[l.get("sample") for l in data_mei],
                                                    header=False)
        matrices.get(mei_type).saveas(os.path.join(options.out_path, mei_type + ".tsv"))

    return matrices


def create_pa_line(l):
    l = l.split("\t")
    samples = set(l[4].split(",")) if "," in l[4] else str(l[4])
    line = []
    for sample in sample_list:  # iterate of available (=present) samples
        if sample in samples:
            line.append(1)  # return 1 if present
        else:
            line.append(0)  # return 0 if absent
    return line


def inverse_samples(samples_present):
    inverted_samples = sample_list - set(samples_present)
    inverted_sample_str = ",".join(inverted_samples)
    return inverted_sample_str


def create_pa_matrix(bt_matrix):
    matrix_lines = str(bt_matrix).split("\n")  # create lines from BedTool object
    loci_list = [line.split("\t")[0:3] for line in matrix_lines]
    pa_matrix = np.matrix([create_pa_line(line) for line in matrix_lines])
    return loci_list, pa_matrix


def create_nexus(bt_matrix, sample_list):
    '''
    create NexusWriter object from BedTool matrix
    BedTool object apparently can contain empty lines, a simple length-check skips these

    requires pybedtools, nexus
    '''

    n = NexusWriter()
    if isinstance(bt_matrix, BedTool):
        matrix_lines = str(bt_matrix).split("\n")  # create lines from BedTool object
    else:
        matrix_lines = bt_matrix
    i = 0
    current_chromosome = ""
    for line in matrix_lines:
        line = line.split("\t")

        if len(line) < 4:  #check for empty line, why 4? it works
            print "skipping: line empty:  %s" % ",".join(line)
            continue

        chrom, start, end, samples_present = line[0], line[1], line[2], line[3].split(
            ",")  # create locus and presnt_ssample objects
        locus = "%s:%i - %i" % (chrom, int(start), int(end))  # define locus
        i += 1

        if chrom != current_chromosome:
            print chrom
            current_chromosome = chrom
        else:
            pass

        for taxon in sample_list:  # iterate over taxons
            presence = 1 if taxon in samples_present else 0
            n.add(taxon, locus, presence)

    return n


def parse_args(args):
    parser = argparse.ArgumentParser(description='Arguments passed to the program')
    parser.add_argument('-i', '--input', required=True, help='input path')
    parser.add_argument('-d', '--distance', default=50, required=False, help='distance for merge', type=int)
    parser.add_argument('-n', '--nexus', action='store_true', default=False, required=False,
                        help='create nexus output?')
    parser.add_argument('-v', '--inverse_deletions', action='store_true', default=False, required=False,
                        help='inverse deletions')
    parser.add_argument('-p', '--file_pattern', required=True, help='file name suffix')
    parser.add_argument('-o', '--out_path', required=True, help='output path')

    return parser.parse_args(args)


def main(args):



    print """
    This is merge_and_intersect.py

    """



    options = parse_args(args)
    if options.inverse_deletions:
        print "You used the -v option, are you really working on deletion, you want to invert?"
    merge_distance = options.distance


    global sample_list
    print "create outdir %s" % create_out_path(options.out_path)
    prep_analyses = [bt_sort, lambda x: bt_merge(x, merge_distance)]  # order of intital preparation steps done for each supplied VCF file
    # prep_analyses = [bt_sort] # order of intital preparation steps done for each supplied VCF file

    matrices = estimate_intervals(options)  # create MEI seperated matrices

    sample_list = sample_list | set(["UmaBG"])
    # sample_list = sample_list # test without reference
    for mei, matrix in matrices.iteritems():

        matrix_lines = str(matrix).split("\n")
        matrix_list = [l.split("\t")[0:3] + [l.split("\t")[4]] for l in matrix_lines if len(l) > 1]
        #print len(matrix_list)
        #if len(matrix_list) < 1:
        #    # empty bed object, skipping
        #    continue # filter out empty BedTool object and DNA insetions

        bt_mei = BedTool(matrix_list).saveas(os.path.join(options.out_path, "matrix_premerge_%s.tsv" % mei))
        #bt_mei = matrix.saveas(os.path.join(options.out_path, "matrix_premerge_%s.tsv" % mei))

        bt_mei = bt_mei.merge(d=100, c=4, o="collapse")
        bt_mei = bt_mei.saveas(os.path.join(options.out_path, "matrix_merge_%s.tsv" % mei))


        if options.nexus:
            print "create NEXUS object for %s" % mei
            nex_object = create_nexus(bt_mei, sample_list)
            print "saving to file %s" % os.path.join(options.out_path, mei + ".nex")
            nex_object.write_to_file(filename=os.path.join(options.out_path, mei + ".nex"), interleave=True,
                                     charblock=True)


        if options.inverse_deletions:
            print "invert",
            matrix_list = []
            for line in str(bt_mei).split("\n"):
                line = line.split("\t")
                if len(line) < 4: continue  # filter empty linesq
                matrix_list.append(tuple(line[0:3] + [inverse_samples(line[3].split(","))]))
            bt_mei = BedTool(matrix_list).saveas(os.path.join(options.out_path, "matrix_inverted_%s.tsv" % mei))
            nex_object = create_nexus(bt_mei, sample_list)
            print "saving to file %s" % os.path.join(options.out_path, mei + ".inversed.nex")
            nex_object.write_to_file(filename=os.path.join(options.out_path, mei + ".inversed.nex"), interleave=True,
                                     charblock=True)


        else:
            print "no NEXUS output will be created, continue"

    print "done, exiting."

if __name__ == '__main__':
    main(sys.argv[1:])
#EOF
