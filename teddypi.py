#!/usr/bin/env python2.7

"""
teddypi.py

(c) Fritjof Lammers

"""

import os
import yaml
import argparse
import tpi_filter
import tpi_helpers
from tpi_svintegration import cluster_calls, nonredundant_2_sets
from pybedtools import BedTool
from collections import defaultdict
from sys import argv


def parse_args(args):
    parser = argparse.ArgumentParser(description='Arguments passed to the program')
    parser.add_argument('-c', '--config', required=True, help=' config file path')
    return parser.parse_args(args)


def main():
    """
    Start the TeddyPi pipeline, loads main configuration file, collects input files and parses TE/SV caller specific configuration.

    This module returns filtered and integrated datasets for tpi_ortho.py.
    """
    options = parse_args(argv[1:])
    modulename = "TeddyPi"
    print u"TeddyPi - Transposable Element detection and discovery for Phylogenetic Inference"
    print u"---------------------------------------------------------------------------------\n"

    print u"[ {} ] Initialize configuration from {}...".format(modulename, options.config),

    # Load main configuration
    with open(options.config) as fin:
        config = yaml.load(fin)
        programs = config['programs']
    print u"done."

    tpi_helpers.create_out_path(config['out_dir'])  # Create output directory

    transposons = config['refte']  # Load reference TE file

    # 1. Filter operations for each program and species

    filtered_files = defaultdict(dict)
    for samplename in config['samples']:
        print u"[ {} ] Loading data for sample {}; ".format(modulename, samplename)
        print u"[ {} ] Config has info on these TE/SV callers: {}".format(modulename,
                                                                          ",".join([elem['name'] for elem in programs]))

        per_sample_files = (fname for fname in os.listdir(config['data_dir']) if
                            fname.startswith(samplename) and fname.endswith(
                                ".vcf"))  # TODO avoid reloading processed files

        for sample_file in per_sample_files:
            # print "%s, " % sample_file
            per_sample_vcf = tpi_filter.LoadVCF(data_dir=config['data_dir'],
                                                out_dir=config['out_dir'],
                                                fname=sample_file,
                                                sname=samplename)
            simple_source = per_sample_vcf.vcf_source.split(" ")[0].lower()

            if config['programs'] == "auto" or simple_source in [elem['name'] for elem in programs]:
                per_sample_vcf.skip = False  # flag to skip filtering
                per_sample_vcf.filter_variants()

                print u"[ {} ] Filtered variants written to: {}\n".format(modulename, per_sample_vcf.out_fname)
                filtered_files[samplename][simple_source] = per_sample_vcf.out_fname
            else:
                print u"[ {} ] Error: Auto-detection of TE/SV callers disabled and VCF-source {} not mentioned in " \
                      u"config.\nskipping...".format(modulename, simple_source)

    # 2. Integrate SV-deletions and convert to Ref+ TE calls
    # tpi_svintegration.py

    if 'call_operations' in config.keys():
        print u"[ {} ] Call operations found in configfile".format(modulename)
        for op, sources in config['call_operations'].iteritems():
            try:
                assert (set([elem['name'] for elem in programs]) >= (set(sources)))
            except AssertionError:
                print u"VCF sources for operations have not been parsed."
                print u"[ {} ] For operation {}, sources {} were not parsed. Check \' programs \' parameter in {}" \
                    .format(modulename, op, ",".join(sources), options.config)

                continue

            if op == "non_redundant":
                print u"[ {} ] Starting operation {} on sources {} over all samples,".format(modulename, op,
                                                                                             ",".join(sources))
                for sample in filtered_files.keys():
                    print u"[ %s ] %s " % (op, sample)
                    sets = (BedTool(os.path.join(config['out_dir'], filtered_files[sample][src])) for src in sources)

                    nr = nonredundant_2_sets(sets)
                    nr_set_outfile = "{s}.{t}.nr.bed".format(s=sample, t="DEL")
                    nr_set_outfile = os.path.join(config['out_dir'], nr_set_outfile)
                    nr.saveas(nr_set_outfile)
                    print u"[ {} ] non_redundant set saved to {}".format(op, nr_set_outfile)

                    te_isect_outfile = "{s}.{t}.bed".format(s=sample, t="TE")
                    sv_set = nr.window(transposons, w=config['ortho_merge_distance']).saveas(
                        os.path.join(config['out_dir'], te_isect_outfile))
                    print u"[ {} ] TE intersected set saved to {}".format(op, os.path.join(config['out_dir'],
                                                                                           te_isect_outfile))

                    te_cls_outfile = "{s}.{t}.cls.bed".format(s=sample, t="TE")
                    sv_set = BedTool(cluster_calls(sv_set)).saveas(os.path.join(config['out_dir'], te_cls_outfile))
                    print u"[ {} ] clustered set saved to {}".format(op,
                                                                     os.path.join(config['out_dir'], te_cls_outfile))

            elif op == "intersection":
                print u"[ {} ] Starting operation {} on sources {} over all samples,".format(modulename, op,
                                                                                             ",".join(sources))
                for sample in filtered_files.keys():
                    print u"[ %s ] %s " % (op, sample)
                    sets = (BedTool(os.path.join(config['out_dir'], filtered_files[sample][src])) for src in sources)

                    isect = sets.next().window(sets.next(), w=100, u=True).sort()
                    isect_set_outfile = "{s}.{t}.is.vcf".format(s=sample, t="NONREF_ISEC")
                    isect.saveas(os.path.join(config['out_dir'], isect_set_outfile))

                    print u"[ {} ] intersected set saved to {}".format(op, os.path.join(config['out_dir'],
                                                                                        isect_set_outfile))
            elif op == "te_intersect":
                print u"[ {} ] Starting operation {} on sources {} over all samples,".format(modulename, op,
                                                                                             ",".join(sources))
                for sample in filtered_files.keys():
                    print u"[ %s ] %s " % (op, sample)
                    assert len(sources) == 1
                    src = sources[0]
                    bt_set = tpi_helpers.make_BED_fromVCF(os.path.join(config['out_dir'], filtered_files[sample][src]))

                    te_isect_outfile = "{s}.{t}.bed".format(s=sample, t="TE")
                    sv_set = bt_set.window(transposons, w=50).saveas(
                        os.path.join(config['out_dir'], te_isect_outfile))
                    print u"[ {} ] TE intersected set saved to {}".format(op, os.path.join(config['out_dir'],
                                                                                           te_isect_outfile))

                    te_cls_outfile = "{s}.{t}.cls.bed".format(s=sample, t="TE")
                    sv_set = BedTool(cluster_calls(sv_set)).saveas(os.path.join(config['out_dir'], te_cls_outfile))
                    print u"[ {} ] clustered set saved to {}".format(op,
                                                                     os.path.join(config['out_dir'], te_cls_outfile))
            else:
                print u"[ {} ] Operation '{}' not known. Nothing will be done. Check the configuration file.".format(
                    modulename, op)

    return 1


if __name__ == '__main__':
    main()
