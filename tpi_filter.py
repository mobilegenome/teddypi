#!/usr/bin/env python2.7

import sys
import os
import tempfile
from operator import itemgetter

import yaml
import vcf
from pybedtools import BedTool



class LoadVCF():
    def __init__(self, data_dir, out_dir, fname, sname):
        """
        initialize VCF object by reading an VCF file usinf pyvcf
        :param data_dir (str) Input directory of files
        :param out_dir  (str) ouput directory of files
        :param fname  (str) Filename to load
        :param sname  (str) Sample name (identifier)
        """
        self.sample = sname
        self.modulename = "Filter"
        print "[ %s ] Load VCF file %s..." % (self.modulename, fname),
        self.data_dir = data_dir
        self.out_dir = out_dir
        self.fname = fname
        self.vcf_obj = vcf.Reader(open(os.path.join(data_dir, fname), "r"))
        self.vcf_template = self.vcf_obj
        print "done."

        self.vcf_source = str(self.vcf_obj.metadata["source"][0])  # identify creator of VCF file

        print "[ %s ] Load YAML config..." % self.modulename,
        self.config = self.__load_config()

    def __load_config(self):
        """
        Load config file for TE/SV caller
        :return:
        """
        yaml_filename = os.path.join(self.data_dir,
                                     "%s.yaml" % self.vcf_source.split(" ")[0].lower())
        try:
            config = yaml.load(file(yaml_filename, "r"))
            print "loaded YAML config for %s (%s)" % (self.vcf_source, yaml_filename)
        except:
            print "No configuration file provided for %s (%s)" % (self.vcf_source, yaml_filename)
            sys.exit(1)

        return config

    def __iterate_filter(self, entries, assigned_filters):
        """
        Iterate through filters defined in YAML file.
        """
        print "[ %s ] Loading filters..." % self.modulename,
        entries = list(entries)  # convert entries generator to list
        print "done."
        print "[ %s ] Found %i entries in VCF. " % (self.modulename, len(entries))
        print "[ %s ] applying filter" % self.modulename
        print "[ %s ] " % self.modulename,
        passed_filter = {}
        for fstep in sorted(assigned_filters, key=itemgetter('order')):  # start iteration
            if fstep.has_key("deactivate") and fstep['deactivate']:  # if filter has deactive-flag in config: skip
                continue
            else:
                print "%s\t" % fstep['name'],
                entries = self.__apply_filter_step(fstep, entries)
                passed_filter[fstep['order']] = len(entries)
        print "\n[ %s ] " % self.modulename,
        print "\t".join([str(i) for i in passed_filter.values()])

        return entries

    def __apply_filter_step(self, step_dict, entries_list):
        """
        Apply filter method on entries in VCF/BED as defined in class VCFfilter
        """
        self.out_fname = self.fname.replace(".vcf", ".%02d_%s.vcf" % (step_dict['order'], step_dict['name']))
        try:
            vcfout_filtered = vcf.VCFWriter(open(os.path.join(self.out_dir, self.out_fname), "w"), self.vcf_template)
            entries_list = list(
                getattr(VCFfilters(), step_dict['method'])(entries_list, template=self.vcf_template, sample=self.sample, **step_dict))

            for record in entries_list:
                vcfout_filtered.write_record(record)
            vcfout_filtered.close()
        except AttributeError:
            print u"[Error] Method %s not defined." % step_dict['method']
            raise
        return entries_list

    def filter_variants(self):
        self.__iterate_filter(self.vcf_obj, self.config["filters"])
        return 1


class VCFfilters():
    @staticmethod
    def retroseq_gt_hom(entries, **kwargs):
        passed = [record for record in entries if record.genotype(record.samples[0].sample)["GT"] == "1/1"]
        return passed

    @staticmethod
    def retroseq_qual(entries, **kwargs):
        passed = [record for record in entries if (
            record.genotype(record.samples[0].sample)["FL"] == 6 and record.genotype(record.samples[0].sample)[
                "GQ"] >= 28) or (
                      record.genotype(record.samples[0].sample)["FL"] >= 7 and
                      record.genotype(record.samples[0].sample)[
                          "GQ"] >= 20)]
        return passed

    @staticmethod
    def exclude(entries, distance, exfile, **kwargs):
        """
        Performs exclusion operations with bedtools window -v
        :param entries:
        :param distance:
        :param exfile:
        :param kwargs:
        :return:
        """
        # create temporary VCF files
        vcf_temp_in = tempfile.NamedTemporaryFile()  # write IN-VCF to disk temporarily
        vcf_temp_out = tempfile.NamedTemporaryFile()  # write OUT-VCF to disk temporarily
        vcfin = vcf.VCFWriter(vcf_temp_in, kwargs["template"])
        vcfout = vcf.VCFWriter(vcf_temp_out, kwargs["template"])

        for record in entries:  # write IN -VCF
            vcfin.write_record(record)
        vcfin.flush()  # Inbetween flushing (avoid clogging)

        entries = BedTool(vcf_temp_in.name)  # generate BedTool object from VCF we just wrote
        # format string (e.g. for sample-specific exclude files)
        if 'sample' in kwargs.keys():
            exfile = exfile.format(SAMPLE=kwargs['sample'])

        entries = entries.window(exfile, w=distance, v=True, output=vcf_temp_out.name)  # apply window operation
        if len(entries) > 0:
            passed_vcf = vcf.VCFReader(vcf_temp_out)  # return VCF object
        else:
            passed_vcf = [] # return empty list
        vcfin.close()  # close writer, delete temporary file
        # TODO implement closing of vcfout without breaking reading
        return passed_vcf

    @staticmethod
    def svtype(entries, type, **kwargs):
        passed = [record for record in entries if record.INFO["SVTYPE"] == type]
        return passed

    @staticmethod
    def meinfo_type(entries, type, **kwargs):
        passed = [record for record in entries if record.INFO["MEINFO"][0] == type]
        return passed

    @staticmethod
    def svlen(entries, min_size=0, max_size=10000, **kwargs):
        passed = [record for record in entries if
                  abs(record.INFO["SVLEN"]) >= min_size and abs(record.INFO["SVLEN"]) <= max_size]
        return passed


def main(args):
   # inputvcf = LoadVCF("test/testB_retroseq.vcf")
    inputvcf = LoadVCF("test/", "test/output/", "testC_pindel.vcf", "testC")
    inputvcf.filter_variants()

#    inputvcf = LoadVCF("test/", "testC_breakdancer.vcf", "testC")



if __name__ == '__main__':
    main(sys.argv[1:])
