#!/usr/bin/env python2.7

from collections import defaultdict
import sys
import tempfile
from pybedtools import BedTool
from tpi_filter import LoadVCF
from tpi_helpers import make_BED_fromVCF


class BEDFunctions():
    def __init__(self, bt_instance, ref):
        self.ref = ref
        self.exclude = BedTool(options.coverage_exclude) if options.coverage_exclude else False
        self.bt_object = bt_instance
        self._bt_fname = bt_instance.fn
        self._filter_sequence = []
        self._fout = options.output

    def length(self):
        return len(self.bt_object)

    def __fname(self):
        fname_ending = self._bt_fname[-4:]
        return self._bt_fname.replace(fname_ending, "." + ".".join(self._filter_sequence) + fname_ending)

    def __saveas(self):
        out_fname = self.__fname()
        self.bt_object.saveas(out_fname)

        return 1

    def merge(self):
        __name__ = "merge"

        self._filter_sequence.append(__name__)
        self.bt_object = self.bt_object.merge()
        self.__saveas()
        return self.bt_object

    def filter_gaps(self):
        __name__ = "filter_gaps"

        self._filter_sequence.append(__name__)
        self.bt_object = self.bt_object.window(self.ref.gaps, v=True, w=200)
        self.__saveas()
        return self.bt_object

    def exclude_repeats(self):
        __name__ = "exclude_repeats"

        self._filter_sequence.append(__name__)
        self.bt_object = self.bt_object.window(self.ref.exclude_regions, v=True, w=50)
        self.__saveas()
        return self.bt_object

    def exclude_coverage(self):
        __name__ = "exclude_coverage"

        self._filter_sequence.append(__name__)
        self.bt_object = self.bt_object.intersect(self.exclude, v=True)
        self.__saveas()
        return self.bt_object

    def sort(self):
        __name__ = "sort"

        self._filter_sequence.append(__name__)
        self.bt_object = self.bt_object.sort()
        self.__saveas()
        return self.bt_object

    def intersect_TE(self):
        __name__ = "intersect_TE"

        self._filter_sequence.append(__name__)
        self.bt_object = self.bt_object.window(ref.transposons, w=50)
        self.__saveas()

        return self.bt_object


with open("activeTEs.list") as fin:  # load list of active TEs
    active_TEs = [elem for elem in fin.read().split("\n")]


def refTEs(fname):
    bedtool = BedTool(fname)
    return bedtool


def nonredundant_2_sets(sets):
    # setgen: generator of BedTools obejects
    set_a = sets.next()
    set_b = sets.next()
    intersection = set_a.intersect(set_b, header=True)  # calculate intersection of sets A and B (AuB)
    set_a_complement = set_a.subtract(intersection)  # get complement of A and AB (A\AuB) = exclusive to A
    set_b_complement = set_b.subtract(intersection)  # get complement of B and AB (B\AuB) = exclusive to B
    symdiff_ab = set_a_complement.cat(set_b_complement, postmerge=False, header=True).sort(
        header=True)  # get symmetric difference of A and B
    # combine intersection and
    tmp = tempfile.NamedTemporaryFile()
    nr_set = intersection.cat(symdiff_ab,
                              postmerge=False, header=True).sort(header=True).saveas(tmp.name)
    nr_set = make_BED_fromVCF(tmp.name)  # create BED file
    nr_set_merged = BedTool(nr_set).merge(header=True)  # combine symmetric difference of A and B and intersection = non redundant set
    #os.remove("test/nr_set.sorted.vcf")
    return nr_set_merged


def cluster_calls(variants):
    '''
    variants: BedTool object?
    '''
    clustered_variants = []
    cluster_dic = defaultdict(list)

    # iterate DEL calls and aggregate hits on RepeatMasker track
    for variant in variants:
        variant = str(variant).split("\t")
        RM_hit = variant[3:]
        call = ".".join(variant[0:3])
        cluster_dic[call].append(RM_hit)

    # iterate aggregated calls and ...
    for call, rm_list in cluster_dic.iteritems():
        types = [l[3] for l in rm_list]
        classes = [l[11] for l in rm_list]
        chrom_values = [l[0] for l in rm_list]
        start_values = [l[1] for l in rm_list]
        end_values = [l[2] for l in rm_list]

        # get coordinates of clusters
        start = int(min(start_values))
        end = int(max(end_values))
        size = end - start
        chrom = str(list(set(chrom_values))[0]) if not len(set(chrom_values)) != 1 else False

        # check if call intersects with more than one repeat-type
        if len(set(classes)) != 1:
            match_active_L1 = [f for f in types if f in active_TEs] # are they active TEs
            if "SINE" in classes and not match_active_L1:
                rm_class = "SINE_potential"
                rm_type = ",".join(set(types))
            else:
                rm_class = "complex"
                rm_type = ",".join(set(types))
        else:
            rm_class = str(classes[0])
            if rm_class == "LINE" and size < 5000:  # if LINE smaller 5 kb store as LINE_fragment (Nellaker et al. (2010) Genome Biology)
                rm_class = "LINE_frag"
            rm_type = ",".join(set(types))
        clustered_variants.append(call.split(".") + [chrom, start, end, rm_class, rm_type, size])

    return clustered_variants


def merge_callsets(dsets):
    """

    :param dsets: List of BedTool objects, each representing a dataset to be merged
    :return:
    """
    if len(dsets) == 1:
        print "Only 1 SV set provided. Skip merging."
        return dsets
    elif len(dsets) == 2:
        print "Two SV sets provided. Merging sets,..."
        merged_set = intersect_2_sets(dsets)
        return merged_set
    elif len(dsets) > 2:
        print "More than 2 SV set provided. Merging of three sets not yet implemented. Exiting"
        sys.exit()


def main():
    # 1. Load datasets

    setdict = {"pindel": "test/testA_pindel.vcf",
               "bd": "test/testA_breakdancer.vcf"}

    for setname, setdata in setdict.iteritems():
        set_vcf = LoadVCF(set)
        set_bed = make_BED_fromVCF(set_vcf)
        setdict[setname] = set_bed.sort().merge()

    sets = setdict.itervalues()

    # 2. Merging
    sv_set = merge_callsets(sets)

    # 3. TE intersection
    transposons = refTE("../../data/reference_genome/UrsMar1.rmsk.SINE.bed")
    sv_set = sv_set.window(transposons, w=50).saveas("test/final_del_te.bed")

    # 4. Cluster

    sv_set = cluster_calls(sv_set)


if __name__ == '__main__':
    main()
