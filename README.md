# TeddyPi: Transposable Element detection and discovery for Phylogentic inference

TeddyPi is a modular pipeline to collect predicted transposable element (TE) and structural variation (SV) calls from a variety of programs and taxa.  TeddiPi filters these variants for those of high-quality. Orthologous TE loci among a set of taxa are utilized to create a presence / absence matrix in NEXUS format for easy phylogenetic inference in phylogeny programs such as PAUP or SplitsTree.

By using the simple-but-powerful YAML syntax for configuration, TeddyPi can easily be adapted for other TE callers. Support for RetroSeq, Pindel and Breakdancer is included.


## Dependencies
TeddyPi runs on Python 2.7. All Python dependencies can be installed using `pip install`.
- BioPython (http://biopython.org)
- pybedtools
- pyVCF
- python-nexus
- pyYAML
- Bedtools (http://bedtools.readthedocs.io/en/latest/)

## Installation

```
git clone https://github.com/mobilegenome/teddypi
```


## Run 
```
./teddypi.py -i conifig.yaml 
```

##  Modules


### teddypi.py
(not existing yet)
Run everything.

### tpi_core.py
(former `merge_and_intersect.py`)

The core functionality in of TeddyPi is to load a set of transposable element variants (TE) from VCF files.
VCF files have to be created per sample, and are intersected by TeddyPi to create a presence / absence matrix.


### tpi_filter.py

Loads VCF files and applies filter operations as defined in YAML-formatted config file.
### tpi_svintegration.py
(former `postprocess_pindel_breakdancer.py`)

This module loads filtered deletion calls files two SV callers and unifies it. Currently it is tailored towards Pindel and Breakdancer, but other SV can be implemented in the future. The module filters out regions with putative bad assembly quality.

### tpi_unite.py
(former `merged_matrix.py`)

The module `tpi_unite.py` combines presence / absence matrices from Ref+ and Ref- calls to final matrix and generates a NEXUS file that can be utilized by phylogeny-software such as PAUP, SplitsTree4 or others.

### utils

Several utilies are included to create VCF files from SV and TE callers that are often proprietary.

## Configuration

TeddyPi is configured with configuration files in the YAML format ( [yaml.org]() ). The main configuration is in `teddypi.yaml` stores the list of samples,list of programs and some parameters.

Example:
```
# teddypi.yaml
samples: testA, testB, testC
ortho_merge_distance: 50

```

For each TE/SV callers an additional configuration file must be provided, that stores general information about the data supplied by the caller and the filtering pipeline to be utilized.

Filters are supplied as list after `filters:`, each list element can contain muliple lines with key:value assignments. The filters are applied by TeddyPi in ascending order given by the `order` key. The `method` field refers to method defined in the `VCFfilters` class in `tpi_filters.py`. By implementing or changing methods, TeddyPi can easily be adapted with new filters or meet requirements for additional or novel TE and SV callers.

Example:
```
name: Pindel
rel: ref # ref or nonref (Ref+ / Ref-)
type: DEL # TE or DEL
filters:
- name: FILTER_1 # Name to identify filter
  order: 01 # This is the first filter step
  method: retroseq_gt_hom # call this method from tpi_filter.py
- name: FILTER_2
  order: 02
  method: method_xyz # another method defined in tpi_filter.py
  additional_key: value # additional parameters for the method
```

## Old Modules

### teddypi.py
(not existing yet)
Run everything.

### tpi_core.py
(former `merge_and_intersect.py`)

The core functionality in of TeddyPi is to load a set of transposable element variants (TE) from VCF files.
VCF files have to be created per sample, and are intersected by TeddyPi to create a presence / absence matrix.

### tpi_svintegration.py
(former `postprocess_pindel_breakdancer.py`)
This module loads VCF files containing deletion calls from two SV callers and unifies it. Currently it is tailored towards Pindel and Breakdancer, but other SV can be implemented in the future. The module filters out regions with putative bad assembly quality.

### tpi_unite.py
(former `merged_matrix.py`)
The module `tpi_unite.py` combines presence / absence matrices from Ref+ and Ref- calls to final matrix and generates a NEXUS file that can be utilized by phylogeny-software such as PAUP, SplitsTree4 or others.

### utils

Several utilies are included to create VCF files from SV and TE callers that are often proprietary.
