# `teatime` 
A collection of python script for TE age assignment using TEA-TIME (Transposable Element Age - Time of Integration Mapping in Evolution) approach.

# Quick links
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#documentation)


# Dependencies
- `python` 3.10
- `biopython` 1.83
- `compress-pickle` 2.1.0
- `numpy` 1.21.5
- `pandas` 1.3.5

# Installation

No installation is needed, the scripts can be downloaded directly from `/pipeline` and use as is

# Usage

# Documentation

### Quick links for docmentation
- [Input preparation](#input-preparation)
- [Annotation mending](#annotation-mending)
- [Sequence similarity analysis](#sequence-similarity-analysis)
- [TE age assignment](#te-age-assignment)
- [LTR pair correction](#ltr-pair-correction)

## List of scripts
scripts were named according to the step in the pipeline
- `01progressive_annotation_mending.py`: Annotation mending 
- `02e_value_calculation.py`: Sequence similarity analysis
- `03age_assignment.py`: TE age assignment
- `04ltr_patch.py`: LTR pair correction

## Input preparation

The required inputs for the pipeline are TE coordinates and multispecies alignment files. 

### TE coordinates

The TE coordinates are also extracted from `RepeatMasker` output (publicly available at: https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz). Instead of fomatting the coordinate table as BED file containing `chromosome`, `start`, `end`, and `strand`, the full table was subsetted and extracted as is, the pipeline also requires additional fields such as `id`, `repStart`, `repEnd` and `repLeft` from the `RepeatMasker` output.

### Multispecies alignment

The  updated Zoonomia dataset (primary chromosome only) were obtained from UCSC database (publicly available at: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/maf/). As `MafIO.MafIndex()` from `biopython` does not support a compressed MAF file, all files were unzipped before indexing for fast data retrieval. For this framework, `.mafindex` files are required to have exactly the same name as their `.maf` counterparts as shown below.
```bash
	chr1.maf
	chr1.mafindex
	...
	chrY.maf
	chrY.mafindex
```
The MAF files were indexed using the python code snippet below:
```python
	# build mafindex file 
	from bio import MafIO
	maf_file = '/path/to/maf/'
	target_species = 'hg38'
	maf_id = f'{target_species}.{chrom}'
	mafindex_filedir = '.'.join(str.split(maf_file, sep='.')[:-1]) #remove .maf 
	mafindex_filepath = f'{mafindex_filedir}.mafindex'
	index_maf = MafIO.MafIndex(mafindex_filepath, maf_file, maf_id)
```

## TEA-TIME Pipeline

### Annotation mending

Assign new internal identifiers for TE entries using `01progressive_annotation_mending.py` 

the script can be executed from the command line as follows:
```bash
	python3 /path/to/annotation/mending/script/01progressive_annotation_mending.py \
	-r /path/to/repeatmasker/output.tsv \
	-o /path/to/output/internal_id/directory \
	-t 0.1 \
	-s subfamily1 subfamily2
	#OR
	-S /path/to/subfamily_list.txt 
```

As the consensus position fields are necessary for detecting and merging fragmented TE annotations, a full or a partial RepeatMasker output table is required for this step. The required inputs for this step is a RepeatMasker output file in `.tsv` format, with headers modified to be compatible with `pandas.read_csv()`.

This step is implemented as a command-line script and accepts the following arguments:
- `-r`, `--repeatmasker`: Path to the RepeatMasker output file `.tsv` format.
- `-o`, `--output`: Output directory for the resulting internal ID tables.
- `-s`, `--subfamily_list`: Optional list of TE subfamilies to process. If not provided, all subfamilies in the input file will be processed.
- `-S`, `--subfamily_file`: Optional path to a plain-text file containing a list of TE subfamilies to process (one per line). This option is mutually exclusive with `--subfamily_list`.
- `-t`, `--overlap_threshold`: A float value (default: 0.1) specifying the maximum allowed overlap between fragments for merging.

TE fragments are considered mergeable only if their overlap is below the specified threshold. This reflects the assumption that small overlaps are more likely caused by annotation artifacts, whereas large overlaps are more likely to indicate distinct or nested elements. If the subfamily list is not specified, the script will automatically process all subfamilies present in the RepeatMasker table file.

The outputs of this step are tab-delimited internal ID tables named after corresponding TE subfamily: `{subfamily}.internal_id.txt`. These tables contain the row indices from the original RepeatMasker output and internal IDs for each TE fragment. The internal ID tables act as reference tables for linking information throughout the pipeline, and thus allowing consistent referencing and integration of data from other tables.

### Sequence similarity analysis

Compare TE alignment sequences in multispecies alignment to sequences from the reference species using `02e_value_calculation.py`.

the script can be executed from the command line as follows:
```bash
	python3 /path/to/sequence/similarity/analysis/script/02e_value_calculation.py \
	-r /path/to/repeatmasker/output.tsv \
	-d /path/to/internal_id/directory \
	-x /path/to/external/data/table.tsv \
	-m /path/to/maf/directory \
	-o path/to/e_tables/directory \
	-t hg38 \
	-l 5000 \
	-e 1e-3 \
	-s subfamily1 subfamily2 
	#OR
	-S /path/to/subfamily_list.txt 
```
This step takes four inputs: the path to a directory of internal ID tables from the previous step, a `RepeatMasker` output (for TE coordinates), a directory of multispecies alignment files, and an external data table containing metadata such as `species divergence`, `genome ID`, and `aligned genome size`.

This step is implemented as a command-line script and accepts the following arguments:

- `-r`, `--repeatmasker`: Path to the RepeatMasker output file (`.tsv` format).
- `-d`, `--internal_id_dir`: Directory containing internal ID tables
- `-m`, `--maf_dir`: Directory of multispecies alignment files (`.maf` format)
- `-x`, `--external_data_table`: Path to the external data/divergence table
- `-o`, `--output`: Output directory for the resulting  E-value tables.
- `-s`, `--subfamily_list`: Optional list of TE subfamilies to process. If not provided, all subfamilies in the input file will be processed.
- `-S`, `--subfamily_file`: Optional path to a plain-text file containing a list of TE subfamilies to process (one per line). This option is mutually exclusive with `--subfamily_list`
- `-t`, `--target_species`: Target/reference genome ID (default: `hg38`)
- `-l`, `--extension_length`: Length of flanking sequence around each TE (default: `5000` bp)
- `-e`, `--e_value_cutoff`: E-value cutoff threshold for filtering (default: `1e-3`)

The `target_species` parameter allows users to specify reference genome ID (default: `hg38` for Zoonomia dataset). Th `extension_length` parameter control the length of the window for the surrounding regions (default:`5000` bp). The `e_value_cutoff` sets the threshold for filtering alignment based on statistical significance (`E-value`). If neither `-s` nor `-S is specified, the script will automatically process all subfamilies present in the RepeatMasker table file.

The outputs are tab-delimited E-value tables named after the corresponding TE subfamily: \texttt{{subfamily}.e\_table.txt}. These tables contain alignment statistics for each TE, including E-values, scores, matched bases, gap percentages, divergence times, and flanking/surrounding regions. 

### TE age assignment

Assign TE age estimation to TE internal identifiers using `03age_assignment.py`.

the script can be executed from the command line as follows:
```bash
	python3 /path/to/sequence/similarity/analysis/script/03age_assignment.py \
	-d /path/to/internal_id/directory \
	-v /path/to/E-value/directory
	-r /path/to/repeatmasker/output.tsv \
	-m /path/to/maf/directory \
	-x /path/to/external/data/table.tsv \
	-o path/to/age/tables/directory \
	-t hg38 \
	-l 5000 \
	-e 1e-3 \
	-a 1e-3 \
	-s subfamily1 subfamily2 
	#OR
	-S /path/to/subfamily_list.txt 
```
This step requires the following inputs: a directory of `internal ID` tables, `E-value` tables, along with a `RepeatMasker` output file (for TE coordinates), multispecies alignment files, and external data table used for resolving 0 MYA cases.

This step is implemented as a command-line script and accepts the following arguments:
- `-d`, `--internal_id_dir`: Directory containing internal ID tables
- `-v`, `--e_table_dir`: Directory containing E-value ID tables
- `-r`, `--repeatmasker`: Path to the RepeatMasker output file (`.tsv` format).
- `-m`, `--maf_dir`: 0 MYA case:	Directory of multispecies alignment files (`.maf` format)
- `-x`, `--external_data_table`: 0 MYA case - Path to the external data/divergence table
- `-o`, `--output`: Output directory for the resulting TE age tables.
- `-s`, `--subfamily_list`: (Optional) A list of TE subfamilies to process. If not provided, all subfamilies in the RepeatMasker output will be processed.
- `-S`, `--subfamily_file`: (Optional) Path to a plain-text file containing a list of TE subfamilies to process (one per line). This option is mutually exclusive with `--subfamily_list`
- `-c`, `--additional_evalue_cutoff`: Optional stringent cutoff applied to already filtered E-value tables.
- `-t`, `--target_species`: 0 MYA case:	Target/reference genome ID (default: `hg38`)
- `-l`, `--extension_length`: 0 MYA case:	Length of flanking sequence around each TE (default: `5000` bp)
- `-e`, `--e_value_cutoff`: 0 MYA case:	E-value cutoff threshold for filtering alignment in 0 MYA cases (default: `1e-3`)

An optional E-value filter for additional stringency can be applied using `additional_evalue_cutoff`. If neither `-s` nor `-S` is provided, the script will process all available TE subfamilies in the RepeatMasker table.

The outputs of this step are tab-delimited TE age tables named after the corresponding TE subfamily: `{subfamily}.teatime.txt`, containing internal IDs and their assigned TE ages. For 0 MYA cases, the script also produces additional tables summarizing statistics from surrounding regions: `subfamily.insertion.txt` and `subfamily.segmental.txt` based on categorization from the script.

### LTR pair correction

Readjust LTR age estimates of the intact LTR pairs using `04ltr_patch.py`

the script can be executed like any standard Python script:
```bash
	python3 /path/to/sequence/similarity/analysis/script/04ltr\_patch.py \
	-d /path/to/internal_id/directory \
	-a /path/to/age/table/directory \
	-r /path/to/repeatmasker/output.tsv \
	-o /path/to/corrected/age/tables/directory \
	-m greedy \
	-s subfamily1 subfamily2 
	#OR
	- S /path/to/subfamily_list.txt 
```
This post-assignment step is specific to LTR TEs. The inputs are path to a directory of `TE age` tables, `internal ID` tables, and a `RepeatMasker` output file. 

This step is implemented as a command-line script and accepts the following arguments:
- `-d`, `--internal_id_dir`: Path to the directory containing internal ID tables.
- `-a`, `--age_table_dir`: Path to the directory containing the initial TE age tables.
- `-r`, `--repeatmasker`: Path to the RepeatMasker output file (in `.tsv` format).
- `-o`, `--output`: Output directory for the corrected TE age tables.
- `-s`, `--subfamily_list`: (Optional) A list of TE subfamilies to process. If not provided, all subfamilies in the RepeatMasker output will be processed.
- `-S`, `--subfamily_file`: (Optional) Path to a plain-text file containing a list of TE subfamilies to process (one per line). This option is mutually exclusive with `--subfamily_list`
- `-m`, `--correction_mode`: (Optional) Strategy for age correction. Available options:
	- `strict` (default): Uses only the age of the LTR.
	- `ltr-int`: If a matching internal element (e.g., `LTR-int`) exists, its TE age is included in the 	correction.
	- `greedy`: Attempts to identify any likely matching internal fragment based on naming similarity, even across subfamilies.

As with previous scripts, this script is designed to iterate through all LTR subfamilies. However, users can specify a custom list of TE subfamilies either by passing names directly with `-s` or providing a file with `-S`. If neither is provided, all subfamilies present in the RepeatMasker output will be processed.

The outputs of this step are tab-delimited corrected TE age table named after the corresponding TE subfamily: `subfamily.ltr_fix.txt`, containing internal IDs and their corrected TE ages.



## License
This repository is licensed under the MIT License. 