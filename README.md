# TEATIME 0.5.1
Transposable Element Age - Time of Integration Mapping in Evolution

**This repo**: a collection of python script for TE age assignment using TEATIME approach

## Quick links
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Older versions](#older-versions)
- [Debug](#debug)
- [Script](#script)
## Dependencies
- ma_mapper(https://github.com/pkncsk/ma_mapper)
- pandas
- numpy
- biopython
- matplotlib

## Usage 
### (optional) annotation mending  
```
/core_engine/annotation_mending.py
```
In this step, fragmented annotations in repeatmasker table were assigned an ID to be calculated as one unit in the next step, by using RepLeft column in the table as a determiner

**!** filtering is limited to subfamily level only

**!** requires repeatmasker output table
### blast and e-value calculation
```
/core_engine/e-value_calculation.py
```
Extract multi-species alingment (MAF) file and calculate BLAST-score (inspired by NCBI-BLAST) and convert them into E-value
### age assignment
```
/core_engine/age_assignment.py
```
After generating E-value table, E-value will be used in age assignment. For each TE, the oldest diveregence from species whose e-value pass the threshold will be assigned as the age of TE

## Older versions
are actually new scripts that replicate logic from older versions (for demonstration)
```
/older_version/teatime_0.1_specieslist.py
```
the simplest logic: scan for the oldest divergence, no filter
```
/older_version/teatime_0.2_blastimplement.py
```
implement NCBI-BLAST scoring
```
/older_version/teatime_0.3_annotationmending.py
```
implement annotation mending, use internal id as proxy
```
/older_version/teatime_0.4_twopass.py
```
calculate and filter from flanked regions to eliminate random alignment

## Debug
Non fucntioned scripts. These should be easier to debug and apply to the main scripts
```
/debug/annotation_mending_dev.py
/debug/eval_debug.py

```
## Script
python script for post-teatime table processing, mostly for visualization and formatting
```
/script/extract_div_from_repeatmasker.py
```
repeatmasker already calculated kimura distance (sort of), so it could be extracted directly from the output table
```
/script/internal_id_age_div_merge_BATCH.py
```
combine internal_id, teatime, and kimura distance tables into one table, convinient for further analysis
```
/script/merge_data_to_repeatmasker.py
```
merging combined table into repeatmasker table (easier for other people to work with)
### teatime analyses
script dedicated for figures in the thesis, generally the analyses on teatime data distribution
```
/script/teatime_analyses/genome_level.py
```
most of the figures came from here
```
/script/teatime_analyses/subfamily_level.py
```
supposedly at subfamily level but might not be used 
```
/script/teatime_analyses/phyloPvsAge_subfamily_level.py
```
plotting mean phyloP vs median teatime of te subfamilies

## License
This repo is licensed under the MIT License. 