#%%
import pandas as pd
from Bio.AlignIO import MafIO
import numpy as np
import os
import re
import math
from itertools import repeat
from concurrent.futures import ProcessPoolExecutor
from Bio.Seq import Seq
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
#%%
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table
#%%
repeatmasker_path = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
repeatmasker_table=repeatmasker_prep(repeatmasker_path)
#%%
repeatmasker_mwe=repeatmasker_table[repeatmasker_table.index.isin(range(7633,12476))]
#%%
repeatmasker_mwe.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/test/repeatmasker.mwe.txt', sep='\t')
# %%
repeatmasker_mwe
# %%
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table
#%%
repeatmasker_path = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
repeatmasker_table=repeatmasker_prep(repeatmasker_path)


# %%
from Bio import SeqIO

from Bio import SeqIO

fasta_path = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa"
output_path = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/test/hg38.mwe.fasta"

chrom = "chr1"
start = 1
end = 7798441

# Load the fasta file indexed dictionary
fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

# Extract sequence (remember Python slicing is zero-based, end-exclusive)
seq = fasta_dict[chrom].seq[start - 1 : end]

# Save to new fasta file with header
with open(output_path, "w") as out_f:
    out_f.write(f">{chrom}:{start}-{end}\n")
    out_f.write(str(seq) + "\n")

#%%
def block_overlaps(query_start, query_end, block_start, block_length):
    block_end = block_start + block_length
    return not (block_end < query_start or block_start > query_end)

def filter_maf_blocks(maf_in_path, maf_out_path, chrom, start, end):
    with open(maf_in_path, 'r') as fin, open(maf_out_path, 'w') as fout:
        # Write header lines starting with ##
        while True:
            pos = fin.tell()
            line = fin.readline()
            if not line:
                break
            if line.startswith('##'):
                fout.write(line)
            else:
                fin.seek(pos)
                break

        block_lines = []
        for line in fin:
            if line.strip() == '':
                # End of a block
                if block_lines:
                    header = block_lines[0].strip().split()
                    if len(header) >= 6:
                        block_chrom = header[0]
                        try:
                            block_start = int(header[1])
                            block_length = int(header[2])
                        except ValueError:
                            block_start = None
                            block_length = None

                        if block_start is not None and block_chrom == chrom:
                            if block_overlaps(start, end, block_start, block_length):
                                fout.writelines(block_lines)
                                fout.write('\n')
                    # Reset for next block
                    block_lines = []
                else:
                    # Empty block, just continue
                    continue
            else:
                block_lines.append(line)

        # Check last block if file doesn't end with newline
        if block_lines:
            header = block_lines[0].strip().split()
            if len(header) >= 6:
                block_chrom = header[0]
                try:
                    block_start = int(header[1])
                    block_length = int(header[2])
                except ValueError:
                    block_start = None
                    block_length = None
                if block_start is not None and block_chrom == chrom:
                    if block_overlaps(start, end, block_start, block_length):
                        fout.writelines(block_lines)
                        fout.write('\n')

# Example usage:
maf_input = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chr1.maf"
maf_output = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/test/chr1.maf"
chromosome = "hg38.chr1"  # use the chromosome format as in your MAF header lines
region_start = 5152456
region_end = 7798441

filter_maf_blocks(maf_input, maf_output, chromosome, region_start, region_end)
print(f"Filtered MAF written to: {maf_output}")

# %%
def block_overlaps(query_start, query_end, block_start, block_length):
    block_end = block_start + block_length
    return not (block_end < query_start or block_start > query_end)

def filter_maf_blocks(maf_in_path, maf_out_path, chrom, start, end, verbose=False):
    with open(maf_in_path, 'r') as fin, open(maf_out_path, 'w') as fout:
        # Write header lines starting with ##
        while True:
            pos = fin.tell()
            line = fin.readline()
            if not line:
                break
            if line.startswith('##'):
                fout.write(line)
            else:
                fin.seek(pos)
                break

        block_lines = []
        block_num = 0
        for line in fin:
            print(line)
            if line.strip() == '':
                # End of a block
                if block_lines:
                    block_num += 1
                    header = block_lines[0].strip().split()
        
                    if len(header) >= 6:
                        block_chrom = header[0]
                        try:
                            block_start = int(header[1])
                            block_length = int(header[2])
                        except ValueError:
                            block_start = None
                            block_length = None

                        keep = False
                        if block_start is not None and block_chrom == chrom:
                            if block_overlaps(start, end, block_start, block_length):
                                keep = True
                                fout.writelines(block_lines)
                                fout.write('\n')
                        if verbose:
                            status = "KEPT" if keep else "SKIPPED"
                            print(f"Block {block_num}: {block_chrom}:{block_start}-{block_start+block_length} {status}")
                    else:
                        if verbose:
                            print(f"Block {block_num}: Malformed header, SKIPPED")
                    block_lines = []
                else:
                    # Empty block, continue
                    continue
            else:
                block_lines.append(line)

        # Check last block if file doesn't end with newline
        if block_lines:
            block_num += 1
            header = block_lines[0].strip().split()
            if len(header) >= 6:
                block_chrom = header[0]
                try:
                    block_start = int(header[1])
                    block_length = int(header[2])
                except ValueError:
                    block_start = None
                    block_length = None

                keep = False
                if block_start is not None and block_chrom == chrom:
                    if block_overlaps(start, end, block_start, block_length):
                        keep = True
                        fout.writelines(block_lines)
                        fout.write('\n')
                if verbose:
                    status = "KEPT" if keep else "SKIPPED"
                    print(f"Block {block_num}: {block_chrom}:{block_start}-{block_start+block_length} {status}")
            else:
                if verbose:
                    print(f"Block {block_num}: Malformed header, SKIPPED")

# Example usage:
maf_input = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chr1.maf"
maf_output = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/test/chr1.maf"
chromosome = "hg38.chr1"
region_start = 5152456
region_end = 7798441

filter_maf_blocks(maf_input, maf_output, chromosome, region_start, region_end, verbose=True)
print(f"Filtered MAF written to: {maf_output}")

# %%
