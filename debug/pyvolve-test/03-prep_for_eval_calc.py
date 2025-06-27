#%%
from ma_mapper import sequence_alignment
#%%
sequence_alignment.mafft_align(
    input_filepath='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simAlign.fasta',
    nthread=6,
    output_filepath='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simAlign.aligned'
)
# %%
from ma_mapper import mapper
alignment_matrix  = mapper.parse_alignment(alignment_file='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simAlign.aligned')

# %%
from ma_mapper import plots
plots.plot(alignment=alignment_matrix,
    show_alignment=True,
    alignment_col = 'dna',
    heatmap_color=["viridis"], 
    vlim =[[-0.5,0.5]],
    )
# %%
import hashlib
import random
import string
from Bio import AlignIO
from Bio.AlignIO.MafIO import MafWriter
from Bio.Align import MultipleSeqAlignment
def count_ungapped(seq):
    return sum(1 for base in seq if base not in {'-', 'N', 'n'})
def generate_seed_from_sequence(seq):
    # Convert the sequence to a string (just in case it’s not already)
    seq_str = str(seq)
    
    # Use hashlib to create a hash of the sequence string
    hash_value = hashlib.md5(seq_str.encode()).hexdigest()
    
    # Convert the hash to an integer and use it as the seed
    seed = int(hash_value, 16) % (2**32)  # Modulo to keep it within a valid range
    
    return seed

def random_chr(seq):
    seed = generate_seed_from_sequence(seq)
    random.seed(seed)
    return ''.join(random.choices('ABCDEFGHIJKLMNOPQRSTUVWXYZ', k=5))

def fix_alignment(fasta_in, maf_out, chr_table=None):
    # Load aligned FASTA
    aln = AlignIO.read(fasta_in, "fasta")
    
    # Optional: reorder so sp11 is last
    aln.sort(key=lambda r: r.id != "hg38")

    new_aln = MultipleSeqAlignment([])

    for rec in aln:
        name = rec.id

        # Get chromosome name from table or randomly generate
        if chr_table and name in chr_table:
            chrom = chr_table[name]
        else:
            chrom = random_chr(rec.seq)

        rec.id = f"{name}.{chrom}"

        
        rec.annotations["start"] = 0
        rec.annotations["size"] = count_ungapped(rec.seq)
        rec.annotations["strand"] = "+"
        rec.annotations["srcSize"] = 1_000_000  # Mock size, replace if needed

        new_aln.append(rec)

    with open(maf_out, "w") as f:
        writer = MafWriter(f)
        writer.write_header()
        writer.write_alignment(new_aln)

# Example usage
chr_table = {
    "hg38": "chr1"
}

fix_alignment('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simAlign.aligned', '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.maf', chr_table)

# %%
# make mafindex
import os
import pickle
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from Bio.AlignIO import MafIO

# Initialize
strand = 1
chrom = 'chr1'
maf_filepath = f'/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.maf'
mafindex_filepath = f'/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
# %%

