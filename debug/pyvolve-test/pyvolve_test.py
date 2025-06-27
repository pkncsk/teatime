#%%
import pyvolve
from io import StringIO
from Bio import Phylo
#%%
# Settings
divergence_times_mya = list(range(10, 110, 10))  # 10, 20, ..., 100 MYA
substitution_rate = 2e-9  # 2e-9 subs/site/year

# Convert MYA to substitutions/site (total distance from root to tip)
subs_per_site = [t * 1e6 * substitution_rate for t in divergence_times_mya]

# Build Newick tree manually: a ladder-like topology where each taxon branches off at increasing times
# Tip names: sp1 to sp11
tips = [f"sp{i+1}" for i in range(len(subs_per_site) + 1)]

# Start from the most ancient split
tree_string = f"{tips[-1]}:{subs_per_site[-1]}"
for i in reversed(range(len(subs_per_site))):
    left = f"{tips[i]}:{subs_per_site[i] - (subs_per_site[i-1] if i > 0 else 0)}"
    right = f"({left},{tree_string}):{(subs_per_site[i-1] if i > 0 else 0)}"
    tree_string = right

# Finalize Newick string
newick_tree = f"{tree_string};"

# Parse and show the tree
tree = Phylo.read(StringIO(newick_tree), "newick")

# Display the tree as ASCII
Phylo.draw_ascii(tree)

# Output the Newick string
newick_tree
#%%
phylogeny = pyvolve.read_tree(tree = newick_tree)
#%%
pyvolve.print_tree(phylogeny)
# %%
# Define a default nucleotide model
my_model = pyvolve.Model("nucleotide")

# Define a Partition object with a specified ancestral sequence
my_partition = pyvolve.Partition(models = my_model, root_sequence = "tgttgtgggaagtcagggacccctgaatggagggacccgctgaagccgcagcacaggaacataaattgtgaagatttcatggacatttatcagttcccaaataatacttttataatttcttatgcttgtctttactttaatctcttaatcctgttatcatcctaagctgaggatgtacgtcacctcaggaccactgtgataattgtgttaactgtacaaactgattgtaaaacatgtgtgtttgaacaatatgaaatcagtgcaccttgaaaaagaacagaataacagcattttttagggagcaagggaagacaaccataagatcttactacctgtggggtcgggcaaaaagagccatatttttcttcttgcagagagcctataaatggacgtgcaagtaggagagatatccctaaattcttttcctagcaatgaataataaaatattaataccctgggaaaggaatgtgttcctcgggggaggtctataaatggccactctgggaatgtctgtcttatgcagttgagataaggactgatacaccctggtctcctgcagtaccctcaggcttattagagtgggggaaatctctgccctggtaaatttgtggtcagaccagttggttctctgctctcaaaccctgttttctgttgtttaagatgtttatcaagacaatgcgtgcaccgctgaacatagacccttatcagtagttctgcttttgccctttgccttgtgatctttgctggacccttatcagtagttctgcttttgccctttgccttgtgatctttgctggatccttatcagtagttctgctttttgccatttgaagcatgtgatctttgtacctactccctgtgcttacatcccctccccttttcaaacccttaataaaaacgtgctggtttgaggctcaggtgggcatcatggtcctaccaatatgtgatgtcacccccagcggcccagctgtaaaattcctctctttttaccctctctctttatttctcaactggccaacacttatggaaaatagaaagaacctacattgaaatattggggatgggttcccccaata".upper())
#%%
# Define an Evolver instance to evolve a single partition
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
# %%
my_evolver(seqfile = 'sim_te.fasta')
# %%
import pandas as pd
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
#%%
repeatmasker_table
# %%
