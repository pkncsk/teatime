#%%
import pyvolve
from io import StringIO
from Bio import Phylo
#%%
newick_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/phylogenetic_tree/hg38.447way.nh.txt'
phylotree=Phylo.read(newick_path, "newick")
#%%
from ete3 import Tree

# Load the Newick tree with named internal nodes
ete_tree = Tree(newick_path, format=1)
ete_tree.write()
#%%
for i, n in enumerate(ete_tree.traverse("postorder")):
    if not n.is_leaf() and not n.name:
        n.name = f"ac{i}"
#%%
# Define the internal node name where the TE is inserted
insertion_root_name = "ac208"  # example
insertion_root = ete_tree.search_nodes(name=insertion_root_name)[0]
# %%
# Copy the full tree first
tree1 = ete_tree.copy()
#%%
for child in insertion_root.get_children():
    print(child.name)
    child_ = tree1.search_nodes(name=child.name)[0]
    removed_node=child_.detach()
#%%
tree1_root=tree1.get_tree_root()
# %%
# Newick string of tree1 (pre-insertion evolution)
print("Tree1 (pre-insertion):")
print(tree1.write(format=1, format_root_node=True))  # format=1 includes branch lengths
print(tree1.get_ascii(show_internal=True))
#%%
# If you want to double-check the full tree as well
print("\nTree2 (post-insertion):")
print(insertion_root.write(format=1))
print(insertion_root.get_ascii(show_internal=True))

# %%
import random

def generate_random_sequence(length=10000, seed=42):
    random.seed(seed)
    return ''.join(random.choices("ACGT", k=length))

random_seq = generate_random_sequence()
print(f"Length: {len(random_seq)}\nFirst 100 bases: {random_seq[:100]}")
#%%
# setup tree
newick_tree=tree1.write(format=1, format_root_node=True)
#newick_tree=insertion_root.write(format=1)
#%%
phylogeny = pyvolve.read_tree(tree = newick_tree)
# Define a basic model (e.g., HKY85)
my_model = pyvolve.Model("nucleotide")
#%% extract root seq from hg38

#%%
#non overlap region
#chr1:2,680,000-2,690,000
import pandas as pd

# Define BED fields
data = {
    "chrom": ["chr1"],
    "chromStart": [2679999],
    "chromEnd": [2690000],
    "name": ["region1"],       # Optional, can be customized
    "score": [0],              # Optional, default to 0
    "strand": ["+"]            # Optional, assume '+' strand
}

# Create DataFrame
bed_df = pd.DataFrame(data)

# Display DataFrame
print(bed_df)

#%%
from ma_mapper import sequence_alignment
sequence_alignment.sequence_io(
    coordinate_table=bed_df,
    source_fasta='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa',
    save_to_file='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.root.fa')
#%%
root_seqeunces=sequence_alignment.sequence_io(
    coordinate_table=bed_df,
    source_fasta='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa',
    save_to_file=False)
#%%
rootseq=str(root_seqeunces[0].seq)
#%%
# Partition with root sequence
my_partition = pyvolve.Partition(models=my_model, root_sequence=rootseq)
#%%
# Then evolve using tree1
evolver = pyvolve.Evolver(partition=my_partition, tree=phylogeny)
#%%
evolver(seqfile="/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.preinsertion.fasta", ratefile=None, infofile=None)
# %%
import pandas as pd

#TE
#chr1:7790171-7792441
chrom = 'chr1'
start = 7790171
end = 7792441
strand = '+'
# Define BED fields
data = {
    "chrom": [chrom],
    "chromStart": [start],
    "chromEnd": [end],
    "name": ["simTHE1C_COMPLETE_0"],       # Optional, can be customized
    "score": [10],              # Optional, default to 0
    "strand": [strand]            # Optional, assume '+' strand
}

# Create DataFrame
bed_df = pd.DataFrame(data)

from ma_mapper import sequence_alignment
sequence_alignment.sequence_io(
    coordinate_table=bed_df,
    source_fasta='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa',
    save_to_file='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.fa')
#%%
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Load ref table
ref_table_path =  "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.ref.txt"
ref_df = pd.read_csv(ref_table_path, sep="\t", index_col=0)

# Load genome sequences (Pyvolve output)
genome_fasta = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.preinsertion.fasta"
genome_seqs = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# Load TE sequences
te_fasta = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.fa"
te_seqs = {record.id.split("::")[0]: str(record.seq).upper() for record in SeqIO.parse(te_fasta, "fasta")}
#%%

target_genome_seq = str(genome_seqs[insertion_root_name].seq)
#%%
# Process insertions
for idx, row in ref_df.iterrows():
    internal_id = row["internal_id"]
    insert_pos = row["genoStart"]
    strand = row["strand"]
    te_seq = te_seqs[internal_id]
    if strand == "-":
            te_seq = str(Seq(te_seq).reverse_complement())
    
    target_genome_seq = target_genome_seq[:insert_pos] + te_seq + target_genome_seq[insert_pos:]
    print(f"Inserted {internal_id} to simchr at {insert_pos} on strand {strand}")

    # Write the modified root sequence
    out_record = SeqRecord(Seq(target_genome_seq), id=insertion_root_name, description="with_TE_insertions")
    SeqIO.write(out_record, f"/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/{insertion_root_name}_with_te.fasta", "fasta")
    print(f"✅ Written {insertion_root_name}_with_te.fasta")
#%%
# setup tree
newick_tree=insertion_root.write(format=1, format_root_node=True)
my_tree = Phylo.read(StringIO(newick_tree), "newick")
Phylo.draw_ascii(my_tree)
phylogeny = pyvolve.read_tree(tree = newick_tree)
# Define a basic model (e.g., HKY85)
my_model = pyvolve.Model("nucleotide")
#%%
# Partition with root sequence
my_partition = pyvolve.Partition(models=my_model, root_sequence=target_genome_seq)
#%%
# Then evolve using tree1
evolver = pyvolve.Evolver(partition=my_partition, tree=phylogeny)
#%%
evolver(seqfile="/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.postinsertion.fasta", ratefile=None, infofile=None)
# %%
from Bio import SeqIO

# Input paths
pre_insertion_fasta = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.preinsertion.fasta"
post_insertion_fasta = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simTHE1C.postinsertion.fasta"
output_fasta = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simAlign.fasta"
records = []
# Read all sequences from tree1_output.fasta
records.extend(record for record in SeqIO.parse(pre_insertion_fasta, "fasta") if record.id != insertion_root_name)

# Read the modified ac3 sequence
records.extend(record for record in SeqIO.parse(post_insertion_fasta, "fasta") if record.id != insertion_root_name)

# Write merged FASTA
SeqIO.write(records, output_fasta, "fasta")
print(f"✅ Merged FASTA written to: {output_fasta}")

# %%
