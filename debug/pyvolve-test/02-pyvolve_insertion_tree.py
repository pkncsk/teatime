#%%
import pyvolve
from io import StringIO
from Bio import Phylo
#%%
# Settings
work_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C'
divergence_times_mya = list(range(10, 110, 10))  # 10, 20, ..., 100 MYA
substitution_rate = 2e-9  # 2e-9 subs/site/year

# Convert MYA to substitutions/site (total distance from root to tip)
subs_per_site = [t * 1e6 * substitution_rate for t in divergence_times_mya]

# Build Newick tree manually: a ladder-like topology where each taxon branches off at increasing times
# Tip names: sp1 to sp11
tips = [f"sp{i+1}" for i in range(len(subs_per_site) + 1)]

# Start from the most ancient split
tree_string = f"{tips[-1]}:{subs_per_site[-1]}"
ac_counter=0
for i in reversed(range(len(subs_per_site))):
    branch_len = subs_per_site[i] - (subs_per_site[i-1] if i > 0 else 0)
    left = f"{tips[i]}:{branch_len}"
    ac_name = f"ac{ac_counter}"
    right = f"({left},{tree_string}){ac_name}:{(subs_per_site[i-1] if i > 0 else 0)}"
    tree_string = right
    ac_counter += 1

newick_tree = f"{tree_string};"

# Parse and show the tree
tree = Phylo.read(StringIO(newick_tree), "newick")

# Display the tree as ASCII
Phylo.draw_ascii(tree)

# Output the Newick string
newick_tree
#%%
nw447 = f'/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/phylogenetic_tree/hg38.447way.nh.txt'
tree447=Phylo.read(nw447, "newick")
#%%
from ete3 import Tree

# Load the Newick tree with named internal nodes
ete_tree = Tree(nw447, format=1)
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
#%%
# Partition with root sequence
my_partition = pyvolve.Partition(models=my_model, root_sequence=random_seq)
#%%
# Then evolve using tree1
evolver = pyvolve.Evolver(partition=my_partition, tree=phylogeny)
#%%
tree1_filepath = f"{work_dir}/simTHE1C.tree1_output.fasta"
evolver(seqfile=tree1_filepath, ratefile=None, infofile=None)
# %%
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Parameters
target_length = 357
output_fasta = "simTHE1C.fasta"
fasta_path = f"/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta"
ref_table_path = f"{work_dir}/simTHE1C.ref.txt"

# Load ref_table
ref_table = pd.read_csv(ref_table_path, sep="\t", index_col=0)

# Create an iterator over the matching ref_table rows
ref_iter = ref_table.iterrows()

# To collect final records
records_to_write = []

# Iterate through the original FASTA file
for record in SeqIO.parse(fasta_path, "fasta"):
    if len(record.seq) == target_length:
        try:
            _, row = next(ref_iter)
        except StopIteration:
            break  # We've used up all the matching ref entries

        internal_id = row["internal_id"]
        genoName = row["genoName"]
        genoStart = row["genoStart"]
        genoEnd = row["genoEnd"]
        strand = row["strand"]
        
        new_header = f"{internal_id}::{genoName}:{genoStart}-{genoEnd}({strand})"
        print(f"✔ Matched: {record.id} → {new_header}")
        
        new_record = SeqRecord(Seq(str(record.seq).upper()), id=new_header, description="")
        records_to_write.append(new_record)

# Write the result
SeqIO.write(records_to_write, output_fasta, "fasta")
print(f"✅ Wrote {len(records_to_write)} sequences to {output_fasta}")
#%%
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Load ref table
ref_table_path =  f"{work_dir}/simTHE1C.ref.txt"
ref_df = pd.read_csv(ref_table_path, sep="\t", index_col=0)

# Load genome sequences (Pyvolve output)
genome_fasta = f"{work_dir}/simTHE1C.tree1_output.fasta"
genome_seqs = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# Load TE sequences
te_fasta = f"{work_dir}/simTHE1C.fasta"
te_seqs = {record.id.split("::")[0]: str(record.seq) for record in SeqIO.parse(te_fasta, "fasta")}
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
    SeqIO.write(out_record, f"{work_dir}/{insertion_root_name}_with_te.fasta", "fasta")
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
evolver(seqfile=f"{work_dir}/simTHE1C.tree2_output.fasta", ratefile=None, infofile=None)
# %%
from Bio import SeqIO

# Input paths
pre_insertion_fasta = f"{work_dir}/simTHE1C.tree1_output.fasta"
post_insertion_fasta = f"{work_dir}/simTHE1C.tree2_output.fasta"
output_fasta = f"{work_dir}/simTHE1C.simTEinsertion.fasta"
records = []
# Read all sequences from tree1_output.fasta
records.extend(record for record in SeqIO.parse(pre_insertion_fasta, "fasta") if record.id != "ac208")

# Read the modified ac3 sequence
records.extend(record for record in SeqIO.parse(post_insertion_fasta, "fasta") if record.id != "ac208")

# Write merged FASTA
SeqIO.write(records, output_fasta, "fasta")
print(f"✅ Merged FASTA written to: {output_fasta}")

# %%
