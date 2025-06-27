#%%
import pyvolve
from io import StringIO
from Bio import Phylo
from ete3 import Tree
import random
from ma_mapper import sequence_alignment
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def generate_random_sequence(seed,length=10000):
    random.seed(int(seed))
    return ''.join(random.choices("ACGT", k=length))

#%%
work_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-run/simTHE1C'
newick_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/phylogenetic_tree/hg38.447way.nh.txt'
#%%
phylotree=Phylo.read(newick_path, "newick")
# Load the Newick tree with named internal nodes
ete_tree = Tree(newick_path, format=1)
for i, n in enumerate(ete_tree.traverse("postorder")):
    if not n.is_leaf() and not n.name:
        n.name = f"ac{i}"
# Define the internal node name where the TE is inserted
insertion_root_name = "ac484"  # example
insertion_root = ete_tree.search_nodes(name=insertion_root_name)[0]
# Copy the full tree first
tree1 = ete_tree.copy()
for child in insertion_root.get_children():
    child_ = tree1.search_nodes(name=child.name)[0]
    removed_node=child_.detach()
tree1_root=tree1.get_tree_root()
#%%
# Define BED fields
data = {#non overlap region #chr1:2,680,000-2,690,000
    "chrom": ["chr1"],
    "chromStart": [2679999],
    "chromEnd": [2690000],
    "name": ["region1"],       # Optional, can be customized
    "score": [0],              # Optional, default to 0
    "strand": ["+"]            # Optional, assume '+' strand
}

# Create DataFrame
bed_df = pd.DataFrame(data)
sequence_alignment.sequence_io(
    coordinate_table=bed_df,
    source_fasta='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa',
    save_to_file='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-run/simTHE1C/simTHE1C.root.fa')

root_seqeunces=sequence_alignment.sequence_io(
    coordinate_table=bed_df,
    source_fasta='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa',
    save_to_file=False)

rootseq=str(root_seqeunces[0].seq)

#Define TE
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
    save_to_file='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-run/simTHE1C/simTHE1C.fa')

# Load TE sequences
te_fasta = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-run/simTHE1C/simTHE1C.fa"
te_seqs = {record.id.split("::")[0]: str(record.seq).upper() for record in SeqIO.parse(te_fasta, "fasta")}
# %%
np.random.seed(575)
seed_list = np.random.randint(0, 2**32, size=100, dtype=np.uint32) 
for idx, seed_n in enumerate(seed_list):
    work_subfolder=f'{work_dir}/{idx}/'
    os.makedirs(work_subfolder, exist_ok=True)
    random_seq = generate_random_sequence(seed = seed_n)
    # setup tree
    newick_tree=tree1.write(format=1, format_root_node=True)
    phylogeny = pyvolve.read_tree(tree = newick_tree)
    # Define a basic model (e.g., HKY85)
    my_model = pyvolve.Model("nucleotide")
    # Partition with root sequence
    my_partition = pyvolve.Partition(models=my_model, root_sequence=rootseq)
    # Then evolve using tree1
    evolver = pyvolve.Evolver(partition=my_partition, tree=phylogeny)
    preinsertion_filepath = f"{work_subfolder}/simTHE1C_{idx}.preinsertion.fasta"
    evolver(seqfile=preinsertion_filepath, ratefile=None, infofile=None)
    # Load ref table
    ref_table_path =  f"{work_subfolder}/simTHE1C_{idx}.ref.txt"
    ref_df = pd.read_csv(ref_table_path, sep="\t", index_col=0)
    # Load genome sequences (Pyvolve output)
    genome_fasta = preinsertion_filepath
    genome_seqs = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    target_genome_seq = str(genome_seqs[insertion_root_name].seq)
    # Process insertions
    for idx2, row in ref_df.iterrows():
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

        SeqIO.write(out_record, f"{work_subfolder}/{insertion_root_name}_with_te.fasta", "fasta")
        print(f"✅ Written {insertion_root_name}_with_te.fasta")
    
    # setup tree
    newick_tree2=insertion_root.write(format=1, format_root_node=True)
    #Phylo.draw_ascii(my_tree2)
    phylogeny2 = pyvolve.read_tree(tree = newick_tree2)
    # Define a basic model (e.g., HKY85)
    # Partition with root sequence
    my_partition2 = pyvolve.Partition(models=my_model, root_sequence=target_genome_seq)
    # Then evolve using tree1
    evolver = pyvolve.Evolver(partition=my_partition2, tree=phylogeny2)
    evolver(seqfile=f"{work_subfolder}/simTHE1C_{idx}.postinsertion.fasta", ratefile=None, infofile=None)
    # Input paths
    print(idx)
    pre_insertion_fasta = f"{work_subfolder}/simTHE1C_{idx}.preinsertion.fasta"
    print(pre_insertion_fasta)
    post_insertion_fasta = f"{work_subfolder}/simTHE1C_{idx}.postinsertion.fasta"
    print(post_insertion_fasta)
    output_fasta = f"{work_subfolder}/simTHE1C_{idx}.prealign.fasta"
    records = []
    # Read all sequences from tree1_output.fasta
    records.extend(record for record in SeqIO.parse(pre_insertion_fasta, "fasta") if record.id != insertion_root_name)

    # Read the modified ac3 sequence
    records.extend(record for record in SeqIO.parse(post_insertion_fasta, "fasta") if record.id != insertion_root_name)

    # Write merged FASTA
    SeqIO.write(records, output_fasta, "fasta")
    print(f"✅ Merged FASTA written to: {output_fasta}")

# %%
