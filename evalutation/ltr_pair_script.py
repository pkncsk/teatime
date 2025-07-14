import os
import pandas as pd
from Bio import SeqIO
from ma_mapper import sequence_alignment

def load_filtered_blocks_table(filepath):
    """
    Load the filtered block list (contains block_id, internal_id, coords, etc).
    """
    df = pd.read_csv(filepath, sep='\t')
    print(f"Loaded {len(df)} rows from {filepath}")
    return df

def batch_extract_all_seq_records(df_filtered, genome_fasta):
    """
    Run sequence_io once to extract all LTR sequences as SeqRecords.
    """
    # Sort so aINT/bINT pairs are grouped together
    df_sorted = df_filtered.sort_values(['block_id', 'genoName', 'genoStart'])

    block_coords = df_sorted[['genoName', 'genoStart', 'genoEnd', 'internal_id', 'strand']].copy()
    block_coords['score'] = 10
    block_coords = block_coords[['genoName', 'genoStart', 'genoEnd', 'internal_id', 'score', 'strand']]
    
    # Run sequence_io once for all coordinates
    seq_records = sequence_alignment.sequence_io(block_coords, source_fasta=genome_fasta, save_to_file=False)
    print(f"Extracted {len(seq_records)} sequence records.")
    return seq_records

def write_seq_records_per_block(seq_records, output_dir="./tmp"):
    """
    Write a separate FASTA file per block_id using the SeqRecords.
    """
    os.makedirs(output_dir, exist_ok=True)
    record_dict = {}

    for rec in seq_records:
        parts = rec.id.split("_", 2)
        block_id = parts[0] + "_" + parts[1]  # e.g., THE1C_4250279
        record_dict.setdefault(block_id, []).append(rec)

    for block_id, records in record_dict.items():
        fasta_path = os.path.join(output_dir, f"{block_id}.fasta")
        if not os.path.isfile(fasta_path):
            with open(fasta_path, "w") as f:
                SeqIO.write(records, f, "fasta")
        else:
            print(f"FASTA already exists: {fasta_path}")

    return list(record_dict.keys())  # return list of processed block_ids

def run_mafft_alignment(fasta_path, aligned_path, nthread=6):
    """
    Run MAFFT alignment if aligned file is not present.
    """
    if not os.path.isfile(aligned_path):
        sequence_alignment.mafft_align(fasta_path, output_filepath=aligned_path, nthread=nthread)
    else:
        print(f"Alignment already exists: {aligned_path}")

def run_batch_alignment(block_ids, tmp_dir="./tmp", nthread=6):
    """
    Align all block FASTA files in tmp_dir using MAFFT.
    """
    for block_id in block_ids:
        print(f'align: {block_id}')
        fasta_path = os.path.join(tmp_dir, f"{block_id}.fasta")
        aligned_path = f"{fasta_path}.aligned"
        run_mafft_alignment(fasta_path, aligned_path, nthread=nthread)

# ------------------ MAIN ------------------

if __name__ == "__main__":
    # Required input
    filtered_blocks_path = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/final_LTR_sample_3000.txt"
    genome_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa'
    tmp_dir = "./tmp"
    nthread = 6

    # Step 1: Load filtered block data
    df_filtered = load_filtered_blocks_table(filtered_blocks_path)

    # Step 2: Extract all sequences in one go
    seq_records = batch_extract_all_seq_records(df_filtered, genome_fasta)

    # Step 3: Write per-block FASTA files
    block_ids = write_seq_records_per_block(seq_records, output_dir=tmp_dir)

    # Step 4: Run MAFFT alignments
    run_batch_alignment(block_ids, tmp_dir=tmp_dir, nthread=nthread)

    print("✅ All alignments done.")

