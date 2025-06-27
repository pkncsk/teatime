#%%
from Bio.AlignIO import MafIO
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.SeqRecord import SeqRecord
import os

#%%
def get_spliced_mod_unalignmask(self, starts, ends, strand=1):
    # Dictionary for IUPAC ambiguity codes for 2-base combinations
    iupac_code = {
        frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S', frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }

    def convert_to_iupac(sequence):
        unique_bases = frozenset(sequence)
        if len(unique_bases) == 1:
            return sequence[0].upper()  
        return iupac_code.get(unique_bases, 'N')  # Default to 'N' for any unhandled cases
    
    def process_sequence_localized(sequence):
        sequence = sequence.upper()
        filtered_sequence = [base for base in sequence if base != '-']

        #base_counts = Counter(sequence)
        #most_common_bases = base_counts.most_common()
        #max_count = most_common_bases[0][1]
        #consensus_bases = [base for base, count in most_common_bases if count == max_count]
        new_base = convert_to_iupac(filtered_sequence)
        return new_base

    if strand not in (1, -1): 
        raise ValueError("Strand must be 1 or -1, got %s" % strand)
    fetched = list(self.search(starts, ends))
    #return fetched
    expected_letters = sum(end - start for start, end in zip(starts, ends))
    if len(fetched) == 0:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})
    if len(fetched) == 1:
        multiseq = fetched[0]
        species = set(seqrec.id for seqrec in multiseq)
        if len(species) == 1 and self._target_seqname in species:
            return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})
    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    split_by_position
    total_rec_length = 0
    ref_first_strand = None
    for multiseq in fetched:
        if len(multiseq) == 1 and multiseq[0].id == self._target_seqname:
            # Fast-path for human-only alignment block
            rec = multiseq[0]
            rec_start = rec.annotations["start"]
            real_pos = rec_start
            for base in rec.seq:
                if base != "-":
                    split_by_position[rec.id][real_pos] = 'N'
                    real_pos += 1
            total_rec_length += rec.annotations["size"]
            continue  # skip rest of loop for this block

        for seqrec in multiseq:
            if seqrec.id == self._target_seqname:
                try:
                    if ref_first_strand is None:
                        ref_first_strand = seqrec.annotations["strand"]

                        if ref_first_strand not in (1, -1):
                            raise ValueError("Strand must be 1 or -1")
                    elif ref_first_strand != seqrec.annotations["strand"]:
                        raise ValueError(
                            "Encountered strand='%s' on target seqname, "
                            "expected '%s'"
                            % (seqrec.annotations["strand"], ref_first_strand)
                        )
                except KeyError:
                    raise ValueError(
                        "No strand information for target seqname (%s)"
                        % self._target_seqname
                    ) from None

                rec_length = len(seqrec)
                rec_start = seqrec.annotations["start"]
                ungapped_length = seqrec.annotations["size"]
                rec_end = rec_start + ungapped_length - 1
                total_rec_length += ungapped_length
                
                for seqrec in multiseq:
                    for pos in range(rec_start, rec_end + 1):
                        split_by_position[seqrec.id][pos] = ""

                break 
            else:
                raise ValueError(
                    "Did not find %s in alignment bundle" % (self._target_seqname,)
                )
        real_pos = rec_start
        edit_id = []
        edit_pos = []
        for gapped_pos in range(rec_length):
            previous_id = ''
            for seqrec in multiseq:
                
                if seqrec.id == self._target_seqname:
                    track_val = seqrec.seq[gapped_pos]
                
                
                split_by_position[seqrec.id][real_pos] += seqrec.seq[gapped_pos]
                if previous_id == seqrec.id:
                        edit_id.append(seqrec.id)
                        edit_pos.append(real_pos)
                previous_id = seqrec.id
            if track_val != "-" and real_pos < rec_end:
                real_pos += 1
        # Debugging: Print lengths of sequences in split_by_position
        for i in range(len(edit_id)):
            _sequence=split_by_position[edit_id[i]][edit_pos[i]]
            new_sequence=process_sequence_localized(_sequence)
            split_by_position[edit_id[i]][edit_pos[i]] = new_sequence
        
        if len(split_by_position[self._target_seqname]) != total_rec_length:
            raise ValueError(
                "Target seqname (%s) has %s records, expected %s"
                % (
                    self._target_seqname,
                    len(split_by_position[self._target_seqname]),
                    total_rec_length,
                )
            )

    realpos_to_len = {
        pos: len(gapped_fragment)
        for pos, gapped_fragment in split_by_position[self._target_seqname].items()
        if len(gapped_fragment) > 1
    }

    seqid_list = []
    seq_list = []
    
    for seqid in all_seqnames:
        seq_split = split_by_position[seqid]
        seq_splice = []
        filler_char = "N" if seqid == self._target_seqname else "-"
        append = seq_splice.append
        for exonstart, exonend in zip(starts, ends):
            for real_pos in range(exonstart, exonend):
                if real_pos in seq_split:
                    append(seq_split[real_pos])
                elif real_pos in realpos_to_len:
                    append(filler_char * realpos_to_len[real_pos])
                else:
                    append(filler_char)

        seqid_list.append(seqid)
        seq_list.append(Seq("".join(seq_splice))) 
    target_index = seqid_list.index(self._target_seqname)
    if len(seq_list[target_index].replace("-", "")) != expected_letters:
        raise ValueError(
            "Returning %s letters for target seqname (%s), expected %s"
            % (
                len(seq_list[target_index].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(seq_list[target_index])
    for seqid, seq in zip(seqid_list, seq_list):
        if len(seq) != ref_subseq_len:
            raise ValueError(
                "Returning length %s for %s, expected %s"
                % (len(seq), seqid, ref_subseq_len)
            )

    # Create a DataFrame
    df = pd.DataFrame({
        'seqid': seqid_list,
        'seq': [seq.reverse_complement() if strand != ref_first_strand else seq for seq in seq_list]
    })
    return df

MafIO.MafIndex.get_spliced = get_spliced_mod_unalignmask
#%%
# Initialize
strand = 1
chrom = 'chr1'
maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
#%%
chromosome_ranges = {
    "chr1": {"start_alignment": 9000, "end_alignment": 248947000},
    "chr2": {"start_alignment": 9000, "end_alignment": 242184000},
    "chr3": {"start_alignment": 1000, "end_alignment": 198236000},
    "chr4": {"start_alignment": 10000, "end_alignment": 190205000},
    "chr5": {"start_alignment": 13000, "end_alignment": 181479000},
    "chr6": {"start_alignment": 60000, "end_alignment": 170746000},
    "chr7": {"start_alignment": 7000, "end_alignment": 159336000},
    "chr8": {"start_alignment": 59000, "end_alignment": 145087000},
    "chr9": {"start_alignment": 10000, "end_alignment": 138335000},
    "chr10": {"start_alignment": 10000, "end_alignment": 133788000},
    "chr11": {"start_alignment": 60000, "end_alignment": 135086000},
    "chr12": {"start_alignment": 9000, "end_alignment": 133265000},
    "chr13": {"start_alignment": 18171000, "end_alignment": 114355000},
    "chr14": {"start_alignment": 16019000, "end_alignment": 106884000},
    "chr15": {"start_alignment": 17723000, "end_alignment": 101982000},
    "chr16": {"start_alignment": 9000, "end_alignment": 90229000},
    "chr17": {"start_alignment": 59000, "end_alignment": 83252000},
    "chr18": {"start_alignment": 10000, "end_alignment": 80271000},
    "chr19": {"start_alignment": 60000, "end_alignment": 58608000},
    "chr20": {"start_alignment": 60000, "end_alignment": 64335000},
    "chr21": {"start_alignment": 5005000, "end_alignment": 46700000},
    "chr22": {"start_alignment": 10510000, "end_alignment": 50817000},
    "chrX": {"start_alignment": 9000, "end_alignment": 156032000},
    "chrY": {"start_alignment": 2781000, "end_alignment": 56888000},
}

start_alignment = chromosome_ranges[chrom]["start_alignment"]
end_alignment = chromosome_ranges[chrom]["end_alignment"]
step_size = 1000  # Process in 1k blocks
batch_size = 100 
max_workers = 8
#%%
chunk_ranges = [(pos, min(pos + step_size, end_alignment)) for pos in range(start_alignment, end_alignment, step_size)]
#%%
# Directory for temp results
TEMP_DIR = f"/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/aligned_count_eff_update/{chrom}_temp"
os.makedirs(TEMP_DIR, exist_ok=True)

# Function to count ungapped bases
def count_ungapped(seq, mask):
    return sum(1 for base, valid in zip(seq, mask) if valid and base not in {'-', 'N', 'n'})

# Worker function for parallel processing
def process_chunk(start_pos,end_pos):
    
    temp_filepath = os.path.join(TEMP_DIR, f"chunk_{start_pos}_{end_pos}.pkl")

    # Skip if temp file exists
    if os.path.exists(temp_filepath):
        #print(f"Skipping {start_pos}-{end_pos}, already processed.", flush=True)
        return None  

    #print(f"Processing {start_pos}-{end_pos}", flush=True)

    # Fetch MAF data
    index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
    spliced_maf_full = index_maf.get_spliced([start_pos], [end_pos], strand)
    
    if not isinstance(spliced_maf_full, pd.DataFrame):
        #print(f"Empty or invalid chunk: {start_pos}-{end_pos}", flush=True)
        return None

    hg38_row = spliced_maf_full[spliced_maf_full['seqid'].str.startswith("hg38")]
    if hg38_row.empty:
        #print(f"No hg38 data at {start_pos}-{end_pos}", flush=True)
        return None
    
    hg38_seq = list(hg38_row.iloc[0]['seq'])  # Extract hg38 sequence
    hg38_mask = [base not in {'-', 'N', 'n'} for base in hg38_seq]  # Masking
    # Count ungapped bases for all species
    species_counts = {}
    for _, row in spliced_maf_full.iterrows():
        species = row['seqid'].split('.')[0]  # Extract species name
        ungapped_count = count_ungapped(row['seq'], hg38_mask)
        species_counts[species] = species_counts.get(species, 0) + ungapped_count #more than one chromosome can align to hg38

    # Save result as pickle
    with open(temp_filepath, "wb") as f:
        pickle.dump(species_counts, f)
        #print(f'DONE: {start_pos}-{end_pos}', flush=True)

    return species_counts
#%%
from tqdm import tqdm
import gc
# Process in batches
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    with tqdm(total=len(chunk_ranges), desc="Processing", unit="chunk") as pbar:
        for i in range(0, len(chunk_ranges), batch_size):
            batch = chunk_ranges[i : i + batch_size]  # Get a batch of tasks
            
            # Submit batch to the executor
            future_to_range = {executor.submit(process_chunk, start, end): (start, end) for start, end in batch}
            
            # Wait for the batch to finish before queuing more
            for future in as_completed(future_to_range):
                try:
                    _ = future.result()  # Ensures failed tasks raise errors
                    #print(_)
                except Exception as e:
                    print(f"Task failed: {e}", flush=True)    # Get the result (if any)
                del future_to_range[future]
                #print(result)  # Print the completed task
                pbar.update(1)

            if i % (batch_size * 5) == 0:  # Every 5 batches, force cleanup
                gc.collect()
#%%
#%%
import os
import pickle
import sys


TEMP_DIR = f"/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/aligned_count_eff_update/{chrom}_temp"
OUTPUT_FILE = f"/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/aligned_count_eff_update/{chrom}_merged.pkl"

# Generate expected filenames
expected_files = {f"chunk_{pos}_{min(pos + step_size, end_alignment)}.pkl" 
                  for pos in range(start_alignment, end_alignment, step_size)}

# Check if all expected files exist
existing_files = set(os.listdir(TEMP_DIR))
missing_files = expected_files - existing_files
#%%
if missing_files:
    print(f"Error: Missing files detected! {len(missing_files)} chunks missing.", file=sys.stderr)
    print("\n".join(missing_files), file=sys.stderr)
    sys.exit(1)
#%%
import re
# Extract numerical start positions and sort filenames
def extract_start_pos(filename):
    match = re.match(r"chunk_(\d+)_\d+\.pkl", filename)
    return int(match.group(1)) if match else float('inf')

sorted_files = sorted(existing_files, key=extract_start_pos)
#%%
# Aggregate counts
total_counts = {}
for idx,file in enumerate(sorted_files):
    file_path = os.path.join(TEMP_DIR, file)
    print(f'{idx}/{len(sorted_files)}:{file_path}')
    try:
        with open(file_path, "rb") as f:
            chunk_counts = pickle.load(f)
            for species, count in chunk_counts.items():
                total_counts[species] = total_counts.get(species, 0) + count
    except (EOFError, pickle.UnpicklingError) as e:
        print(f"Error: Corrupt file detected -> {file}. Terminating.", file=sys.stderr)
        sys.exit(1)

# Save the merged result
try:
    with open(OUTPUT_FILE, "wb") as f:
        pickle.dump(total_counts, f)
    print(f"Merged result saved to {OUTPUT_FILE}")

    # Delete temp files after successful save
    for file in expected_files:
        os.remove(os.path.join(TEMP_DIR, file))
    os.rmdir(TEMP_DIR)
    print("Temporary files deleted.")

except Exception as e:
    print(f"Error saving merged result: {e}", file=sys.stderr)
    sys.exit(1)

# %%
