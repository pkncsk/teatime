#%%
from Bio.AlignIO import MafIO
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ma_mapper import gzmaf
#%%
def get_spliced_mod(self, starts, ends, strand=1):
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
    #exit early if no alignment
    #print(expected_letters)
    print(fetched)
    if len(fetched) == 0:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})
    #if len(fetched) == 1:
    #    return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})
    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}    
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    split_by_position
    total_rec_length = 0
    ref_first_strand = None
    for multiseq in fetched:
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

    if len(seq_list[seqid_list.index(self._target_seqname)].replace("-", "")) != expected_letters:
        raise ValueError(
            "Returning %s letters for target seqname (%s), expected %s"
            % (
                len(seq_list[seqid_list.index(self._target_seqname)].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(seq_list[seqid_list.index(self._target_seqname)])
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

MafIO.MafIndex.get_spliced = get_spliced_mod

#%%
def search_mod(self, starts, ends):
    # verify the provided exon coordinates
    if len(starts) != len(ends):
        raise ValueError("Every position in starts must have a match in ends")

    # Could it be safer to sort the (exonstart, exonend) pairs?
    for exonstart, exonend in zip(starts, ends):
        exonlen = exonend - exonstart
        if exonlen < 1:
            raise ValueError(
                "Exon coordinates (%d, %d) invalid: exon length (%d) < 1"
                % (exonstart, exonend, exonlen)
            )
    con = self._con

    yielded_rec_coords = set()
    # search for every exon
    for exonstart, exonend in zip(starts, ends):
        try:
            possible_bins = ", ".join(
                map(str, self._region2bin(exonstart, exonend))
            )
        except TypeError:
            raise TypeError(
                "Exon coordinates must be integers "
                "(start=%d, end=%d)" % (exonstart, exonend)
            ) from None
        result = con.execute(
            "SELECT DISTINCT start, end, offset FROM offset_data "
            "WHERE bin IN (%s) "
            "AND (end BETWEEN %s AND %s OR %s BETWEEN start AND end) "
            "ORDER BY start, end, offset ASC;"
            % (possible_bins, exonstart, exonend - 1, exonend - 1)
        )
        rows = result.fetchall()
        for rec_start, rec_end, offset in rows:
            if (rec_start, rec_end) in yielded_rec_coords:
                continue
            else:
                yielded_rec_coords.add((rec_start, rec_end))


            fetched = self._get_record(int(offset))

            for record in fetched:
                if record.id == self._target_seqname:
    
                    start = record.annotations["start"]
                    end = start + record.annotations["size"] - 1

                    if not (start == rec_start and end == rec_end):
                        raise ValueError(
                            "Expected %s-%s @ offset %s, found %s-%s"
                            % (rec_start, rec_end, offset, start, end)
                        )

            yield fetched
#%%
def search_debug(self, starts, ends):
    """Search index database for MAF records overlapping ranges provided."""
    if len(starts) != len(ends):
        raise ValueError("Every position in starts must have a match in ends")
    
    for exonstart, exonend in zip(starts, ends):
        exonlen = exonend - exonstart
        if exonlen < 1:
            raise ValueError(f"Exon coordinates ({exonstart}, {exonend}) invalid: exon length ({exonlen}) < 1")
        
        print(f"Processing exon from {exonstart} to {exonend}")
        con = self._con
        
        yielded_rec_coords = set()
        
        try:
            possible_bins = ", ".join(map(str, self._region2bin(exonstart, exonend)))
            print(f"Possible bins for start, end: {possible_bins}")
        except TypeError:
            raise TypeError(f"Exon coordinates must be integers (start={exonstart}, end={exonend})") from None
        
        sql_query = f"SELECT DISTINCT start, end, offset FROM offset_data WHERE bin IN ({possible_bins}) AND (end BETWEEN {exonstart} AND {exonend - 1} OR {exonend - 1} BETWEEN start AND end) ORDER BY start, end, offset ASC;"
        print(f"SQL query: {sql_query}")
        
        result = con.execute(sql_query)
        rows = result.fetchall()
        print(f"Rows retrieved: {rows}")
        
        for rec_start, rec_end, offset in rows:
            print(f"Processing record from {rec_start} to {rec_end} with offset {offset}")
            if (rec_start, rec_end) in yielded_rec_coords:
                continue
            else:
                yielded_rec_coords.add((rec_start, rec_end))
            
            fetched = self._get_record(int(offset))
            
            for record in fetched:
                if record.id == self._target_seqname:
                    start = record.annotations["start"]
                    end = start + record.annotations["size"] - 1
                    
                    print(f"Fetched record start: {start}, end: {end}")
                    
                    if not (start == rec_start and end == rec_end):
                        raise ValueError(f"Expected {rec_start}-{rec_end} @ offset {offset}, found {start}-{end}")
            
            print(f"Yielded record coordinates: {yielded_rec_coords}")
            yield fetched


MafIO.MafIndex.search = search_debug

#%%
def search_debug2(self, starts, ends):
    """Search index database for MAF records overlapping ranges provided."""
    # verify the provided exon coordinates
    if len(starts) != len(ends):
        raise ValueError("Every position in starts must have a match in ends")

    # search for every exon
    for exonstart, exonend in zip(starts, ends):
        exonlen = exonend - exonstart
        if exonlen < 1:
            raise ValueError(
                "Exon coordinates (%d, %d) invalid: exon length (%d) < 1"
                % (exonstart, exonend, exonlen)
            )
        con = self._con

        # Construct the SQL query without using binning
        sql_query = f"""
            SELECT DISTINCT start, end, offset FROM offset_data
WHERE end BETWEEN 28805000 AND 28806000 OR 28806000 BETWEEN start AND end
ORDER BY start, end, offset ASC;
        """

        # Print the query for debugging
        print(f"SQL query: {sql_query}")

        # Execute the query
        result = con.execute(sql_query)

        # Fetch and print the rows retrieved
        rows = result.fetchall()
        print(f"Rows retrieved: {rows}")

        # Check if rows are empty
        if not rows:
            print(f"No rows found for coordinates {exonstart} to {exonend}")
            continue

        # Process the rows
        for rec_start, rec_end, offset in rows:
            print(f"Processing record from {rec_start} to {rec_end} with offset {offset}")
            fetched = self._get_record(int(offset))

            for record in fetched:
                if record.id == self._target_seqname:
                    start = record.annotations["start"]
                    end = start + record.annotations["size"] - 1

                    if not (start == rec_start and end == rec_end):
                        raise ValueError(
                            "Expected %s-%s @ offset %s, found %s-%s"
                            % (rec_start, rec_end, offset, start, end)
                        )

            # Debug output for fetched results
            print(f"Fetched record start: {rec_start}, end: {rec_end}")
            yield fetched
#%%
def search_debug3(self, starts, ends):
    # verify the provided exon coordinates
    if len(starts) != len(ends):
        raise ValueError("Every position in starts must have a match in ends")

    # Could it be safer to sort the (exonstart, exonend) pairs?
    for exonstart, exonend in zip(starts, ends):
        exonlen = exonend - exonstart
        if exonlen < 1:
            raise ValueError(
                "Exon coordinates (%d, %d) invalid: exon length (%d) < 1"
                % (exonstart, exonend, exonlen)
            )
    con = self._con

    yielded_rec_coords = set()
    # search for every exon
    for exonstart, exonend in zip(starts, ends):
        try:
            possible_bins = ", ".join(
                map(str, self._region2bin(exonstart, exonend))
            )
        except TypeError:
            raise TypeError(
                "Exon coordinates must be integers "
                "(start=%d, end=%d)" % (exonstart, exonend)
            ) from None
        result = con.execute(
            "SELECT DISTINCT start, end, offset FROM offset_data "
            "WHERE bin IN (%s) "
            "AND (end BETWEEN %s AND %s OR %s BETWEEN start AND end) "
            "ORDER BY start, end, offset ASC;"
            % (possible_bins, exonstart, exonend - 1, exonend - 1)
        )
        rows = result.fetchall()
        for rec_start, rec_end, offset in rows:
            print(f"Processing record from {rec_start} to {rec_end} with offset {offset}")
            if (rec_start, rec_end) in yielded_rec_coords:
                continue
            else:
                yielded_rec_coords.add((rec_start, rec_end))


            fetched = self._get_record(int(offset))

            for record in fetched:
                if record.id == self._target_seqname:
    
                    start = record.annotations["start"]
                    end = start + record.annotations["size"] - 1

                    if not (start == rec_start and end == rec_end):
                        raise ValueError(
                            "Expected %s-%s @ offset %s, found %s-%s"
                            % (rec_start, rec_end, offset, start, end)
                        )

            yield fetched
#%%
MafIO.MafIndex.search = search_debug3
#%%
strand = 1
#chrom = 'chr1'
#start_list = [119563] #[156425872]
#end_list = [119944] #[156425923]
#start_list = [8000]
#end_list = [9000]
start_list = [28804000]
end_list = [28805000]
chrom = 'chrY'
maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
#%%
spliced_maf_full = index_maf.get_spliced(start_list, end_list, strand)
spliced_maf_full
#%%
import sqlite3

# Connect to the .mafindex file
conn = sqlite3.connect(mafindex_filepath)  # Replace with your actual path

# Create a cursor to execute queries
cursor = conn.cursor()
#%%
cursor.execute("SELECT * FROM offset_data;")
rows = cursor.fetchall()
for row in rows:
    print(row)
#%%
start_value = 28804000
end_value = 28806010
cursor.execute("SELECT * FROM offset_data WHERE start >= ? AND end <= ?;", (start_value, end_value))
rows = cursor.fetchall()
for row in rows:
    print(row)
#%%
fetched = list(index_maf.search(start_list, end_list))
print(len(fetched))  # How many alignments are returned?
for alignment in fetched:
    print(alignment)  # Print alignment details
# %%
def count_ungapped(seq, mask):
    return sum(1 for base, valid in zip(seq, mask) if valid and base not in {'-', 'N', 'n'})

#find hg38 row
hg38_row = spliced_maf_full[spliced_maf_full['seqid'].str.startswith("hg38")]
if not hg38_row.empty:
    hg38_seq = list(hg38_row.iloc[0]['seq'])  # Extract the first matching row
else:
    raise ValueError("No hg38 reference found in seqid column!")
#masking
hg38_mask = [base not in {'-', 'N', 'n'} for base in hg38_seq]
#apply counter
spliced_maf_full['ungapped_count'] = spliced_maf_full['seq'].apply(lambda seq: count_ungapped(seq, hg38_mask))
# %%
import pandas as pd
from Bio.AlignIO import MafIO

# Initialize
strand = 1
chrom = 'chr6'
maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
#%%
# count ungapped alignment
def count_ungapped(seq, mask):
    return sum(1 for base, valid in zip(seq, mask) if valid and base not in {'-', 'N', 'n'})

# Initialize tracking
ungapped_counts = {}  # Store results as {species: cumulative_count}
start_pos = 0
step_size = 1000
found_valid_alignment = False  # Track when we start counting
iteration_count = 0  # Counter for iterations
#%%
# Iterate until we hit empty alignments or reach 100 iterations
while iteration_count < 1000:
    print(f"Iteration {iteration_count}: {start_pos}")
    end_pos = start_pos + step_size
    spliced_maf_full = index_maf.get_spliced([start_pos], [end_pos], strand)

    # **Check if spliced_maf_full is a DataFrame**
    if isinstance(spliced_maf_full, pd.DataFrame):
        df = spliced_maf_full  # Use directly
    else:
        # If not a DataFrame, check if it contains only 'NNN...N'
        if all(record.seq == 'N' * len(record.seq) for record in spliced_maf_full):
            if found_valid_alignment:
                break  # Stop if we've started counting and now hit an empty region
            else:
                start_pos += step_size  # Otherwise, keep moving forward
                print("Empty alignment found, shifting position to", start_pos)
                continue

    found_valid_alignment = True  # Found real alignment → Start counting

    # Find hg38 row dynamically
    hg38_row = df[df['seqid'].str.startswith("hg38")]
    if hg38_row.empty:
        start_pos += step_size
        continue  # Skip if no hg38 alignment is found

    hg38_seq = list(hg38_row.iloc[0]['seq'])  # Extract hg38 sequence
    hg38_mask = [base not in {'-', 'N', 'n'} for base in hg38_seq]  # Masking

    # Compute ungapped counts per species
    for _, row in df.iterrows():
        species = row['seqid']
        ungapped_count = count_ungapped(row['seq'], hg38_mask)

        # Accumulate counts
        ungapped_counts[species] = ungapped_counts.get(species, 0) + ungapped_count

    # Move forward
    start_pos += step_size
    iteration_count += 1  # Increment iteration counter
# %%
grouped_counts = pd.Series(ungapped_counts).groupby(lambda x: x.split('.')[0]).sum()

# Print final grouped results
for species, total_count in grouped_counts.items():
    print(f"{species}: {total_count}")
#%%
import pandas as pd
from Bio.AlignIO import MafIO

# Initialize
strand = 1
chrom = 'chrY'
maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)

# Function to count ungapped bases
def count_ungapped(seq, mask):
    return sum(1 for base, valid in zip(seq, mask) if valid and base not in {'-', 'N', 'n'})

# Function to check if the region is empty (only 'N' bases)
def is_empty_region(spliced_maf_full):
    return all(record.seq == 'N' * len(record.seq) for record in spliced_maf_full)

# Initialize tracking variables
ungapped_counts = {}  # Store results as {species: cumulative_count}
start_pos = 0
step_size = 1000  # Window size for processing
skip_size = 10000000  # Large skip size to scan chunks faster
found_valid_alignment = False  # Track when we find the first valid alignment

# Step 1: Scan the first chunk normally, starting from position 0
while not found_valid_alignment:
    print(f"Checking chunk: {start_pos} - {start_pos + step_size}")
    spliced_maf_full = index_maf.get_spliced([start_pos], [start_pos + step_size], strand)

    # Check if the result is a DataFrame (indicating valid data)
    if isinstance(spliced_maf_full, pd.DataFrame):
        # Look for hg38 data
        hg38_row = spliced_maf_full[spliced_maf_full['seqid'].str.startswith('hg38')]

        # Process if hg38 data exists
        if not hg38_row.empty:
            found_valid_alignment = True
        else:
            print(f"Found NNNNNN at position {start_pos}. Moving to next chunk.")
            start_pos += step_size  # Move to the next chunk of 1000 bases
    
    else:
        print(f"Empty region at position {start_pos}. Moving to next chunk.")
        start_pos += step_size  # Move to the next chunk of 1000 bases

# Step 2: Once valid region is found, skip by 1 million bases
print(f"Found valid alignment at position {start_pos}. Now skipping by 10M bases.")

while found_valid_alignment:
    print(f"checking chunk: {start_pos} - {start_pos + step_size}")
    spliced_maf_full = index_maf.get_spliced([start_pos], [start_pos + step_size], strand)
    # Check if the result is a DataFrame (indicating valid data)
    if isinstance(spliced_maf_full, pd.DataFrame):
        # Look for hg38 data again
        hg38_row = spliced_maf_full[spliced_maf_full['seqid'].str.startswith('hg38')]

        # Check for empty regions
        if not hg38_row.empty:
            print(f"Valid data found at position {start_pos} - {start_pos + step_size}")
            start_pos += skip_size  # Break when we find valid data, else continue skipping by 1M bases
        else:
            print(f"Empty region found at position {start_pos} - {start_pos + step_size}")
            found_valid_alignment = False

    else:
        if skip_size > step_size:
            print(f"Empty region at position {start_pos}. reset the loop.")
            start_pos=int(start_pos-skip_size)
            skip_size=int(skip_size/10)
        else:
            print(f"Empty region at position {start_pos}. terminate the loop.")
            found_valid_alignment = False 
# %%
#note:
#chr1: 9000-248947000 
#chr2: 9000-242184000 
#chr3: 1000-198236000
#chr4: 10000-190205000
#chr5: 13000-181479000
#chr6: 60000-170746000
#chr7: 7000-159336000
#chr8: 59000-145087000
#chr9: 10000-138335000
#chr10: 10000-133788000
#chr11: 60000-135086000
#chr12: 9000-133265000
#chr13: 18171000-114355000
#chr14: 16019000-106884000
#chr15: 17723000-101982000
#chr16: 9000-90229000
#chr17: 59000-83252000
#chr18: 10000-80271000
#chr19: 60000-58608000
#chr20: 60000-64335000
#chr21: 5005000-46700000 
#chr22: 10510000-50817000
#chrX: 9000-156032000
#chrY: 2781000-56888000
# %%
import os
import pickle
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from Bio.AlignIO import MafIO

# Initialize
strand = 1
chrom = 'chr1'
maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)

#%%
# Directory for temp results
TEMP_DIR = "/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/aligned_count/temp_results"
os.makedirs(TEMP_DIR, exist_ok=True)

# Function to count ungapped bases
def count_ungapped(seq, mask):
    return sum(1 for base, valid in zip(seq, mask) if valid and base not in {'-', 'N', 'n'})

# Worker function for parallel processing
def process_chunk(chunk_range):
    start_pos, end_pos = chunk_range
    temp_filepath = os.path.join(TEMP_DIR, f"chunk_{start_pos}_{end_pos}.pkl")

    # Skip if temp file exists
    if os.path.exists(temp_filepath):
        print(f"Skipping {start_pos}-{end_pos}, already processed.")
        return None  

    print(f"Processing {start_pos}-{end_pos}")

    # Fetch MAF data
    spliced_maf_full = index_maf.get_spliced([start_pos], [end_pos], strand)
    
    if not isinstance(spliced_maf_full, pd.DataFrame):
        print(f"Empty or invalid chunk: {start_pos}-{end_pos}")
        return None

    hg38_row = spliced_maf_full[spliced_maf_full['seqid'].str.startswith("hg38")]
    if hg38_row.empty:
        print(f"No hg38 data at {start_pos}-{end_pos}")
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
        print(f'DONE: {start_pos}-{end_pos}')

    return species_counts

# Main execution
def process_all_chunks(start, end, step_size):
    chunk_ranges = [(pos, min(pos + step_size, end)) for pos in range(start, end, step_size)]

    with ProcessPoolExecutor() as executor:
        results = executor.map(process_chunk, chunk_ranges)


#%%
start_alignment = 9000  # Example start
end_alignment = 248947000  # Example end
step_size = 1000  # Process in 1k blocks

process_all_chunks(start_alignment, end_alignment, step_size)

#%%
# Aggregate counts from temp files
total_counts = {}
for file in os.listdir(TEMP_DIR):
    if file.endswith(".pkl"):
        with open(os.path.join(TEMP_DIR, file), "rb") as f:
            chunk_counts = pickle.load(f)
            for species, count in chunk_counts.items():
                total_counts[species] = total_counts.get(species, 0) + count
#%%
# %%
