# %%
from tracemalloc import start
from unittest import result
from numpy import interp, repeat
import pandas as pd

from dev.installation.biopython.Bio.Wise import align

# Load RepeatMasker table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility

repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)

# Select LTR subfamily
subfamily = 'THE1C'
subfamily_table = repeatmasker_table[repeatmasker_table['repName'] == subfamily]

# Extract row indices
subfamily_indices = subfamily_table.index.tolist()

# Store blocks separately instead of merging
categorized_blocks = {
    'complete': [],
    'one-sided': [],
    'standalone': []
}

# Track used indices to prevent duplication
used_indices = set()

# Check for patterns in a ±10 row window
for idx in subfamily_indices:
    if idx in used_indices:  # Skip if already used
        continue

    # Define the ±10 row window
    start_idx = max(0, idx - 20)
    end_idx = min(len(repeatmasker_table) - 1, idx + 20)
    window = repeatmasker_table.loc[start_idx:end_idx]

    # Check strand consistency
    target_strand = repeatmasker_table.loc[idx, 'strand']
    window = window[window['strand'] == target_strand]

    # Get unique repNames in the window
    rep_names = list(window['repName'])

    # Locate occurrences of subfamily elements
    ltr_positions = [i for i, name in enumerate(rep_names) if name == subfamily]
    int_positions = [i for i, name in enumerate(rep_names) if name == f"{subfamily}-int"]

    # Check for complete pattern: 'xxx' - 'xxx-int' - 'xxx'
    if len(ltr_positions) >= 2 and len(int_positions) >= 1:
        first_ltr = ltr_positions[0]
        last_ltr = ltr_positions[-1]
        has_int_between = any(first_ltr < pos < last_ltr for pos in int_positions)

        if has_int_between:
            extracted_table = window.iloc[first_ltr:last_ltr + 1]  # Include everything between first and last LTR
            categorized_blocks['complete'].append(extracted_table)

            # Mark these indices as used
            used_indices.update(extracted_table.index)
            continue  # No need to check further, already categorized

    # Check for partial cases and extract all rows in between
    if len(ltr_positions) >= 1 and len(int_positions) >= 1:
        first = min(ltr_positions + int_positions)
        last = max(ltr_positions + int_positions)
        extracted_table = window.iloc[first:last + 1]  # Include everything in between
        categorized_blocks['one-sided'].append(extracted_table)

        # Mark these indices as used
        used_indices.update(extracted_table.index)
        continue  # No need to check standalone

    elif len(ltr_positions) > 1:  # Case: 'xxx' - 'xxx'
        first = ltr_positions[0]
        last = ltr_positions[-1]
        extracted_table = window.iloc[first:last + 1]  # Include everything in between
        categorized_blocks['one-sided'].append(extracted_table)

        # Mark these indices as used
        used_indices.update(extracted_table.index)
        continue  # No need to check standalone

    # If none of the above, classify as standalone
    extracted_table = window[window['repName'] == subfamily]

    # Make sure rows aren’t already used
    extracted_table = extracted_table[~extracted_table.index.isin(used_indices)]
    if not extracted_table.empty:
        categorized_blocks['standalone'].append(extracted_table)
        used_indices.update(extracted_table.index)

# Print summary
for category, blocks in categorized_blocks.items():
    print(f"{category}: {len(blocks)} blocks")

# Each block is a separate DataFrame in the list
complete_blocks = categorized_blocks['complete']
one_sided_blocks = categorized_blocks['one-sided']
standalone_blocks = categorized_blocks['standalone']

# %%
global_block_count = 1  # Global running counter

# Add block_id column to each block in categorized blocks
for category, blocks in categorized_blocks.items():
    for i, block in enumerate(blocks):
        block_id = f"{subfamily}_{category}_{global_block_count}"
        block['block_id'] = block_id  # Add block_id column
        global_block_count += 1  # Increment counter

# After adding the block IDs, you can concatenate the blocks as before
all_blocks = pd.concat([block for blocks in categorized_blocks.values() for block in blocks])

# Now `all_blocks` will contain a `block_id` column for each individual block
#%%
all_blocks.to_csv(f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/{subfamily}.blocks.txt', sep='\t', index=False)
# %%
internal_id_table = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/{subfamily}.internal_id.txt'
age_table = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/{subfamily}.txt'
# %%
internal_id_df = pd.read_csv(internal_id_table, sep='\t')
age_df = pd.read_csv(age_table, sep='\t')
# %%
age_internal_id_df = pd.merge(age_df, internal_id_df, on='internal_id')
#%%
internal_id_table_int = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/{subfamily}-int.internal_id.txt'
age_table_int = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/{subfamily}-int.txt'
# %%
internal_id_df_int = pd.read_csv(internal_id_table_int, sep='\t')
age_df_int = pd.read_csv(age_table_int, sep='\t')
# %%
age_internal_id_df_int = pd.merge(age_df_int, internal_id_df_int, on='internal_id')
#%%
age_internal_id = pd.concat([age_internal_id_df, age_internal_id_df_int])
# %%
all_blocks_age = pd.merge(all_blocks, age_internal_id,left_index=True, right_on='rmsk_index', how='left')
#%%
complete_block_id=all_blocks_age[all_blocks_age['block_id'].str.contains('complete')]['block_id'].unique()
#%%

# %%
from ma_mapper import sequence_alignment
from Bio import SeqIO
import sys
sys.path.append('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime')
from core_engine import e_value_calculation
source_fasta = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa'
records = SeqIO.to_dict(SeqIO.parse(open(source_fasta), 'fasta'))
#%%
def align_and_compare(row1, row2):
    """
    Aligns two sequences onto their consensus positions and compares them.
    """
    print(f'original row1: {row1["repStart"]} {row1["repEnd"]}')
    print(f'original row2: {row2["repStart"]} {row2["repEnd"]}')
    start1, end1 = row1["repStart"]-1, row1["repEnd"]-1
    start2, end2 = row2["repStart"]-1, row2["repEnd"]-1
    seq1, seq2 = row1["sequence"], row2["sequence"]
    print(f'row1: {start1} {end1} {len(seq1)}')
    print(f'row2: {start2} {end2} {len(seq2)}')
    # Determine alignment range
    align_start = min(start1, start2)
    align_end = max(end1, end2)
    align_length = align_end - align_start+1
    print(align_start, align_end, align_length)
    # Create blank aligned sequences
    aligned_seq1 = ["-"] * align_length
    aligned_seq2 = ["-"] * align_length

    # Insert sequences at correct positions
    aligned_seq1[start1:end1] = seq1
    aligned_seq2[start2:end2] = seq2

    # Convert to strings
    aligned_seq1 = "".join(aligned_seq1)
    aligned_seq2 = "".join(aligned_seq2)
    print(aligned_seq1)
    print(aligned_seq2)
    # Compute similarity
    return e_value_calculation.affine_count_simple(aligned_seq1, aligned_seq2)
#%%
from Bio.Align import PairwiseAligner
import pandas as pd
import numpy as np
from itertools import combinations

def calculate_div(df):
    """Calculate weighted divergence for merged fragments."""
    frag_length = df['genoEnd'] - df['genoStart']
    frag_div = frag_length * df['milliDiv']
    return frag_div.sum() / frag_length.sum() if frag_length.sum() != 0 else 0  # Avoid division by zero

aligner = PairwiseAligner()
aligner.mode = 'global'  
aligner.match_score = 1  
aligner.mismatch_score = -1  
aligner.open_gap_score = -4  
aligner.extend_gap_score = -1  

mutation_rate = 2.2e-9  # Mutation rate per site per year

alignment_results = []

for block_id in all_blocks_age['block_id'].unique():
    block_df = all_blocks_age[all_blocks_age['block_id'] == block_id].copy()

    # Get unique internal IDs, filtering out NaNs and ensuring they're strings
    valid_internal_ids = block_df['internal_id'].dropna()
    valid_internal_ids = valid_internal_ids[valid_internal_ids.apply(lambda x: isinstance(x, str))].unique()

    # Merge sequences for valid internal IDs
    merged_sequences = {}
    div_values = {}
    age_values = {}

    for internal_id in valid_internal_ids:
        sub_df = block_df[block_df['internal_id'] == internal_id]

        # Merge sequences
        merged_seq = "".join(
            sequence_alignment.extract_sequence(row['genoName'], row['genoStart'], row['genoEnd'], row['strand'], records)
            for _, row in sub_df.iterrows()
        )
        merged_sequences[internal_id] = merged_seq
        
        # Calculate divergence for merged fragments
        div_values[internal_id] = calculate_div(sub_df)
        age_values[internal_id] = sub_df['te_age'].values[0] 
    # Get LTR sequences only
    ltr_sequences = {k: v for k, v in merged_sequences.items() if "int" not in k.lower()}
    int_age = [round(v,2) for k, v in age_values.items() if "int" in k.lower()]
    # Pairwise align all LTR sequences
    for (id1, seq1), (id2, seq2) in combinations(ltr_sequences.items(), 2):
        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]
        aligned_seq1, aligned_seq2 = best_alignment

        #LTR DIVERGENCE CALCULATION
        # Initialize mismatch counts
        match_count = 0
        mismatch_count = 0
        gap_count = 0
        transition_count = 0
        transversion_count = 0

        # Define transition and transversion sets
        transitions = {("a", "g"), ("g", "a"), ("c", "t"), ("t", "c")}
        transversions = {("a", "c"), ("c", "a"), ("a", "t"), ("t", "a"),
                        ("g", "c"), ("c", "g"), ("g", "t"), ("t", "g")}
        
        # Count matches, mismatches, and gaps
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == b and a != "-":
                match_count += 1
            elif a == "-" or b == "-":
                gap_count += 1
            else:
                mismatch_count += 1
                if (a, b) in transitions:
                    transition_count += 1
                elif (a, b) in transversions:
                    transversion_count += 1
        # Total sites for K2P (excluding gaps)
        total_sites = match_count + mismatch_count
        p = transition_count / total_sites
        q = transversion_count / total_sites

        # Standard K2P (equal Ti/Tv rates)
        K = -0.5 * np.log((1 - 2*(p+q)) * np.sqrt(1 - 2*q)) if (1 - 2*(p+q)) > 0 and (1 - 2*q) > 0 else np.nan
        
        # Refined K2P (actual Ti/Tv ratio)
        refined_K = -0.5 * np.log((1 - 2*p - q) * np.sqrt(1 - 2*q)) if (1 - 2*p - q) > 0 and (1 - 2*q) > 0 else np.nan

        # Convert K2P distances into age estimates
        ldAge = K / (2 * mutation_rate) if not np.isnan(K) else np.nan
        refined_ldAge = refined_K / (2 * mutation_rate) if not np.isnan(refined_K) else np.nan

        seq_length = max(len(aligned_seq1), len(aligned_seq2))
        identity = (match_count / seq_length) * 100

        # Calculate kAge
        div1 = div_values[id1]
        age1 = age_values[id1]
        div2 = div_values[id2]
        age2 = age_values[id2]
        kAge1 = round(div1 / (100 * mutation_rate * 1e6), 2)
        kAge2 = round(div2 / (100 * mutation_rate * 1e6), 2)

        alignment_results.append({
            "Block ID": block_id,
            "LTR1 ID": id1,
            "LTR2 ID": id2,
            "LTR ldAge": round(ldAge/1e+6, 2),
            "LTR refined ldAge": round(refined_ldAge/1e+6, 2),
            "LTR1 age": age1,
            "LTR2 age": age2,
            "int Age": int_age,
            "LTR1 kAge": kAge1,
            "LTR2 kAge": kAge2,
            "LTR1 div": round(div1,2),
            "LTR2 div": round(div2,2),
            "LTR Similarity (%)": round(identity, 2),
            "LTR Matches": match_count,
            "LTR Mismatches": mismatch_count,
            "LTR Gaps": gap_count,
            "LTR Transitions": transition_count,
            "LTR Transversions": transversion_count,
            "LTR K2P": round(K,2),
            "LTR refined K2P": round(refined_K,2),
        })

alignment_df = pd.DataFrame(alignment_results)
#%%
alignment_df.to_csv(f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/{subfamily}.ltr.txt', sep='\t', index=False)
#%%
alignment_df
#%%
complete_block_df=alignment_df[alignment_df['Block ID'].str.contains('complete')]
block_num=len(complete_block_df['Block ID'].unique())
#%%
filtered_df = complete_block_df[complete_block_df['Block ID'].duplicated(keep=False) == False]
twoLTR_num=len(filtered_df['Block ID'].unique())
equal_teatime_df=filtered_df[filtered_df['LTR1 age']==filtered_df['LTR2 age']]
equalTEATIME_num = len(equal_teatime_df['Block ID'].unique())
#%% ref1
print(f'subfamily: {subfamily}')
print(f'all complete LTR blocks:{block_num}')
print(f'blocks with two LTRs:{twoLTR_num}')
perc_hit=equalTEATIME_num/twoLTR_num
print(f'Two LTR with equal teatime age:{equalTEATIME_num}')
print(f'simple case with equal TEATIME/LTR block with two LTRs: {perc_hit}')

#%%

#%%

#%%
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# Extract values for plotting
x = alignment_df["LTR1 age"]
y = alignment_df["LTR2 age"]
z = alignment_df["LTR refined ldAge"]

# Create 3D scatter plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with color mapping to similarity
sc = ax.scatter(x, y, z, c=z, cmap='viridis', edgecolor='k', s=50)

# Labels and title
ax.set_xlabel("LTR1 age")
ax.set_ylabel("LTR2 age")
ax.set_zlabel("LTR ldAge")
ax.set_title("3D Plot of LTR Age estimated by teatime and pairwise distance")

# Color bar
cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
cbar.set_label("LTR ldAge")

plt.show()
# %%
for complete_block in complete_block_id:
    block_df = all_blocks_age[all_blocks_age['block_id'] == complete_block]

    # Find the expected full length of the TE
    te_repEnd = block_df['repEnd'].median()
    te_repStart = 1
    te_length = te_repEnd - te_repStart + 1  
    print(f"Solving the simplest block: {complete_block}")

    ltr_df = block_df[block_df['repName'] == subfamily].copy().reset_index(drop=True)
    int_df = block_df[block_df['repName'].str.contains("-int", case=False, na=False)].copy().reset_index(drop=True)

    if len(ltr_df) == 2 and len(int_df) == 1:
        # Assume LTRs are at index 0 (head) and 1 (tail), while INT is index 0 in int_df
        head_ltr = ltr_df.iloc[0]
        tail_ltr = ltr_df.iloc[1]
        int_element = int_df.iloc[0]

        # Distance from LTR to INT (Head)
        distance1 = abs(int_element['genoStart'] - head_ltr['genoEnd'])

        # Distance from INT to LTR (Tail)
        distance2 = abs(tail_ltr['genoStart'] - int_element['genoEnd'])

        print(f"H: {distance1} ({head_ltr['repStart'],head_ltr['repEnd']}), T: {distance2} ({tail_ltr['repStart'],tail_ltr['repEnd']})")
    else:
        print(f"Skipping block {complete_block}: complex structure")
#%%
import plotly.express as px
import pandas as pd

# Convert alignment results to DataFrame (if not already)
alignment_df = pd.DataFrame(alignment_results)

# Create interactive 3D scatter plot
fig = px.scatter_3d(
    alignment_df,
    x="LTR1 kAge",
    y="LTR2 kAge",
    z="LTR refined ldAge",
    color="LTR Similarity (%)",  # Color by estimated age
    hover_data=["Block ID", "LTR1 age", "LTR2 age"],  # Extra info on hover
    title="LTR kAge vs. ldAge",
)
fig.update_traces(marker=dict(size=3))
fig.update_layout(
    scene=dict(
        xaxis_title="LTR1 kAge",
        yaxis_title="LTR2 kAge",
        zaxis_title="LTR ldAge",
    ),
    coloraxis_colorbar=dict(title="LTR Similarity (%)"),
)

fig.show()
#%%
import plotly.io as pio
pio.write_html(fig, "the1c_kAgekAge_ldAge.html")
#%%
import plotly.express as px
import pandas as pd

# Convert alignment results to DataFrame (if not already)
alignment_df = pd.DataFrame(alignment_results)

# Create interactive 3D scatter plot
fig = px.scatter_3d(
    alignment_df,
    x="LTR1 age",
    y="LTR2 age",
    z="LTR refined ldAge",
    color="LTR Similarity (%)",  # Color by estimated age
    hover_data=["Block ID", "LTR1 age", "LTR2 age"],  # Extra info on hover
    title="LTR kAge vs. ldAge",
)
fig.update_traces(marker=dict(size=3))
fig.update_layout(
    scene=dict(
        xaxis_title="LTR1 age",
        yaxis_title="LTR2 age",
        zaxis_title="LTR ldAge",
    ),
    coloraxis_colorbar=dict(title="LTR Similarity (%)"),
)

fig.show()
#%%
import plotly.io as pio
pio.write_html(fig, "the1c_ageage_ldAge.html")
#%%
import matplotlib.pyplot as plt
import pandas as pd

# Load data (if not already in a DataFrame)
alignment_df = pd.DataFrame(alignment_results)

# Create figure
fig, ax = plt.subplots(figsize=(7, 6))

# Scatter plot
sc = ax.scatter(
    alignment_df["LTR1 kAge"],
    alignment_df["LTR2 kAge"],
    c=alignment_df["LTR refined ldAge"],  # Color by LTR ldAge
    cmap="viridis",  # Colormap (change to "plasma", "coolwarm", etc. if needed)
    marker='o',
    rasterized=True,
    s=1,  # Marker size
    alpha=0.75  # Transparency for better overlap visibility
)

# Add colorbar
cbar = plt.colorbar(sc)
cbar.set_label("LTR ldAge (MYA)")

# Labels and title
ax.set_xlabel("LTR1 kAge (MYA)")
ax.set_ylabel("LTR2 kAge (MYA)")
ax.set_title("LTR pairwise age comparison (Kimura vs. LTR Divergence)")

# Show plot
plt.show()
#%%
import seaborn as sns
import matplotlib.pyplot as plt

# Joint Plot for LTR1 Age vs LTR2 Age, colored by ldAge
sns.jointplot(
    x="LTR1 age", 
    y="LTR2 age", 
    data=alignment_df, 
    kind="hex",  # hexbin for density of points
    cmap="viridis",
    hue="LTR refined ldAge"
)

plt.suptitle("Joint Plot: LTR1 vs LTR2 Age, Colored by ldAge", y=1.02)
plt.show()
#%%
import seaborn as sns
import matplotlib.pyplot as plt

# Create pairplot for LTR1, LTR2, and ldAge
sns.pairplot(alignment_df[['LTR1 age', 'LTR2 age', 'LTR refined ldAge']], hue='LTR refined ldAge',
             palette="viridis", plot_kws={'alpha':0.6, 's':10})

plt.suptitle("Pairwise Relationships of LTR1 Age, LTR2 Age, and ldAge", y=1.02)
plt.show()
# %%
ideal_cases = 0

for complete_block in complete_block_id:
    block_df = all_blocks_age[all_blocks_age['block_id'] == complete_block]

    # Find expected TE length
    te_repEnd = block_df['repEnd'].median()
    te_length = te_repEnd - 1 + 1  # Ensuring correct base-1 to base-0 handling

    if block_df.shape[0] == 3:
        ltr_df = block_df[block_df['repName'] == subfamily].copy().reset_index(drop=True)
        int_df = block_df[block_df['repName'].str.contains("-int", case=False)].copy().reset_index(drop=True)

        if ltr_df.shape[0] == 2 and int_df.shape[0] == 1:
            # Check LTR completeness criteria
            if all((ltr_df['repStart'] == 1) & (ltr_df['repEnd'] == te_length)):
                ideal_cases += 1
                ltr_lengths = (ltr_df['genoEnd'] - ltr_df['genoStart']).tolist()
                int_length = (int_df['genoEnd'] - int_df['genoStart']).values[0]
                print(f"Perfect case: {complete_block} (LTRs: {ltr_lengths}, INT: {int_length})")

print(f"Total ideal cases: {ideal_cases}")

# %%
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.mode = 'global'  # Use global alignment
aligner.match_score = 1  # Matches get +1
aligner.mismatch_score = -1  # Mismatches get -1
aligner.open_gap_score = -4  # Gap opening penalty
aligner.extend_gap_score = -1  # Gap extension penalty
alignment_results = []
for complete_block in complete_block_id:
    block_df = all_blocks_age[all_blocks_age['block_id'] == complete_block]

    te_repEnd = block_df['repEnd'].median()
    te_length = te_repEnd - 1 + 1  

    if block_df.shape[0] == 3:
        ltr_df = block_df[block_df['repName'] == subfamily].copy().reset_index(drop=True)
        int_df = block_df[block_df['repName'].str.contains("-int", case=False)].copy().reset_index(drop=True)
        ltr_df['sequence'] = ltr_df.apply(
            lambda row: sequence_alignment.extract_sequence(
                row['genoName'], row['genoStart'], row['genoEnd'], row['strand'], records
            ), axis=1
        )
        int_age = int_df['te_age'].values[0]
        int_div = int_df['milliDiv'].values[0]
        if ltr_df.shape[0] == 2 and int_df.shape[0] == 1:
            if all((ltr_df['repStart'] == 1) & (ltr_df['repEnd'] == te_length)):
                ltr1_seq = ltr_df.iloc[0]['sequence']
                ltr2_seq = ltr_df.iloc[1]['sequence']
                ltr1_age = ltr_df.iloc[0]['te_age']
                ltr2_age = ltr_df.iloc[1]['te_age']
                ltr1_div = ltr_df.iloc[0]['milliDiv']
                ltr2_div = ltr_df.iloc[1]['milliDiv']
                # Perform global alignment
                alignments = aligner.align(ltr1_seq, ltr2_seq)
                best_alignment = alignments[0]  # Get best alignment

                # Extract aligned sequences
                aligned_seq1, aligned_seq2 = best_alignment

                # Count match/mismatch/gap
                match_count = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != "-")
                mismatch_count = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != "-" and b != "-")
                gap_count = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == "-" or b == "-")

                seq_length = max(len(aligned_seq1), len(aligned_seq2))
                identity = (match_count / seq_length) * 100

                #print(f"Perfect case: {complete_block} Age: {ltr1_age},{ltr2_age} (Identity: {identity:.2f}%)")
                #print(f"Matches: {match_count}, Mismatches: {mismatch_count}, Gaps: {gap_count}")
                #print(best_alignment)
                alignment_results.append({
                    "Block ID": complete_block,
                    "LTR1 div": ltr1_div,
                    "int div": int_div,
                    "LTR2 div": ltr2_div,
                    "LTR1 kAge": round(ltr1_div/(100*(2.2e-9)*1e+6),2), 
                    "int kAge": round(int_div/(100*(2.2e-9)*1e+6),2),
                    "LTR2 kAge": round(ltr2_div/(100*(2.2e-9)*1e+6),2), 
                    "LTR1 Age": ltr1_age,
                    "int Age": int_age,
                    "LTR2 Age": ltr2_age,
                    "LTR Similarity (%)": round(identity, 2),
                    "LTR Matches": match_count,
                    "LTR Mismatches": mismatch_count,
                    "LTR Gaps": gap_count
                })
alignment_df = pd.DataFrame(alignment_results)
# %%
alignment_df
# %%
from Bio.Align import PairwiseAligner
import pandas as pd
import numpy as np
from itertools import combinations

def calculate_div(df):
    """Calculate weighted divergence for merged fragments."""
    frag_length = df['genoEnd'] - df['genoStart']
    frag_div = frag_length * df['milliDiv']
    return frag_div.sum() / frag_length.sum() if frag_length.sum() != 0 else 0  # Avoid division by zero

aligner = PairwiseAligner()
aligner.mode = 'global'  
aligner.match_score = 1  
aligner.mismatch_score = -1  
aligner.open_gap_score = -4  
aligner.extend_gap_score = -1  

mutation_rate = 2.2e-9  # Mutation rate per site per year

alignment_results = []

for block_id in all_blocks_age['block_id'].unique():
    block_df = all_blocks_age[all_blocks_age['block_id'] == block_id].copy()

    # Get unique internal IDs, filtering out NaNs and ensuring they're strings
    valid_internal_ids = block_df['internal_id'].dropna()
    valid_internal_ids = valid_internal_ids[valid_internal_ids.apply(lambda x: isinstance(x, str))].unique()

    # Merge sequences for valid internal IDs
    merged_sequences = {}
    div_values = {}

    for internal_id in valid_internal_ids:
        sub_df = block_df[block_df['internal_id'] == internal_id]

        # Merge sequences
        merged_seq = "".join(
            sequence_alignment.extract_sequence(row['genoName'], row['genoStart'], row['genoEnd'], row['strand'], records)
            for _, row in sub_df.iterrows()
        )
        merged_sequences[internal_id] = merged_seq
        
        # Calculate divergence for merged fragments
        div_values[internal_id] = calculate_div(sub_df)

    # Get LTR sequences only
    ltr_sequences = {k: v for k, v in merged_sequences.items() if "int" not in k.lower()}
    ltr_divs = {k: v for k, v in div_values.items() if "int" not in k.lower()}

    # Pairwise align all LTR sequences
    for (id1, seq1), (id2, seq2) in combinations(ltr_sequences.items(), 2):
        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]
        aligned_seq1, aligned_seq2 = best_alignment

        # Count matches, mismatches, and gaps
        match_count = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != "-")
        mismatch_count = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != "-" and b != "-")
        gap_count = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == "-" or b == "-")

        seq_length = max(len(aligned_seq1), len(aligned_seq2))
        identity = (match_count / seq_length) * 100

        # Calculate kAge
        div1 = ltr_divs[id1]
        div2 = ltr_divs[id2]
        kAge1 = round(div1 / (100 * mutation_rate * 1e6), 2)
        kAge2 = round(div2 / (100 * mutation_rate * 1e6), 2)

        alignment_results.append({
            "Block ID": block_id,
            "LTR1 ID": id1,
            "LTR2 ID": id2,
            "LTR1 div": div1,
            "LTR2 div": div2,
            "LTR1 kAge": kAge1,
            "LTR2 kAge": kAge2,
            "LTR Similarity (%)": round(identity, 2),
            "LTR Matches": match_count,
            "LTR Mismatches": mismatch_count,
            "LTR Gaps": gap_count
        })

alignment_df = pd.DataFrame(alignment_results)

# %% investigation
subfamily = "THE1A"
e_table_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/e_value/{subfamily}.txt'
e_table = pd.read_csv(e_table_filepath, sep='\t')
# %%
e_table[e_table['internal_id']=='THE1A_SINGLE_3164']
# %%
e_table[e_table['internal_id']=='THE1A_COMPLETE_3165']
# %%
internal_id_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/{subfamily}.internal_id.txt'
internal_id_df = pd.read_csv(internal_id_filepath, sep='\t')
# %%
internal_id_df[internal_id_df['internal_id']=='THE1A_SINGLE_3164']
# %%
repeatmasker_table[repeatmasker_table.index==862705]
# %%
teatime_unequal=filtered_df[filtered_df['LTR1 age']!=filtered_df['LTR2 age']]
# %%
teatime_unequal['age_diff'] = teatime_unequal['LTR1 age'] - teatime_unequal['LTR2 age']
# %%
