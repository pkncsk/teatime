#%%
import numpy as np
import pandas as pd
#%%
ltr_age_table_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/THE1C.ltr_fix.txt'

ltr_age_table=pd.read_csv(ltr_age_table_filepath, sep='\t')
#%%
# pattern looks for one of the three intTypes surrounded by underscores
pattern = r'_(aINT|bINT|nINT)_'

ltr_age_table['int_type'] = ltr_age_table['internal_id'].str.extract(pattern)

print(ltr_age_table[['internal_id', 'int_type']].head())

#%%
internal_id_tbl_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline/THE1C.internal_id.txt'
internal_id_tbl=pd.read_csv(internal_id_tbl_filepath, sep='\t')
#%%x
ltr_age_id_table=ltr_age_table.merge(internal_id_tbl, on='internal_id')
#%%
repeatmasker_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'

repeatmasker_table = pd.read_csv(repeatmasker_filepath, sep='\t', index_col=0)
#%%
ltr_age_id_coord_tbl=ltr_age_id_table.merge(repeatmasker_table[['genoName', 'genoStart', 'genoEnd', 'strand']], left_on='rmsk_index', right_index=True)
#%%
summary = (
    ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['int_type'].isin(['aINT', 'bINT'])]
    .groupby('block_id')
    .agg(num_rows=('internal_id', 'count'),
         unique_ids=('internal_id', 'nunique'),
         unique_types=('int_type', 'nunique'))
)

# Focus only on candidates with both aINT and bINT
candidates = summary[summary['unique_types'] == 2]

# Now check how many have more than 2 internal_ids (e.g. due to fragments)
problem_blocks = candidates[candidates['unique_ids'] > 2]
good_blocks = candidates[candidates['unique_ids'] == 2]
#%%
ltr_age_id_coord_tbl['length'] = ltr_age_id_coord_tbl['genoEnd'] - ltr_age_id_coord_tbl['genoStart']
valid_blocks = ltr_age_id_coord_tbl.groupby(['block_id', 'int_type'])['internal_id'].nunique().unstack(fill_value=0)
#%%
valid_id_list = valid_blocks[(valid_blocks['aINT'] == 1) & (valid_blocks['bINT'] == 1)].index
df_filtered = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(valid_id_list)]
#%%
block_lengths_ok  = (
    df_filtered.groupby('block_id')['length']
    .apply(lambda x: (x >= 300).sum() == 2)
)
#%%
valid_id_list2 = block_lengths_ok[block_lengths_ok].index
df_filtered = df_filtered[df_filtered['block_id'].isin(valid_id_list2)]
#%%
import os
from ma_mapper import sequence_alignment
genome_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa'
for block_id in valid_id_list:
    print(block_id)
    fasta_path = f"./tmp/{block_id}.fasta"
    aligned_path = f"{fasta_path}.aligned"
    if (os.path.isfile(aligned_path) == False):
        block_df = df_filtered[df_filtered['block_id'] == block_id]
        block_coords = block_df[block_df['int_type'].isin(['aINT', 'bINT'])][
            ['genoName', 'genoStart', 'genoEnd', 'strand', 'internal_id']
        ].copy()
        block_coords = block_coords.rename(columns={
        'genoName': 'chrom',
        'genoStart': 'start',
        'genoEnd': 'end',
        'internal_id': 'name',
        'strand': 'strand'
        })
        block_coords['score'] = 10

        # Reorder columns to BED standard
        block_coords = block_coords[['chrom', 'start', 'end', 'name', 'score', 'strand']]

        seq_records = sequence_alignment.sequence_io(block_coords, source_fasta=genome_fasta, save_to_file=False)

        
        with open(fasta_path, "w") as f:
            from Bio import SeqIO
            SeqIO.write(seq_records, f, "fasta")
        
        # Align with MAFFT
        
        sequence_alignment.mafft_align(fasta_path, output_filepath=aligned_path, nthread=6)
    else:
        print('already done')

#%%
from Bio import SeqIO

def calculate_identity_from_alignment(aligned_fasta_path):
    records = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    if len(records) != 2:
        raise ValueError(f"Expected 2 sequences in alignment, got {len(records)}")

    seq1 = str(records[0].seq).upper()
    seq2 = str(records[1].seq).upper()

    if len(seq1) != len(seq2):
        raise ValueError("Aligned sequences are not the same length")

    aligned_length = len(seq1)
    matches = sum((a == b) and (a != '-') for a, b in zip(seq1, seq2))
    gaps = sum((a == '-' or b == '-') for a, b in zip(seq1, seq2))

    identity_pct = matches / aligned_length * 100
    gap_pct = gaps / aligned_length * 100

    return {
        'aligned_length': aligned_length,
        'matches': matches,
        'identity_pct': identity_pct,
        'gaps': gaps,
        'gap_pct': gap_pct
    }
#%%
result = {}
for block_id in valid_id_list2:
    fasta_path = f"./tmp/{block_id}.fasta"
    aligned_path = f"{fasta_path}.aligned"
    ltr_id_stats=calculate_identity_from_alignment(aligned_path)
    result[block_id] = ltr_id_stats
#%%
identity_df = pd.DataFrame.from_dict(result, orient='index')
identity_df.index.name = 'block_id'
identity_df.reset_index(inplace=True)
# %%
age_diff_df = (
    ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['int_type'].isin(['aINT', 'bINT'])]
    .groupby('block_id')[['te_age', 'corrected_te_age']]
    .agg(['min', 'max'])
)

# Flatten MultiIndex
age_diff_df.columns = ['_'.join(col).strip() for col in age_diff_df.columns.values]
age_diff_df.reset_index(inplace=True)

# Calculate age difference (pre- and post-correction)
age_diff_df['age_diff_raw'] = age_diff_df['te_age_max'] - age_diff_df['te_age_min']
age_diff_df['age_diff_corrected'] = (
    age_diff_df['corrected_te_age_max'] - age_diff_df['corrected_te_age_min']
)
# %%
full_stats = pd.merge(identity_df, age_diff_df, on='block_id')
# %%
print(full_stats[['identity_pct', 'age_diff_raw','corrected_te_age_max']].corr())

#%%
bins = [-0.1]
genomesize = 3049315783
for i in range(0,100):
    bins.append(i+0.9)
full_stats['identity_pct_binned'] = pd.cut(full_stats['identity_pct'], bins=bins, labels=list(range(0,100)))
#%%
age_count_overall=full_stats.groupby(['corrected_te_age_max', 'identity_pct_binned']).count().reset_index()[['corrected_te_age_max','identity_pct_binned','block_id']].rename(columns={'block_id':'count'})
#%%
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_overall["corrected_te_age_max"].unique())  # Original x-labels
y_labels = sorted(age_count_overall["identity_pct_binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_overall.iterrows():
    x_idx = x_labels.index(row["corrected_te_age_max"])
    y_idx = y_labels.index(row["identity_pct_binned"])
    heatmap_array[y_idx, x_idx] = row["count"]
# Mask zero or negative values
masked_array = np.ma.masked_less_equal(heatmap_array, 0)

# Set color for masked values
cmap = plt.cm.viridis.copy()
cmap.set_bad(color='lightgray')  # Replace with desired hash-like color

# Use LogNorm
norm = mcolors.LogNorm(vmin=1e-4, vmax=np.nanmax(heatmap_array))

# Define log-normalization for colors (ignoring zero values)
norm = mcolors.Normalize(vmin=0, vmax=np.nanmax(heatmap_array))

# 🔹 Step 5: Plot using `pcolormesh()` with custom bin edges
fig, ax = plt.subplots(figsize=(8, 6))
cax = ax.pcolormesh(x_edges, y_edges, masked_array, cmap=cmap, norm=norm, shading="flat")

for x in x_edges:
    ax.axvline(x, color="white", linewidth=0.5, alpha=0.5)
for y in y_edges:
    ax.axhline(y, color="white", linewidth=0.5, alpha=0.5)

# Add colorbar
cbar = fig.colorbar(cax)
cbar.set_label("Count")

# 🔹 Step 6: Set x-ticks at bin edges
ax.set_xticks(x_edges)
ax.set_xticklabels([str(x) for x in x_labels] + ["43+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in y_tick_positions])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120
ax.set_ylim(50, 90)
ax.set_xlim(8.6, None)
ax.set_xlabel("tAge (MYA)")
ax.set_ylabel("%identity of LTR pair")
ax.set_title("THE1C LTR pair counts grouped by TEATIME and %identity of LTR pair")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig("tage_iden.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
x_raw = full_stats['corrected_te_age_max']
y_raw = full_stats['identity_pct']
#%%
def percentile_bounds(series, lower_pct=0.01, upper_pct=0.99):
    lower = series.quantile(lower_pct)
    upper = series.quantile(upper_pct)
    return lower, upper
x_min, x_max = percentile_bounds(x_raw, 0.01, 0.99)
y_min, y_max = percentile_bounds(y_raw, 0.01, 0.99)
print("x:",x_min,x_max,"y:",y_min,y_max)
# %%
x_min, x_max = percentile_bounds(x_raw, 0.10, 0.90)
y_min, y_max = percentile_bounds(y_raw, 0.10, 0.90)
print("x:",x_min,x_max,"y:",y_min,y_max)
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_overall["corrected_te_age_max"].unique())  # Original x-labels
y_labels = sorted(age_count_overall["identity_pct_binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_overall.iterrows():
    x_idx = x_labels.index(row["corrected_te_age_max"])
    y_idx = y_labels.index(row["identity_pct_binned"])
    heatmap_array[y_idx, x_idx] = row["count"]
# Calculate per-column mean and std ignoring NaNs
col_means = np.nanmean(heatmap_array, axis=0)
col_stds = np.nanstd(heatmap_array, axis=0)
col_stds[col_stds == 0] = 1  # avoid division by zero

# Z-center columns
heatmap_z = (heatmap_array - col_means) / col_stds

vmax = np.nanmax(np.abs(heatmap_z))
norm = mcolors.Normalize(vmin=-4, vmax=4)

# Mask NaNs and actual zeros
masked_z = np.ma.masked_where((heatmap_array == 0) | np.isnan(heatmap_array), heatmap_z)

# Set color for masked values
cmap = plt.cm.seismic.copy()
cmap.set_bad(color='lightgray') 

# 🔹 Step 5: Plot using `pcolormesh()` with custom bin edges
fig, ax = plt.subplots(figsize=(8, 6))
cax = ax.pcolormesh(x_edges, y_edges, masked_z, cmap=cmap, norm=norm, shading="flat")

for x in x_edges:
    ax.axvline(x, color="white", linewidth=0.5, alpha=0.5)
for y in y_edges:
    ax.axhline(y, color="white", linewidth=0.5, alpha=0.5)

# Add colorbar
cbar = fig.colorbar(cax)
cbar.set_label("Z-score")

# 🔹 Step 6: Set x-ticks at bin edges
ax.set_xticks(x_edges)
ax.set_xticklabels([str(x) for x in x_labels] + ["43+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in y_tick_positions])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120
ax.set_ylim(50, 90)
ax.set_xlim(8.6, None)
ax.set_xlabel("TEATIME age (MYA)")
ax.set_ylabel("%identity of LTR pair")
ax.set_title("Per column Z-scores of THE1C LTR pair counts grouped by TEATIME and %identity of LTR pair")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig("overall_te_heatmap_z.png", dpi=300, bbox_inches="tight")
plt.show()
#%%
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_overall["corrected_te_age_max"].unique())  # Original x-labels
y_labels = sorted(age_count_overall["identity_pct_binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_overall.iterrows():
    x_idx = x_labels.index(row["corrected_te_age_max"])
    y_idx = y_labels.index(row["identity_pct_binned"])
    heatmap_array[y_idx, x_idx] = row["count"]
# Calculate per-column mean and std ignoring NaNs
column_sums = np.nansum(heatmap_array, axis=0)
column_sums[column_sums == 0] = np.nan  # Avoid division by zero
heatmap_frac = heatmap_array / column_sums  # Broadcasting

# Normalize color range between 0 and 1 (fractions)
norm = mcolors.Normalize(vmin=0, vmax=0.15)

#norm = mcolors.Normalize(vmin=0, vmax=np.nanmax(heatmap_frac))  # or fixed vmax like 0.4

# Mask NaNs and actual zeros
masked_z = np.ma.masked_where((heatmap_frac == 0) | np.isnan(heatmap_frac), heatmap_frac)

# Set color for masked values
cmap = plt.cm.plasma.copy()
cmap.set_bad(color='lightgray') 

# 🔹 Step 5: Plot using `pcolormesh()` with custom bin edges
fig, ax = plt.subplots(figsize=(8, 6))
cax = ax.pcolormesh(x_edges, y_edges, masked_z, cmap=cmap, norm=norm, shading="flat")


for x in x_edges:
    ax.axvline(x, color="white", linewidth=0.5, alpha=0.5)
for y in y_edges:
    ax.axhline(y, color="white", linewidth=0.5, alpha=0.5)

# Add colorbar
cbar = fig.colorbar(cax)
cbar.set_label("Fraction of column total")

# 🔹 Step 6: Set x-ticks at bin edges
ax.set_xticks(x_edges)
ax.set_xticklabels([str(x) for x in x_labels] + ["43+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in y_tick_positions])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120
ax.set_ylim(50, 90)
ax.set_xlim(8.6, None)
ax.set_xlabel("TEATIME age (MYA)")
ax.set_ylabel("%identity of LTR pair")
ax.set_title("Per column Fraction of column total of THE1C LTR pair counts\ngrouped by TEATIME and %identity of LTR pair")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig("overall_te_heatmap_f.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
# Flatten the matrix and get indices of top 10 values
flat_indices = np.argpartition(heatmap_frac.flatten(), -10)[-10:]
top10_flat = flat_indices[np.argsort(heatmap_array.flatten()[flat_indices])[::-1]]

# Convert flat indices back to 2D coordinates
top10_coords = [np.unravel_index(idx, heatmap_frac.shape) for idx in top10_flat]

# Show value + coordinates
top10_values = [(coord, heatmap_frac[coord]) for coord in top10_coords]

# Print nicely
for (i, j), val in top10_values:
    print(f"Row {i} (y: {y_labels[i]}), Col {j} (x: {x_labels[j]}): Count = {val}")
# %%
n_unique_internal_id = ltr_age_id_coord_tbl['internal_id'].nunique()
block_counts = ltr_age_id_coord_tbl['block_id'].value_counts()
n_block_gt2 = (block_counts > 2).sum()
valid_blocks = ltr_age_id_coord_tbl.groupby(['block_id', 'int_type'])['internal_id'].nunique().unstack(fill_value=0)
valid_id_list = valid_blocks[(valid_blocks['aINT'] == 1) & (valid_blocks['bINT'] == 1)].index
n_eligible_blocks = len(valid_id_list)
print("Unique internal_id:", n_unique_internal_id)
print("block_id with >2 rows:", n_block_gt2)
print("Eligible block_id:", n_eligible_blocks)

# %%
n_total_id = ltr_age_id_coord_tbl['internal_id'].nunique()
id_counts = ltr_age_id_coord_tbl['block_id'].value_counts()
singleton_blocks = id_counts[id_counts == 1].index
n_singleton_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(singleton_blocks)]['internal_id'].nunique()
multi_blocks = id_counts[id_counts > 1].index
n_multi_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(multi_blocks)]['internal_id'].nunique()
# You already have:
valid_blocks = ltr_age_id_coord_tbl.groupby(['block_id', 'int_type'])['internal_id'].nunique().unstack(fill_value=0)
valid_id_list = valid_blocks[(valid_blocks['aINT'] == 1) & (valid_blocks['bINT'] == 1)].index
eligible_df = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(valid_id_list)]
n_eligible_id = eligible_df['internal_id'].nunique()
def has_mismatch(x):
    return x['te_age'].nunique() > 1

mismatch_blocks = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(multi_blocks)] \
    .groupby('block_id').filter(has_mismatch)['block_id'].unique()

n_mismatch_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(mismatch_blocks)]['internal_id'].nunique()
match_blocks = list(set(multi_blocks) - set(mismatch_blocks))
n_match_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(match_blocks)]['internal_id'].nunique()
eligible_mismatch_blocks = list(set(valid_id_list) & set(mismatch_blocks))
n_eligible_mismatch_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(eligible_mismatch_blocks)]['internal_id'].nunique()
eligible_match_blocks = list(set(valid_id_list) & set(match_blocks))
n_eligible_match_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(eligible_match_blocks)]['internal_id'].nunique()
summary = {
    'Total unique internal_id': n_total_id,
    'Singleton internal_id': f"{n_singleton_id} ({n_singleton_id/n_total_id:.1%})",
    'Multi-block internal_id': f"{n_multi_id} ({n_multi_id/n_total_id:.1%})",
    'Eligible internal_id': f"{n_eligible_id} ({n_eligible_id/n_total_id:.1%})",
    'Multi age mismatch': f"{n_mismatch_id} ({n_mismatch_id/n_total_id:.1%})",
    'Multi age match': f"{n_match_id} ({n_match_id/n_total_id:.1%})",
    'Eligible mismatch': f"{n_eligible_mismatch_id} ({n_eligible_mismatch_id/n_total_id:.1%})",
    'Eligible match': f"{n_eligible_match_id} ({n_eligible_match_id/n_total_id:.1%})",
}

import pandas as pd
summary_df = pd.DataFrame.from_dict(summary, orient='index', columns=['Count (Percentage)'])
print(summary_df)

# %%
# Length column (already done, but keeping for completeness)
ltr_age_id_coord_tbl['length'] = ltr_age_id_coord_tbl['genoEnd'] - ltr_age_id_coord_tbl['genoStart']

# Unique internal_id total
n_total_id = ltr_age_id_coord_tbl['internal_id'].nunique()

# Singleton and multi-block internal_id
block_counts = ltr_age_id_coord_tbl['block_id'].value_counts()
singleton_blocks = block_counts[block_counts == 1].index
multi_blocks = block_counts[block_counts > 1].index

n_singleton_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(singleton_blocks)]['internal_id'].nunique()
n_multi_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(multi_blocks)]['internal_id'].nunique()
valid_blocks = ltr_age_id_coord_tbl.groupby(['block_id', 'int_type'])['internal_id'].nunique().unstack(fill_value=0)
valid_id_list = valid_blocks[(valid_blocks['aINT'] == 1) & (valid_blocks['bINT'] == 1)].index
eligible_df = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(valid_id_list)]
n_eligible_id = eligible_df['internal_id'].nunique()
# Any block with >1 corrected_te_age value
def corrected_age_mismatch(g):
    return g['corrected_te_age'].nunique() > 1

mismatch_blocks = (
    ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(multi_blocks)]
    .groupby('block_id')
    .filter(corrected_age_mismatch)['block_id']
    .unique()
)

n_mismatch_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(mismatch_blocks)]['internal_id'].nunique()
match_blocks = list(set(multi_blocks) - set(mismatch_blocks))
n_match_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(match_blocks)]['internal_id'].nunique()
eligible_mismatch_blocks = list(set(valid_id_list) & set(mismatch_blocks))
n_eligible_mismatch_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(eligible_mismatch_blocks)]['internal_id'].nunique()
eligible_match_blocks = list(set(valid_id_list) & set(match_blocks))
n_eligible_match_id = ltr_age_id_coord_tbl[ltr_age_id_coord_tbl['block_id'].isin(eligible_match_blocks)]['internal_id'].nunique()
summary = {
    'Total unique internal_id': n_total_id,
    'Singleton internal_id': f"{n_singleton_id} ({n_singleton_id/n_total_id:.1%})",
    'Multi-block internal_id': f"{n_multi_id} ({n_multi_id/n_total_id:.1%})",
    'Eligible internal_id': f"{n_eligible_id} ({n_eligible_id/n_total_id:.1%})",
    'Multi corrected_age mismatch': f"{n_mismatch_id} ({n_mismatch_id/n_total_id:.1%})",
    'Multi corrected_age match': f"{n_match_id} ({n_match_id/n_total_id:.1%})",
    'Eligible corrected_age mismatch': f"{n_eligible_mismatch_id} ({n_eligible_mismatch_id/n_total_id:.1%})",
    'Eligible corrected_age match': f"{n_eligible_match_id} ({n_eligible_match_id/n_total_id:.1%})",
}

import pandas as pd
summary_df = pd.DataFrame.from_dict(summary, orient='index', columns=['Count (Percentage)'])
print(summary_df)

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
# Aggregate (mean y per x), also collect count as weight
agg = full_stats.groupby('corrected_te_age_max').agg(
    mean_idt=('identity_pct', 'mean'),
    count=('identity_pct', 'count')
).reset_index()

X_mean = agg['corrected_te_age_max'].values.reshape(-1, 1)
y_mean = agg['mean_idt'].values
weights = agg['count'].values


# Weighted linear regression
model = LinearRegression()
model.fit(X_mean, y_mean, sample_weight=weights)
y_pred = model.predict(X_mean)

# Weighted R²
def weighted_r2(y_true, y_pred, weights):
    ss_res = np.sum(weights * (y_true - y_pred) ** 2)
    ss_tot = np.sum(weights * (y_true - np.average(y_true, weights=weights)) ** 2)
    return 1 - ss_res / ss_tot

r2 = weighted_r2(y_mean, y_pred, weights)

# Effective sample size for weighted data
n_eff = (np.sum(weights))**2 / np.sum(weights**2)
p = 1  # number of predictors

# Adjusted R²
adj_r2 = 1 - (1 - r2) * (n_eff - 1) / (n_eff - p - 1)

slope = model.coef_[0]
intercept = model.intercept_

# Correlations (Spearman/Pearson)
pearson_corr, pearson_p = pearsonr(agg['corrected_te_age_max'], agg['mean_idt'])
spearman_corr, spearman_p = spearmanr(agg['corrected_te_age_max'], agg['mean_idt'])

# Plot
plt.figure(figsize=(6,6))
plt.plot(X_mean, y_pred, color='red', linewidth=2, label='Weighted fit on means')
plt.scatter(X_mean, y_mean, color='black', s=10, label='Mean identity perc per corrected tAge')
N_true = len(agg)
# Annotation text
textstr = (
    f"Weighted regression:\ny = {slope:.3f}x + {intercept:.3f}\n"
    f"Weighted R² = {r2:.3f}\n"
    f"Adjusted R² = {adj_r2:.3f}\n"
    f"N (raw TEs) = {N_true}, n = {len(agg)}, effective n = {n_eff:.1f}\n\n"
    f"Pearson r = {pearson_corr:.3f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_corr:.3f} (p = {spearman_p:.2e})"
)

plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
               verticalalignment='top', bbox=dict(boxstyle="square, pad=0.5", facecolor='white', alpha=0.7))

plt.ylim(50, 90)
plt.xlim(8.6, None)
plt.xlabel('tAge (MYA)', fontsize=12)
plt.ylabel('kAge (MYA)', fontsize=12)
plt.title('Aggregated Regression Fit and Raw Data Correlation')
plt.legend()
plt.tight_layout()
plt.savefig(f"tage_iden_fitw.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression

# Filter raw data for repName == 'LINE'
x_raw=age_count_overall['identity_pct_binned'].astype(int)
y_raw=age_count_overall['corrected_te_age_max'].astype(float)

# Correlations
pearson_r, pearson_p = pearsonr(x_raw, y_raw)
spearman_rho, spearman_p = spearmanr(x_raw, y_raw)

# Linear regression fit
X = x_raw.values.reshape(-1, 1)
model = LinearRegression()
model.fit(X, y_raw)
y_pred = model.predict(X)

# R² and adjusted R²
r2 = model.score(X, y_raw)
n = len(y_raw)
p = 1  # number of predictors
adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p - 1)

# Regression equation text
slope = model.coef_[0]
intercept = model.intercept_
equation_text = f"y = {slope:.3f}x + {intercept:.3f}"

# Plot
plt.figure(figsize=(6, 6))
plt.scatter(x_raw, y_raw, s=5, alpha=0.3, label='Data points')
plt.plot(x_raw, y_pred, color='red', linewidth=2, label='Linear fit')

# Stats annotation
textstr = (
    f"Regression equation:\n{equation_text}\n"
    f"R² = {r2:.3f}\n"
    f"Adjusted R² = {adj_r2:.3f}\n\n"
    f"Pearson r = {pearson_r:.3f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_rho:.3f} (p = {spearman_p:.2e})"
)
plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
               verticalalignment='top', bbox=dict(boxstyle="square,pad=0.5", facecolor='white', alpha=0.7))

plt.ylim(50, 90)
plt.xlim(8.6, None)
plt.xlabel('tAge (MYA)', fontdict={'fontsize': 12})
plt.ylabel('kAge (MYA)', fontdict={'fontsize': 12})
plt.title('Raw Data and Linear Regression Fit')
plt.legend()
plt.tight_layout()
plt.savefig(f"tage_fit.png", dpi=300, bbox_inches="tight")
plt.show()

# %%
