#%%
import pandas as pd
import matplotlib.pyplot as plt
from ma_mapper import utility

#%%
repeatmasker_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'

repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
#%%
LTR_list=repeatmasker_table[repeatmasker_table['repClass'].str.contains('LTR', na=False)]['repName'].unique()
# %%
# Step 1: Convert to set for fast lookup
LTR_set = set(LTR_list)
#%%
ltr_names = [name for name in LTR_list if not name.endswith('-int')]
#%%
ltr_only = repeatmasker_table[repeatmasker_table['repClass'].str.contains('LTR', na=False)].copy()
ltr_only = ltr_only[ltr_only['genoName'].str.match(r'chr[\dXY]+$')]
cols_to_keep = ['genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass']
ltr_only = ltr_only[cols_to_keep]
#%%
ltr_by_chr_strand = {
    (chrom, strand): df.sort_values('genoStart' if strand == '+' else 'genoEnd')
    for (chrom, strand), df in ltr_only.groupby(['genoName', 'strand'])
}
# %%
max_gap = 15000
max_span = 15000
min_span = 1000
intact_ltrs = []

for subfam in ltr_names:
    print(subfam)
    subfam_table=repeatmasker_table[repeatmasker_table['repName']==subfam]
    subfam_table = subfam_table[subfam_table['genoName'].str.match(r'chr[\dXY]+$')]
    cols_to_keep = ['genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass']
    subfam_table = subfam_table[cols_to_keep]
    consumed_indices = set()
    for run_id in range(len(subfam_table)):
        print(run_id)
        row_idx = subfam_table.iloc[run_id].name
        print(row_idx)
        current_row = ltr_only[ltr_only.index==int(row_idx)]

        # Skip if already part of an earlier trio
        if row_idx in consumed_indices:
            continue

        chrom = current_row['genoName'].values[0]
        strand = current_row['strand'].values[0]
        anchor_pos = current_row['genoStart'].values[0] if strand == '+' else current_row['genoEnd'].values[0]
        ltr_class = current_row['repClass'].values[0]

        try:
            candidates = ltr_by_chr_strand[(chrom, strand)]
        except KeyError:
            continue  # no LTRs on this chrom/strand

        # Get nearby entries within 10kb window
        if strand == '+':
            window_df = candidates[
                (candidates['repClass'] == ltr_class) &
                (candidates['genoStart'] >= anchor_pos) &
                (candidates['genoStart'] <= anchor_pos + max_gap)
            ].copy()
            window_df = window_df.sort_values('genoStart')
        else:
            window_df = candidates[
                (candidates['repClass'] == ltr_class) &
                (candidates['genoEnd'] <= anchor_pos) &
                (candidates['genoEnd'] >= anchor_pos - max_gap)
            ].copy()
            window_df = window_df.sort_values('genoEnd', ascending=False)
        if len(window_df) <3:
            continue

        # Try to find LTR -int LTR pattern
        for i in range(len(window_df) - 2):
            a, b, c = window_df.iloc[i], window_df.iloc[i+1], window_df.iloc[i+2]
            idx_a, idx_b, idx_c = a.name, b.name, c.name

            if idx_a in consumed_indices or idx_b in consumed_indices or idx_c in consumed_indices:
                continue

            # Check repName pattern: LTR - -int - LTR
            if not (
                not a['repName'].endswith('-int') and
                b['repName'].endswith('-int') and
                not c['repName'].endswith('-int')
            ):
                continue

            # Strand-aware span
            start = min(a['genoStart'], b['genoStart'], c['genoStart'])
            end   = max(a['genoEnd'],   b['genoEnd'],   c['genoEnd'])
            span  = end - start

            if span > max_span or span < min_span:
                continue

            # Passed all: record and mark as used
            trio = window_df.loc[[idx_a, idx_b, idx_c]]
            ltr_pair_id = f"{a['repName']}_{run_id}"
            trio['is_int'] = trio['repName'].str.endswith('-int')
            trio['ltr_pair_id'] = ltr_pair_id
            trio['rmsk_index'] = trio.index
            intact_ltrs.append(trio)
            consumed_indices.update([idx_a, idx_b, idx_c])
            break  # move to next anchor LTR
# Final output
intact_ltrs_df = pd.concat(intact_ltrs).sort_values(['genoName', 'genoStart']).reset_index(drop=True)

# %%
#intact_ltrs_df.to_csv('./ltr_pair_10_lookup_nonredundant.txt', sep='\t')
intact_ltrs_df=pd.read_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/ltr_pair_10_lookup_nonredundant.txt', sep='\t', index_col=0)
#%%
LTR_list=intact_ltrs_df['repName'].unique()
#%%
with open("LTRs_subfam_list.txt", "w") as f:
    for ltr in LTR_list:
        f.write(f"{ltr}\n")
# %%
internal_id_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/internal_id_simple'
internal_id_list = []
total=0
for ltr in LTR_list:
    print(ltr)
    rmskindex_list=intact_ltrs_df[intact_ltrs_df["repName"]==ltr]['rmsk_index'].values
    input_filepath = f'{internal_id_dir}/{ltr}.internal_id.txt'
    internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
    internal_id_subset=internal_id_tbl[internal_id_tbl['rmsk_index'].isin(rmskindex_list)]
    internal_id_list.append(internal_id_subset)
    total += internal_id_subset.shape[0]
    print("total:",total)
intact_ltr_internal=pd.concat(internal_id_list)

#%%
intact_ltr_internal.to_csv('./ltr_pair_10_lookup_nonredundant.internal_id.txt', sep='\t')
# %%
merged_intact_ltr=intact_ltr_internal.merge(intact_ltrs_df, on='rmsk_index')
# %%
te_age_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/teatime_simple/ltr_pair_10_lookup_nonredundant.teatime.txt'
te_age_df=pd.read_csv(te_age_path, sep='\t')
#%%
intact_ltr_idt=merged_intact_ltr.merge(te_age_df, on = 'internal_id')
# %%
intact_ltr_idt
# %%
import pandas as pd

# Assume df is your DataFrame

# Step 1: Remove rows where is_int is True
df_filtered = intact_ltr_idt[~intact_ltr_idt['is_int']]

# Step 2: Remove rows where te_age is NaN
df_filtered = df_filtered.dropna(subset=['te_age'])

# Step 3–5: Group and filter
grouped = df_filtered.groupby('ltr_pair_id')

# Get groups with 2 rows and equal te_age
valid_ids = [
    name for name, group in grouped
    if group['internal_id'].nunique() == 2 and #group['te_age'].nunique() == 1 and 
    group['repName'].nunique() == 1
]

# Step 6: Filter original df_filtered by valid ltr_pair_ids
df_filtered = df_filtered[df_filtered['ltr_pair_id'].isin(valid_ids)]

# %%
df_filtered.to_csv('./ltr_pair_valid.txt', sep='\t')
# %%
df_filtered
# %%
import os
from Bio import SeqIO
tmp_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/tmp_intactpair'
def calculate_identity_from_alignment(aligned_fasta_path):
    records = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    seq1 = str(records[0].seq).upper()
    seq2 = str(records[1].seq).upper()

    aligned_length = len(seq1)
    matches = sum((a == b) and (a != '-') for a, b in zip(seq1, seq2))
    gaps = sum((a == '-' or b == '-') for a, b in zip(seq1, seq2))

    return {
        'aligned_length': aligned_length,
        'matches': matches,
        'identity_pct': matches / aligned_length * 100,
        'gaps': gaps,
        'gap_pct': gaps / aligned_length * 100
    }

# Step 3: Run on all
result = {}
for block_id in valid_ids:
    
    aligned_path = os.path.join(tmp_dir, f"{block_id}.fasta.aligned")
    print(aligned_path)
    try:
        result[block_id] = calculate_identity_from_alignment(aligned_path)
    except Exception as e:
        print(f"Error in {block_id}: {e}")

# Step 4: Wrap into a DataFrame
identity_df = pd.DataFrame.from_dict(result, orient='index')
identity_df.index.name = 'ltr_pair_id'
identity_df.reset_index(inplace=True)

# Now identity_df is ready!
identity_df.head()
# %%
age_diff_df = (
    df_filtered.groupby('ltr_pair_id')['te_age']
    .agg(['min', 'max'])
    .assign(age_diff=lambda x: x['max'] - x['min'])
    .reset_index()
)
#%%
# Step 1: Keep only groups of size 2
grouped = df_filtered.groupby('ltr_pair_id')
valid_groups = grouped.filter(lambda g: len(g) == 2)

# Step 2: Build the summary with internal_ids and age_diff
age_diff_df = (
    valid_groups
    .groupby('ltr_pair_id')
    .apply(lambda g: pd.Series({
        'te_age_min': g['te_age'].min(),
        'te_age_max': g['te_age'].max(),
        'age_diff': g['te_age'].max() - g['te_age'].min(),
        'internal_id1': g['internal_id'].iloc[0],
        'internal_id2': g['internal_id'].iloc[1],
    }))
    .reset_index()
)
full_stats=age_diff_df.merge(identity_df, on='ltr_pair_id')
# %%
full_stats['identity_diff'] = 100 - full_stats['identity_pct']
# %%
# %%
import matplotlib.pyplot as plt

# Calculate % divergence
full_stats['identity_diff'] = 100 - full_stats['identity_pct']
#%%
same_age_pairs = full_stats[full_stats['age_diff'] == 0].copy()
same_age_pairs = same_age_pairs[same_age_pairs['aligned_length']>=300]
same_age_pairs['divergence'] = 100 - same_age_pairs['identity_pct']
#%%
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Sanity: assert min == max for this filtered group
assert (same_age_pairs['te_age_min'] == same_age_pairs['te_age_max']).all()

# Pearson correlation
r, p = pearsonr(same_age_pairs['te_age_min'], same_age_pairs['divergence'])

plt.figure(figsize=(6,6))
plt.scatter(
    same_age_pairs['te_age_min'],
    same_age_pairs['divergence'],
    s=10, alpha=0.3, edgecolors='k'
)

# Regression line (optional but nice)
import numpy as np
from sklearn.linear_model import LinearRegression

X = same_age_pairs['te_age_min'].values.reshape(-1,1)
y = same_age_pairs['divergence'].values


model = LinearRegression().fit(X, y)
x_line = np.linspace(X.min(), X.max(), 100).reshape(-1,1)
y_line = model.predict(x_line)
plt.plot(x_line, y_line, color='red', label='Regression line')

# Annotate correlation
textstr = f'Pearson r = {r:.2f}\nP = {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props)

plt.xlabel("tAge (MYA)")
plt.ylabel("Percent Divergence (100 - Identity%)")
plt.title("LTR Pair Divergence vs. Shared tAge")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.show()

# %%
same_age_pairs
# %%
# Compute average divergence per tAge
binned = (
    same_age_pairs.groupby('te_age_min')['divergence']
    .mean()
    .reset_index()
)

# Plot everything as before...
plt.figure(figsize=(6,6))

# Average dots
plt.scatter(
    binned['te_age_min'],
    binned['divergence'],
    c='black', s=20, label='Average per tAge'
)

# Regression line
plt.plot(x_line, y_line, color='red', label='Regression line')

# Pearson r annotation
textstr = f'Pearson r = {r:.2f} p-value= {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props)

# Labels and formatting
plt.xlabel("tAge (MYA)", fontdict={'fontsize':10})
plt.ylabel("mean Percent Divergence (100 - Identity%)", fontdict={'fontsize':10})
plt.title("LTR Pair Divergence vs. Shared tAge", fontdict={'fontsize':12})
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("tAge_div_avg_corr.png", dpi=300, bbox_inches="tight")
plt.show()
# %%


#%%
# %%
# Filter rows where age == 0
age_0_pairs = same_age_pairs[same_age_pairs['min'] == 0]

# Optional: sort by a relevant column, e.g. 'score'
age_0_pairs_sorted = age_0_pairs.sort_values(by='divergence', ascending=False)

# Show the top N (e.g., 10) entries
print(age_0_pairs_sorted.head(10))

# %%
print(age_0_pairs_sorted[['block_id','identity_diff']])

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
bins = [-1]
for i in range(0,100,10):
    bins.append(i+10)
same_age_pairs['binned_div'] = pd.cut(same_age_pairs['divergence'], bins=bins, labels=list(range(0,100,10)))

counted = same_age_pairs.groupby(['te_age_min', 'binned_div']).count().reset_index()[['te_age_min','binned_div','ltr_pair_id']].rename(columns={'ltr_pair_id':'count'})
#%%
# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(counted["te_age_min"].unique())  # Original x-labels
y_labels = sorted(counted["binned_div"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in counted.iterrows():
    x_idx = x_labels.index(row["te_age_min"])
    y_idx = y_labels.index(row["binned_div"])
    heatmap_array[y_idx, x_idx] = row["count"]
# Mask zero or negative values
masked_array = np.ma.masked_less_equal(heatmap_array, 0)

# Set color for masked values
cmap = plt.cm.viridis.copy()
cmap.set_bad(color='lightgray')  # Replace with desired hash-like color

# Use LogNorm
norm = mcolors.LogNorm(vmin=1e-4, vmax=np.nanmax(heatmap_array))

# Define log-normalization for colors (ignoring zero values)
norm = mcolors.LogNorm(vmin=max(1, np.nanmin(heatmap_array)), vmax=1000)

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
ax.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, 100, 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in np.arange(len(y_tick_positions))])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120


ax.set_xlabel("tAge (MYA)", fontdict={'fontsize':12})
ax.set_ylabel("% LTR divergence",fontdict={'fontsize':12})
ax.set_title(f"LTR Pair Divergence vs. Shared tAge")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig(f"tAge_div_te_heatmap10-2.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
bins = [-1]
for i in range(0,100):
    bins.append(i+0.99)
same_age_pairs['binned_div'] = pd.cut(same_age_pairs['divergence'], bins=bins, labels=list(range(0,100)))

counted = same_age_pairs.groupby(['min', 'binned_div']).count().reset_index()[['min','binned_div','ltr_pair_id']].rename(columns={'ltr_pair_id':'count'})
#%%
# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(counted["min"].unique())  # Original x-labels
y_labels = sorted(counted["binned_div"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in counted.iterrows():
    x_idx = x_labels.index(row["min"])
    y_idx = y_labels.index(row["binned_div"])
    heatmap_array[y_idx, x_idx] = row["count"]
# Mask zero or negative values
masked_array = np.ma.masked_less_equal(heatmap_array, 0)

# Set color for masked values
cmap = plt.cm.viridis.copy()
cmap.set_bad(color='lightgray')  # Replace with desired hash-like color

# Use LogNorm
norm = mcolors.LogNorm(vmin=1e-4, vmax=np.nanmax(heatmap_array))

# Define log-normalization for colors (ignoring zero values)
norm = mcolors.LogNorm(vmin=max(1, np.nanmin(heatmap_array)), vmax=100)

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
ax.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in np.arange(len(y_tick_positions))])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120


ax.set_xlabel("tAge (MYA)", fontdict={'fontsize':12})
ax.set_ylabel("% LTR divergence",fontdict={'fontsize':12})
ax.set_title(f"LTR Pair Divergence vs. Shared tAge")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig(f"tAge_div_te_heatmap100-2.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
intact_ltrs_df
# %%
intact_ltrs_df[intact_ltrs_df['repName']=='THE1C']
# %%
repeatmasker_table[repeatmasker_table['repName']=='LTR5_Hs']
# %%
nofilter_iden=df_filtered.merge(identity_df, on='ltr_pair_id')
#%%
nofilter_iden['divergence'] = 100- nofilter_iden['identity_pct']
#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming your DataFrame is named df and contains columns:
# 'repName', 'te_age', 'divergence'

# 1. Calculate mean te_age for sorting
rep_counts = nofilter_iden['repName'].value_counts()
valid_repnames = rep_counts[rep_counts >= 5].index
df_filtered = nofilter_iden[nofilter_iden['repName'].isin(valid_repnames)]
mean_age = df_filtered.groupby('repName')['te_age'].mean().sort_values()

# 2. Set the category order for repName
df_filtered['repName'] = pd.Categorical(df_filtered['repName'], categories=mean_age.index, ordered=True)

# 3. Create the boxplot
plt.figure(figsize=(12, 6))
sns.boxplot(data=df_filtered, x='repName', y='divergence',  showfliers=False)

# 4. Beautify
plt.xticks(rotation=90)
plt.title('Divergence by repName (ordered by mean TE age)')
plt.tight_layout()
plt.show()


# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
nofilter_iden2=nofilter_iden[nofilter_iden['repName']=='LTR5_Hs']
# Step 1: Filter repNames with ≥ 50 entries
rep_counts = nofilter_iden2['repName'].value_counts()
valid_repnames = rep_counts[rep_counts >= 5].index
df_filtered = nofilter_iden2[nofilter_iden2['repName'].isin(valid_repnames)].copy()

# Step 2: Compute median te_age for each repName
rep_median_age = df_filtered.groupby('repName')['te_age'].median()

# Step 3: Map the median te_age back to dataframe
df_filtered['rep_median_age'] = df_filtered['repName'].map(rep_median_age)

# Step 4: Remove outliers (IQR method) within each repName
def remove_outliers_iqr(group):
    q1 = group['divergence'].quantile(0.25)
    q3 = group['divergence'].quantile(0.75)
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    return group[(group['divergence'] >= lower) & (group['divergence'] <= upper)]

df_no_outliers = df_filtered.groupby('repName', group_keys=False).apply(remove_outliers_iqr)

# Step 5: Sort repName categories by median te_age
rep_order = df_no_outliers.groupby('repName')['rep_median_age'].first().sort_values().index
df_no_outliers['repName'] = pd.Categorical(df_no_outliers['repName'], categories=rep_order, ordered=True)

# Step 6: Plot
plt.figure(figsize=(14, 6))
sns.stripplot(data=df_no_outliers, x='repName', y='divergence', jitter=True)
plt.xticks(rotation=90)
plt.title('Divergence by Subfamily (repName), ordered by median TE age')
plt.tight_layout()
plt.show()


# %%
# Replace te_age of each ltr_pair_id with its maximum te_age
df_filtered['te_age_max'] = df_filtered.groupby('ltr_pair_id')['te_age'].transform('max')

# %%
identity_df.merge(df_filtered, on='ltr_pair_id')

# %%
mean_age_by_repName = nofilter_iden.groupby('repName')['te_age'].mean().reset_index()
# %%
nofilter_iden
# %%
modes = df_filtered.groupby('repName')['te_age'].agg(pd.Series.mode)
# %%
nofilter_iden[nofilter_iden['repName']=='LTR5_Hs']
# %%
