#%%
import numpy as np
import pandas as pd
#%%
repeatmasker_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'

repeatmasker_table = pd.read_csv(repeatmasker_filepath, sep='\t', index_col=0)
#%%
LTR_list=repeatmasker_table[repeatmasker_table['repClass'].str.contains('LTR', na=False)]['repName'].unique()
# %%
# Step 1: Convert to set for fast lookup
LTR_set = set(LTR_list)

# Step 2: Check which names have a corresponding -int entry
with_int = [ltr for ltr in LTR_set if f"{ltr}-int" in LTR_set]

# %%
with open("LTRs_with_internal.txt", "w") as f:
    for ltr in with_int:
        f.write(f"{ltr}\n")

# %%
# List of -int counterparts
int_list = [f"{ltr}-int" for ltr in with_int]

# Save to text file
with open("LTRs_internal_counterparts.txt", "w") as f:
    for item in int_list:
        f.write(f"{item}\n")

# %%
import numpy as np
import pandas as pd

import numpy as np
import pandas as pd
import os
import math

def process_internal_id_file(internal_id_filepath, repeatmasker_filepath,
                             ltr_name, per_ltr_sample_size=50, length_cutoff=300):
    internal_id_tbl = pd.read_csv(internal_id_filepath, sep='\t')

    # Filter for this LTR only
    internal_id_tbl = internal_id_tbl[internal_id_tbl['internal_id'].str.startswith(f"{ltr_name}_")].copy()

    # Extract block_id
    def extract_block_id(x):
        suffix = x[len(ltr_name) + 1:]
        rmsk_id = suffix.split('_')[0]
        return f"{ltr_name}_{rmsk_id}"

    internal_id_tbl['block_id'] = internal_id_tbl['internal_id'].apply(extract_block_id)

    # INT type
    pattern = r'_(aINT|bINT|nINT)_'
    internal_id_tbl['int_type'] = internal_id_tbl['internal_id'].str.extract(pattern)

    # RepeatMasker
    repeatmasker_table = pd.read_csv(repeatmasker_filepath, sep='\t', index_col=0)

    merged_tbl = internal_id_tbl.merge(
        repeatmasker_table[['genoName', 'genoStart', 'genoEnd', 'strand', 'repName']],
        left_on='rmsk_index', right_index=True
    )

    merged_tbl['length'] = merged_tbl['genoEnd'] - merged_tbl['genoStart']

    # One aINT and one bINT per block
    valid_blocks = (
        merged_tbl.groupby(['block_id', 'int_type'])['internal_id']
        .nunique().unstack(fill_value=0)
    )
    # Ensure 'aINT' and 'bINT' columns exist
    for col in ['aINT', 'bINT']:
        if col not in valid_blocks.columns:
            valid_blocks[col] = 0

    valid_id_list = valid_blocks[
        (valid_blocks['aINT'] == 1) & (valid_blocks['bINT'] == 1)
    ].index

    df_filtered = merged_tbl[merged_tbl['block_id'].isin(valid_id_list)]

    # Length check
    block_lengths_ok = df_filtered.groupby('block_id')['length'].apply(
        lambda x: (x >= length_cutoff).sum() == 2
    )
    valid_id_list2 = block_lengths_ok[block_lengths_ok].index
    df_filtered = df_filtered[df_filtered['block_id'].isin(valid_id_list2)]

    # Sample per-ltr limit
    valid_block_ids = df_filtered['block_id'].unique()
    n_block=len(valid_block_ids)
    print(f'found: {n_block}')
    n = min(per_ltr_sample_size, len(valid_block_ids))
    if n == 0:
        return pd.DataFrame()
    sampled_ids = np.random.choice(valid_block_ids, size=n, replace=False)
    df_sample = df_filtered[df_filtered['block_id'].isin(sampled_ids)].copy()
    df_sample['element_family'] = ltr_name  # Add family tag

    return df_sample


#%%
from math import ceil

total_sample_target = 3000
ltr_list = with_int
#per_ltr_sample_size = ceil(total_sample_target / len(ltr_list))
per_ltr_sample_size = 100
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
internal_id_base = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/internal_id/'

results = []
for ltr in ltr_list:
    internal_id_path = os.path.join(internal_id_base, f"{ltr}.internal_id.txt")
    if not os.path.exists(internal_id_path):
        print(f"⚠️ Missing: {internal_id_path}")
        continue

    print(f"Processing {ltr}...")
    df_sample = process_internal_id_file(
        internal_id_filepath=internal_id_path,
        repeatmasker_filepath=repeatmasker_filepath,
        ltr_name=ltr,
        per_ltr_sample_size=per_ltr_sample_size,
        length_cutoff=300
    )
    if not df_sample.empty:
        results.append(df_sample)
        print(f" → Got {df_sample['block_id'].nunique()} blocks")
    else:
        print(f" → No valid blocks")

# Combine and trim to exactly 1000 blocks if needed
combined = pd.concat(results, ignore_index=True)

# Keep only 1000 unique block_ids if oversampled
final_blocks = combined['block_id'].drop_duplicates().sample(n=min(3000, combined['block_id'].nunique()), random_state=42)
final_sample = combined[combined['block_id'].isin(final_blocks)]

print(f"\n✅ Final sample: {final_sample['block_id'].nunique()} blocks across {final_sample['element_family'].nunique()} families")

# Save if needed
final_sample.to_csv("final_LTR_sample_2000.txt", sep='\t', index=False)
internal_id_df=final_sample[['rmsk_index','tag','internal_id']]
internal_id_df.to_csv("final_LTR_sample_2000.internal_id.txt", sep='\t', index=False)
#%%
import os
from Bio import SeqIO
import pandas as pd

# Define directory where aligned files are
tmp_dir = "./tmp"

# Step 1: Gather valid aligned files (2 sequences, equal length)
valid_id_list2 = []
for fname in os.listdir(tmp_dir):
    if fname.endswith(".fasta.aligned"):
        path = os.path.join(tmp_dir, fname)
        try:
            records = list(SeqIO.parse(path, "fasta"))
            if len(records) == 2 and len(records[0].seq) == len(records[1].seq):
                block_id = fname.replace(".fasta.aligned", "")
                valid_id_list2.append(block_id)
        except Exception as e:
            print(f"Skipping {fname}: {e}")

print(f"✅ {len(valid_id_list2)} valid aligned blocks found")

# Step 2: Calculate % identity for each valid block
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
for block_id in valid_id_list2:
    aligned_path = os.path.join(tmp_dir, f"{block_id}.fasta.aligned")
    try:
        result[block_id] = calculate_identity_from_alignment(aligned_path)
    except Exception as e:
        print(f"Error in {block_id}: {e}")

# Step 4: Wrap into a DataFrame
identity_df = pd.DataFrame.from_dict(result, orient='index')
identity_df.index.name = 'block_id'
identity_df.reset_index(inplace=True)

# Now identity_df is ready!
identity_df.head()

# %%
import pandas as pd

# Load your TE age table
age_table_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/teatime_simple/final_LTR_sample_3000.teatime.txt'
age_df = pd.read_csv(age_table_path, sep='\t')  # adjust sep if needed
#%%
# Extract block_id from internal_id
def extract_block_id(internal_id):
    return "_".join(internal_id.split("_")[:2])

age_df['block_id'] = age_df['internal_id'].apply(extract_block_id)
# Identify block_ids with ANY NaN in te_age
bad_blocks = age_df[age_df['te_age'].isna()]['block_id'].unique()

# Drop entire blocks with NaN te_age
filtered_age_df = age_df[~age_df['block_id'].isin(bad_blocks)].copy()

# Group by block_id and calculate age difference
# We assume there are exactly 2 rows per block (aINT + bINT)
age_diff_df = (
    filtered_age_df.groupby('block_id')['te_age']
    .agg(['min', 'max'])
    .assign(age_diff=lambda x: x['max'] - x['min'])
    .reset_index()
)

# Show the result
age_diff_df.head()

# %%
full_stats=age_diff_df.merge(identity_df, on='block_id')
# %%
full_stats['identity_diff'] = 100 - full_stats['identity_pct']

# %%
import matplotlib.pyplot as plt

# Calculate % divergence
full_stats['identity_diff'] = 100 - full_stats['identity_pct']

plt.figure(figsize=(6, 6))
plt.scatter(
    full_stats['age_diff'],
    full_stats['identity_diff'],
    s=10, alpha=0.6, edgecolors='k'
)

plt.xlabel('Age Difference (LTR pair; MYA)')
plt.ylabel('Percent Divergence (100 - Identity%)')
plt.title('LTR Pair Divergence vs Age Difference')
plt.grid(True)
plt.tight_layout()
plt.savefig("tAge_diff_div.png", dpi=300, bbox_inches="tight")
plt.show()


# %%
# %%
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

# Pearson correlation
r, p = pearsonr(full_stats['age_diff'], full_stats['identity_diff'])

# Group by raw age_diff and compute mean identity
binned_df = (
    full_stats.groupby('age_diff')['identity_diff']
    .mean()
    .reset_index()
)

# Linear regression (on full data)
X = full_stats['age_diff'].values.reshape(-1, 1)
y = full_stats['identity_diff'].values
model = LinearRegression().fit(X, y)
x_line = np.linspace(X.min(), X.max(), 200).reshape(-1, 1)
y_line = model.predict(x_line)

# Plot
plt.figure(figsize=(6,6))
plt.scatter(
    full_stats['age_diff'],
    full_stats['identity_diff'],
    s=10, alpha=0.2, edgecolors='k'
)
plt.scatter(
    binned_df['age_diff'],
    binned_df['identity_diff'],
    s=10, c='red', label='Average divergence'
)
plt.plot(x_line, y_line, color='blue', linewidth=1.5, label='Regression line')

# Add text box with Pearson r and p
textstr = f'Pearson r = {r:.2f}; p-value = {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(
    0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props
)

plt.xlabel('Age Difference (LTR pair; MYA)')
plt.ylabel('Percent Divergence (100 - Identity%)')
plt.title('LTR Pair Divergence vs Age Difference')
plt.grid(True)
plt.legend()
plt.savefig("tAge_diff_div_avg_corr.png", dpi=300, bbox_inches="tight")
plt.tight_layout()
plt.show()

# %%
# Example block_id: "MLT1E1A_907097"
full_stats['subfamily'] = full_stats['block_id'].str.extract(r'^([^_]+(?:_[^_]+)*)')
from scipy.stats import pearsonr

def compute_subfam_corr(df):
    subfam_results = []
    for subfam, group in df.groupby('subfamily'):
        if len(group) >= 10:  # minimum data points for meaningful stats
            r, p = pearsonr(group['age_diff'], group['identity_diff'])
            subfam_results.append({'subfamily': subfam, 'r': r, 'p': p, 'n': len(group)})
    return pd.DataFrame(subfam_results)

subfam_corr_df = compute_subfam_corr(full_stats)
subfam_corr_df.sort_values('r', ascending=False, inplace=True)
#%%
import matplotlib.pyplot as plt

plt.figure(figsize=(10,5))
plt.bar(subfam_corr_df['subfamily'], subfam_corr_df['r'], color='skyblue')
plt.xticks(rotation=90)
plt.axhline(0, color='gray', linestyle='--')
plt.ylabel("Pearson r")
plt.title("Correlation between Age Difference and Divergence per Subfamily")
plt.tight_layout()
plt.show()


#%% REAL TRUE CORRECT
same_age_pairs = full_stats[full_stats['age_diff'] == 0].copy()
same_age_pairs['divergence'] = 100 - same_age_pairs['identity_pct']
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Sanity: assert min == max for this filtered group
assert (same_age_pairs['min'] == same_age_pairs['max']).all()

# Pearson correlation
r, p = pearsonr(same_age_pairs['min'], same_age_pairs['divergence'])

plt.figure(figsize=(6,6))
plt.scatter(
    same_age_pairs['min'],
    same_age_pairs['divergence'],
    s=10, alpha=0.3, edgecolors='k'
)

# Regression line (optional but nice)
import numpy as np
from sklearn.linear_model import LinearRegression

X = same_age_pairs['min'].values.reshape(-1,1)
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
# Compute average divergence per tAge
binned = (
    same_age_pairs.groupby('min')['divergence']
    .mean()
    .reset_index()
)

# Plot everything as before...
plt.figure(figsize=(6,6))

# Raw scatter
plt.scatter(
    same_age_pairs['min'],
    same_age_pairs['divergence'],
    s=10, alpha=0.3, edgecolors='k', label='Individual pairs'
)

# Average dots
plt.scatter(
    binned['min'],
    binned['divergence'],
    c='red', s=20, label='Average per tAge'
)

# Regression line
plt.plot(x_line, y_line, color='blue', label='Regression line')

# Pearson r annotation
textstr = f'Pearson r = {r:.2f}\nP = {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props)

# Labels and formatting
plt.xlabel("tAge (MYA)")
plt.ylabel("Percent Divergence (100 - Identity%)")
plt.title("LTR Pair Divergence vs. Shared tAge")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("tAge_div_avg_corr.png", dpi=300, bbox_inches="tight")
plt.show()

# %%
import pandas as pd
import matplotlib.pyplot as plt

age_table_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/LTR5_Hs.txt'
age_df = pd.read_csv(age_table_path, sep='\t')
print(age_df[['internal_id', 'te_age']].head())
filtered_age_df = age_df[age_df['te_age'].notna()]
plt.figure(figsize=(6,4))
plt.hist(filtered_age_df['te_age'], bins=50, color='salmon', edgecolor='black')
plt.xlabel('TE Age (MYA)')
plt.ylabel('Count')
plt.title('Histogram of TE Ages (segmental duplications filtered out)')
plt.tight_layout()
plt.show()

# %%
# Filter rows where age == 0
age_0_pairs = same_age_pairs[same_age_pairs['min'] == 0]

# Optional: sort by a relevant column, e.g. 'score'
age_0_pairs_sorted = age_0_pairs.sort_values(by='divergence', ascending=False)

# Show the top N (e.g., 10) entries
print(age_0_pairs_sorted.head(10))

# %%
print(age_0_pairs_sorted[['block_id','identity_diff']])
# %%
