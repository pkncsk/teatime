#%%
import pandas as pd
import matplotlib.pyplot as plt
from ma_mapper import utility

#%%
repeatmasker_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'

repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
LTR_list=repeatmasker_table[repeatmasker_table['repClass'].str.contains('LTR', na=False)]['repName'].unique()
# Step 1: Convert to set for fast lookup
LTR_set = set(LTR_list)
ltr_names = [name for name in LTR_list if not name.endswith('-int')]
#%%
ltr_only = repeatmasker_table[repeatmasker_table['repClass'].str.contains('LTR', na=False)].copy()
ltr_only = ltr_only[ltr_only['genoName'].str.match(r'chr[\dXY]+$')]
cols_to_keep = ['genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass']
ltr_only = ltr_only[cols_to_keep]
ltr_by_chr_strand = {
    (chrom, strand): df.sort_values('genoStart' if strand == '+' else 'genoEnd')
    for (chrom, strand), df in ltr_only.groupby(['genoName', 'strand'])
}
#%%
intact_ltrs_df=pd.read_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/ltr_pair_10_lookup_nonredundant.txt', sep='\t', index_col=0)
#%%
LTR_list=intact_ltrs_df['repName'].unique()
#%%
intact_ltr_internal=pd.read_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/ltr_pair_10_lookup_nonredundant.internal_id.txt', sep='\t', index_col=0)
merged_intact_ltr=intact_ltr_internal.merge(intact_ltrs_df, on='rmsk_index')
te_age_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/teatime_simple/ltr_pair_10_lookup_nonredundant.teatime.txt'
te_age_df=pd.read_csv(te_age_path, sep='\t')
intact_ltr_idt=merged_intact_ltr.merge(te_age_df, on = 'internal_id')
#%%
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
#%%
valid_ids = [
    name for name, group in grouped
    if group['internal_id'].nunique() == 2 and #group['te_age'].nunique() == 1 and 
    group['repName'].nunique() == 1
]
# Step 6: Filter original df_filtered by valid ltr_pair_ids
df_filtered = df_filtered[df_filtered['ltr_pair_id'].isin(valid_ids)]
#%%
df_filtered.to_csv('./ltr_pair_valid.txt', sep='\t')
#%%
import os
from Bio import SeqIO
tmp_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/tmp_intactpair'
def calculate_identity_from_alignment(aligned_fasta_path):
    records = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    seq1 = str(records[0].seq).upper()
    seq2 = str(records[1].seq).upper()

    aligned_length = min([len(seq1),len(seq2)])
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
#%%# 
ageiden_df=df_filtered.merge(identity_df, on='ltr_pair_id')
# %%
ageiden_df['pct_divergence'] = 100 -ageiden_df['identity_pct']
#%%
ageiden_df['te_length'] = ageiden_df['genoEnd']-ageiden_df['genoStart']
#%%
ageiden_df["min_te_length_in_pair"] = ageiden_df.groupby("ltr_pair_id")["te_length"].transform("min")
#%%
ageiden_df['identity_pct_by_min_length'] = ageiden_df['matches']/ageiden_df['min_te_length_in_pair']*100
#%%
ageiden_df['pct_divergence_by_min_length']=100-ageiden_df['identity_pct_by_min_length']
#%% 
#%% 
ageiden_df_long=ageiden_df[ageiden_df['min_te_length_in_pair']>=300]
# %% 
# Assuming your DataFrame is called `df`
ageiden_df_long["max_te_age_in_pair"] = ageiden_df_long.groupby("ltr_pair_id")["te_age"].transform("max")
#%% 
dups = ageiden_df_long[ageiden_df_long.duplicated(subset=['ltr_pair_id', 'te_age'], keep=False)]
#%% 
ageiden_slim=ageiden_df_long[['ltr_pair_id','max_te_age_in_pair','pct_divergence_by_min_length']]
ageiden_slim=ageiden_slim.drop_duplicates()
# %%
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Sanity: assert min == max for this filtered group
ageiden_slim = ageiden_slim[
    ageiden_slim.groupby('max_te_age_in_pair')['pct_divergence_by_min_length']
    .transform('count') >= 3
]

# Pearson correlation
r, p = pearsonr(ageiden_slim['max_te_age_in_pair'], ageiden_slim['pct_divergence_by_min_length'])

plt.figure(figsize=(6,6))


binned = (
    ageiden_slim.groupby('max_te_age_in_pair')['pct_divergence_by_min_length']
    .mean()
    .reset_index()
)
# Average dots
plt.scatter(
    binned['max_te_age_in_pair'],
    binned['pct_divergence_by_min_length'],
    c='black', s=10, label='Average per tAge'
)
# Regression line (optional but nice)
import numpy as np
from sklearn.linear_model import LinearRegression

X = ageiden_slim['max_te_age_in_pair'].values.reshape(-1,1)
y = ageiden_slim['pct_divergence_by_min_length'].values


model = LinearRegression().fit(X, y)
x_line = np.linspace(X.min(), X.max(), 100).reshape(-1,1)
y_line = model.predict(x_line)
plt.plot(x_line, y_line, color='red', label='Regression line')

# Annotate correlation
textstr = f'Pearson r = {r:.2f}; p-value = {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props)

plt.xlabel("tAge (MYA)",fontdict={'fontsize':12})
plt.ylabel("Percent Divergence (100 - Identity%)",fontdict={'fontsize':12})
plt.title("LTR Pair Divergence vs. tAge",fontdict={'fontsize':12})
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('./tAge_div_nofilter.png',dpi=300)
plt.show()
#%%
import numpy as np
import matplotlib.pyplot as plt

def remove_outliers(arr):
    q1 = np.percentile(arr, 25)
    q3 = np.percentile(arr, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return [x for x in arr if lower_bound <= x <= upper_bound]

# Group and remove outliers
grouped = ageiden_slim.groupby('max_te_age_in_pair')['pct_divergence_by_min_length'].apply(list)

ages = []
cleaned_data = []

for age, values in grouped.items():
    filtered_values = remove_outliers(values)
    if len(filtered_values) > 0:
        ages.append(age)
        cleaned_data.append(filtered_values)

# Sort by age
sorted_pairs = sorted(zip(ages, cleaned_data))
ages_sorted, data_sorted = zip(*sorted_pairs)

# Plot boxplot without outliers
plt.figure(figsize=(6, 6))
box=plt.boxplot(data_sorted, positions=ages_sorted, widths=1.5, patch_artist=True, showfliers=False)
for b in box['boxes']:
    b.set(facecolor='white', edgecolor='black')
plt.xlabel('tAge (MYA)',fontdict={'fontsize':12})
plt.ylabel('Percent Divergence (100 - Identity%)',fontdict={'fontsize':12})
plt.title('Boxplot of % Divergence per tAge',fontdict={'fontsize':12})
plt.grid(True, axis='y')
plt.xticks(rotation=45) 
plt.tight_layout()
plt.savefig('./box_tAge_div_nofilter.png',dpi=300)
plt.show()

# %%
ageiden_fam=ageiden_df_long[['ltr_pair_id','max_te_age_in_pair','pct_divergence_by_min_length','repName']]
ageiden_fam=ageiden_fam.drop_duplicates()
# %%
# Step 1: Group by repName and te_age, then count occurrences
mode_counts = ageiden_fam.groupby(['repName', 'max_te_age_in_pair']).size().reset_index(name='count')
# Step 2: For each repName, find the te_age with the highest count
mode_per_fam = mode_counts.sort_values('count', ascending=False).drop_duplicates('repName')
# Step 3: Rename for clarity
mode_per_fam = mode_per_fam.rename(columns={'max_te_age_in_pair': 'common_te_age'})
# Step 4 (optional): Reset index
mode_per_fam = mode_per_fam[['repName', 'common_te_age']]

# %%
# Merge df with mode_per_rep on 'repName'
df_merged = ageiden_fam.merge(mode_per_fam, on='repName', how='left')

# Filter where max_te_age_in_pair matches the mode (common_te_age)
ageiden_common_age = df_merged[df_merged['max_te_age_in_pair'] == df_merged['common_te_age']]
# %% 
# %%
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Sanity: assert min == max for this filtered group
# Step 1: Filter groups with fewer than 3 data points

# Pearson correlation
r, p = pearsonr(ageiden_common_age['max_te_age_in_pair'], ageiden_common_age['pct_divergence_by_min_length'])

plt.figure(figsize=(6,6))


binned = (
    ageiden_common_age.groupby('max_te_age_in_pair')['pct_divergence_by_min_length']
    .mean()
    .reset_index()
)
# Average dots
plt.scatter(
    binned['max_te_age_in_pair'],
    binned['pct_divergence_by_min_length'],
    c='black', s=10, label='Average per tAge'
)
# Regression line (optional but nice)
import numpy as np
from sklearn.linear_model import LinearRegression

X = ageiden_common_age['max_te_age_in_pair'].values.reshape(-1,1)
y = ageiden_common_age['pct_divergence_by_min_length'].values


model = LinearRegression().fit(X, y)
x_line = np.linspace(X.min(), X.max(), 100).reshape(-1,1)
y_line = model.predict(x_line)
plt.plot(x_line, y_line, color='red', label='Regression line')

# Annotate correlation
textstr = f'Pearson r = {r:.2f}; p-value = {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props)

plt.xlabel("tAge (MYA)",fontdict={'fontsize':12})
plt.ylabel("Percent Divergence (100 - Identity%)",fontdict={'fontsize':12})
plt.title("LTR Pair Divergence vs. Shared tAge",fontdict={'fontsize':12})
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('./tAge_div_common.png',dpi=300)
plt.show()
# %%

ageiden_common_age = ageiden_common_age[
    ageiden_common_age.groupby('max_te_age_in_pair')['pct_divergence_by_min_length']
    .transform('count') >= 5
]
#%%
# %%
import numpy as np
import matplotlib.pyplot as plt

def remove_outliers(arr):
    q1 = np.percentile(arr, 25)
    q3 = np.percentile(arr, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return [x for x in arr if lower_bound <= x <= upper_bound]

# Group and remove outliers
grouped = ageiden_common_age.groupby('common_te_age')['pct_divergence_by_min_length'].apply(list)

ages = []
cleaned_data = []

for age, values in grouped.items():
    filtered_values = remove_outliers(values)
    if len(filtered_values) > 0:
        ages.append(age)
        cleaned_data.append(filtered_values)

# Sort by age
sorted_pairs = sorted(zip(ages, cleaned_data))
ages_sorted, data_sorted = zip(*sorted_pairs)

# Plot boxplot without outliers
plt.figure(figsize=(6, 6))
plt.boxplot(data_sorted, positions=ages_sorted, widths=1.5, patch_artist=True, showfliers=False)

plt.xlabel('tAge (MYA)',fontdict={'fontsize':12})
plt.ylabel('Percent Divergence',fontdict={'fontsize':12})
plt.title('Boxplot of % Divergence per tAge',fontdict={'fontsize':12})
plt.grid(True, axis='y')

plt.xticks(rotation=45) 
plt.tight_layout()
plt.savefig('./box_tAge_div_common.png',dpi=300)
plt.show()

# %%
ageiden_common_age[ageiden_common_age['common_te_age']==15.2]
# %%
ageiden_df_long[ageiden_df_long['ltr_pair_id']=='MER50_2125']
# %%

#%%
ageiden_common_age = ageiden_common_age[
    ageiden_common_age.groupby('max_te_age_in_pair')['pct_divergence_by_min_length']
    .transform('count') >= 5
]

# %%
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Sanity: assert min == max for this filtered group
# Step 1: Filter groups with fewer than 3 data points

# Pearson correlation
r, p = pearsonr(ageiden_common_age['max_te_age_in_pair'], ageiden_common_age['pct_divergence_by_min_length'])

plt.figure(figsize=(6,6))


binned = (
    ageiden_common_age.groupby('max_te_age_in_pair')['pct_divergence_by_min_length']
    .mean()
    .reset_index()
)
# Average dots
plt.scatter(
    binned['max_te_age_in_pair'],
    binned['pct_divergence_by_min_length'],
    c='black', s=10, label='Average per tAge'
)
# Regression line (optional but nice)
import numpy as np
from sklearn.linear_model import LinearRegression

X = ageiden_common_age['max_te_age_in_pair'].values.reshape(-1,1)
y = ageiden_common_age['pct_divergence_by_min_length'].values


model = LinearRegression().fit(X, y)
x_line = np.linspace(X.min(), X.max(), 100).reshape(-1,1)
y_line = model.predict(x_line)
plt.plot(x_line, y_line, color='red', label='Regression line')

# Annotate correlation
textstr = f'Pearson r = {r:.2f}; p-value = {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props)

plt.xlabel("tAge (MYA)",fontdict={'fontsize':12})
plt.ylabel("Percent Divergence (100 - Identity%)",fontdict={'fontsize':12})
plt.title("LTR Pair Divergence vs. Shared tAge",fontdict={'fontsize':12})
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('./tAge_div_common_filtered.png',dpi=300)
plt.show()
# %%

ageiden_common_age = ageiden_common_age[
    ageiden_common_age.groupby('max_te_age_in_pair')['pct_divergence_by_min_length']
    .transform('count') >= 5
]
#%%
# %%
import numpy as np
import matplotlib.pyplot as plt

def remove_outliers(arr):
    q1 = np.percentile(arr, 25)
    q3 = np.percentile(arr, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return [x for x in arr if lower_bound <= x <= upper_bound]

# Group and remove outliers
grouped = ageiden_common_age.groupby('common_te_age')['pct_divergence_by_min_length'].apply(list)

ages = []
cleaned_data = []

for age, values in grouped.items():
    filtered_values = remove_outliers(values)
    if len(filtered_values) > 0:
        ages.append(age)
        cleaned_data.append(filtered_values)

# Sort by age
sorted_pairs = sorted(zip(ages, cleaned_data))
ages_sorted, data_sorted = zip(*sorted_pairs)

# Plot boxplot without outliers
plt.figure(figsize=(6, 6))
box=plt.boxplot(data_sorted, positions=ages_sorted, widths=1.5, patch_artist=True, showfliers=False)
for b in box['boxes']:
    b.set(facecolor='white', edgecolor='black')
plt.plot(x_line, y_line, color='red', label='Regression line')

# Annotate correlation
textstr = f'Pearson r = {r:.2f}; p-value = {p:.1e}'
props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
    fontsize=10, verticalalignment='top', bbox=props)

plt.xlabel('tAge (MYA)',fontdict={'fontsize':12})
plt.ylabel('Percent Divergence',fontdict={'fontsize':12})
plt.title('Boxplot of % Divergence per tAge',fontdict={'fontsize':12})
plt.grid(True, axis='y')

plt.xticks(rotation=45) 
plt.tight_layout()
plt.savefig('./box_tAge_div_common_filtered.png',dpi=300)
plt.show()
# %%
#Total (with identity): 16210
#Total (with alignment>300 bp): 15,552 total 7776 pairs non shorts alignment
# LTR with same ages: 11714 (5857 pairs) - 1919 pairs rcovered
#*second run removed TE with ungapped length < 300 bp:11790 rows -5895 pairs
#*second run same age 9236 (4618 pairs) - 1277 pairs recovered
#%%
import matplotlib.pyplot as plt
import numpy as np

# Raw counts
same_age = [5857, 7776]
diff_age = [1919, 0]
total = [sum(x) for x in zip(same_age, diff_age)]

# Calculate percentages
same_pct = [s / t * 100 for s, t in zip(same_age, total)]
diff_pct = [d / t * 100 for d, t in zip(diff_age, total)]

labels = ['Before', 'After']
x = np.arange(len(labels))  # [0, 1] positions
width = 0.35  # bar width

# Plot
fig, ax = plt.subplots(figsize=(7, 6))

bars1 = ax.bar(x - width/2, same_pct, width, color='white', edgecolor='black', label='Same Age')
bars2 = ax.bar(x + width/2, diff_pct, width, color='black', label='Different Age')

# Annotate with raw counts
for i, (b1, b2) in enumerate(zip(bars1, bars2)):
    ax.text(b1.get_x() + b1.get_width()/2, b1.get_height() + 1, f"{same_age[i]}", 
            ha='center', va='bottom', fontsize=10)
    ax.text(b2.get_x() + b2.get_width()/2, b2.get_height() + 1, f"{diff_age[i]}", 
            ha='center', va='bottom', fontsize=10)

# Axis labels and formatting
ax.set_ylabel('Percentage of LTR Pairs',fontdict={'fontsize':12})
ax.set_title('LTR Pair Age Match Before and After correction',fontdict={'fontsize':12})
ax.set_xticks(x)
ax.set_xticklabels(labels,fontdict={'fontsize':12})
ax.set_ylim(0, 110)
ax.legend()
#ax.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig('./before_after.png',dpi=300)
plt.show()


# %%
