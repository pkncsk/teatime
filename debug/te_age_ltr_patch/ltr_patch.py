#%%
import pandas as pd
#%%
subfamily = 'THE1C'
internal_id_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/'
old_internal_id_filepath = f'{internal_id_dir}/{subfamily}.internal_id.txt'
old_internal_id=pd.read_csv(old_internal_id_filepath, sep='\t')
#%%
new_internal_id_dir =  '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/annotation_mending_revise/workflow_test/'
new_internal_id_filepath = f'{new_internal_id_dir}/{subfamily}.internal_id.txt'
new_internal_id=pd.read_csv(new_internal_id_filepath, sep='\t', index_col = 0)
#%%
internal_ids=old_internal_id.merge(new_internal_id, on='rmsk_index')
# %%
age_table_dir ='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient'
age_table_filepath = f'{age_table_dir}/{subfamily}.txt'
age_table = pd.read_csv(age_table_filepath, sep='\t')
# %%
age_table
# %%
internal_ids_age = internal_ids.merge(age_table, right_on='internal_id', left_on='internal_id_x') 
# %%
# Extract INT status
internal_ids_age['int_type'] = internal_ids_age['internal_id_y'].str.extract(r'_(aINT|bINT|nINT|singleton)')

# Extract block_id from internal_id_y
internal_ids_age['block_id'] = internal_ids_age['internal_id_y'].str.extract(r'^(.*?_\d+)_')[0]

# Get counts of INT types per block_id
int_counts = internal_ids_age.pivot_table(index='block_id',
                                      columns='int_type',
                                      aggfunc='size',
                                      fill_value=0)

# Identify eligible blocks: both aINT and bINT must be present
eligible_blocks = int_counts[(int_counts.get('aINT', 0) > 0) & (int_counts.get('bINT', 0) > 0)].index

# Compute corrected_te_age: max age per eligible block
block_max_age = internal_ids_age[internal_ids_age['block_id'].isin(eligible_blocks) &
                             internal_ids_age['int_type'].isin(['aINT', 'bINT'])] \
    .groupby('block_id')['te_age'].max()

# Add corrected_te_age column
internal_ids_age['corrected_te_age'] = internal_ids_age.apply(
    lambda row: block_max_age[row['block_id']] if row['block_id'] in eligible_blocks and row['int_type'] in ['aINT', 'bINT']
    else row['te_age'], axis=1
)


# %%
internal_ids_age[internal_ids_age['int_type'].isin(['aINT','bINT'])]
# %%
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(8, 5))
plt.hist(internal_ids_age['te_age'], bins=30, alpha=0.5, label='Original te_age')
plt.hist(internal_ids_age['corrected_te_age'], bins=30, alpha=0.5, label='Corrected te_age')
plt.xlabel('TE Age')
plt.ylabel('Count')
plt.title('Histogram: TE Age Comparison')
plt.legend()
plt.tight_layout()
plt.show()

# %%
plt.figure(figsize=(6, 6))
plt.scatter(internal_ids_age['te_age'], internal_ids_age['corrected_te_age'], alpha=0.6)
plt.plot([internal_ids_age['te_age'].min(), internal_ids_age['te_age'].max()],
         [internal_ids_age['te_age'].min(), internal_ids_age['te_age'].max()],
         color='red', linestyle='--', label='No change line')

plt.xlabel('Original te_age')
plt.ylabel('Corrected te_age')
plt.title('Scatter: Original vs Corrected TE Age')
plt.legend()
plt.tight_layout()
plt.show()

# %%
total = len(internal_ids_age)
changed = (internal_ids_age['te_age'] != internal_ids_age['corrected_te_age']).sum()
unchanged = total - changed

print(f"Total entries: {total}")
print(f"Affected (age changed): {changed}")
print(f"Unchanged: {unchanged}")

# %%
summary = (
    internal_ids_age
    .loc[internal_ids_age['te_age'] != internal_ids_age['corrected_te_age']]
    .groupby(['te_age', 'corrected_te_age'])
    .size()
    .reset_index(name='count')
    .sort_values('count', ascending=False)
)

print("Age Correction Summary (Original ➜ Corrected):")
print(summary)

# %%
print("Original te_age summary:")
print(internal_ids_age['te_age'].describe())

print("\nCorrected te_age summary:")
print(internal_ids_age['corrected_te_age'].describe())
#%%
internal_ids_age.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/te_age_ltr_patch/the1c_age_bandaid_fix.txt', sep='\t')
# %%
from Bio.AlignIO import MafIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
maf_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species/241-mammalian-2020v2b.maf.chr11'
chrom = 'chr11'
target_species='Homo_sapiens'
from collections import Counter, defaultdict
maf_id = f'{target_species}.{chrom}'
mafindex_filepath = f'{maf_file}.mafindex'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_file, maf_id) 
#%%