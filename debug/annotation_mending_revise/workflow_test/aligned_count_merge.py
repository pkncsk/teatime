#%%
import pandas as pd
import pickle
#%%
chrom = 'chr2'
input_filepath=f'/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/aligned_count_eff_update/{chrom}_merged.pkl'
#%%
with open(input_filepath, 'rb') as f:
    imported_table=pickle.load(f)
#%%
import pandas as pd
import pickle
import os

# List of chromosomes — adapt as needed
chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# Path template
base_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/aligned_count_eff_update'

# Dictionary to hold data per chromosome
all_data = {}

# Load each chromosome's pickle
for chrom in chromosomes:
    file_path = os.path.join(base_path, f'{chrom}_merged.pkl')
    with open(file_path, 'rb') as f:
        chrom_data = pickle.load(f)  # This is a dict {species: count}
        all_data[chrom] = chrom_data

# Convert to DataFrame
df = pd.DataFrame(all_data)

# Add a column for total aligned bases per species
df['total'] = df.sum(axis=1)

# Optional: sort by total
df = df.sort_values('total', ascending=False)

# Show the final DataFrame
print(df)

# %%
"""
species_info_filepath  = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_b_partial_fix_a.txt'
sp_df=pd.read_csv(species_info_filepath, sep='\t', index_col = 0)

emergency fix species table track name is supposed to have two hylobates pileastus just like in maf but in tree one of them was changed to lar instead
sp_df.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_b_partial_fix_b.txt', sep = '\t')
manual fix 
"""
#%%
species_info_filepath  = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_b_partial_fix_b.txt'
sp_df=pd.read_csv(species_info_filepath, sep='\t')
# %%
sp_df_aligned_count=sp_df.merge(df['total'], left_on = 'track_name', right_index = True)
# %%
sp_df_aligned_count.rename(columns={'total': 'ungapped_length'}, inplace=True)
#%%
sp_df_aligned_count.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_c.txt', sep = '\t', index = False)
# %%
