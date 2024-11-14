#%%
import numpy as np
import pandas as pd
from Bio.AlignIO import MafIO
from concurrent.futures import ProcessPoolExecutor
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_early_test as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_maf
MafIO.MafIndex.get_spliced = extract_maf.get_spliced_mod
#import config_mm39_dfam as config
#sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/varyzer/stable')
#import config
#import config_baseline as config
#%% 
#%%
subfamily = 'THE1C'
filtered_table=config.filtered_table
subfam_table=filtered_table[filtered_table['repName']==subfamily]
# %%
species_table_list = []
splice_maf_age_list = []
i = 0
for idx, row in subfam_table.iterrows():
    i = i+1
    print(f'process TE repeatmasker_index:\t{idx}\t{i}/{subfam_table.shape[0]}')
    start = row['genoStart']
    end = row['genoEnd']
    chrom = row['genoName']
    strand = row['strand']
    if strand=='-':
        strand = -1
    else:
        strand = 1
    target_species = config.target_species
    species_table = config.species_table
    if target_species == 'Homo_sapiens':
        mafpath = f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.{chrom}'

    target_chrom = f'{target_species}.{chrom}'
    index_maf = MafIO.MafIndex(f'{mafpath}.mafindex', mafpath, target_chrom)
    spliced_maf_full =index_maf.get_spliced([start],[end],strand)
    spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
    spliced_maf_age=spliced_maf_full.merge(species_table[['meta_name','ungapped_length','Estimated Time (MYA)']], how='left',on='meta_name')
    splice_maf_age_list.append(spliced_maf_age)
    species_table = spliced_maf_age[['meta_name', 'Estimated Time (MYA)']]
    species_table_list.append(species_table)
#%%
import compress_pickle
output_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age/{subfamily}_raw_v0.1.lzma'
compress_pickle.dump(species_table_list, output_filepath, compression="lzma")
#%%
import compress_pickle
output_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age/{subfamily}_raw_v0.1.lzma'
species_table_list=compress_pickle.load(output_filepath, compression="lzma")
# %%
te_name_list = []
te_age_list = []
for idx, species_table in enumerate(species_table_list):
    te_age = species_table['Estimated Time (MYA)'].max()
    repeatmasker_id=subfam_table.iloc[idx].name
    te_name = f'{subfamily}_{repeatmasker_id}'
    te_name_list.append(te_name)
    te_age_list.append(te_age)

te_age_df=pd.DataFrame({
    'te_name': te_name_list,
    'te_age': te_age_list
})
# %%
age_ref_table = pd.DataFrame(data=config.age_ref_table_template)
age_set=te_age_df['te_age'].sort_values().unique()

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(5,2))
grid = fig.add_gridspec(nrows = 10, ncols = 10, hspace=0)

bar_axes = fig.add_subplot(grid[0:10,0:10])
bar_axes.tick_params(axis='x', which='major', labelsize=8, labelrotation=90)
bar_axes.tick_params(axis='y', which='major', labelsize=8)
bar_container = bar_axes.bar(age_ref_table[age_ref_table['age'].isin(age_set)]['representative_age'], te_age_df.groupby('te_age', dropna=False).count()['te_name'], color = 'black')
bar_axes.set_title(f'Distribution of TE age in {subfamily}')
bar_axes.set_ylabel('counts')
bar_axes.bar_label(bar_container, fmt='{:,.0f}', fontsize=6)
#%%
#%%
seqid_list = []
seq_list = []
perc_gap_list = []
for idx, row in splice_maf_age_list[0].iterrows():
    seqid=row['seqid']
    seqid_list.append(seqid)
    seq=''.join(row['seq'])
    seq_list.append(seq)
    perc_gap_list.append(seq.count('-')/len(seq)*100)
#%%
def calculate_perc_gap(row):
    seq = ''.join(row['seq'])  
    perc_gap = seq.count('-') / len(seq) * 100  
    return perc_gap
bins = [-0.1,0,10,20,30,40,50,60,70,80,90,100]
perc_gap_count_list = []
for sub_df in splice_maf_age_list:
    sub_df['perc_gap'] = sub_df.apply(calculate_perc_gap, axis=1)
    sub_df['perc_binned'] = pd.cut(sub_df['perc_gap'], bins=bins, labels=list(range(1, len(bins))))
    perc_binned_count=sub_df.groupby(['perc_binned']).count()
    perc_gap_count=perc_binned_count['seqid'].values
    perc_gap_count_list.append(perc_gap_count)
#%%
import compress_pickle
output_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age/{subfamily}_v0.1.lzma'
compress_pickle.dump(splice_maf_age_list, output_filepath, compression="lzma")
#%%
import compress_pickle
output_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age/{subfamily}_v0.1.lzma'
splice_maf_age_list=compress_pickle.load(output_filepath, compression="lzma")
#%%
perc_gap_freq_array = np.array([arr / np.sum(arr) for arr in perc_gap_count_list])
#%%
seq_tbl=pd.DataFrame({
    'species': seqid_list,
    'alignment': seq_list
})
# %%
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(6,4))
grid = fig.add_gridspec(nrows = 5, ncols = 5, hspace=0)
labels = ['0','10','20','30','40','50','60','70','80','90','100']
box_axes = fig.add_subplot(grid[0:10,0:10])
box_axes.boxplot(perc_gap_freq_array, 1,'', labels=labels)
box_axes.set_ylabel('ratio')
box_axes.set_xlabel('gap percentage')
box_axes.set_title(f'THE1C gap percentage of alignment in MSA')
# %%
