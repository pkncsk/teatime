#%%
import pandas as pd
import os
from ma_mapper import utility
import numpy as np
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
#%% INPUT PARAMETERS
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
mean_phyloP_filepath  = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/phyloP.txt'
combined_table_dir = None
subfamily_list = ['MER11A']
#%% INITIATION
age_canon = utility.age_canon()[:-1]
mean_phyloP = pd.read_csv(mean_phyloP_filepath, sep='\t', index_col = 0).reset_index()
if combined_table_dir is None:
    combined_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
combined_te_age_filepath = f'{combined_table_dir}/all_subfamilies.itd.txt'
combined_te_age_df = pd.read_csv(combined_te_age_filepath, sep='\t')
#%%
collector = []
collector_ = []
_collector=[]
__collector=[]
_collector__=[]
for subfamily in subfamily_list:
    print(subfamily)
    subfamily_filename = subfamily.replace('/','_')
    #subfamily = 'THE1C'
    _collector.append(subfamily_filename)
    combined_filepath = f'{combined_table_dir}/{subfamily_filename}.itd.txt'
    subfam_combined=pd.read_csv(combined_filepath, sep='\t')
    median_raw = subfam_combined['te_age'].median()
    collector.append(median_raw)
    median_adjusted=find_nearest(age_canon,median_raw)
    _collector__.append(median_adjusted)
    collector_.append(subfam_combined['te_age'].mode()[0])
    __collector.append(subfam_combined['te_age'].mean())
prep_dict = {'subfamily':_collector, 'median_te_age':collector, 'median_adj_te_age':_collector__,'mode_te_age':collector_,'mean_te_age':__collector}
subfam_age=pd.DataFrame(prep_dict)


# %%
subfam_age_phyloP=subfam_age.merge(mean_phyloP, on='subfamily')
#%%
subset = subfam_age_phyloP[(subfam_age_phyloP['median_adj_te_age'] > 0) & (subfam_age_phyloP['median_adj_te_age'] <= 96.0)]
#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10,4))
flierprops = dict(marker='.', markersize=3, linestyle='none', markeredgecolor='black', alpha=0.6)
boxprops = dict(facecolor = 'white',alpha=0.5)
grouped_data = [subfam_age_phyloP[subfam_age_phyloP['median_adj_te_age'] == age]['mean_phyloP'] for age in age_canon]

for i, d in enumerate(grouped_data):
        y = d
        x = np.random.normal(age_canon[i], 0.4, len(y))  # Jitter added to the x-values
        ax.plot(x, y, 'bo', alpha=0.05, markersize=1)
ax.boxplot(grouped_data, positions=age_canon, widths=2, flierprops=flierprops,patch_artist=True, boxprops=boxprops)

from scipy import stats
# First regression line (0 to 43.2)
res = stats.linregress(subset['median_adj_te_age'], subset['mean_phyloP'])
ax.plot(subset['median_adj_te_age'], res.intercept + res.slope * subset['median_adj_te_age'], color='red', label='Regression Line (0-96)', linewidth=1, alpha=0.5)

set_size=subfam_age_phyloP[(subfam_age_phyloP['median_adj_te_age'] > 0)].shape[0]
ax.text(0.99, 0.01, f'TE subfamilies n={set_size}', 
        transform=ax.transAxes, 
        fontsize=12, 
        verticalalignment='bottom', 
        horizontalalignment='right')
ax.axhline(y=0, color='grey', linewidth=1, alpha=0.5)
ax.axhline(y=3.342447251181431, color='blue',linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlim(1,110)
ax.set_xlabel('median TE age (MYA) estimated by TEA-TIME', fontsize=12)
ax.set_ylabel('mean phyloP', fontsize=12)
ax.set_title('Mean phyloP grouped by TEA-TIME', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.show()
# %%
