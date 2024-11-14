#%%
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
#import config_mm39_dfam as config
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
age_canon = config.age_canon
age_ref_table = pd.DataFrame(data=config.age_ref_table_reorg)
age_ref_table = pd.DataFrame(data=config.age_ref_table_template).iloc[0:-1]
#%%
repeatmasker_table  = config.rmskout_table
repeatmasker_table['rmsk_index'] = repeatmasker_table.index
combined_te_age_div_folder = config.combined_age_div_folder
combined_te_age_tbl = f'{combined_te_age_div_folder}/all_subfamilies.txt'
combined_te_age_df = pd.read_csv(combined_te_age_tbl, sep='\t')
#%%
repeatmasker_update=repeatmasker_table.merge(combined_te_age_df, how = 'left', on='rmsk_index')
#%%
age_df = repeatmasker_update[~repeatmasker_update['te_age'].isna()].copy()
#%% group correction
#age_df.loc[age_df['te_age'] == 82, 'te_age'] = 76
#age_df.loc[age_df['te_age'] == 90, 'te_age'] = 96
#%%
age_df['div_percent'] = age_df['te_div']
age_df['div_fraction'] = age_df['div_percent']/100
age_df['div_age'] = age_df['div_fraction']/(2*2e-9)
age_df['div_age_mya'] = age_df['div_age']/1e+6
age_df['length'] =  age_df.genoEnd - age_df.genoStart
bins = [-0.1]
genomesize = 3049315783
for i in range(0,62):
    bins.append(i+0.9)
age_df['binned'] = pd.cut(age_df['div_percent'], bins=bins, labels=list(range(0,62)))
bins = [-1]
for i in range(0,158):
    bins.append(i+0.9)
age_df['mya_binned'] = pd.cut(age_df['div_age_mya'], bins=bins, labels=list(range(0,158)))
#%%
age_df['percent_coverage'] = age_df.length/genomesize*100
comparable_cat = ['SINE/MIR','SINE/tRNA-Deu','SINE/tRNA-RTE','SINE/tRNA','SINE/Alu','SINE/5S-Deu-L2','Retroposon/SVA','LINE/Penelope','LINE/Dong-R4','LINE/Jockey','LINE/L2','LINE/CR1','LINE/RTE-X','LINE/RTE-BovB','LINE/L1','LINE/L1-Tx1', 'LTR/ERVK','LTR/ERV1','LTR','LTR/ERVL','LTR/ERVL-MaLR','LTR/Gypsy','RC/Helitron','DNA/TcMar','DNA/TcMar-Mariner','DNA/TcMar-Pogo','DNA/TcMar-Tc1','DNA/TcMar-Tc2','DNA/TcMar-Tigger','DNA/PiggyBac','DNA/MULE-MuDR','DNA/Merlin','DNA','DNA/Kolobok','DNA/hAT','DNA/hAT-Ac','DNA/hAT-Blackjack','DNA/hAT-Charlie','DNA/hAT-Tag1','DNA/hAT-Tip100','DNA/PIF-Harbinger','Unknown']
comparable_cat_color = ['#D7B4F8','#CE9BF7','#C481F5','#B966F4','#B358F3','#A637F1','#FF4D4D','#ACD8E5','#99B3D7','#8FA1CF','#625CB1','#483AA2','#38299A','#38299A','#00008B','#00008B','#90ED90','#73CD70','#65BD61','#57AE51','#57AE51','#489E42','#FF00FF','#FF936C','#FF936C','#FF936C','#FF936C','#FF936C','#FF936C','#FF865E','#FF7850','#FF7149','#FF6A42','#FF5A34','#FF512D','#FF512D','#FF512D','#FF512D','#FF512D','#FF512D','#FF4825','#999999']
comparable_cat.reverse()
comparable_cat_color.reverse()
# %%
#age_df.to_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/repeatmasker_update/age_hg38.txt', sep='\t')
# %%

# Create a custom mapping for specific subclasses
custom_groups = {
    'SINE/Alu': 'SINE/Alu',
    'SINE/MIR': 'SINE/MIR',
    'Retroposon/SVA': 'Retroposon/SVA',
    'LINE/L1': 'LINE/L1',
    'LINE/L2': 'LINE/L2',
    'LINE/CR1': 'LINE/CR1',
    'LTR/ERV1': 'LTR/ERV1',
    'LTR/ERVL': 'LTR/ERVL',
    'LTR/ERVL-MaLR': 'LTR/ERVL',
    'RC/Helitron':'RC/Helitron',
    'DNA/TcMar':'DNA/TcMar',
    'DNA/TcMar-Mariner':'DNA/TcMar',
    'DNA/TcMar-Pogo':'DNA/TcMar',
    'DNA/TcMar-Tc1':'DNA/TcMar',
    'DNA/TcMar-Tc2':'DNA/TcMar',
    'DNA/TcMar-Tigger':'DNA/TcMar',
    'DNA/hAT': 'DNA/hAT',
    'DNA/hAT-Ac': 'DNA/hAT',
    'DNA/hAT-Blackjack': 'DNA/hAT',
    'DNA/hAT-Charlie': 'DNA/hAT',
    'DNA/hAT-Tag1': 'DNA/hAT',
    'DNA/hAT-Tip100': 'DNA/hAT',
    'Unknown': 'Unknown'
    # Add more mappings as needed
}
# Group other subclasses into their respective major class followed by "/others"
for subclass in comparable_cat:
    if subclass not in custom_groups:
        major_class = subclass.split('/')[0]
        custom_groups[subclass] = f'{major_class}/others'

# Apply the custom group to the repClass_ column
age_df['custom_group'] = age_df['repClass'].map(custom_groups).fillna(age_df['repClass'])
# %%
age_count = age_df.groupby(['te_age', 'custom_group']).count().reset_index()[['te_age','custom_group','genoName']].rename(columns={'genoName':'count'})
#%%
repclass = ['Unknown', 'DNA/others','DNA/TcMar','DNA/hAT','RC/Helitron','LTR/others','LTR/ERVL','LTR/ERV1','LINE/others','LINE/CR1','LINE/L1','LINE/L2','Retroposon/SVA','SINE/others','SINE/Alu','SINE/MIR']
age_count_pivot=age_count.pivot(index = 'te_age',columns='custom_group', values='count')[repclass]
#%%
age_count_pivot.fillna(0, inplace=True)
age_count_by_class = [age_count_pivot[col].tolist() for col in repclass]


col_dict={
    'Unknown': '#999999',
    'DNA/others': '#FF6A42',
    'DNA/TcMar':'#FF512D',
    'DNA/hAT': '#FF936C',
    'RC/Helitron': '#FF00FF',
    'LTR/others': '#90ED90',
    'LTR/ERVL': '#57AE51',
    'LTR/ERV1': '#73CD70',
    'LINE/others': '#ACD8E5',
    'LINE/CR1':'#99B3D7',
    'LINE/L1': '#00008B',
    'LINE/L2': '#625CB1',
    'Retroposon/SVA': '#FF4D4D',
    'SINE/others': '#C481F5',
    'SINE/Alu': '#B358F3',
    'SINE/MIR': '#D7B4F8'
}

# %%
age_ref_table=age_ref_table[~age_ref_table['age'].isna()]
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,10))
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
main = fig.add_subplot(grid[0:50,20:70])
legend = fig.add_subplot(grid[0:30,72:75])
main.tick_params(axis='both', which='major', labelsize=8, labelrotation=90)
bottom = np.zeros(len(age_ref_table.representative_age))

# Create the stacked bar chart
for values, label in zip(age_count_by_class, repclass):
    main.bar(age_ref_table.representative_age, values, label=label, bottom=bottom, color=col_dict[label])
    bottom += values
import matplotlib.ticker as ticker
def format_func(value, tick_number):
    if value >= 1e6:
            return f'{value/1e6:.1f}M'
    elif value >= 1e3:
        return f'{value/1e3:.0f}k'
    else:
        return f'{value:.0f}'
formatter = ticker.FuncFormatter(format_func)
main.yaxis.set_major_formatter(formatter)
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib
cmap = ListedColormap(colors=col_dict.values())
bounds = range(cmap.N+1)
norm = BoundaryNorm(bounds,cmap.N)
data_legend=fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm), cax=legend, orientation='vertical')
data_legend.set_ticks([(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds) - 1)])
data_legend.set_ticklabels(repclass, fontsize = 'small')
data_legend.ax.set_title('TE Class', fontsize ='small')
main.set_title('Distribution of TE age estimated by TEA-TIME\ngrouped by class')
main.set_ylabel('counts')
# %%
#%%

#%%
age_count = age_df.groupby(['mya_binned', 'custom_group']).count().reset_index()[['mya_binned','custom_group','genoName']].rename(columns={'genoName':'count'})
#%%
repclass = ['Unknown', 'DNA/others','DNA/TcMar','DNA/hAT','RC/Helitron','LTR/others','LTR/ERVL','LTR/ERV1','LINE/others','LINE/CR1','LINE/L1','LINE/L2','Retroposon/SVA','SINE/others','SINE/Alu','SINE/MIR']
age_count_pivot=age_count.pivot(index = 'mya_binned',columns='custom_group', values='count')[repclass]
age_count_pivot.fillna(0, inplace=True)
age_count_by_class = [age_count_pivot[col].tolist() for col in repclass]


col_dict={
    'Unknown': '#999999',
    'DNA/others': '#FF6A42',
    'DNA/TcMar':'#FF512D',
    'DNA/hAT': '#FF936C',
    'RC/Helitron': '#FF00FF',
    'LTR/others': '#90ED90',
    'LTR/ERVL': '#57AE51',
    'LTR/ERV1': '#73CD70',
    'LINE/others': '#ACD8E5',
    'LINE/CR1':'#99B3D7',
    'LINE/L1': '#00008B',
    'LINE/L2': '#625CB1',
    'Retroposon/SVA': '#FF4D4D',
    'SINE/others': '#C481F5',
    'SINE/Alu': '#B358F3',
    'SINE/MIR': '#D7B4F8'
}

# %%
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,10))
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
main = fig.add_subplot(grid[0:50,20:70])
legend = fig.add_subplot(grid[0:30,72:75])
main.tick_params(axis='both', which='major', labelsize=8, labelrotation=90)
bottom = np.zeros(len(list(range(0,158))))
main.set_xlim(0,110)
# Create the stacked bar chart
for values, label in zip(age_count_by_class, repclass):
    main.bar(list(range(0,158)), values, label=label, bottom=bottom, color=col_dict[label])
    bottom += values
import matplotlib.ticker as ticker
def format_func(value, tick_number):
    if value >= 1e6:
            return f'{value/1e6:.1f}M'
    elif value >= 1e3:
        return f'{value/1e3:.0f}k'
    else:
        return f'{value:.0f}'
formatter = ticker.FuncFormatter(format_func)
main.yaxis.set_major_formatter(formatter)
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib
cmap = ListedColormap(colors=col_dict.values())
bounds = range(cmap.N+1)
norm = BoundaryNorm(bounds,cmap.N)
data_legend=fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm), cax=legend, orientation='vertical')
data_legend.set_ticks([(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds) - 1)])
data_legend.set_ticklabels(repclass, fontsize = 'small')
data_legend.ax.set_title('TE Class', fontsize ='small')
main.set_title('Distribution of TE age, estimated using\nKimura distance, grouped by class')
main.set_ylabel('counts')
main.set_xlabel('TE age (MYA)')
#%%
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,30))
grid = fig.add_gridspec(nrows = 300, ncols = 100, hspace=0)
legend = fig.add_subplot(grid[0:30,72:75])
bar_axes=[]
for idx, age in enumerate(age_canon):
    bottom = np.zeros(len(list(range(0,158))))
    age_df_sub=age_df[age_df['te_age']==age]
    age_count = age_df_sub.groupby(['mya_binned', 'custom_group'], observed=False).count().reset_index()[['mya_binned','custom_group','genoName']].rename(columns={'genoName':'count'})
    repclass = ['Unknown', 'DNA/others','DNA/TcMar','DNA/hAT','RC/Helitron','LTR/others','LTR/ERVL','LTR/ERV1','LINE/others','LINE/CR1','LINE/L1','LINE/L2','Retroposon/SVA','SINE/others','SINE/Alu','SINE/MIR']
    age_count_pivot=age_count.pivot(index = 'mya_binned',columns='custom_group', values='count')
    age_count_pivot.fillna(0, inplace=True)
    available_columns = [col for col in repclass if col in age_count_pivot.columns]
    age_count_pivot = age_count_pivot[available_columns]
    age_count_by_class = [age_count_pivot[col].tolist() for col in available_columns]
    bar_axes.append(fig.add_subplot(grid[0+20*idx:20+20*idx,0:70]))
    for values, label in zip(age_count_by_class, available_columns):
        bar_axes[idx].bar(list(range(0,158)), values, label=label, bottom=bottom, color=col_dict[label])
        bottom += values
    age_anno=age_ref_table[age_ref_table['age']==age]['representative'].values[0]
    bar_axes[idx].text(0.95, 0.85, f'Age: {age} MYA\n{age_anno}', 
            transform=bar_axes[idx].transAxes, 
            fontsize=20, 
            verticalalignment='top', 
            horizontalalignment='right')
    bar_axes[idx].set_xlim(0,110)
cmap = ListedColormap(colors=col_dict.values())
bounds = range(cmap.N+1)
norm = BoundaryNorm(bounds,cmap.N)
data_legend=fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm), cax=legend, orientation='vertical')
data_legend.set_ticks([(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds) - 1)])
data_legend.set_ticklabels(repclass, fontsize = 'small')
data_legend.ax.set_title('TE Class', fontsize ='small')

# %%
age_df[age_df['te_age']==105]['binned']
# %%
xlim_min = int(min(age_df[age_df['te_age']==105]['binned']))
# %%
custom_groups_list=['SINE/Alu', 'SINE/MIR', 'SINE/others','LINE/CR1', 'LINE/L1', 'LINE/L2','LINE/others','LTR/ERV1', 'LTR/ERVL', 'LTR/others', 'DNA/TcMar',  'DNA/hAT','DNA/others', 'RC/Helitron','Retroposon/SVA',  'Unknown']

# %%
all_sample_counts = [
    len(age_df[age_df['custom_group'] == group].groupby('te_age')['div_age_mya'].apply(list).iloc[i])
    for group in custom_groups_list
    for i in range(len(age_df[age_df['custom_group'] == group].groupby('te_age')['div_age_mya']))
]
#%%
set_pick=range(12,16)
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import seaborn as sns
fig = plt.figure(figsize=(20,60))
grid = fig.add_gridspec(nrows = 400, ncols = 110, hspace=0)

flierprops = dict(marker='o', markersize=1, linestyle='none', markeredgecolor='black', alpha=0.6)
boxprops = dict(facecolor = 'white',alpha=0.5)
box_axes = []
dot_axes = []
base_dot_size = 3000  

max_log_sample_size = np.log10(max(all_sample_counts))

for idx, custom_group in enumerate(custom_groups_list[i] for i in set_pick):
    custom_df=age_df[age_df['custom_group']==custom_group]
    grouped_data = custom_df.groupby('te_age')['div_age_mya'].apply(list)
    grouped_values = grouped_data.tolist()
    unique_ages = grouped_data.index.tolist()
    print(unique_ages)
    box_axes.append(fig.add_subplot(grid[0+35*idx:20+35*idx,0:100]))
    dot_axes.append(fig.add_subplot(grid[20+35*idx:32+35*idx,0:100]))
    dot_axes[idx].set_xticks(age_canon[:-1])
    box_axes[idx].set_xticks(age_canon[:-1])
    for i, d in enumerate(grouped_data):
        y = d
        x = np.random.normal(unique_ages[i], 0.2, len(y))  # Jitter added to the x-values
        box_axes[idx].plot(x, y, 'o', alpha=0.15, markersize=0.1, color = col_dict[custom_group])

    sample_counts = [len(values) for values in grouped_values]  # Counts of samples in each category

    dot_axes[idx].bar(unique_ages, sample_counts, color=col_dict[custom_group])
    # Plot a dot for each unique age showing sample size
    
    box_axes[idx].boxplot(grouped_values, positions=unique_ages, widths=2, flierprops=flierprops, patch_artist=True, boxprops=boxprops)
    #ax.set_xticks(range(len(unique_ages)))

    if idx < len(set_pick)-1:
        dot_axes[idx].set_xticklabels([])
        
    else:
        
        dot_axes[idx].set_xticklabels(unique_ages, fontsize=16, rotation=45)
        dot_axes[idx].set_xlabel('TE age estimated by TEA-TIME (MYA)', fontsize=22)
    dot_axes[idx].set_yscale('log')
    dot_axes[idx].set_ylim([1,1000000])
    box_axes[idx].set_yticks([0,25,50,75,100,125])
    box_axes[idx].set_yticklabels(box_axes[idx].get_yticks(),fontsize=16)
    box_axes[idx].set_ylim([-10,130])
    box_axes[idx].set_xlim([-4,110])
    dot_axes[idx].set_xlim([-4,110])
    dot_axes[idx].yaxis.set_label_position("right")
    dot_axes[idx].yaxis.tick_right()
    dot_axes[idx].set_yticks([1,10,100,1000,10000,100000,1000000])
    dot_axes[idx].yaxis.set_tick_params(labelsize=16)
    #ax.set_ylabel()
    if idx == 0:
        box_axes[idx].set_title(f'Distribution of TE age estimated using Kimura distance, grouped by TEA-TIME', fontsize = 28)
    box_axes[idx].text(0.99, 0.89, f'{custom_group} n={custom_df.shape[0]}', 
            transform=box_axes[idx].transAxes, 
            fontsize=16, 
            verticalalignment='bottom', 
            horizontalalignment='right')



fig.text(0.08, 0.76, 'TE age estimated by Kimura distance (MYA)', va='center', rotation='vertical', fontsize=22)
fig.text(0.87, 0.75, 'Barplots: counts', va='center', rotation='vertical', fontsize=22)
plt.show()
# %%
custom_groups_list=['L1HS', 'L1PA2', 'L1PA3', 'L1PA4',  'L1PA5','L1PA6', 'L1PA7','L1PA8']

all_sample_counts = [
    len(age_df[age_df['repName'] == group].groupby('te_age')['div_age_mya'].apply(list).iloc[i])
    for group in custom_groups_list
    for i in range(len(age_df[age_df['repName'] == group].groupby('te_age')['div_age_mya']))
]
# %%
set_pick=range(4,8)
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import seaborn as sns
fig = plt.figure(figsize=(20,60))
grid = fig.add_gridspec(nrows = 400, ncols = 110, hspace=0)

flierprops = dict(marker='o', markersize=1, linestyle='none', markeredgecolor='black', alpha=0.6)
boxprops = dict(facecolor = 'white',alpha=0.5)
box_axes = []
dot_axes = []
base_dot_size = 3000  

max_log_sample_size = np.log10(max(all_sample_counts))
for idx, custom_group in enumerate(custom_groups_list[i] for i in set_pick):
    custom_df=age_df[age_df['repName']==custom_group]
    grouped_data = custom_df.groupby('te_age')['div_age_mya'].apply(list)
    grouped_values = grouped_data.tolist()
    unique_ages = grouped_data.index.tolist()

    box_axes.append(fig.add_subplot(grid[0+35*idx:20+35*idx,0:100]))
    dot_axes.append(fig.add_subplot(grid[20+35*idx:32+35*idx,0:100]))
    dot_axes[idx].set_xticks(age_canon[:-1])

    for i, d in enumerate(grouped_data):
        y = d
        x = np.random.normal(unique_ages[i], 0.2, len(y))  # Jitter added to the x-values
        box_axes[idx].plot(x, y, 'bo', alpha=0.2, markersize=0.1)

    sample_counts = [len(values) for values in grouped_values]  # Counts of samples in each category

    dot_axes[idx].bar(unique_ages, sample_counts, color='black')
    # Plot a dot for each unique age showing sample size
    box_axes[idx].boxplot(grouped_values, positions=unique_ages, widths=2, flierprops=flierprops, patch_artist=True, boxprops=boxprops)
    #ax.set_xticks(range(len(unique_ages)))

    if idx < len(set_pick)-1:
        dot_axes[idx].set_xticklabels([])
        
    else:
        
        dot_axes[idx].set_xticklabels(age_canon[:-1], fontsize=16, rotation=45)
        dot_axes[idx].set_xlabel('TE age estimated by TEA-TIME (MYA)', fontsize=22)
    dot_axes[idx].set_yscale('log')
    dot_axes[idx].set_ylim([1,10000])
    box_axes[idx].set_yticks([0,25,50,75])
    box_axes[idx].set_yticklabels(box_axes[idx].get_yticks(),fontsize=16)
    box_axes[idx].set_ylim([-10,100])
    box_axes[idx].set_xlim([-4,110])
    dot_axes[idx].set_xlim([-4,110])
    dot_axes[idx].yaxis.set_label_position("right")
    dot_axes[idx].yaxis.tick_right()
    dot_axes[idx].set_yticks([1,10,100,1000,10000])
    dot_axes[idx].yaxis.set_tick_params(labelsize=16)

    #ax.set_ylabel()
    if idx == 0:
        box_axes[idx].set_title(f'Distribution of TE age of primate-specific L1\n estimated using Kimura distance, grouped by TEA-TIME', fontsize=28)
    if idx == len(custom_groups_list)-1:
        box_axes[idx].set_xlabel('TE age estimated by TEA-TIME (MYA)',fontsize=12)
    box_axes[idx].text(0.99, 0.85, f'{custom_group} n={custom_df.shape[0]}', 
            transform=box_axes[idx].transAxes, 
            fontsize=20, 
            verticalalignment='bottom', 
            horizontalalignment='right')
fig.text(0.08, 0.75, 'Boxplots: TE age estimated by Kimura distance (MYA)', va='center', rotation='vertical', fontsize=22)
fig.text(0.87, 0.75, 'Barplots: counts', va='center', rotation='vertical', fontsize=22)
plt.tight_layout()
plt.show()
# %%
