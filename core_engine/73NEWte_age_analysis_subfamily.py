#%%
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import dev.ma_mapper.script.config_main as config
#%%
#%%
subfamily = 'MER11C'
age_ref_table = pd.DataFrame(data=config.age_ref_table_template)
stringent_file = f'{config.main_folder}/age_stringent/{subfamily}.txt'
lenient_file = f'{config.main_folder}/age_lenient/{subfamily}.txt'
age_stringent=pd.read_csv(stringent_file, sep='\t')
age_lenient=pd.read_csv(lenient_file, sep='\t')
#%%
count_stringent=age_stringent.groupby('te_age', dropna=False).size().reset_index(name='te_age_STRINGENT')
count_lenient=age_lenient.groupby('te_age', dropna=False).size().reset_index(name='te_age_LENIENT')
count_dist = pd.merge(count_lenient, count_stringent, on='te_age', how='outer')
count_dist['te_age_LENIENT'] = count_dist['te_age_LENIENT'].fillna(0)
count_dist['te_age_STRINGENT'] = count_dist['te_age_STRINGENT'].fillna(0)
# %%
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig, ax = plt.subplots()
#grid = fig.add_gridspec(nrows = 300, ncols = 100, hspace=0)
# Define bar height
bar_height = 0.25
y = range(len(count_dist['te_age']))

# Plot each table's counts as horizontal bars with an offset on the y-axis
bars2=ax.barh([i + bar_height for i in y], count_dist['te_age_LENIENT'], height=bar_height, label='2nd pass lenient', align='center', color = 'grey')
bars3=ax.barh(y, count_dist['te_age_STRINGENT'], height=bar_height, label='2nd pass stringent', align='center', color='black')

# Set y-ticks to match the categories
ax.set_yticks([i + bar_height for i in y])
ax.set_yticklabels(count_dist['te_age'])

# Labels and title
ax.set_xlabel('Count')
ax.set_ylabel('te_age')
#ax.set_xscale('log')
ax.tick_params(axis='both', which='major', labelsize=8)
ax.set_title(f'Distribution of TE age in {subfamily}')
#ax.set_title('Category Count Comparison across 3 Tables (Horizontal Bars)')
ax.legend()
# Add annotations to each bar
# For bars from Table 2
for bar in bars2:
    ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
            f'{int(bar.get_width())}', va='center', fontdict={'fontsize':6})

# For bars from Table 3
for bar in bars3:
    ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
            f'{int(bar.get_width())}', va='center', fontdict={'fontsize':6})

plt.show()
#%%