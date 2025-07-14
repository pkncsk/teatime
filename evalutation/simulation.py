#%%
import os
import pandas as pd
import numpy as np
#%%
work_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-run/simTHE1C'
np.random.seed(575)
seed_list = np.random.randint(0, 2**32, size=100, dtype=np.uint32) 
#%%
age_tables = []
for idx, seed in enumerate(seed_list):
    age_table_filepath = f'{work_dir}/{idx}/simTHE1C.teatime.txt'
    age_table=pd.read_csv(age_table_filepath, sep='\t')
    age_tables.append(age_table)
#%%
pd.concat(age_tables)
# %%
