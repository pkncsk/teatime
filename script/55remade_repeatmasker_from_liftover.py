#%%
import pandas as pd
import os
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
#%%
repeatmasker_update=pd.read_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/repeatmasker_update/repeatmasker_update.txt', sep='\t', low_memory=False)
#%%