#%%
import pandas as pd
import os
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
#%%
repeatmasker_table  = config.rmskout_table
repeatmasker_table['rmsk_index'] = repeatmasker_table.index
combined_te_age_div_folder = config.combined_age_div_folder
combined_te_age_tbl = f'{combined_te_age_div_folder}/all_subfamilies.txt'
combined_te_age_df = pd.read_csv(combined_te_age_tbl, sep='\t')
#%%
repeatmasker_update=repeatmasker_table.merge(combined_te_age_df, how = 'left', on='rmsk_index')
# %%
filtered_table = repeatmasker_update[~repeatmasker_update.tag.isna()]
lift_over_repeatmasker=pd.DataFrame()
lift_over_repeatmasker[['chrom','start','end','name']] = filtered_table[['genoName','genoStart','genoEnd','#entry']]
lift_over_repeatmasker['score'] = 10
lift_over_repeatmasker['strand'] = filtered_table.strand
#%%
lift_over_repeatmasker.to_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/repeatmasker_update/for_liftover.txt', sep='\t', index=False, header=False)
# %%
liftover=pd.read_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/repeatmasker_update/liftover_hg19.bed',sep = '\t',header=None,dtype={'0': str, '1': int, '2': int,'3':str,'4':int,'5':str})
liftover.columns = ['chrom','start','end','#entry','score','strand_liftover']
# %%
repeatmasker_df=repeatmasker_update.merge(liftover, how = 'left', on = '#entry')
repeatmasker_df=repeatmasker_df.drop(columns=['score'])

# %%
repeatmasker_df.to_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/repeatmasker_update/repeatmasker_update.txt', sep='\t',index=False)
# %%
