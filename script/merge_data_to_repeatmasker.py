#%%
import pandas as pd
import os
import sys
from ma_mapper import utility
#%% INPUT PARAMETERS
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
combined_table_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/combined_age_div_lenient/'
output_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/updated_repeatmasker_table_for_diana_and_michael/'
#%% INITIATION
repeatmasker_table=utility.repeatmasker_prep(repeatmasker_filepath)
repeatmasker_table['rmsk_index'] = repeatmasker_table.index
if combined_table_dir is None:
    combined_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
combined_te_age_filepath = f'{combined_table_dir}/all_subfamilies.itd.txt'
combined_te_age_df = pd.read_csv(combined_te_age_filepath, sep='\t')
#%%
repeatmasker_update=repeatmasker_table.merge(combined_te_age_df, how = 'left', on='rmsk_index')
# %%
lift_over_repeatmasker=pd.DataFrame()
lift_over_repeatmasker[['chrom','start','end','name']] = repeatmasker_update[['genoName','genoStart','genoEnd','#entry']]
lift_over_repeatmasker['score'] = 10
lift_over_repeatmasker['strand'] = repeatmasker_update.strand
#%%
if output_dir is None:
    output_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
lift_over_repeatmasker.to_csv('{output_dir}/for_liftover.txt', sep='\t', index=False, header=False)
#%% MANUALLY SUBMITTING the table to UCSC LIFTOVER
#%% IMPORT LIFTOVER TABLE
if output_dir is None:
    output_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
liftover=pd.read_csv(f'{output_dir}/liftover_hg19.bed',sep = '\t',header=None,dtype={'0': str, '1': int, '2': int,'3':str,'4':int,'5':str})
liftover.columns = ['chrom','start','end','#entry','score','strand_liftover']
# %%
repeatmasker_df=repeatmasker_update.merge(liftover, how = 'left', on = '#entry')
repeatmasker_df=repeatmasker_df.drop(columns=['score'])
# %%
repeatmasker_df.to_csv(f'{output_dir}/repeatmasker_update.txt', sep='\t',index=False)
# %%
