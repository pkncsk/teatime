#%%
import pandas as pd
#%%
#load age table
age_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/combined_age_div_lenient/all_subfamilies.itd.txt'
age_table=pd.read_csv(age_table_filepath, sep='\t')
#%%
# filter for subfamily (repName)
filtered_group=age_table[age_table['#entry'].str.contains('LTR104')]
#%% filter for extremes
Q1 = filtered_group['te_age'].quantile(0.25)
Q3 = filtered_group['te_age'].quantile(0.75)
IQR = Q3 - Q1
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

outliers = filtered_group[(filtered_group['te_age'] < lower_bound) | (filtered_group['te_age'] > upper_bound)]
outliers
# from the table, we can see that #entry LTR104_Mam_0, LTR104_Mam_101, LTR104_Mam_807, LTR104_Mam_1928 and LTR104_Mam_1929  are the one with 6.7 MYA
# we can check merge pattern from internal ID
# internal_id_structure 
# <subfamily>_<mode of merging>_<id>
# mode of merging: SINGLE = TE fragment can't find neighboring TE fragment, COMPLETE = TE fragment is a complete TE, no merge, JOIN = TE fragment is merged with neighboring TE fragment with the same internal ID 
# for example, the internal_id LTR104_Mam_JOIN_1248 has 2 #entries LTR104_Mam_1301 and LTR104_Mam_1302
# LTR104_Mam_0, LTR104_Mam_101, LTR104_Mam_807, LTR104_Mam_1928 and LTR104_Mam_1929's internal_id are LTR104_Mam_SINGLE_0, LTR104_Mam_SINGLE_101, LTR104_Mam_SINGLE_807, and LTR104_Mam_JOIN_1527 respectively
# this means the first three entries are single TEs, while the last two are joined.
# check e-value next to see if it is the case where alignments barely pass the 1e-3 threshold 
# %%
e_value_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/e_value/LTR104_Mam.txt'
e_value_table = pd.read_csv(e_value_table_filepath, sep='\t')
# 
# %% filter for internal_id
e_value_table[e_value_table['internal_id']=='LTR104_Mam_SINGLE_0']
# %%
e_value_table[e_value_table['internal_id']=='LTR104_Mam_SINGLE_101']
# %%
e_value_table[e_value_table['internal_id']=='LTR104_Mam_SINGLE_807']
# %%
e_value_table[e_value_table['internal_id']=='LTR104_Mam_JOIN_1527']
# %%
# all of them seem  to pass the threshold, the next step is to check the alignment or the e-value calculation step-by-step 
# %% 
# extract