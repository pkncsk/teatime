#%%
import pandas as pd
import sys
import os
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
#%%
rmsk_table=config.rmskout_table
filtered_table =config.filtered_table
# %%
subfamily = 'THE1C'
#for subfamily in config.subfamily_list:
overlap_cov_threshold = 0.25
output_folder = config.internal_id_folder
subfamily_filename = subfamily.replace('/', '_')
output_filepath = f'{output_folder}/{subfamily_filename}.txt'

# Prepare TE data
subfam_table = filtered_table[filtered_table['repName'] == subfamily].copy()
subfam_table['length'] = subfam_table['genoEnd'] - subfam_table['genoStart']
subfam_table['repLength'] = subfam_table['repEnd'] - subfam_table['repStart'] + 1

# Calculate common length
if subfam_table[subfam_table.repLeft==0].shape[0] != 0:
    common_length = subfam_table[subfam_table.repLeft==0].repEnd.mode()[0]
else: 
    common_length = (subfam_table.repLeft + subfam_table.repEnd).mode()[0]

# Group fragments
singlet_id = []
singlet_index = []
all_complete_frag_id = []
all_complete_frag_index = []
partial_frag_id = []

for idx, row in subfam_table.groupby('id'):
    if all(row.repLength >= common_length):
        all_complete_frag_index.extend(row.index.tolist())
        all_complete_frag_id.extend(row.id.unique())
    elif row.shape[0] == 1:
        singlet_id.extend(row.id.tolist())
        singlet_index.extend(row.index.tolist())
    else:
        partial_frag_id.extend(row.id.unique())
#%%
# assign running numbers 
internal_id_list = []
finalised_index_list = []
internal_id = 0
for idx in singlet_index:
    finalised_index_list.append(idx)
    internal_id_list.append(f'{subfamily}_SINGLE_{internal_id}')
    internal_id +=1
for idx in all_complete_frag_index:
    finalised_index_list.append(idx)
    internal_id_list.append(f'{subfamily}_COMPLETE_{internal_id}')
    internal_id +=1
#%%
import numpy as np
finalised_index_list_test = []
internal_id_list_test = []
#sometimes the task will be simpler if the code imitates human behavior
for idx in [1493635]:
    table_by_id = filtered_table[(filtered_table.id == idx) & (filtered_table.repName == subfamily)].copy()
    if table_by_id.strand.unique()[0] =='-':
        table_by_id  = table_by_id.sort_index(ascending=False)
    table_by_id['repLength'] = table_by_id.repEnd - table_by_id.repStart + 1
    table_by_id['next_repLength'] = table_by_id['repLength'].shift(-1)
    table_by_id['next_repStart'] = table_by_id['repStart'].shift(-1)
    table_by_id['next_repEnd'] = table_by_id['repEnd'].shift(-1)
    table_by_id['next_repLeft'] = table_by_id['repLeft'].shift(-1)
    table_by_id['tag'] = None
    #1 look for rows that align nicely first
    indext_list=table_by_id.index.to_list()
    for idx, row in table_by_id.iterrows():
        if row.tag is None:
            if row.repLeft == 0 and row.repLength >= common_length:
                print(0)
                table_by_id.at[idx,'tag'] = 'COMPLETE'
            elif idx == indext_list[-1]:
                print(1)
                table_by_id.at[idx,'tag'] = 'CLOSE'
            elif row.next_repStart > row.repStart:
                if row.repLeft >= row.next_repLength:
                    print(2)
                    table_by_id.at[idx,'tag'] = 'OPEN'
                else:
                    if row.repLeft > row.next_repLeft:
                        print(3)
                        table_by_id.at[idx,'tag'] = 'OPEN'
                    else:
                        print(4)
                        table_by_id.at[idx,'tag'] = 'CLOSE'
            else:
                if row.next_repLength +row.repEnd >= common_length:
                    print(5)
                    table_by_id.at[idx,'tag'] = 'CLOSE'
                else:
                    print(99)
                    table_by_id.at[idx,'tag'] = 'POSTHOC'
    table_by_id['next_tag'] = table_by_id['tag'].shift(-1)
    for idx, row in table_by_id.iterrows():
        if row.tag == 'POSTHOC':
            if row.next_tag == 'OPEN':
                print(6)
                table_by_id.at[idx,'tag'] = 'CLOSE'
            elif row.next_repStart == row.repStart or row.next_repEnd == row.repEnd:
                print(7)
                table_by_id.at[idx,'tag'] = 'OPEN'
    current_label = None
    for index, row in table_by_id.iterrows():
        if row['tag'] == 'OPEN':
            if current_label is None:
                current_label = f'{subfamily}_JOIN_{internal_id}'
                internal_id += 1
            finalised_index_list_test.append(index)
            internal_id_list_test.append(current_label)
            print(current_label)
        elif row['tag'] == 'CLOSE':
            if current_label is None:
                current_label = f'{subfamily}_SINGLE_{internal_id}'
                internal_id += 1
            finalised_index_list_test.append(index)
            internal_id_list_test.append(current_label)
            print(current_label)
            current_label = None
        elif row['tag'] == 'COMPLETE':
            finalised_index_list_test.append(index)
            internal_id_list_test.append(f'{subfamily}_COMPLETE_{internal_id}')
            print(f'{subfamily}_COMPLETE_{internal_id}')
            internal_id += 1
            current_label = None
table_by_id[['repStart','repEnd','repLeft','repLength','tag']]
#%%
dict_prep = {'rmsk_index': finalised_index_list,'internal_id': internal_id_list,}
internal_id_tbl=pd.DataFrame(dict_prep)

#table_by_id[['repStart','repEnd','repLeft','repLength','tag']]
#%%
#1493635
#1945352
#2783791
#for idx in partial_frag_id:
for idx in [874833]:
    table_by_id = filtered_table[(filtered_table.id == idx) & (filtered_table.repName == subfamily)].copy()
    if table_by_id.strand.unique()[0] =='-':
        table_by_id  = table_by_id.sort_index(ascending=False)
    table_by_id['repLength'] = table_by_id.repEnd - table_by_id.repStart + 1
    table_by_id['next_repLength'] = table_by_id['repLength'].shift(-1)
    table_by_id['next_repStart'] = table_by_id['repStart'].shift(-1)
    table_by_id['next_repEnd'] = table_by_id['repEnd'].shift(-1)
    table_by_id['tag'] = None
    table_by_id.loc[(table_by_id['repLength'] > common_length * 0.95) or table_by_id, 'tag'] = 'COMPLETE'
    for idx, row in table_by_id.iterrows():
        if row.tag is None:
            if ((row.repEnd == row.next_repEnd) or (row.repStart == row.next_repStart)) and (row.next_repLength >1):
                table_by_id.loc[idx,'tag'] = 'CLOSE'
            #elif row.repLeft >= row.next_repLength:
            #    table_by_id.loc[idx,'tag'] = 'OPEN'
            else:
                if (row.repEnd - row.next_repStart)/common_length <= overlap_cov_threshold: 
                    table_by_id.loc[_idx,'tag'] = 'OPEN'
                else:
                    table_by_id.loc[_idx,'tag'] = 'CLOSE'
table_by_id
# %%
for idx in partial_frag_id:
#for idx in [3889163]:
    table_by_id = filtered_table[(filtered_table.id == idx) & (filtered_table.repName == subfamily)].copy()
    if table_by_id.shape[0]>3:
        print(idx)
# %%
    for idx in partial_frag_id:
    #for idx in [439967]:
        table_by_id = filtered_table[filtered_table.id == idx]
        if table_by_id.strand.unique()[0] =='-':
            table_by_id  = table_by_id.sort_index(ascending=False)
        table_by_id=table_by_id.reset_index()
        merge_switch = False
        before = []
        for row_num, (table_idx,table_row) in enumerate(table_by_id.iterrows()):
            
            prev_subfam = table_by_id.loc[row_num - 1, 'repName'] if row_num > 0 else None
            next_subfam = table_by_id.loc[row_num + 1, 'repName'] if row_num < len(table_by_id) - 1 else None
            repLength  = table_row.repEnd - table_row.repStart
        # Check if previous or next value equals the target value
            #print(prev_subfam,next_subfam)
            if ((prev_subfam is None and next_subfam) or (prev_subfam and next_subfam is None) or (prev_subfam and next_subfam)) and (table_row.repName == subfamily):
                #print('check1', next_subfam)
                
                if next_subfam != subfamily and next_subfam:
                    before.append([table_row['index'],idx,next_subfam,repLength,table_row.repLeft])
                #print('check2', prev_subfam)
                
                if before:
                    if prev_subfam != subfamily and prev_subfam:
                        
                        index_b,idx_b,next_subfam_b,repLength_b,refleft_b=before[0] 
                        if refleft_b > common_length * 0.05:
                            print(index_b,idx_b,'before',next_subfam_b,repLength_b,refleft_b) 
                            before = []
                            print(table_row['index'],idx,'after',prev_subfam,repLength,table_row.repLeft)