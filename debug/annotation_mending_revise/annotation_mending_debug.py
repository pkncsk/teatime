#%% 
import pandas as pd
import itertools
from concurrent.futures import ProcessPoolExecutor
from ma_mapper import utility
#%% INPUT PARAMETERS
output_dir = None #None will save the output file in the same dir as repeatmasker table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility
repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
subfamily = 'THE1C'
#%% table prep
subfamily_filename = subfamily.replace('/', '_')
subfam_table = repeatmasker_table[repeatmasker_table['repName'] == subfamily].copy()
subfam_table['length'] = subfam_table['genoEnd'] - subfam_table['genoStart']
subfam_table['repLength'] = subfam_table['repEnd'] - subfam_table['repStart'] + 1
# Calculate common length
if subfam_table[subfam_table.repLeft==0].shape[0] != 0:
    common_length = subfam_table[subfam_table.repLeft==0].repEnd.mode()[0]
else: 
    common_length = (subfam_table.repLeft + subfam_table.repEnd).mode()[0]
#%% grouping by completion
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
#sometimes the task will be simpler if the code imitates human behavior
#%%
#%%for idx in partial_frag_id:
idx = 186368
table_by_id = repeatmasker_table[(repeatmasker_table.id == idx) & (repeatmasker_table.repName == subfamily)].copy()
if table_by_id.strand.unique()[0] =='-':
    table_by_id  = table_by_id.sort_index(ascending=False)
table_by_id['repLength'] = table_by_id.repEnd - table_by_id.repStart + 1
table_by_id['next_repLength'] = table_by_id['repLength'].shift(-1)
table_by_id['next_repStart'] = table_by_id['repStart'].shift(-1)
table_by_id['next_repEnd'] = table_by_id['repEnd'].shift(-1)
table_by_id['next_repLeft'] = table_by_id['repLeft'].shift(-1)
table_by_id['tag'] = None
#%%
#1 look for rows that align nicely first
indext_list=table_by_id.index.to_list()
for idx, row in table_by_id.iterrows():
    if row.tag is None:
        print(f'check for potential complete sequnce that somehow avoid the filter: repLeft {row.repLeft}: repLength {row.repLength}>= commonlength {common_length}?')
        if row.repLeft == 0 and row.repLength >= common_length:
            print(f'COMPLETE sequnce repLeft == 0 and >= {common_length}')
            table_by_id.at[idx,'tag'] = 'COMPLETE'
        elif idx == indext_list[-1]:
            print(f'not complete but is the last row CLOSE')
            table_by_id.at[idx,'tag'] = 'CLOSE'
        elif row.next_repStart > row.repStart:
            #NOTE: weakest logic gate here:
            print(f'not COMPLETE not last row')
            print(f'fragment are potentially next to eachother with repStart current:{row.repStart} next:{row.next_repStart}')
            print(f'can the next fragment fit into the space left repLeft {row.repLeft} >= next row repLength {row.next_repLength}')
            if row.repLeft >= row.next_repLength:
                print('can fit OPEN')
                table_by_id.at[idx,'tag'] = 'OPEN'
            else:
                print(f'cannot fit because repLeft {row.repLeft} < next row repLength {row.next_repLength}') 
                #NOTE: this somehow doesn't make senses, 
                if row.repLeft > row.next_repLeft:
                    print(f'does not make sense but repLeft {row.repLeft} > next row repLeft {row.next_repLeft} OPEN')
                    table_by_id.at[idx,'tag'] = 'OPEN'
                else:
                    print(f'default')
                    table_by_id.at[idx,'tag'] = 'CLOSE'
        else:
            print(f'this is also not making any sense kind of desperate attemp')
            print(f'fragments are not next to each other/possible overlap?')
            if row.next_repLength +row.repEnd >= common_length:
                print(f'the current repEnd {row.repEnd}+next fragment length {row.next_repLength} >= {common_length}  CLOSE')
                table_by_id.at[idx,'tag'] = 'CLOSE'
            else:
                print(f'the current repEnd {row.repEnd}+next fragment length {row.next_repLength} < {common_length}  POSTHOC')
                table_by_id.at[idx,'tag'] = 'POSTHOC'
#%%
table_by_id['next_tag'] = table_by_id['tag'].shift(-1)
for idx, row in table_by_id.iterrows():
    if row.tag == 'POSTHOC':
        if row.next_tag == 'OPEN':
            #print(6)
            table_by_id.at[idx,'tag'] = 'CLOSE'
        elif row.next_repStart == row.repStart or row.next_repEnd == row.repEnd:
            #print(7)
            table_by_id.at[idx,'tag'] = 'OPEN'
current_label = None
# assign running numbers 
for index, row in table_by_id.iterrows():
    if row['tag'] == 'OPEN':
        if current_label is None:
            current_label = f'{subfamily}_JOIN_{internal_id}'
            internal_id += 1
        finalised_index_list.append(index)
        internal_id_list.append(current_label)
        #print(current_label)
    elif row['tag'] == 'CLOSE':
        if current_label is None:
            current_label = f'{subfamily}_SINGLE_{internal_id}'
            internal_id += 1
        finalised_index_list.append(index)
        internal_id_list.append(current_label)
        current_label = None
        #print(current_label)
    elif row['tag'] == 'COMPLETE':
        finalised_index_list.append(index)
        internal_id_list.append(f'{subfamily}_COMPLETE_{internal_id}')
        internal_id += 1
        current_label = None
#%%
dict_prep = {'rmsk_index': finalised_index_list,'internal_id': internal_id_list,}
internal_id_tbl=pd.DataFrame(dict_prep)
internal_id_tbl.to_csv(output_filepath, sep='\t', index=False)

#%%
