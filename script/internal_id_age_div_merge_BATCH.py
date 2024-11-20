#%%
import pandas as pd
import sys
#%% INPUT PARAMETERS
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
internal_id_dir = None
kimura_distance_dir = None
age_table_dir = None
combined_table_dir = None
subfamily_list = ['MER11A']
#%%
def merge_tables(subfamily, internal_id_dir, age_table_dir, kimura_distance_dir,combined_table_dir):
    subfamily_filename = subfamily.replace('/','_')
    print(f'process: {subfamily}')
    if internal_id_dir is None:
        internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    internal_id_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
    internal_id_df = pd.read_csv(internal_id_filepath, sep='\t')

    if age_table_dir is None:
        age_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    teatime_filepath = f'{age_table_dir}/{subfamily_filename}.teatime.txt'
    te_age_df = pd.read_csv(teatime_filepath, sep='\t')

    if kimura_distance_dir is None:
        kimura_distance_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    te_div_tbl = f'{kimura_distance_dir}/{subfamily_filename}.div_table.txt'
    te_div_df = pd.read_csv(te_div_tbl, sep='\t')

    internal_id_div = internal_id_df.merge(te_div_df, on = 'internal_id', how='left')
    internal_id_age=internal_id_div.merge(te_age_df, on = 'internal_id', how='left')
    internal_id_age['#entry'] = subfamily+'_'+ internal_id_age.index.astype(str)
    
    output_filepath = f'{combined_table_dir}/{subfamily_filename}.itd.txt'
    internal_id_age.to_csv(output_filepath, sep='\t', index=False)
    print(f'done: {subfamily}')
    return internal_id_age
# %%
def main():
    global combined_table_dir
    if combined_table_dir is None:
        combined_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    collector = []
    for subfamily in subfamily_list:
        collector.append(merge_tables(subfamily, internal_id_dir, age_table_dir, kimura_distance_dir, combined_table_dir))
    merged_table=pd.concat(collector)
    output_filepath = f'{combined_table_dir}/all_subfamilies.itd.txt'
    merged_table.to_csv(output_filepath, sep='\t', index=False)
# %%
if __name__ == '__main__':
    main()
#%%