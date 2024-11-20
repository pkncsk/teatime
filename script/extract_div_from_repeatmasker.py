#%%
import pandas as pd
import itertools
from ma_mapper import utility
from concurrent.futures import ProcessPoolExecutor
#%% INPUT PARAMETERS
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
internal_id_dir = None
kimura_distance_dir = None
subfamily_list = ['MER11A']
#%%
def calculate_div(group):
            frag_length = group['genoEnd'] - group['genoStart']
            frag_div = frag_length * group['milliDiv']
            return frag_div.sum() / frag_length.sum()

def extract_millDiv(subfamily, internal_id_dir, kimura_distance_dir):
    subfamily_filename = subfamily.replace('/','_')
    print(f'process: {subfamily}')
    if internal_id_dir is None:
        internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    input_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
    internal_id_df = pd.read_csv(input_filepath, sep='\t', index_col= 0)
    intersect_df = internal_id_df.merge(repeatmasker_table, left_on='rmsk_index', right_index=True)
    div_tbl = intersect_df.groupby('internal_id').apply(calculate_div).reset_index(name='te_div')
    if kimura_distance_dir is None:
        kimura_distance_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    output_filepath = f'{kimura_distance_dir}/{subfamily_filename}.div_table.txt'
    div_tbl.to_csv(output_filepath, sep='\t', index = False)
    print(f'done: {subfamily}')
#%%
def main():
    global repeatmasker_table
    repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(extract_millDiv,subfamily_list, itertools.repeat(internal_id_dir), itertools.repeat(kimura_distance_dir))
#%%
if __name__ == '__main__':
    main()
#%%