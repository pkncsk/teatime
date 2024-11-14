#%%
import pandas as pd
import os
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import dev.ma_mapper.script.config_main as config
#%%
#%%
def merge_tables(subfamily):
    #subfamily = 'THE1C'
    subfamily_filename = subfamily.replace('/','_')
    print(f'process: {subfamily}')
    #subfamily = 'MER5C'
    #coord_internal_id_folder = config.coord_internal_id_folder
    #te_coord_tbl = f'{coord_internal_id_folder}/{subfamily_filename}.txt'
    #te_coord_df = pd.read_csv(te_coord_tbl, sep='\t')
    internal_id_folder = config.internal_id_folder
    internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
    te_age_folder = config.te_age_folder
    te_age_tbl = f'{te_age_folder}/{subfamily_filename}.txt'
    te_age_df = pd.read_csv(te_age_tbl, sep='\t')
    te_div_folder = config.kimura_distance_folder
    te_div_tbl = f'{te_div_folder}/{subfamily_filename}.txt'
    te_div_df = pd.read_csv(te_div_tbl, sep='\t')
    #te_tag_folder = config.te_tag_folder
    #te_tag_tbl = f'{te_tag_folder}/{subfamily_filename}.txt'
    #te_tag_df = pd.read_csv(te_tag_tbl, sep='\t')
    internal_id_div = internal_id_df.merge(te_div_df, on = 'internal_id', how='left')
    internal_id_age=internal_id_div.merge(te_age_df, on = 'internal_id', how='left')
    #internal_id_tag=internal_id_age.merge(te_tag_df, on = 'internal_id', how = 'left')
    internal_id_age['#entry'] = subfamily+'_'+ internal_id_age.index.astype(str)
    combined_te_age_div_folder = config.combined_age_div_folder
    output_filepath = f'{combined_te_age_div_folder}/{subfamily_filename}.txt'
    internal_id_age.to_csv(output_filepath, sep='\t', index=False)
    print(f'done: {subfamily}')
    return internal_id_age
# %%
def main():
    #for subfamily in ['THE1C']:
    collector = []
    for subfamily in config.subfamily_list:
        collector.append(merge_tables(subfamily))
    merged_table=pd.concat(collector)
    combined_te_age_div_folder = config.combined_age_div_folder
    output_filepath = f'{combined_te_age_div_folder}/all_subfamilies.txt'
    merged_table.to_csv(output_filepath, sep='\t', index=False)
# %%
if __name__ == '__main__':
    main()
#%%