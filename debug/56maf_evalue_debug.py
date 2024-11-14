#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import custom_cmap
#%%
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_table)
# %%
coord_file = f'{config.coord_internal_id_folder}/{subfamily}.txt'
maf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf'
e_value_table = f'{config.e_value_folder}/{subfamily}.txt'
#%%
separated_maf = True
meta_id = metadata_age.meta_id.unique()
grouped = metadata_age.groupby('meta_id', sort=False)
chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
start_list = grouped.apply(lambda x: x.iloc[:,1].astype(int).tolist(), include_groups=False).tolist()
end_list = grouped.apply(lambda x: x.iloc[:,2].astype(int).tolist(), include_groups=False).tolist()
strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0], include_groups=False).tolist()
maf_call_list = []
for chrom in chrom_list:
    if separated_maf == True:
        maf_file = f'{maf_dir}.{chrom}'
    else:
        maf_file = maf_dir
    maf_call_list.append(maf_file)
#%%
from Bio.AlignIO import MafIO
import numpy as np
target_species = 'Homo_sapiens'
import time
count_arg='base_count'
start_time = time.time()
for idx, internal_id in enumerate([meta_id[0]]):
    print(idx, internal_id)
    start=start_list[idx]
    end=end_list[idx]
    strand=strand_list[idx]
    e_value_df = pd.read_csv(e_value_table, sep='\t')
    maf_id = f'{target_species}.{chrom}'
    
    index_maf = MafIO.MafIndex(f'{maf_file}.mafindex', maf_file, maf_id) 
    n_strand = -1 if strand == '-' else 1
    results =index_maf.get_spliced(start,end,n_strand)
    if e_value_df is None:
        collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results}
    else:
        e_value_internal_id=e_value_df[e_value_df.internal_id == internal_id]
        if e_value_internal_id.shape[0] <= 1:
            collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results if seqrec.id.split('.')[0] == 'Homo_sapiens'}
        else:
            try:
                collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results if seqrec.id.split('.')[0] in e_value_internal_id.species.values}
            except KeyError:
                print(f'source: {internal_id}')

    
    array_transposed=np.array(list(collector.values())).transpose()
    output_array=[]
    try: 
        ref_alleles=collector['Homo_sapiens']
    except KeyError:
                print(f'source: {internal_id}')
    for idx,ref_pos in enumerate(array_transposed):
        (unique, counts) = np.unique(ref_pos, return_counts=True)
        frequencies = dict(zip(unique, counts))
        ref_allele=ref_alleles[idx]
        #frequencies.pop('-', None) #->count deletion/insertion
        total = sum(frequencies.values())
        if count_arg == 'human_ref':
            alt_count = total - frequencies[ref_allele]
            alt_freq=alt_count/total
            if total == 1:
                alt_freq = np.nan
            output_array.append(alt_freq)
        elif count_arg in ['coverage', 'total_raw']:
            output_array.append(total)
        elif count_arg in ['common', 'common_raw']:
            common_allele = max(frequencies, key=frequencies.get)
            common_count = frequencies.get(common_allele)
            if count_arg == 'common':
                common_freq = common_count/total
                if total == 1:
                    alt_freq = np.nan
                output_array.append(common_freq)
            else:
                output_array.append(common_count)
        elif count_arg in ['a','t','c','g']:
            nucl_count = frequencies.get(count_arg, 0)
            output_array.append(nucl_count)
        elif count_arg == 'base_count':
            output_array.append(frequencies)
time.time() - start_time
#%%