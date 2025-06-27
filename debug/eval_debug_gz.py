#%%
import sys
import pandas as pd
sys.path.append('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime')
from core_engine import e_value_calculation
#%% INPUT PARAMETERS
target_species = 'Homo_sapiens'
divergence_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species241_info.tsv'
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
maf_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species'
maf_file_prefix = '241-mammalian-2020v2b.maf'
subfamily_list = ['MER11A']
internal_id_dir = None
e_table_dir = None
#%%
from ma_mapper import utility
repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
divergence_table=utility.divergence_table_prep(divergence_table_filepath)
#%%
subfamily_filename= subfamily_list[0]
if internal_id_dir is None:
        internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
input_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
# %%
test_id=internal_id_tbl.iloc[0]['internal_id']
#%% test e_val_calc 
internal_id_tbl_subset = internal_id_tbl[internal_id_tbl.internal_id == test_id]
subset_index=internal_id_tbl_subset.rmsk_index.to_list()
rmsk_subset=repeatmasker_table[repeatmasker_table.index.isin(subset_index)]
chrom=rmsk_subset.genoName.unique()[0]
strand=rmsk_subset.strand.unique()[0]
if strand =='-':
    internal_id_tbl_subset  = internal_id_tbl_subset.sort_index(ascending=False)
start_list=rmsk_subset.genoStart.to_list()
end_list=rmsk_subset.genoEnd.to_list()
if strand=='-':
    strand = -1
else:
    strand = 1
start_flanked=[min(start_list)-5000] + start_list + [max(end_list)]
end_flanked = [min(start_list)] + end_list + [max(end_list)+5000]
maf_filepath = f'{maf_dir}/{maf_file_prefix}.{chrom}.gz'

# %%
E_table = e_value_calculation.e_val_engine_full(chrom, strand, start_flanked, end_flanked,maf_filepath, e_cutoff=1e-3, target_species=target_species)
# %% test e_val_engine_full part1
from ma_mapper import gzmaf
target_chrom = f'{target_species}.{chrom}'
mafindex_filepath = f'{maf_dir}/{maf_file_prefix}.{chrom}.mafindex'
#print(mafindex_filepath)
index_maf = gzmaf.gzMafIndex(mafindex_filepath, maf_filepath, target_chrom)
#%%
print(start_list, end_list, strand)
spliced_maf_full =index_maf.get_spliced(start_list,end_list,strand)
#%%
start_flanked=[min(start_list)-100] + start_list + [max(end_list)]
end_flanked = [min(start_list)] + end_list + [max(end_list)+100]
#%%
spliced_maf_full =index_maf.get_spliced(start_flanked,end_flanked,strand)
# %% test e_val_engine_full part2
import numpy as np
target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains(target_species)]['seq'].to_list())[0]
#%%
# Apply the optimized function to the DataFrame
spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: e_value_calculation.affine_count_simple_optimized(np.array(row['seq'][5000:-5000]), target_seq[5000:-5000]), axis=1)
spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
spliced_maf_age=spliced_maf_full.merge(divergence_table[['meta_name','ungapped_length','Estimated Time (MYA)']], how='left',on='meta_name')
spliced_maf_age['seq_length'] = spliced_maf_full['seq'].apply(len) - 10000
lambda_ = 1.08;K = 0.28;H = 0.54;alpha = 2.0;beta = -2

spliced_maf_age[['p_value', 'E_value']] = spliced_maf_age.apply(lambda row: pd.Series(e_value_calculation.BLAST_StoP(
    alignment_score=row['alignment_score'],
    m=row['seq_length'],
    n=row['ungapped_length'],
    lambda_=lambda_,
    K=K,
    H=H,
    alpha=alpha,
    beta=beta,
    gapped=row['gapped']
)), axis=1)
first_pass = spliced_maf_age[spliced_maf_age['E_value']<=e_cutoff].copy()
#%%
e_cutoff = 1e-3
first_pass[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = first_pass.apply(lambda row: e_value_calculation.affine_count_simple_optimized(np.array(row['seq'][:5000]), target_seq[:5000]), axis=1)
first_pass[['p_value_front', 'E_value_front']] = first_pass.apply(lambda row: pd.Series(e_value_calculation.BLAST_StoP(
    alignment_score=row['alignment_score_front'],
    m=5000,
    n=row['ungapped_length'],
    lambda_=lambda_,
    K=K,
    H=H,
    alpha=alpha,
    beta=beta,
    gapped=row['gapped_front']
)), axis=1)

first_pass[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = first_pass.apply(lambda row: e_value_calculation.affine_count_simple_optimized(np.array(row['seq'][-5000:]), target_seq[-5000:]), axis=1)
first_pass[['p_value_back', 'E_value_back']] = first_pass.apply(lambda row: pd.Series(e_value_calculation.BLAST_StoP(
    alignment_score=row['alignment_score_back'],
    m=5000,
    n=row['ungapped_length'],
    lambda_=lambda_,
    K=K,
    H=H,
    alpha=alpha,
    beta=beta,
    gapped=row['gapped_back']
)), axis=1)
second_pass=first_pass[(first_pass['E_value_front']<=e_cutoff)|(first_pass['E_value_back']<=e_cutoff)].copy()
e_table = second_pass.apply(e_value_calculation.calculate_metrics, axis=1).sort_values('divergence',ascending =True)
if second_pass.shape[0] >1:
    return e_table
else:
    flanked_tbl = spliced_maf_age.copy()
    flanked_tbl[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = flanked_tbl.apply(lambda row: e_value_calculation.affine_count_simple_optimized(np.array(row['seq'][:5000]), target_seq[:5000]), axis=1)
    flanked_tbl[['p_value_front', 'E_value_front']] = flanked_tbl.apply(lambda row: pd.Series(e_value_calculation.BLAST_StoP(
        alignment_score=row['alignment_score_front'],
        m=5000,
        n=row['ungapped_length'],
        lambda_=lambda_,
        K=K,
        H=H,
        alpha=alpha,
        beta=beta,
        gapped=row['gapped_front']
    )), axis=1)
    flanked_tbl[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = flanked_tbl.apply(lambda row: e_value_calculation.affine_count_simple_optimized(np.array(row['seq'][-5000:]), target_seq[-5000:]), axis=1)
    flanked_tbl[['p_value_back', 'E_value_back']] = flanked_tbl.apply(lambda row: pd.Series(BLAST_StoP(
        alignment_score=row['alignment_score_back'],
        m=5000,
        n=row['ungapped_length'],
        lambda_=lambda_,
        K=K,
        H=H,
        alpha=alpha,
        beta=beta,
        gapped=row['gapped_back']
    )), axis=1)
    # Add columns to determine match types
    flanked_tbl['match_front'] = flanked_tbl['E_value_front'] <= e_cutoff
    flanked_tbl['match_back'] = flanked_tbl['E_value_back'] <= e_cutoff

    # Create additional columns for summary
    flanked_tbl['front_only'] = flanked_tbl['match_front'] & ~flanked_tbl['match_back']
    flanked_tbl['back_only'] = ~flanked_tbl['match_front'] & flanked_tbl['match_back']
    flanked_tbl['both'] = flanked_tbl['match_front'] & flanked_tbl['match_back']
    flanked_tbl['nonmatch'] = ~flanked_tbl['match_front'] & ~flanked_tbl['match_back']

    # Summarize the counts
    summary = {
        'total_count': len(flanked_tbl),
        'match_count': flanked_tbl['match_front'].sum() + flanked_tbl['match_back'].sum() -  flanked_tbl['both'].sum(),
        'front_only_count': flanked_tbl['front_only'].sum(),
        'back_only_count': flanked_tbl['back_only'].sum(),
        'both_count': flanked_tbl['both'].sum(),
        'nonmatch_count': flanked_tbl['nonmatch'].sum()
    }

    e_table = e_table[['chr_code','divergence','%iden','%gap','BLAST','E_value','%iden_flanks']]
    e_table['match_total'] = [[summary['match_count'],summary['total_count']]]
    e_table['front_back'] = [[summary['front_only_count'],summary['back_only_count']]]
    e_table['both_non'] = [[summary['both_count'],summary['nonmatch_count']]]
    e_table.columns = ['species','chr_code','divergence','%iden','%gap','BLAST','E_value','%iden_flanks','%gap_flanks','E_val_flanks']