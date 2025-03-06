#%%
#load package
from numpy import repeat
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
# all of them have the 6.7MYA alingment pass the threshold, and not barely (near 1e-3), the next step is to check the alignment or the e-value calculation step-by-step 
import sys
import pandas as pd
sys.path.append('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime')
from core_engine import e_value_calculation
#%%
# find genomic addresses of these TEs
# first extract repeatmasker table indices
outliers[outliers['te_age'] == 6.7]
# %%
rmsk_index=outliers[outliers['te_age'] == 6.7]['rmsk_index'].to_list()
#%% load repeatmasker table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
repeatmasker_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col=0)
#%% filter repeatmasker table
filtered_repeatmasker=repeatmasker_table[repeatmasker_table.index.isin(rmsk_index)]
#%% extract alignment from maf file using genomic addresses from the repeatmasker table (genoName=chrom, genoStart=start, genoEnd=end)
from ma_mapper import mafio_patch
from ma_mapper import utility
import numpy as np
maf_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species'
divergence_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species241_info.tsv'
divergence_table=utility.divergence_table_prep(divergence_table_filepath)
maf_file_prefix = '241-mammalian-2020v2b.maf'
target_species = 'Homo_sapiens'
chrom='chr1'
start_list=[224008920]
end_list=[224009176]
strand=-1
e_cutoff=1e-3
extend_length=5000
target_chrom = f'{target_species}.{chrom}'
maf_filepath = f'{maf_dir}/{maf_file_prefix}.{chrom}.gz'
mafindex_filedir = '.'.join(str.split(maf_filepath, sep ='.')[:-1])
mafindex_filepath = f'{mafindex_filedir}.mafindex'
#print(mafindex_filepath)
index_maf = mafio_patch.gzMafIndex(mafindex_filepath, maf_filepath, target_chrom)
start_flanked=[min(start_list)-extend_length] + start_list + [max(end_list)]
end_flanked = [min(start_list)] + end_list + [max(end_list)+extend_length]
spliced_maf_full =index_maf.get_spliced(start_flanked,end_flanked,strand)
target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains(target_species)]['seq'].to_list())[0]
# Apply the optimized function to the DataFrame
spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: e_value_calculation.affine_count_simple_optimized(np.array(row['seq'][extend_length:-extend_length]), target_seq[extend_length:-extend_length]), axis=1)
spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
spliced_maf_age=spliced_maf_full.merge(divergence_table[['meta_name','ungapped_length','Estimated Time (MYA)']], how='left',on='meta_name')
spliced_maf_age['seq_length'] = spliced_maf_full['seq'].apply(len) - (2*extend_length)
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
#%% full e-table
spliced_maf_age.to_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/LTR104_Mam_101.spliced_maf.csv', sep='\t')
#%% TE alignment with 1e-3 cutoff
first_pass.to_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/LTR104_Mam_807.first_pass.csv', sep='\t')
#%%
spliced_maf_age
#%% 
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
#%% 
first_pass.to_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/LTR104_Mam_0.first_pass_with_flank.csv', sep='\t')
# %%
#from the pre-filtered e-value table, they seem normal so i want to check the annotation 
# first load internal_id table
internal_id_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/LTR104_Mam.txt'
internal_id_table = pd.read_csv(internal_id_table_filepath, sep='\t')
# merge with repeatmasker table
repeatmasker_filtered_internal_id=repeatmasker_table.merge(internal_id_table, how='right', left_index=True, right_on='rmsk_index')[['genoName','genoStart','genoEnd','strand','internal_id','repStart','repEnd','repLeft']]
# %%
str(spliced_maf_age.iloc[67].seq)[5000:-5000]
# %%
spliced_maf_age['seq_substring'] = spliced_maf_age['seq'].apply(lambda seq: str(seq)[5000:-5000])
# %%
#0,101 = flanking region check
#807 = only partial TE matches
#1928,1929 = alignment aligned to either but not both
#%%
pd.read_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/LTR104_Mam_0.first_pass_with_flank.csv', sep='\t')
# %%
