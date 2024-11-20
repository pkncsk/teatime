#%%
import itertools
import numpy as np
import pandas as pd
from Bio.AlignIO import MafIO
from concurrent.futures import ProcessPoolExecutor
from ma_mapper import extract_maf, utility
MafIO.MafIndex.get_spliced = extract_maf.get_spliced_mod
from datetime import datetime
import math
#%% INPUT PARAMETERS
target_species = 'Homo_sapiens'
divergence_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species241_info.tsv'
subfamily = 'MER11A'
output_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/list_of_e_value_tables_for_figures'
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
maf_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species'
maf_file_prefix = '241-mammalian-2020v2b.maf'
internal_id_dir = None
#%% math function for blast score calculation
def affine_count_simple(str1,str2,
    matchWeight = 1,
    mismatchWeight = 11,
    gapWeight = 10,
    extendedgapWeight = 7,
    ):
    gapped = False
    str1 = str1.lower()
    str2 = str2.lower()
    length1 = len(str1)
    if str1 == str2:
        #skip process for perfectmatch
        return [matchWeight*length1, length1,gapped, 0]
    alignment_score = 0 
    matched = 0
    gap_count = 0
    in_gap = False
    for idx in range(length1):
        #check match
        if str1[idx] == str2[idx]:
            alignment_score = alignment_score + matchWeight
            matched = matched +1
            in_gap = False
        elif str2[idx] != '-':
            alignment_score = alignment_score - mismatchWeight
            in_gap = False
        elif str2[idx] == '-':
            gap_count = gap_count+1
            if in_gap == False:
                alignment_score = alignment_score - gapWeight
                in_gap = True
                gapped = True
            if in_gap == True:
                alignment_score = alignment_score -  extendedgapWeight
    #print('gapped: ', str(gapped))
    return [alignment_score, matched, gapped, gap_count]

def BLAST_Expm1(x):
    absx = abs(x);
    #print(x)
    #print(absx)
    #print(np.exp(x))
    if (absx > .33):
        return np.exp(x, dtype = 'float128') - 1.;
    elif (absx < 1.e-16):
        return absx
    else:
        return x * (1. + x *
             (1./2. + x * 
             (1./6. + x *
             (1./24. + x * 
             (1./120. + x *
             (1./720. + x * 
             (1./5040. + x *
             (1./40320. + x * 
             (1./362880. + x *
             (1./3628800. + x * 
             (1./39916800. + x *
             (1./479001600. + 
              x/6227020800.))))))))))));

def BLAST_Expm2(x):
    return np.exp(x, dtype = 'float128') - 1.

def BLAST_StoP(alignment_score, m,n ,lambda_ ,K, H, alpha, beta, gapped):
    if gapped == False:
        eff_l = (np.log (K* m * n)) / H
    else:
        
        #N aka db_num_seq will always be 1 since we search the entire genome
        N=1
        a = N
        mb = m * N + n
        c = n * m - max([m, n])/K
        kMaxIterations = 20
        ell_next = 0
        ell_min = 0
        converged = False
        #mb_power_2 = mb*mb
        #print('a: ',a,' mb: ',mb,' c: ',c,' test1:',mb_power_2, ' check1:', (mb * mb), ' check2:', (-4 * a * c))
        if (c < 0):
            eff_l = 0
        else:
            ell_max = 2 * c / (mb + math.sqrt(mb*mb - 4 * a * c)) 

            for i in range(kMaxIterations):
                ell = ell_next
                ss = (m - ell) * (n - N * ell)
                ell_bar = alpha/lambda_ * (np.log(K) + np.log(ss)) + beta
                if (ell_bar >= ell):
                    ell_min = ell
                    if(ell_bar - ell_min <= 1.0):
                        converged = True
                        break
                    if (ell_min == ell_max):
                        break
                else:
                    ell_max = ell
                
                if  ((ell_min <= ell_bar) & (ell_bar <= ell_max)):
                    ell_next = ell_bar
                else:
                    if i == 1:
                        ell_next = ell_max
                    else:
                        (ell_min + ell_max) / 2
            
            if converged == True:
                eff_l = ell_min
                ell = np.ceil(ell_min)
                if (ell <= ell_max):
                    ss = (m - ell) * (n - N * ell)
                    if  (alpha/lambda_ * (np.log(K) + np.log(ss)) + beta >= ell):
                        eff_l = ell
            else:
                eff_l = ell_min
    # In a large search space, the expected HSP length may be greater than 
    # the length of the query, resulting in a negative effective length, 
    # m´. In practice, if the effective length is less than 1/k, it is set to 1/k, 
    eff_m = m-eff_l
    if eff_m < 1/K:
        eff_m = 1/K    
    eff_n = n-eff_l
    search_sp = eff_m * eff_n
    E = search_sp * np.exp((-(lambda_) * alignment_score)+ np.log(K, dtype = 'float128'), dtype = 'float128')
    #print(E)
    p_value = -BLAST_Expm1(-E)
    #p_value = -BLAST_Expm2(-E)
    #p_value = 1 - math.exp(-E)
    return p_value, E
#%%
def affine_count_simple_optimized(seq_array, target_seq,
    matchWeight=1,
    mismatchWeight=1,
    gapWeight=4,
    extendedgapWeight=1):
    # Define IUPAC ambiguity code mappings in uppercase
    IUPAC_CODES = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'U': {'T'},  # Treat 'U' as 'T'
        'R': {'A', 'G','R'},
        'Y': {'C', 'T','Y'},
        'S': {'G', 'C','S'},
        'W': {'A', 'T','W'},
        'K': {'G', 'T','K'},
        'M': {'A', 'C','M'},
        'B': {'C', 'G', 'T','B'},
        'D': {'A', 'G', 'T','D'},
        'H': {'A', 'C', 'T','H'},
        'V': {'A', 'C', 'G','V'},
        'N': {'A', 'C', 'G', 'T','N'}
    }
    seq_array = np.char.upper(seq_array)
    target_seq = np.char.upper(target_seq)
    #print(''.join(seq_array))
    match = np.array([seq in IUPAC_CODES.get(target, set()) for seq, target in zip(seq_array, target_seq)])
    gap = (seq_array == '-')
    mismatch = ~match & ~gap
    
    
    alignment_score = (match * matchWeight).sum() - (mismatch * mismatchWeight).sum()
    
    gap_diff = np.diff(gap.astype(int))
    gap_starts = np.sum(gap_diff == 1)
    gap_extends = np.sum(gap) - gap_starts
    
    alignment_score -= (gap_starts * gapWeight + gap_extends * extendedgapWeight)
    
    matched = match.sum()
    gapped = gap.any()
    gap_count = gap.sum()
    
    return pd.Series([alignment_score, matched, gapped, gap_count])
def calculate_metrics(row):
    matched = row['matched']
    gap_count = row['gap_count']
    E_value = row['E_value']
    matched_front = row['matched_front']
    gap_count_front = row['gap_count_front']
    E_front = row['E_value_front']
    matched_back = row['matched_back']
    gap_count_back = row['gap_count_back']
    E_back = row['E_value_back']
    seq_length = row['seq_length']
    
    iden = round((matched / seq_length * 100), 2)
    pgap = round((gap_count / seq_length * 100), 2)
    E_score = '{:0.2e}'.format(E_value)
    iden_front = round((matched_front / 5000 * 100), 2)
    pgap_front = round((gap_count_front / 5000 * 100), 2)
    E_score_front = '{:0.2e}'.format(E_front)
    iden_back = round((matched_back / 5000 * 100), 2)
    pgap_back = round((gap_count_back / 5000 * 100), 2)
    E_score_back = '{:0.2e}'.format(E_back)
    
    return pd.Series({
        'species':  row['meta_name'],
        'chr_code': row['chr_code'],
        'divergence': row['Estimated Time (MYA)'],
        '%iden': iden,
        '%gap': pgap,
        'BLAST': row['alignment_score'],
        'E_value': E_score,
        '%iden_flanks': [iden_front, iden_back],
        '%gap_flanks': [pgap_front, pgap_back],
        'E_val_flanks': [E_score_front, E_score_back]
    })

#%% INITIATION

subfamily_filename = subfamily.replace('/','_') 
if internal_id_dir is None:
        internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
input_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
internal_id_list = internal_id_tbl.internal_id.unique()
repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
subfam_table=repeatmasker_table[repeatmasker_table['repName']==subfamily]
divergence_table=utility.divergence_table_prep(divergence_table_filepath)
#%% PARALLELIZE E-VALUE CALCULATION
def twopass_calc(internal_id,e_cutoff = 1e-3, testing_width = 5000, criteria = 'lenient'):
    print(internal_id)
    internal_id_tbl_subset = internal_id_tbl[internal_id_tbl.internal_id == internal_id]
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
    maf_filepath = f'{maf_dir}/{maf_file_prefix}.{chrom}'

    target_chrom = f'{target_species}.{chrom}'
    index_maf = MafIO.MafIndex(f'{maf_filepath}.mafindex', maf_filepath, target_chrom)
    start_flanked=[min(start_list)-testing_width] + start_list + [max(end_list)]
    end_flanked = [min(start_list)] + end_list + [max(end_list)+testing_width]
    
    spliced_maf_full =index_maf.get_spliced(start_flanked,end_flanked,strand)
    
    target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains(target_species)]['seq'].to_list())[0]
 
    # Apply the optimized function to the DataFrame
    spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][testing_width:-testing_width]), target_seq[testing_width:-testing_width]), axis=1)
    #spliced_maf_full['meta_name'] =spliced_maf_full['seqid'].str.split('.').str[0]
    spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
    spliced_maf_age=spliced_maf_full.merge(divergence_table[['meta_name','ungapped_length','Estimated Time (MYA)']], how='left',on='meta_name')
    #splice_maf_age_list.append(spliced_maf_age)
    spliced_maf_age['seq_length'] = spliced_maf_full['seq'].apply(len) - 2*testing_width
    lambda_ = 1.08;K = 0.28;H = 0.54;alpha = 2.0;beta = -2

    spliced_maf_age[['p_value', 'E_value']] = spliced_maf_age.apply(lambda row: pd.Series(BLAST_StoP(
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
    if first_pass.shape[0] < 1:
        e_table = pd.DataFrame()
    else:
        first_pass[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][:testing_width]), target_seq[:testing_width]), axis=1)
        first_pass[['p_value_front', 'E_value_front']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
            alignment_score=row['alignment_score_front'],
            m=testing_width,
            n=row['ungapped_length'],
            lambda_=lambda_,
            K=K,
            H=H,
            alpha=alpha,
            beta=beta,
            gapped=row['gapped_front']
        )), axis=1)
        
        first_pass[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][-testing_width:]), target_seq[-testing_width:]), axis=1)
        first_pass[['p_value_back', 'E_value_back']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
            alignment_score=row['alignment_score_back'],
            m=testing_width,
            n=row['ungapped_length'],
            lambda_=lambda_,
            K=K,
            H=H,
            alpha=alpha,
            beta=beta,
            gapped=row['gapped_back']
        )), axis=1)
        if criteria == 'stringent':
            second_pass=first_pass[(first_pass['E_value_front']<=e_cutoff)&(first_pass['E_value_back']<=e_cutoff)].copy()
        else:
            second_pass=first_pass[(first_pass['E_value_front']<=e_cutoff)|(first_pass['E_value_back']<=e_cutoff)].copy()
        if second_pass.shape[0] < 1:
            e_table = pd.DataFrame()
        else:
            e_table = second_pass
    return [e_table, spliced_maf_age]
#%%
criteria = 'lenient'
with ProcessPoolExecutor(max_workers=80) as executor:
    results = executor.map(twopass_calc, internal_id_list, itertools.repeat(criteria))
e_table_list = []
splice_maf_age_list = []
for result in results:
    e_table_list.append(result[0])
    splice_maf_age_list.append(result[1])
#%% SAVE FILTERED E-VALUE TABLE LIST (LENIENT CRITERIA)
import compress_pickle
output_filepath = f'{output_dir}/{subfamily}_0.4_LENIENT.lzma'
compress_pickle.dump(e_table_list, output_filepath, compression="lzma")
#%% SAVE E-VALUE TABLE LIST (LENIENT CRITERIA)
import compress_pickle
output_filepath = f'{output_dir}/{subfamily}_raw_0.4_LENIENT.lzma'
compress_pickle.dump(splice_maf_age_list, output_filepath, compression="lzma")
#%%
criteria = 'stringent'
with ProcessPoolExecutor(max_workers=80) as executor:
    results = executor.map(twopass_calc, internal_id_list, itertools.repeat(criteria))
e_table_list = []
splice_maf_age_list = []
for result in results:
    e_table_list.append(result[0])
    splice_maf_age_list.append(result[1])
#%% SAVE FILTERED E-VALUE TABLE LIST (STRINGENT CRITERIA)
import compress_pickle
output_filepath = f'{output_dir}/{subfamily}_0.4_LENIENT.lzma'
compress_pickle.dump(e_table_list, output_filepath, compression="lzma")
#%% SAVE E-VALUE TABLE LIST (STRINGENT CRITERIA)
import compress_pickle
output_filepath = f'{output_dir}/{subfamily}_raw_0.4_LENIENT.lzma'
compress_pickle.dump(splice_maf_age_list, output_filepath, compression="lzma")
#%%
def make_teatime_table(e_table_list):
    te_name_list = []
    te_age_list = []
    for idx, e_table in enumerate(e_table_list):
        #print(idx)
        if e_table.shape[0] > 0:
            te_age = e_table['Estimated Time (MYA)'].max()
        else:
            te_age = np.nan
        te_name = internal_id_list[idx]
        te_name_list.append(te_name)
        te_age_list.append(te_age)
        if e_table.shape[0]>0:
            for idx, row  in e_table.iterrows():
                if row['meta_name'] == 'Homo_sapiens':
                    seq = ''.join(row['seq'])  
                    len_list.append(len(seq))
    te_age_df=pd.DataFrame({
    'te_name': te_name_list,
    'te_age': te_age_list,
    })
    return te_age_df
### COMPARE OUTPUT 
#%% LOAD FILTERED E-VALUE TABLE LIST (LENIENT CRITERIA)
import compress_pickle
output_filepath = f'{output_dir}/{subfamily}_0.4_LENIENT.lzma'
e_table_list_lenient=compress_pickle.load(output_filepath, compression="lzma")
te_age_lenient = make_teatime_table(e_table_list_lenient)
#%% LOAD FILTERED E-VALUE TABLE LIST (STRINGENT CRITERIA)
import compress_pickle
output_filepath = f'{output_dir}/{subfamily}_0.4_STRINGENT.lzma'
e_table_list_stringent=compress_pickle.load(output_filepath, compression="lzma")
te_age_stringent = make_teatime_table(e_table_list_stringent)
#%% LOAD FILTERED E-VALUE TABLE LIST (ONE-PASS)
import compress_pickle
output_filepath = f'{output_dir}/{subfamily}_v0.3.lzma'
e_table_list_onepass=compress_pickle.load(output_filepath, compression="lzma")
te_age_onepass = make_teatime_table(e_table_list_onepass)
#%% MERGE TEATIME TO INTERNAL ID
te_age_idv_lenient = internal_id_tbl.merge(te_age_lenient, how='left', right_on='te_name', left_on='internal_id')
te_age_idv_stringent = internal_id_tbl.merge(te_age_stringent, how='left', right_on='te_name', left_on='internal_id')
te_age_idv_onepass = internal_id_tbl.merge(te_age_onepass, how='left', right_on='te_name', left_on='internal_id')
#COMPARISON
compare_te_age = te_age_idv_stringent.merge(te_age_idv_lenient, on = ['rmsk_index','internal_id','te_name'])
#%% COUNT TE IN EACH AGE CATEGORY
countFIRST = te_age_idv_onepass.groupby('te_age').size().reset_index(name='te_age_FIRST')
countOR = te_age_idv_lenient.groupby('te_age').size().reset_index(name='te_age_LENIENT')
countAND = te_age_idv_stringent.groupby('te_age').size().reset_index(name='te_age_STRINGENT')
count_dist = pd.merge(countFIRST, countOR, on='te_age', how='outer').fillna(0)
count_dist = pd.merge(count_dist, countAND, on='te_age', how='outer').fillna(0)
#%% PLOT TEATIME DISTRIBUTION
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig, ax = plt.subplots()
#grid = fig.add_gridspec(nrows = 300, ncols = 100, hspace=0)
# Define bar height
bar_height = 0.25
y = range(len(count_dist['te_age']))
# Plot each table's counts as horizontal bars with an offset on the y-axis
bars1=ax.barh([i + 2 * bar_height for i in y], count_dist['te_age_FIRST'], height=bar_height, label='1st pass', align='center', color='lightgrey')
bars2=ax.barh([i + bar_height for i in y], count_dist['te_age_LENIENT'], height=bar_height, label='2nd pass lenient', align='center', color = 'grey')
bars3=ax.barh(y, count_dist['te_age_STRINGENT'], height=bar_height, label='2nd pass stringent', align='center', color='black')

# Set y-ticks to match the categories
ax.set_yticks([i + bar_height for i in y])
ax.set_yticklabels(count_dist['te_age'])

# Labels and title
ax.set_xlabel('Count')
ax.set_ylabel('te_age')
#ax.set_xscale('log')
ax.tick_params(axis='both', which='major', labelsize=8)
ax.set_title(f'Distribution of TE age in {subfamily}')
#ax.set_title('Category Count Comparison across 3 Tables (Horizontal Bars)')
ax.legend()
# Add annotations to each bar
# For bars from Table 1
for bar in bars1:
    ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
            f'{int(bar.get_width())}', va='center', fontdict={'fontsize':6})

# For bars from Table 2
for bar in bars2:
    ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
            f'{int(bar.get_width())}', va='center', fontdict={'fontsize':6})

# For bars from Table 3
for bar in bars3:
    ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
            f'{int(bar.get_width())}', va='center', fontdict={'fontsize':6})
plt.show()
#%% DEBUG: FOR LOOP 
e_table_list = []
splice_maf_age_list = []
i = 0
e_cutoff=1e-3
test_width = 5000
criteria = 'lenient'
for idx, internal_id in enumerate(['THE1C_COMPLETE_3080']):
#for idx, internal_id in enumerate(internal_id_list):
    i = i+1
    print(f'process TE repeatmasker_index:\t{idx}\t{i}/{len(internal_id_list)}')
    internal_id_tbl_subset = internal_id_tbl[internal_id_tbl.internal_id == internal_id]
    subset_index=internal_id_tbl_subset.rmsk_index.to_list()
    rmsk_subset=repeatmasker_table[repeatmasker_table.index.isin(subset_index)]
    chrom=rmsk_subset.genoName.unique()[0]
    strand=rmsk_subset.strand.unique()[0]
    if strand =='-':
        internal_id_tbl_subset  = internal_id_tbl_subset.sort_index(ascending=False)
    start_list=rmsk_subset.genoStart.to_list()
    end_list=rmsk_subset.genoEnd.to_list()
    print(chrom, start_list,end_list,strand)
    if strand=='-':
        strand = -1
    else:
        strand = 1
    maf_filepath = f'{maf_dir}/{maf_file_prefix}.{chrom}'

    target_chrom = f'{target_species}.{chrom}'
    index_maf = MafIO.MafIndex(f'{maf_filepath}.mafindex', maf_filepath, target_chrom)
    start_flanked=[min(start_list)-test_width] + start_list + [max(end_list)]
    end_flanked = [min(start_list)] + end_list + [max(end_list)+test_width]
    
    spliced_maf_full =index_maf.get_spliced(start_flanked,end_flanked,strand)
    print(start_flanked, end_flanked, strand)
    target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains('Homo_sapiens')]['seq'].to_list())[0]
    print(''.join(target_seq[5000:-5000]))
    print(len(target_seq), len(target_seq[5000:-5000]))
    # Apply the optimized function to the DataFrame
    spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][test_width:-test_width]), target_seq[test_width:-test_width]), axis=1)
    #spliced_maf_full['meta_name'] =spliced_maf_full['seqid'].str.split('.').str[0]
    spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
    spliced_maf_age=spliced_maf_full.merge(divergence_table[['meta_name','ungapped_length','Estimated Time (MYA)']], how='left',on='meta_name')
    #splice_maf_age_list.append(spliced_maf_age)
    spliced_maf_age['seq_length'] = spliced_maf_full['seq'].apply(len) -2*test_width
    lambda_ = 1.08;K = 0.28;H = 0.54;alpha = 2.0;beta = -2

    spliced_maf_age[['p_value', 'E_value']] = spliced_maf_age.apply(lambda row: pd.Series(BLAST_StoP(
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
    if first_pass.shape[0] < 1:
        e_table = pd.DataFrame()
    else:
        first_pass[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][:test_width]), target_seq[:test_width]), axis=1)
        first_pass[['p_value_front', 'E_value_front']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
            alignment_score=row['alignment_score_front'],
            m=test_width,
            n=row['ungapped_length'],
            lambda_=lambda_,
            K=K,
            H=H,
            alpha=alpha,
            beta=beta,
            gapped=row['gapped_front']
        )), axis=1)
        
        first_pass[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][-test_width:]), target_seq[-test_width:]), axis=1)
        first_pass[['p_value_back', 'E_value_back']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
            alignment_score=row['alignment_score_back'],
            m=test_width,
            n=row['ungapped_length'],
            lambda_=lambda_,
            K=K,
            H=H,
            alpha=alpha,
            beta=beta,
            gapped=row['gapped_back']
        )), axis=1)
        if criteria == 'stringent':
            second_pass=first_pass[(first_pass['E_value_front']<=e_cutoff)&(first_pass['E_value_back']<=e_cutoff)].copy()
        else:
            second_pass=first_pass[(first_pass['E_value_front']<=e_cutoff)|(first_pass['E_value_back']<=e_cutoff)].copy()
        if second_pass.shape[0] < 1:
            e_table = pd.DataFrame()
        else:
            e_table = second_pass
#%%
te_name_list = []
te_age_list = []
len_list = []
alignment_row_list = []
alignment_max_list = []
for idx, e_table in enumerate(e_table_list):
    #print(idx)
    if e_table.shape[0] > 0:
        te_age = e_table['Estimated Time (MYA)'].max()
    else:
        te_age = np.nan
    #alignment_row = e_table[e_table['Estimated Time (MYA)'] == te_age].shape[0]
    #alignment_max_list.append(alignment_max)
    #alignment_row_list.append(alignment_row)
    te_name = internal_id_list[idx]
    te_name_list.append(te_name)
    #print(f'{idx}\t{te_name}\t{te_age}\t{alignment_row}/{alignment_max}')
    te_age_list.append(te_age)
    if e_table.shape[0]>0:
        for idx, row  in e_table.iterrows():
            if row['meta_name'] == 'Homo_sapiens':
                seq = ''.join(row['seq'])  
                len_list.append(len(seq))

#%%
te_age_df=pd.DataFrame({
    'te_name': te_name_list,
    'te_age': te_age_list,
    #'te_len': len_list,
    #'aln_row': alignment_row_list,
    #'aln_max':alignment_max_list,
})
#%% CHECK FOR ROWS WITH ONLY 1 ALINGMENT
#te_age_df[(te_age_df['aln_row']==1)&(te_age_df['aln_row']/te_age_df['aln_max']<0.5)]


#%%
# %%
#%%
