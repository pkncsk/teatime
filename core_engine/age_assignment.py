#%%
import pandas as pd
from Bio.AlignIO import MafIO
import numpy as np
import os
import re
import math
import itertools
from concurrent.futures import ProcessPoolExecutor
from ma_mapper import extract_maf, utility
MafIO.MafIndex.get_spliced = extract_maf.get_spliced_mod
#%% INPUT PARAMETERS
subfamily_list = ['MER11A']
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
target_species = 'Homo_sapiens'
divergence_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species241_info.tsv'
maf_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species'
maf_file_prefix = '241-mammalian-2020v2b.maf'

internal_id_dir = None
e_table_dir = None
age_table_dir = None
#%%
def affine_count_simple(str1,str2,
    matchWeight = 1,
    mismatchWeight = 1,
    gapWeight = 4,
    extendedgapWeight = 1,
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
#%%
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
def string_to_float_list(val):
    # First, ensure that val is a list of strings
    if isinstance(val, list):
        return [float(i) if i != 'inf' else float('inf') for i in val]
    return []

#%%
def string_to_list(s):
    try:
        return [int(x) for x in re.findall(r'\d+', s)]
    except (ValueError, TypeError):
        print(f"Warning: Could not parse {s}")
        return None

def filter_e2(e_table,
              subfamily,
              internal_id_dir,
              e_cutoff = 1e-3, 
              test_width = 5000):
    te_age = 0
    segmental = pd.DataFrame()
    unclassified = pd.DataFrame()
    nosig_match = pd.DataFrame()
    internal_id = e_table.internal_id.unique()[0]
    if e_table.shape[0] > 1:
    #test second pass
        e_table=e_table[e_table.E_value.astype('float64')<=e_cutoff].copy()
        #Convert the strings in 'E_val_flanks' to lists
        e_table['E_val_flanks'] = e_table.E_val_flanks.apply(lambda x: string_to_float_list(eval(x)) if isinstance(x, str) else x)

        second_pass_tbl=e_table[e_table['E_val_flanks'].apply(lambda x: any(i <= e_cutoff for i in x))]
        if second_pass_tbl.shape[0] > 1:
            te_age=second_pass_tbl.divergence.max()
        else:
            subfamily_filename = subfamily.replace('/','_') 
            if internal_id_dir is None:
                internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])

            input_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
            internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
            internal_id=second_pass_tbl['internal_id'].unique()[0]
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
            start_flanked=[min(start_list)-test_width] + start_list + [max(end_list)]
            end_flanked = [min(start_list)] + end_list + [max(end_list)+test_width]
            spliced_maf_full =index_maf.get_spliced(start_flanked,end_flanked,strand)
            
            target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains('Homo_sapiens')]['seq'].to_list())[0]
        
            # Apply the optimized function to the DataFrame
            spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][test_width:-test_width]), target_seq[test_width:-test_width]), axis=1)
            spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
            spliced_maf_age=spliced_maf_full.merge(divergence_table[['meta_name','ungapped_length','Estimated Time (MYA)']], how='left',on='meta_name')
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
            first_pass[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][:5000]), target_seq[:5000]), axis=1)
            first_pass[['p_value_front', 'E_value_front']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
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
            
            first_pass[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][-5000:]), target_seq[-5000:]), axis=1)
            first_pass[['p_value_back', 'E_value_back']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
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
            
            e_table = second_pass.apply(calculate_metrics, axis=1).sort_values('divergence',ascending =True)
            if second_pass.shape[0] == 1:
                flanked_tbl = spliced_maf_age.copy()
                flanked_tbl[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = flanked_tbl.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][:5000]), target_seq[:5000]), axis=1)
                flanked_tbl[['p_value_front', 'E_value_front']] = flanked_tbl.apply(lambda row: pd.Series(BLAST_StoP(
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
                flanked_tbl[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = flanked_tbl.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][-5000:]), target_seq[-5000:]), axis=1)
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
                e_table['internal_id'] = internal_id
    
    if e_table.shape[0] == 1:
        te_age = np.nan
        e_table.columns = ['chr_code','divergence','%iden','%gap','BLAST_score','E_value','%iden_flanks','match_total_flank','front_back_matches','both_no_match', 'internal_id']
        e_table.loc[:,'match_total_flank'] = e_table.match_total_flank.astype(str).apply(string_to_list)
        if e_table['match_total_flank'].apply(lambda x: x[0] <= 1).values[0]:
            if e_table['match_total_flank'].apply(lambda x: x[0] == 1).values[0]:
                segmental = e_table
            else:
                unclassified = e_table
            #print(idx, 'possible segmental duplication')
        else:
            nosig_match = e_table
            te_age = 0
    return internal_id, te_age, segmental, unclassified, nosig_match

def filter_e_for_age(subfamily, 
                     e_table_dir,
                     age_table_dir,
                     e_cutoff = 1e-3):
    if e_table_dir is None:
        e_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    e_table = pd.read_csv(f'{e_table_dir}/{subfamily}.e_table.txt',sep = '\t', low_memory=False)
    if age_table_dir is None:
        age_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    output_filepath = f'{age_table_dir}/{subfamily}.teatime.txt'
    if (os.path.isfile(output_filepath) == False):
        print('start',subfamily)
        grouped =e_table.groupby('internal_id', sort=False)
        # Create a list to store the smaller DataFrames
        e_val_table_by_id = [group for _, group in grouped]
       
        global interal_id_dir
        with ProcessPoolExecutor(max_workers=90) as executor:
            results = executor.map(filter_e2, 
                                   e_val_table_by_id,
                                   itertools.repeat(subfamily), 
                                   itertools.repeat(internal_id_dir) ,
                                   itertools.repeat(e_cutoff))           

        id_list = []
        age_list = []
        nosig_match = []
        segmental = []
        unclassified = []
        for idx, result in enumerate(results):
            #print(idx)
            id_list.append(result[0])
            age_list.append(result[1])
            segmental.append(result[2])
            unclassified.append(result[3])
            nosig_match.append(result[4])
        dict_prep = {'internal_id': id_list, 'te_age':age_list,}
        output_table=pd.DataFrame(dict_prep)
        output_table.to_csv(output_filepath, sep='\t', index=False)

        nosig_match_df=pd.concat(nosig_match)
        if nosig_match_df.shape[0] > 0:
            nosig_filepath=f'{age_table_dir}/{subfamily}.insertion.txt'
            nosig_match_df.to_csv(nosig_filepath, sep='\t', index=False)
    
        segmental_df=pd.concat(segmental)
        if segmental_df.shape[0] > 0:
            segmental_filepath=f'{age_table_dir}/{subfamily}.segmental.txt'
            segmental_df.to_csv(segmental_filepath, sep='\t', index=False)
    
        unclassified_df=pd.concat(unclassified)
        if unclassified_df.shape[0] > 0:
            unclassified_filepath=f'{age_table_dir}/{subfamily}.unclassified.txt'
            unclassified_df.to_csv(unclassified_filepath, sep='\t', index=False)
        
        print('done',subfamily)
    else:
        print('already done', subfamily)
# %%
def main():
    global repeatmasker_table, divergence_table, maf_dir, maf_file_prefix
    repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
    divergence_table=utility.divergence_table_prep(divergence_table_filepath)
    #for subfamily in ['THE1C']:
    for subfamily in subfamily_list:
        filter_e_for_age(subfamily,e_table_dir,age_table_dir)
# %%
if __name__ == '__main__':
    main()
#%%