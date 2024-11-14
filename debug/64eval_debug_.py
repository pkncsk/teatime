#%%
import numpy as np
import pandas as pd
from Bio.AlignIO import MafIO
from concurrent.futures import ProcessPoolExecutor
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
#import config_mm39_dfam as config
#sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/varyzer/stable')
#import config
#import config_baseline as config
import os
import time
from datetime import datetime
import shutil
import glob
import logging
import math
import itertools
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
subfamily='MER11A'
input_folder = config.internal_id_folder
output_folder = config.e_value_folder
if os.path.isdir(output_folder) == False:
    os.mkdir(output_folder)
subfamily_filename = subfamily.replace('/','_') 
input_filepath = f'{input_folder}/{subfamily_filename}.txt'
global internal_id_tbl
internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
output_filepath = f'{output_folder}/{subfamily_filename}.txt'
operation_log_path = f'{output_folder}/op_log.log'
#setup logger
#%%
#meta_id = 'MER11A_72'
#internal_id_tbl['meta_id'] = subfamily + '_' + internal_id_tbl.index.astype(str)
#internal_id = internal_id_tbl[internal_id_tbl.meta_id == meta_id]['internal_id'].values[0]
#%%
internal_id = 'MER11A_SINGLE_4'
internal_id_tbl_subset = internal_id_tbl[internal_id_tbl.internal_id == internal_id]
subset_index=internal_id_tbl_subset.rmsk_index.to_list()
rmsk_subset=config.filtered_table[config.filtered_table.index.isin(subset_index)]
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
#%%
e_cutoff=1e-3
#add flank fragments:
start_flanked=[min(start_list)-5000] + start_list + [max(end_list)]
end_flanked = [min(start_list)] + end_list + [max(end_list)+5000]
#%%
#%%
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_maf
MafIO.MafIndex.get_spliced = extract_maf.get_spliced_mod
#%%    
#full_length = [min(genostart_list)-5000,max(genoend_list)+5000]
if strand=='-':
    strand = -1
else:
    strand = 1
target_species = config.target_species
species_table = config.species_table
if target_species == 'Homo_sapiens':
    mafpath = f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.{chrom}'

target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(f'{mafpath}.mafindex', mafpath, target_chrom)
spliced_maf_full =index_maf.get_spliced(start_flanked,end_flanked,strand)
#%%
target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains('Homo_sapiens')]['seq'].to_list())[0]
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
        'T': {'T','U'},
        'U': {'T','U'},  # Treat 'U' as 'T'
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
    print('query')
    print(''.join(seq_array))
    target_seq = np.char.upper(target_seq)
    print('target')
    print(''.join(target_seq))
    match = np.array([seq in IUPAC_CODES.get(target, set()) for seq, target in zip(seq_array, target_seq)])

    #match = np.array([(target in IUPAC_CODES.get(seq, set()) and (seq in IUPAC_CODES[target] or seq == target))for seq, target in zip(seq_array, target_seq)])
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
    print(f'score {alignment_score}')
    return pd.Series([alignment_score, matched, gapped, gap_count])

# Assuming your DataFrame is named df
target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains('Homo_sapiens')]['seq'].to_list())[0]
#%%

# Apply the optimized function to the DataFrame
spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][5000:-5000]), target_seq[5000:-5000]), axis=1)
#spliced_maf_full['meta_name'] =spliced_maf_full['seqid'].str.split('.').str[0]
spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
spliced_maf_age=spliced_maf_full.merge(species_table[['meta_name','ungapped_length','Estimated Time (MYA)']], how='left',on='meta_name')
spliced_maf_age['seq_length'] = spliced_maf_full['seq'].apply(len) - 10000
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
#%%
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
#%%
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
#%%
second_pass=first_pass[(first_pass['E_value_front']<=e_cutoff)|(first_pass['E_value_back']<=e_cutoff)].copy()

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
#%%
e_table = second_pass.apply(calculate_metrics, axis=1).sort_values('divergence',ascending =True)
#%%
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
#%%
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
#%%
e_table['internal_id'] = internal_id
#%%
BLASTNA_TO_NCBI4NA = {
    'A': 1, 'C': 2, 'G': 4, 'T': 8, 'U': 8, 'R': 5, 'Y': 10,
    'S': 6, 'W': 9, 'K': 12, 'M': 3, 'B': 14, 'D': 13,
    'H': 11, 'V': 7, 'N': 15, '-': 0
}

# Function to compute degeneracy score
def calculate_degeneracy(base1, base2):
    if base1 == '-' or base2 == '-':
        return 0
    overlap = BLASTNA_TO_NCBI4NA[base1] & BLASTNA_TO_NCBI4NA[base2]
    if overlap:
        degen = bin(overlap).count('1')  # Count the number of matching bits
        return degen
    return 0

def affine_count_blast_like(seq_array, target_seq,
    matchWeight=1,
    mismatchWeight=1,
    gapWeight=4,
    extendedgapWeight=1):
    
    seq_array = np.char.upper(seq_array)
    target_seq = np.char.upper(target_seq)

    match = np.array([calculate_degeneracy(seq, target) for seq, target in zip(seq_array, target_seq)])
    # Full matches get the full score, partial matches adjusted
    match_scores = np.where(match == 1, matchWeight, (match - 1) * mismatchWeight + matchWeight)
    
    # Gap handling
    gap = (seq_array == '-')
    mismatch = (match == 0) & ~gap
    
    # Initial alignment score: full matches and adjusted partial matches
    alignment_score = match_scores.sum() - mismatch.sum() * mismatchWeight
    
    # Gap penalties
    gap_diff = np.diff(gap.astype(int))
    gap_starts = np.sum(gap_diff == 1)
    gap_extends = np.sum(gap) - gap_starts
    
    # Adjust score for gaps
    alignment_score -= (gap_starts * gapWeight + gap_extends * extendedgapWeight)
    
    matched = match.sum()  # Total matches, including partial
    gapped = gap.any()
    gap_count = gap.sum()
    
    print(f'Score: {alignment_score}')
    return pd.Series([alignment_score, matched, gapped, gap_count])
#%%
# Example usage
seq_array = np.array(list('ATANTGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAGGT'))
target_seq = np.array(list('ATCCTGTTATCTTCGTAAGCTGAGAATATACATCACCTCAGGACCACTATTGTACAAATTGATTGTAGAACATGTGT'))
#affine_count_blast_like(seq_array, target_seq)

# %%
seq_array = np.char.upper(seq_array)
target_seq = np.char.upper(target_seq)
# %%
match = np.array([calculate_degeneracy(seq, target) for seq, target in zip(seq_array, target_seq)])

# %%
reward = 1
penalty = -1
import numpy as np

# Define BLASTNA_TO_NCBI4NA as a dictionary with base to code mappings
BLASTNA_TO_NCBI4NA = {
    'A': 1, 'C': 2, 'G': 4, 'T': 8, 'U': 8, 'R': 5, 'Y': 10,
    'S': 6, 'W': 9, 'K': 12, 'M': 3, 'B': 14, 'D': 13,
    'H': 11, 'V': 7, 'N': 15, '-': 0
}

def calculate_degeneracy(base):
    """Calculate the degeneracy of a given base."""
    if base == '-':
        return 0
    return bin(BLASTNA_TO_NCBI4NA[base]).count('1')

def create_substitution_matrix(reward, penalty):
    bases = list(BLASTNA_TO_NCBI4NA.keys())
    size = len(bases)
    
    matrix = np.zeros((size, size), dtype=np.int32)

    for i, base1 in enumerate(bases):
        for j, base2 in enumerate(bases):
            if base1 == base2:
                matrix[i][j] = reward
            else:
                code1 = BLASTNA_TO_NCBI4NA[base1]
                code2 = BLASTNA_TO_NCBI4NA[base2]
                overlap = code1 & code2

                if overlap:
                    # Calculate degeneracy for the overlap
                    degen = bin(overlap).count('1')
                    # Calculate the score based on degeneracy
                    score = round(((degen - 1) * penalty + reward) / degen)
                    matrix[i][j] = score
                else:
                    matrix[i][j] = penalty
    
    # Handle gaps
    gap_value = np.iinfo(np.int32).min // 2
    matrix[-1, :] = gap_value
    matrix[:, -1] = gap_value

    return matrix

matrix = create_substitution_matrix(reward, penalty)
def calculate_alignment_score(query, subject, matrix):
    m = len(query)
    n = len(subject)
    
    # Initialize score matrix
    dp = np.zeros((m + 1, n + 1), dtype=np.int32)
    
    # Fill the first row and column
    for i in range(1, m + 1):
        dp[i][0] = dp[i-1][0] + matrix[BLASTNA_TO_NCBI4NA['-']][BLASTNA_TO_NCBI4NA[query[i-1]]]
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j-1] + matrix[BLASTNA_TO_NCBI4NA['-']][BLASTNA_TO_NCBI4NA[subject[j-1]]]
    
    # Fill the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + matrix[BLASTNA_TO_NCBI4NA[query[i-1]]][BLASTNA_TO_NCBI4NA[subject[j-1]]]
            delete = dp[i-1][j] + matrix[BLASTNA_TO_NCBI4NA[query[i-1]]][BLASTNA_TO_NCBI4NA['-']]
            insert = dp[i][j-1] + matrix[BLASTNA_TO_NCBI4NA['-']][BLASTNA_TO_NCBI4NA[subject[j-1]]]
            dp[i][j] = max(match, delete, insert)
    
    # The score is in the bottom-right corner
    return dp[m][n]

query = 'ATANTGCTATCTTCGTAAGCTGAGAATATACATCACCTCAGGACCACTATTGTACAAATTGATTGTAGAACATAGGT'
subject = 'ATCCTGTTATCTTCGTAAGCTGAGAATATACATCACCTCAGGACCACTATTGTACAAATTGATTGTAGAACATGTGT'

score = calculate_alignment_score(query, subject, matrix)
print(f"Alignment Score: {score}")
#%%
#import numpy as np

# Define BLASTNA_TO_NCBI4NA as a dictionary with base to code mappings
BLASTNA_TO_NCBI4NA = {
    'A': 1, 'C': 2, 'G': 4, 'T': 8, 'U': 8, 'R': 5, 'Y': 10,
    'S': 6, 'W': 9, 'K': 12, 'M': 3, 'B': 14, 'D': 13,
    'H': 11, 'V': 7, 'N': 15, '-': 0
}

def create_substitution_matrix(reward, penalty):
    bases = list(BLASTNA_TO_NCBI4NA.keys())
    size = len(bases)
    
    matrix = np.full((size, size), penalty, dtype=np.int32)  # Start with all mismatches

    for i, base1 in enumerate(bases):
        for j, base2 in enumerate(bases):
            if base1 == base2:
                matrix[i][j] = reward  # Reward for exact matches

    # Handle gaps
    gap_value = np.iinfo(np.int32).min // 2
    matrix[-1, :] = gap_value
    matrix[:, -1] = gap_value

    return matrix

def calculate_alignment_score(query, subject, matrix):
    m = len(query)
    n = len(subject)
    
    # Initialize score matrix
    dp = np.zeros((m + 1, n + 1), dtype=np.int32)
    
    # Fill the first row and column
    for i in range(1, m + 1):
        dp[i][0] = dp[i-1][0] + matrix[BLASTNA_TO_NCBI4NA[query[i-1]]][BLASTNA_TO_NCBI4NA['-']]
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j-1] + matrix[BLASTNA_TO_NCBI4NA['-']][BLASTNA_TO_NCBI4NA[subject[j-1]]]
    
    # Fill the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + matrix[BLASTNA_TO_NCBI4NA[query[i-1]]][BLASTNA_TO_NCBI4NA[subject[j-1]]]
            delete = dp[i-1][j] + matrix[BLASTNA_TO_NCBI4NA[query[i-1]]][BLASTNA_TO_NCBI4NA['-']]
            insert = dp[i][j-1] + matrix[BLASTNA_TO_NCBI4NA['-']][BLASTNA_TO_NCBI4NA[subject[j-1]]]
            dp[i][j] = max(match, delete, insert)
    
    # The score is in the bottom-right corner
    return dp[m][n]

# Example sequences
query1 = 'ATANTGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAGGT'
subject1 = 'ATCCTGTTATCTTCGTAAGCTGAGAATATACATCACCTCAGGACCACTATTGTACAAATTGATTGTAGAACATGTGT'
query2 = 'ATANTGCTATCTTCGTAAGCTGAGAATATACATCACCTCAGGACCACTATTGTACAAATTGATTGTAGAACATAGGT'
subject2 = 'ATCCTGTTATCTTCGTAAGCTGAGAATATACATCACCTCAGGACCACTATTGTACAAATTGATTGTAGAACATGTGT'

# Create substitution matrix
reward = 1
penalty = -1
matrix = create_substitution_matrix(reward, penalty)

# Calculate alignment scores
score1 = calculate_alignment_score(query1, subject1, matrix)
score2 = calculate_alignment_score(query2, subject2, matrix)

print(f"Alignment Score 1: {score1}")
print(f"Alignment Score 2: {score2}")

# %%
