#%%
import numpy as np
import pandas as pd
from Bio.AlignIO import MafIO
from concurrent.futures import ProcessPoolExecutor
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import dev.ma_mapper.script.config_main as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_maf
MafIO.MafIndex.get_spliced = extract_maf.get_spliced_mod
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
def e_val_engine_full(chrom, strand, start_list, end_list, e_cutoff=1e-3):
 
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
    spliced_maf_full =index_maf.get_spliced(start_list,end_list,strand)

    target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains('Homo_sapiens')]['seq'].to_list())[0]
 
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
    if first_pass.shape[0] < 1:
        return pd.DataFrame()
    else:
        return first_pass

#%%
def e_val_calc(internal_id):
    #print(f'{internal_id}')
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
    start_flanked=[min(start_list)-5000] + start_list + [max(end_list)]
    end_flanked = [min(start_list)] + end_list + [max(end_list)+5000]
    try:
        E_table=e_val_engine_full(chrom, strand, start_flanked, end_flanked)
    except UnboundLocalError:
        print(internal_id)
        # NOTE table join
    E_table['internal_id'] = internal_id
    return E_table
#%% DEBUG
subfamily = 'THE1C'
subfamily_filename = subfamily.replace('/','_') 
input_folder = config.internal_id_folder
input_filepath = f'{input_folder}/{subfamily_filename}.txt'
internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
internal_id_list = internal_id_tbl.internal_id.unique()
#%% DEBUG
results = []
#for internal_id in ['AluYc3_SINGLE_71','AluYc3_SINGLE_72','AluYc3_SINGLE_73']:
for internal_id in internal_id_list:
    print(internal_id)
    results.append(e_val_calc(internal_id))
#%% DEBUG
#output_table=pd.concat(results)
#%%
def e_val_calc_batch(subfamily):
    #load config
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
    logging.root.handlers = []
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    handlers = [
                        logging.FileHandler(operation_log_path, mode = 'a'),
                        logging.StreamHandler()
                        ]
                    )
    if (os.path.isfile(output_filepath) == False):
        # datetime object containing current date and time    
        start_time = time.time()
        now = datetime.now()
        logging.info(f'process: {subfamily}')
        output_file = f'{output_folder}/{subfamily_filename}.txt'
        internal_id_list = internal_id_tbl.internal_id.unique()
        id_count=len(internal_id_list)
        logging.info('entry counts: '+str(id_count))
        #print('checkpoint1')
        if id_count > 100:
            output_temp_folder = f'{output_folder}/{subfamily_filename}'
            if os.path.isdir(output_temp_folder) == False:
                os.mkdir(output_temp_folder)
            logging.info('over 100 count, spitting table')
            pointer = 0
            while pointer < id_count:
                time_start = time.time()
                startID = pointer
                if pointer + 100 < id_count:
                    endID = pointer + 100
                else:
                    endID = id_count
                pointer += 100
                id_sublist = internal_id_list[startID:endID]
                df_list = list()
                output_file_temp = f'{output_temp_folder}/{subfamily_filename}_{startID}_{endID}.txt'
                logging.debug(f'process: {output_file_temp}')
                if (os.path.isfile(output_file_temp) == False) :
                    with ProcessPoolExecutor(max_workers=80) as executor:
                        results = executor.map(e_val_calc, id_sublist)

                    for result in results:
                        df_list.append(result)
                    #for i, df in enumerate(df_list):
                    #    print(f"DataFrame {i}:: internal_id {id_sublist[i]}: Shape = {df.shape}")
                    #for te_id in internal_id_list:
                    #    result=e_val_calc(te_id)
                    #    df_list.append(result)
                    # NOTE: print temporary table out in case of multiple parts to avoid overhead issue
                    output_table=pd.concat(df_list)
                    output_table.to_csv(output_file_temp, sep='\t', index=False)
                    logging.info('done, saving the table temporarily at: '+output_file_temp)
                    logging.debug(" %s seconds ---\n" % (time.time() - time_start))
                else:
                    logging.info('already done with this fragment')
            logging.info('start merging files')
            table_list=glob.glob(f'{output_temp_folder}/{subfamily_filename}*.txt')
            for index, table_filename in enumerate(table_list):
                logging.info(f'merging file {index+1}/{len(table_list)}')
                if index == 0:
                    temp_df = pd.read_csv(table_filename, sep='\t')
                else:
                    currentable = pd.read_csv(table_filename, sep='\t')
                    temp_df=pd.concat([temp_df, currentable], axis=0, ignore_index=True)
            temp_df.to_csv(output_file, sep='\t', index=False)
            logging.info(f'done, saving the table at: {output_file}')
            logging.info('cleaning up temporary files')
            shutil.rmtree(output_temp_folder)

        else:
            df_list = list()
            with ProcessPoolExecutor(max_workers=80) as executor:
                results = executor.map(e_val_calc, internal_id_list)
            #
            for result in results:
                    df_list.append(result)
            #for i, df in enumerate(df_list):
            #    print(f"DataFrame {i}:: internal_id {id_sublist[i]}: Shape = {df.shape}")
            #for te_id in internal_id_list:
            #    result=e_val_calc(te_id)
            #    df_list.append(result)
            output_table=pd.concat(df_list)
            output_table.to_csv(output_file, sep='\t', index=False)
            logging.info('done, saving the table at: '+output_file)
        logging.debug("--- %s seconds ---\n" % (time.time() - start_time))
    else:
        logging.info(subfamily+' already done')
#%%
def main():
    #for subfam in ['MER11C']:
    for subfam in config.subfamily_list:
        e_val_calc_batch(subfam)
#%%
if __name__ == '__main__':
    main()
#%%