#%%
import numpy as np
import pandas as pd
from Bio.AlignIO import MafIO
from concurrent.futures import ProcessPoolExecutor
import os
import time
from datetime import datetime
import shutil
import glob
import logging
import math
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ma_mapper import gzmaf
#%% INPUT PARAMETERS
target_species = 'Homo_sapiens'
divergence_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species241_info.tsv'
subfamily_list = ['THE1C']
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
maf_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species/'
internal_id_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/annotation_mending_revise/workflow_test/'
e_table_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/annotation_mending_revise/workflow_test/'
extended_length = 500
#%%
#%%
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table

def divergence_table_prep(divergence_table_filepath):
    species = pd.read_table(divergence_table_filepath, sep='\t')
    species['meta_name']=species['track_name']

    return species

def get_spliced_mod(self, starts, ends, strand=1):
    # Dictionary for IUPAC ambiguity codes for 2-base combinations
    iupac_code = {
        frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S', frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }

    def convert_to_iupac(sequence):
        unique_bases = frozenset(sequence)
        if len(unique_bases) == 1:
            return sequence[0].upper()  
        return iupac_code.get(unique_bases, 'N')  # Default to 'N' for any unhandled cases
    
    def process_sequence_localized(sequence):
        sequence = sequence.upper()
        filtered_sequence = [base for base in sequence if base != '-']

        #base_counts = Counter(sequence)
        #most_common_bases = base_counts.most_common()
        #max_count = most_common_bases[0][1]
        #consensus_bases = [base for base, count in most_common_bases if count == max_count]
        new_base = convert_to_iupac(filtered_sequence)
        return new_base

    if strand not in (1, -1): 
        raise ValueError("Strand must be 1 or -1, got %s" % strand)
    fetched = list(self.search(starts, ends))
    #return fetched
    expected_letters = sum(end - start for start, end in zip(starts, ends))
    if len(fetched) == 0:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})
    if len(fetched) == 1:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})
    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    split_by_position
    total_rec_length = 0
    ref_first_strand = None
    for multiseq in fetched:
        for seqrec in multiseq:
            if seqrec.id == self._target_seqname:
                try:
                    if ref_first_strand is None:
                        ref_first_strand = seqrec.annotations["strand"]

                        if ref_first_strand not in (1, -1):
                            raise ValueError("Strand must be 1 or -1")
                    elif ref_first_strand != seqrec.annotations["strand"]:
                        raise ValueError(
                            "Encountered strand='%s' on target seqname, "
                            "expected '%s'"
                            % (seqrec.annotations["strand"], ref_first_strand)
                        )
                except KeyError:
                    raise ValueError(
                        "No strand information for target seqname (%s)"
                        % self._target_seqname
                    ) from None

                rec_length = len(seqrec)
                rec_start = seqrec.annotations["start"]
                ungapped_length = seqrec.annotations["size"]
                rec_end = rec_start + ungapped_length - 1
                total_rec_length += ungapped_length
                
                for seqrec in multiseq:
                    for pos in range(rec_start, rec_end + 1):
                        split_by_position[seqrec.id][pos] = ""

                break 
            else:
                raise ValueError(
                    "Did not find %s in alignment bundle" % (self._target_seqname,)
                )
        real_pos = rec_start
        edit_id = []
        edit_pos = []
        for gapped_pos in range(rec_length):
            previous_id = ''
            for seqrec in multiseq:
                
                if seqrec.id == self._target_seqname:
                    track_val = seqrec.seq[gapped_pos]
                
                
                split_by_position[seqrec.id][real_pos] += seqrec.seq[gapped_pos]
                if previous_id == seqrec.id:
                        edit_id.append(seqrec.id)
                        edit_pos.append(real_pos)
                previous_id = seqrec.id
            if track_val != "-" and real_pos < rec_end:
                real_pos += 1
        # Debugging: Print lengths of sequences in split_by_position
        for i in range(len(edit_id)):
            _sequence=split_by_position[edit_id[i]][edit_pos[i]]
            new_sequence=process_sequence_localized(_sequence)
            split_by_position[edit_id[i]][edit_pos[i]] = new_sequence
        
        if len(split_by_position[self._target_seqname]) != total_rec_length:
            raise ValueError(
                "Target seqname (%s) has %s records, expected %s"
                % (
                    self._target_seqname,
                    len(split_by_position[self._target_seqname]),
                    total_rec_length,
                )
            )

    realpos_to_len = {
        pos: len(gapped_fragment)
        for pos, gapped_fragment in split_by_position[self._target_seqname].items()
        if len(gapped_fragment) > 1
    }

    seqid_list = []
    seq_list = []
    
    for seqid in all_seqnames:
        seq_split = split_by_position[seqid]
        seq_splice = []
        filler_char = "N" if seqid == self._target_seqname else "-"
        append = seq_splice.append
        for exonstart, exonend in zip(starts, ends):
            for real_pos in range(exonstart, exonend):
                if real_pos in seq_split:
                    append(seq_split[real_pos])
                elif real_pos in realpos_to_len:
                    append(filler_char * realpos_to_len[real_pos])
                else:
                    append(filler_char)

        seqid_list.append(seqid)
        seq_list.append(Seq("".join(seq_splice))) 

    if len(seq_list[seqid_list.index(self._target_seqname)].replace("-", "")) != expected_letters:
        raise ValueError(
            "Returning %s letters for target seqname (%s), expected %s"
            % (
                len(seq_list[seqid_list.index(self._target_seqname)].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(seq_list[seqid_list.index(self._target_seqname)])
    for seqid, seq in zip(seqid_list, seq_list):
        if len(seq) != ref_subseq_len:
            raise ValueError(
                "Returning length %s for %s, expected %s"
                % (len(seq), seqid, ref_subseq_len)
            )

    # Create a DataFrame
    df = pd.DataFrame({
        'seqid': seqid_list,
        'seq': [seq.reverse_complement() if strand != ref_first_strand else seq for seq in seq_list]
    })
    return df

MafIO.MafIndex.get_spliced = get_spliced_mod
#%% math function for blast score calculation
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
        print('idx: ', idx)
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
    iden_front = round((matched_front / extended_length * 100), 2)
    pgap_front = round((gap_count_front / extended_length * 100), 2)
    E_score_front = '{:0.2e}'.format(E_front)
    iden_back = round((matched_back / extended_length * 100), 2)
    pgap_back = round((gap_count_back / extended_length * 100), 2)
    E_score_back = '{:0.2e}'.format(E_back)
    
    return pd.Series({
        'species':  row['meta_name'],
        'chr_code': row['chr_code'],
        'divergence': row['divergence'],
        '%iden': iden,
        '%gap': pgap,
        'blast': row['alignment_score'],
        'E_value': E_score,
        '%iden_flanks': [iden_front, iden_back],
        '%gap_flanks': [pgap_front, pgap_back],
        'E_val_flanks': [E_score_front, E_score_back]
    })

#%%
def e_val_engine_full(chrom, 
                      strand, 
                      start_list, 
                      end_list, 
                      maf_filepath,
                      e_cutoff=1e-3, 
                      target_species = 'Homo_sapiens',
                      intermidiate_output=None):
    print(f'{chrom}:{start_list}{end_list}({strand})')
    print(maf_filepath)
    target_chrom = f'{target_species}.{chrom}'
    mafindex_filedir = '.'.join(str.split(maf_filepath, sep ='.')[:-1])
    mafindex_filepath = f'{mafindex_filedir}.mafindex'
    #print(mafindex_filepath)
    if '.gz' in mafindex_filepath:
        index_maf = gzmaf.gzMafIndex(mafindex_filepath, maf_filepath, target_chrom)
    else:
        index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
    spliced_maf_full =index_maf.get_spliced(start_list,end_list,strand)

    target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains(target_species)]['seq'].to_list())[0]

    # Apply the optimized function to the DataFrame
    spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][extended_length:-extended_length]), target_seq[extended_length:-extended_length]), axis=1)
    spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
    spliced_maf_age=spliced_maf_full.merge(divergence_table[['meta_name','ungapped_length','divergence']], how='left',on='meta_name')
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
    if intermidiate_output=='nofilter':
        return spliced_maf_age
    else:
        first_pass = spliced_maf_age[spliced_maf_age['E_value']<=e_cutoff].copy()
        if first_pass.shape[0] < 1:
            return pd.DataFrame()
        else:
            first_pass[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][:extended_length]), target_seq[:extended_length]), axis=1)
            first_pass[['p_value_front', 'E_value_front']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
                alignment_score=row['alignment_score_front'],
                m=extended_length,
                n=row['ungapped_length'],
                lambda_=lambda_,
                K=K,
                H=H,
                alpha=alpha,
                beta=beta,
                gapped=row['gapped_front']
            )), axis=1)
            
            first_pass[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = first_pass.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][-extended_length:]), target_seq[-extended_length:]), axis=1)
            first_pass[['p_value_back', 'E_value_back']] = first_pass.apply(lambda row: pd.Series(BLAST_StoP(
                alignment_score=row['alignment_score_back'],
                m=extended_length,
                n=row['ungapped_length'],
                lambda_=lambda_,
                K=K,
                H=H,
                alpha=alpha,
                beta=beta,
                gapped=row['gapped_back']
            )), axis=1)
            if intermidiate_output=='firstpass':
                return first_pass   
            else:
                second_pass=first_pass[(first_pass['E_value_front']<=e_cutoff)|(first_pass['E_value_back']<=e_cutoff)].copy()
                e_table = second_pass.apply(calculate_metrics, axis=1).sort_values('divergence',ascending =True)
                if second_pass.shape[0] >1:
                    return e_table
                else:
                    flanked_tbl = spliced_maf_age.copy()
                    flanked_tbl[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = flanked_tbl.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][:extended_length]), target_seq[:extended_length]), axis=1)
                    flanked_tbl[['p_value_front', 'E_value_front']] = flanked_tbl.apply(lambda row: pd.Series(BLAST_StoP(
                        alignment_score=row['alignment_score_front'],
                        m=extended_length,
                        n=row['ungapped_length'],
                        lambda_=lambda_,
                        K=K,
                        H=H,
                        alpha=alpha,
                        beta=beta,
                        gapped=row['gapped_front']
                    )), axis=1)
                    flanked_tbl[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = flanked_tbl.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][-extended_length:]), target_seq[-extended_length:]), axis=1)
                    flanked_tbl[['p_value_back', 'E_value_back']] = flanked_tbl.apply(lambda row: pd.Series(BLAST_StoP(
                        alignment_score=row['alignment_score_back'],
                        m=extended_length,
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

                    e_table = e_table[['chr_code','divergence','%iden','%gap','blast','E_value','%iden_flanks']]
                    e_table['match_total'] = [[summary['match_count'],summary['total_count']]]
                    e_table['front_back'] = [[summary['front_only_count'],summary['back_only_count']]]
                    e_table['both_non'] = [[summary['both_count'],summary['nonmatch_count']]]
                    e_table.columns = ['species','chr_code','divergence','%iden','%gap','blast','E_value','%iden_flanks','%gap_flanks','E_val_flanks']
                    return e_table
#%%
def get_maf_filepath(maf_dir, chrom):
    files = os.listdir(maf_dir)
    maf_filename = f"241-mammalian-2020v2b.maf.{chrom}.gz" 
    maf_files = [f for f in files if f.endswith(maf_filename)]

    # Determine the appropriate file to use
    maf_filepath = None
    if any(f.endswith('.maf') for f in maf_files):
        maf_filepath = f"{maf_dir}/{maf_filename}" #.maf is more optimal performance wise 
    elif any(f.endswith(f'.gz') for f in maf_files):
        maf_filepath = f"{maf_dir}/{maf_filename}" 
    else:
        raise FileNotFoundError(f"No .maf or .maf.gz file found for chromosome {chrom} in {maf_dir}")

    return maf_filepath
#%%
def e_val_calc(internal_id, target_species = 'Homo_sapiens'):
    internal_id_tbl_subset = internal_id_tbl[internal_id_tbl.internal_id == internal_id]
    subset_index=internal_id_tbl_subset.rmsk_index.to_list()
    rmsk_subset=repeatmasker_table[repeatmasker_table.index.isin(subset_index)]
    chrom=rmsk_subset.genoName.unique()[0]
    print(chrom)
    strand=rmsk_subset.strand.unique()[0]
    if strand =='-':
        internal_id_tbl_subset  = internal_id_tbl_subset.sort_index(ascending=False)
    start_list=rmsk_subset.genoStart.to_list()
    end_list=rmsk_subset.genoEnd.to_list()
    if strand=='-':
        strand = -1
    else:
        strand = 1
    start_flanked=[min(start_list)-extended_length] + start_list + [max(end_list)]
    end_flanked = [min(start_list)] + end_list + [max(end_list)+extended_length]
    print(f'{chrom}:{start_list}{end_list}({strand})')
    try:
        maf_filepath = get_maf_filepath(maf_dir, chrom)
        print(maf_filepath)
        E_table=e_val_engine_full(chrom, strand, start_flanked, end_flanked,maf_filepath, e_cutoff=1e-3, target_species=target_species)
        
    except UnboundLocalError:
        print(internal_id)
    E_table['internal_id'] = internal_id
    return E_table
    
#%%
def e_val_calc_batch(subfamily,
                     internal_id_dir,
                     e_table_dir):
    subfamily_filename = subfamily.replace('/','_') 
    print(repeatmasker_filepath)
    if internal_id_dir is None:
        internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    input_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
    global internal_id_tbl
    internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
    if e_table_dir is None:
        e_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
    if not os.path.exists(e_table_dir):
        os.makedirs(e_table_dir)
    output_filepath = f'{e_table_dir}/{subfamily_filename}.e_table.txt'
    operation_log_path = f'{e_table_dir}/op_log.log'
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
        internal_id_list = internal_id_tbl.internal_id.unique()
        id_count=len(internal_id_list)
        logging.info('entry counts: '+str(id_count))
        #print('checkpoint1')
        if id_count > 100:
            output_temp_folder = f'{e_table_dir}/{subfamily_filename}'
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
                    with ProcessPoolExecutor(max_workers=1) as executor:
                        results = executor.map(e_val_calc, id_sublist)

                    for result in results:
                        df_list.append(result)
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
            temp_df.to_csv(output_filepath, sep='\t', index=False)
            logging.info(f'done, saving the table at: {output_filepath}')
            logging.info('cleaning up temporary files')
            shutil.rmtree(output_temp_folder)

        else:
            df_list = list()
            with ProcessPoolExecutor(max_workers=1) as executor:
                results = executor.map(e_val_calc, internal_id_list)
            #
            for result in results:
                    df_list.append(result)
            output_table=pd.concat(df_list)
            output_table.to_csv(output_filepath, sep='\t', index=False)
            logging.info('done, saving the table at: '+output_filepath)
        logging.debug("--- %s seconds ---\n" % (time.time() - start_time))
    else:
        logging.info(subfamily+' already done')
#%%
def main():
    #initialize 
    global repeatmasker_table, divergence_table,target_species, maf_dir, maf_file_prefix
    repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
    divergence_table=divergence_table_prep(divergence_table_filepath)
    #for subfam in ['THE1C']:
    for subfam in subfamily_list:
        e_val_calc_batch(subfam, internal_id_dir, e_table_dir)
#%%
if __name__ == '__main__':
    main()
#%% manual test
from Bio.AlignIO import MafIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ma_mapper import gzmaf
maf_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species/241-mammalian-2020v2b.maf.chr11'
chrom = 'chr11'
target_species='Homo_sapiens'
from collections import Counter, defaultdict
maf_id = f'{target_species}.{chrom}'
mafindex_filepath = f'{maf_file}.mafindex'
index_maf = gzmaf.gzMafIndex(mafindex_filepath, maf_file+'.gz', maf_id) 
# %%
start =[21327004, 21327353]
end =[21327183, 21327573]
strand = '+'
n_strand = -1 if strand == '-' else 1
results =index_maf.get_spliced(start,end,n_strand)
# %%
results
# %%
