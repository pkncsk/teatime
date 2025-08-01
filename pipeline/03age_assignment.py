#%%
import pandas as pd
from Bio.AlignIO import MafIO
import numpy as np
import os
import re
import math
from itertools import repeat
from concurrent.futures import ProcessPoolExecutor
from Bio.Seq import Seq
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
import argparse
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
    from typing import List, Tuple
    def trim_overlapping_fragments(starts: List[int], ends: List[int]) -> Tuple[List[int], List[int]]:
        """
        Trim overlapping TE fragments by adjusting the start positions to avoid base duplication.

        This is used to **prevent per-base duplication** when processing overlapping TE fragments,
        especially in multi-species alignments. It ensures that each base in the reference genome
        is only used once, even if it is part of multiple fragment calls.

        This function is necessary for **multi-species sequence alignments** where 
            overlapping TE fragments (e.g. from RepeatMasker) can introduce duplicated sequence 
            content if not properly merged. Without this step, fragments with physical overlap 
            in one species can be duplicated in the alignment, causing:

                - False inflation of sequence length
                - Artificially high TE conservation across species
                - Misleading phylogenetic signals
                - Incorrect ancestral reconstruction

            In contrast, if working in **single-genome**, retaining overlaps 
            may preserve internal deletions or tandem repeat events and should be handled differently.

        Parameters:
        -----------
        starts : List[int]
            List of start positions.
        ends : List[int]
            List of end positions (must be same length as `starts`).

        Returns:
        --------
        Tuple[List[int], List[int]]
            Adjusted (trimmed) start and end positions with no overlapping bases.
        """
        if not starts or not ends or len(starts) != len(ends):
            raise ValueError("Starts and ends must be non-empty and of equal length.")

        trimmed_starts = [starts[0]]
        trimmed_ends = [ends[0]]

        for i in range(1, len(starts)):
            prev_end = trimmed_ends[-1]
            curr_start = starts[i]
            curr_end = ends[i]

            # Trim current start to be >= previous end
            adjusted_start = max(curr_start, prev_end)
            if adjusted_start >= curr_end:
                # Skip completely overlapped region
                continue
            trimmed_starts.append(adjusted_start)
            trimmed_ends.append(curr_end)

        return trimmed_starts, trimmed_ends

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
    
    starts, ends = trim_overlapping_fragments(starts, ends)
    fetched = list(self.search(starts, ends))
    #return fetched
    expected_letters = sum(end - start for start, end in zip(starts, ends))
    if len(fetched) == 0:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})

    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    
    exon_positions = []
    for start, end in zip(starts, ends):
        
        exon_positions.extend(range(start, end))
    exon_positions = sorted(set(exon_positions))  # Ensure unique & sorted
    
    total_rec_length = 0
    ref_first_strand = None
    split_by_position = defaultdict(lambda: defaultdict(str))
    for multiseq in fetched:
        if len(multiseq) == 1 and multiseq[0].id == self._target_seqname:
            # Fast-path for human-only alignment block
            rec = multiseq[0]
            rec_start = rec.annotations["start"]
            gapped_seq = rec.seq
            if ref_first_strand is None:
                ref_first_strand = rec.annotations["strand"]
            if "-" not in gapped_seq:
                #print('debug-humanonly: superfast route')
                # Fully ungapped, directly assign bases by offset
                for exonstart, exonend in zip(starts, ends):
                    # Clamp the range to only what's within the alignment block
                    effective_start = max(exonstart, rec_start)
                    effective_end = min(exonend, rec_start + len(gapped_seq))
           
                    if effective_start >= effective_end:
                        # Exon lies completely outside this block — skip
                        continue

                    start_idx = effective_start - rec_start
                    end_idx = effective_end - rec_start
                    for offset, base in zip(range(effective_start, effective_end), gapped_seq[start_idx:end_idx]):
                        split_by_position[rec.id][offset] = base
            else:
                #print('debug-humanonly: fast route')
                # Gapped version: build position map efficiently
                valid_indexes = []
                positions = []
                realpos = rec_start
                for i, base in enumerate(gapped_seq):
                    if base != "-":
                        valid_indexes.append((i, realpos))
                        positions.append(realpos)
                        realpos += 1

                for exonstart, exonend in zip(starts, ends):
                    i_start = bisect_left(positions, exonstart)
                    i_end = bisect_right(positions, exonend - 1)
                 
                    for i in range(i_start, i_end):
                        aln_idx, real_pos = valid_indexes[i]
                    
                        split_by_position[rec.id][real_pos] += gapped_seq[aln_idx]

            continue  # Skip rest of loop for this block
        # Find the human reference sequence
        ref = next((r for r in multiseq if r.id == self._target_seqname), None)
        if ref is None:
            raise ValueError(
                "Did not find %s in alignment bundle" % (self._target_seqname,)
            )
        
        try:
            if ref_first_strand is None:
                ref_first_strand = ref.annotations["strand"]
                if ref_first_strand not in (1, -1):
                    raise ValueError("Strand must be 1 or -1")
            elif ref_first_strand != ref.annotations["strand"]:
                raise ValueError(
                    "Encountered strand='%s' on target seqname, expected '%s'"
                    % (strand, ref_first_strand)
                )
        except KeyError:
            raise ValueError(
                "No strand information for target seqname (%s)" % self._target_seqname
            ) from None

        ref_start = ref.annotations["start"]
        ref_seq = ref.seq

        # Gapped path: build ref_pos_map
        #print("debug-humanref: fast route")
        ref_pos_map = []
        realpos = ref_start
        for i, base in enumerate(ref_seq):
            if base != "-":
                ref_pos_map.append((i, realpos))
                realpos += 1
        
        block_positions  = [pos for _, pos in ref_pos_map]
     
        # Now extract aligned fragments for all sequences using the same alignment coordinates
        for exonstart, exonend in zip(starts, ends):
            i_start = bisect_left(block_positions, exonstart)
            i_end = bisect_right(block_positions, exonend - 1)
   
            if i_start == i_end:
                continue  # This exon not represented in the alignment block

            aln_indexes = [ref_pos_map[i][0] for i in range(i_start, i_end)]
            ref_coords = block_positions[i_start:i_end]  # genomic coords covered in this exon by this block
      
            for rec in multiseq:
                rec_id = rec.id
                gapped_seq = rec.seq
                frag = "".join(gapped_seq[i] for i in aln_indexes)

                for pos, base in zip(ref_coords, frag):
                    split_by_position[rec_id][pos] += base
          
            # Count how many times each seqrec.id appears (detect duplicates)
            id_counts = Counter(seqrec.id for seqrec in multiseq)
            duplicate_ids = {seq_id for seq_id, count in id_counts.items() if count > 1}

            # Process duplicates only
            for seq_id in duplicate_ids:
                for pos, seq_str in split_by_position[seq_id].items():
                    # If sequence length > 1, means concatenation happened (duplication)
                    if len(seq_str) > 1:
                        new_seq = process_sequence_localized(seq_str)
                        split_by_position[seq_id][pos] = new_seq        
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
    target_index = seqid_list.index(self._target_seqname)
    if len(seq_list[target_index].replace("-", "")) != expected_letters:
        raise ValueError(
            "coordinates %s , %s Returning %s letters for target seqname (%s), expected %s"
            % (
                starts,
                ends,
                len(seq_list[target_index].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(seq_list[target_index])
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
#%%
def get_maf_filepath(maf_dir, chrom):
    files = os.listdir(maf_dir)
    maf_filename = f"{chrom}.maf" 
    maf_files = [f for f in files if f.endswith(maf_filename)]

    # Determine the appropriate file to use
    maf_filepath = None
    if any(f.endswith('.maf') for f in maf_files):
        maf_filepath = f"{maf_dir}/{maf_filename}" #.maf is more optimal performance wise 
    elif any(f.endswith('.maf.gz') for f in maf_files):
        maf_filepath = f"{maf_dir}/{maf_filename}.gz" 
    else:
        raise FileNotFoundError(f"No .maf or .maf.gz file found for chromosome {chrom} in {maf_dir}")

    return maf_filepath
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
def calculate_metrics(row, extension_length):
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
    iden_front = round((matched_front / extension_length * 100), 2)
    pgap_front = round((gap_count_front / extension_length * 100), 2)
    E_score_front = '{:0.2e}'.format(E_front)
    iden_back = round((matched_back / extension_length * 100), 2)
    pgap_back = round((gap_count_back / extension_length * 100), 2)
    E_score_back = '{:0.2e}'.format(E_back)
    
    return pd.Series({
        'species':  row['meta_name'],
        'chr_code': row['chr_code'],
        'divergence': row['divergence'],
        'pct_iden': iden,
        'pct_gap': pgap,
        'blast': row['alignment_score'],
        'E_value': E_score,
        'pct_iden_front': iden_front,
        'pct_gap_front':pgap_front,
        'E_val_front':E_score_front,
        'pct_iden_back': iden_back,
        'pct_gap_back':pgap_back,
        'E_val_back': E_score_back
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
              internal_id_dir,
              maf_dir,
              subfamily,
              additional_evalue_cutoff,
              target_species,
              extension_length,
              e_cutoff, 
              ):
    te_age = 0
    segmental = pd.DataFrame()
    unclassified = pd.DataFrame()
    nosig_match = pd.DataFrame()
    internal_id = e_table.internal_id.unique()[0]
    print(internal_id)#debug
    if e_table.shape[0] > 1:
        if additional_evalue_cutoff is not None:
            e_table=e_table[e_table.E_value.astype('float64')<=additional_evalue_cutoff].copy()

        second_pass_tbl=e_table[(e_table['E_val_front'] <= e_cutoff) | (e_table['E_val_back'] <= e_cutoff)]
        if second_pass_tbl.shape[0] > 1:
            te_age=second_pass_tbl.divergence.max()
    
    if e_table.shape[0] == 1:
        subfamily_filename = subfamily.replace('/','_') 
        input_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
        internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
        internal_id=e_table['internal_id'].unique()[0]
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
        maf_filepath = get_maf_filepath(maf_dir, chrom)

        target_chrom = f'{target_species}.{chrom}'
        mafindex_filedir = '.'.join(str.split(maf_filepath, sep ='.')[:-1])
        mafindex_filepath = f'{mafindex_filedir}.mafindex'

        if '.gz' in mafindex_filepath:
            index_maf = mafio_patch.gzMafIndex(mafindex_filepath, maf_filepath, target_chrom)
        else:
            index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
        start_flanked=[min(start_list)-extension_length] + start_list + [max(end_list)]
        end_flanked = [min(start_list)] + end_list + [max(end_list)+extension_length]
        spliced_maf_full =index_maf.get_spliced(start_flanked,end_flanked,strand)
        
        if all(base == 'N' for base in spliced_maf_full['seq'].values[0]):
            unclassified = e_table
            te_age = np.nan
            return internal_id, te_age, segmental, unclassified, nosig_match

        target_seq = np.char.upper(spliced_maf_full[spliced_maf_full.seqid.str.contains(target_species)]['seq'].to_list())[0]
        spliced_maf_full[['alignment_score', 'matched', 'gapped', 'gap_count']] = spliced_maf_full.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][extension_length:-extension_length]), target_seq[extension_length:-extension_length]), axis=1)
        spliced_maf_full[['meta_name', 'chr_code']] = spliced_maf_full['seqid'].str.split('.', n=1, expand=True)
        spliced_maf_age=spliced_maf_full.merge(external_data_table[['meta_name','ungapped_length','divergence']], how='left',on='meta_name')
        spliced_maf_age['seq_length'] = spliced_maf_full['seq'].apply(len) -2*extension_length
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
        
        flanked_tbl = spliced_maf_age.copy()
        flanked_tbl[['alignment_score_front', 'matched_front', 'gapped_front', 'gap_count_front']] = flanked_tbl.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][:extension_length]), target_seq[:extension_length]), axis=1)
        flanked_tbl[['p_value_front', 'E_value_front']] = flanked_tbl.apply(lambda row: pd.Series(BLAST_StoP(
            alignment_score=row['alignment_score_front'],
            m=extension_length,
            n=row['ungapped_length'],
            lambda_=lambda_,
            K=K,
            H=H,
            alpha=alpha,
            beta=beta,
            gapped=row['gapped_front']
        )), axis=1)
        flanked_tbl[['alignment_score_back', 'matched_back', 'gapped_back', 'gap_count_back']] = flanked_tbl.apply(lambda row: affine_count_simple_optimized(np.array(row['seq'][-extension_length:]), target_seq[-extension_length:]), axis=1)
        flanked_tbl[['p_value_back', 'E_value_back']] = flanked_tbl.apply(lambda row: pd.Series(BLAST_StoP(
            alignment_score=row['alignment_score_back'],
            m=extension_length,
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

        summary = {
            'total_count': len(flanked_tbl),
            'match_count': int(flanked_tbl['match_front'].sum() + flanked_tbl['match_back'].sum() -  flanked_tbl['both'].sum()),
            'front_only_count': int(flanked_tbl['front_only'].sum()),
            'back_only_count': int(flanked_tbl['back_only'].sum()),
            'both_count': int(flanked_tbl['both'].sum()),
            'nonmatch_count': int(flanked_tbl['nonmatch'].sum())
                        }
        updated_e_table = e_table.assign(**summary)

        if updated_e_table['match_count'].values[0] > 1:
            te_age = 0
            nosig_match = updated_e_table
        elif updated_e_table['match_count'].values[0] == 1:
            te_age = np.nan
            segmental = updated_e_table
        else:
            te_age = np.nan
            unclassified = updated_e_table
    
    return internal_id, te_age, segmental, unclassified, nosig_match

def filter_e_for_age(subfamily, 
                     internal_id_dir,
                     age_table_dir,
                     e_table_dir,
                     maf_dir,
                     additional_evalue_cutoff,
                     target_species,
                     extension_length,
                     e_cutoff):
    e_table = pd.read_csv(f'{e_table_dir}/{subfamily}.e_table.txt',sep = '\t', low_memory=False)
    output_filepath = f'{age_table_dir}/{subfamily}.teatime.txt'
    if (os.path.isfile(output_filepath) == False):
        print('start',subfamily)
        grouped =e_table.groupby('internal_id', sort=False)
        # Create a list to store the smaller DataFrames
        e_val_table_by_id = [group for _, group in grouped]
       
        with ProcessPoolExecutor(max_workers=80) as executor:
            results = executor.map(filter_e2, 
                                   e_val_table_by_id,
                                   repeat(internal_id_dir), 
                                   repeat(maf_dir),
                                   repeat(subfamily),
                                   repeat(additional_evalue_cutoff),
                                   repeat(target_species), 
                                   repeat(extension_length),
                                   repeat(e_cutoff))           

        id_list = []
        age_list = []
        nosig_match = []
        segmental = []
        unclassified = []
        for idx, result in enumerate(results):
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
def main(internal_id_dir,
         e_table_dir,
         subfamily_list,
         repeatmasker_filepath,
         maf_dir, 
         external_data_table_filepath, 
         age_table_dir,
         additional_evalue_cutoff,
         target_species, 
         extension_length, 
         e_value_cutoff
         ):
    global repeatmasker_table, external_data_table
    repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
    external_data_table=divergence_table_prep(external_data_table_filepath)
    if subfamily_list is None:
        repname_counts = repeatmasker_table['repName'].value_counts().reset_index()
        repname_counts.columns = ['repName', 'count']
        subfamily_list=repname_counts[repname_counts['count']<1000]['repName'].unique()
    for subfamily in subfamily_list:
        filter_e_for_age(
            subfamily=subfamily,
            internal_id_dir=internal_id_dir,
            age_table_dir=age_table_dir, 
            e_table_dir=e_table_dir,
            maf_dir=maf_dir,
            additional_evalue_cutoff=additional_evalue_cutoff,
            target_species=target_species, 
            extension_length=extension_length,
            e_cutoff=e_value_cutoff)
# %%
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign TE ages using internal IDs and E-values.")
    parser.add_argument('-d', '--internal_id_dir', required=True,
                        help="Directory containing internal ID tables")
    parser.add_argument('-e', '--e_table_dir', required=True,
                        help="Directory containing E-value tables")
    parser.add_argument('-r', '--repeatmasker', required=True,
                        help="Path to RepeatMasker output file (.tsv) with modified headers compatible with pandas.read_csv()")
    parser.add_argument('-m', '--maf_dir', required=True,
                        help="Directory containing multispecies alignment (MAF) files")
    parser.add_argument('-x', '--external_data_table', required=True,
                        help="Path to external data table containing genome ID, species divergence, and aligned genome size")
    parser.add_argument('-t', '--target_species', default='hg38',
                        help="Target/reference genome ID (default: hg38)")
    parser.add_argument('-l', '--extension_length', type=int, default=5000,
                        help="Length of flanking sequence to extract (default: 5000 bp)")
    parser.add_argument('-c', '--e_value_cutoff', type=float, default=1e-3,
                        help="E-value threshold for 0 MYA region filtering (default: 1e-3)")
    parser.add_argument('-a', '--additional_evalue_cutoff', type=float, default=None,
                        help="Optional second E-value cutoff for more stringent filtering")
    parser.add_argument('-o', '--output', required=True,
                        help="Output directory for TE age tables")
    parser.add_argument('-s', '--subfamily_list', nargs='+', default=None,
                        help="List of subfamily names to process (optional). If not given, the script will run through all subfamilies available on the RepeatMasker output table")
    parser.add_argument("-S", "--subfamily_file",  type=str,
                        help="Optional path to a text file with a list of subfamilies (one per line).")
    args = parser.parse_args()
    subfamily_list = None
    if args.subfamily_list:
        subfamily_list = args.subfamily_list
    elif args.subfamily_file:
        with open(args.subfamily_file) as f:
            subfamily_list = [line.strip() for line in f if line.strip()]
    main(
        internal_id_dir=args.internal_id_dir,
        e_table_dir=args.e_table_dir,
        subfamily_list=subfamily_list,
        repeatmasker_filepath=args.repeatmasker,
        maf_dir=args.maf_dir,
        external_data_table_filepath=args.external_data_table,
        age_table_dir=args.output,
        additional_evalue_cutoff=args.additional_evalue_cutoff,
        target_species=args.target_species,
        extension_length=args.extension_length,
        e_value_cutoff=args.e_value_cutoff
    )
#%%
#%%