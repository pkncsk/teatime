#%% 
import pandas as pd
import re
import os
import re
from collections import defaultdict
from typing import List
import argparse
#%%
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table

def tag_completeness(df):
    def completeness_mask(group):
        return (group['repStart'] == 1) & (group['repLeft'] == 0)

    # Apply the completeness_mask per group by 'id' and assign back with transform
    df['is_complete'] = df.groupby('id', group_keys=False).apply(completeness_mask,include_groups=False)

    return df



def ltr_handler(subfamily):
    global repeatmasker_table, ltr_by_chr_strand, ltr_only  # use the global precomputed dict

    # Filter subfamily repeats on valid chromosomes & needed columns
    subfam_table = repeatmasker_table[
        (repeatmasker_table['repName'] == subfamily) &
        (repeatmasker_table['genoName'].str.match(r'chr[\dXY]+$'))
    ].copy()

    cols_to_keep = ['id', 'repName', 'repClass', 'genoStart', 'genoEnd', 'repStart', 'repEnd', 'repLeft', 'strand']
    subfam_trimmed = subfam_table[cols_to_keep]

    consumed_indices = set()
    intact_ltrs = []

    max_gap = 15000
    max_span = 15000
    min_span = 1000 
    subfam_table['ltr_pair_id'] = pd.Series(dtype='Int64') 
    subfam_table['is_intact_ltr'] = False
    for run_id in range(len(subfam_trimmed)):
        
        row_idx = subfam_trimmed.iloc[run_id].name
        print(run_id, row_idx)
        if row_idx in consumed_indices:
            continue

        current_row = ltr_only.loc[[row_idx]]  

        chrom = current_row['genoName'].values[0]
        strand = current_row['strand'].values[0]
        anchor_pos = current_row['genoStart'].values[0] if strand == '+' else current_row['genoEnd'].values[0]
        ltr_class = current_row['repClass'].values[0]

        try:
            candidates = ltr_by_chr_strand[(chrom, strand)]
        except KeyError:
            continue

        # Select window of candidates within max_gap around anchor_pos
        if strand == '+':
            window_df = candidates[
                (candidates['repClass'] == ltr_class) &
                (candidates['genoStart'] >= anchor_pos) &
                (candidates['genoStart'] <= anchor_pos + max_gap)
            ].copy().sort_values('genoStart')
        else:
            window_df = candidates[
                (candidates['repClass'] == ltr_class) &
                (candidates['genoEnd'] <= anchor_pos) &
                (candidates['genoEnd'] >= anchor_pos - max_gap)
            ].copy().sort_values('genoEnd', ascending=False)

        if len(window_df) < 3:
            continue

        # Scan window for LTR - int - LTR trio pattern
        for i in range(len(window_df) - 2):
            a, b, c = window_df.iloc[i], window_df.iloc[i + 1], window_df.iloc[i + 2]
            idx_a, idx_b, idx_c = a.name, b.name, c.name

            if idx_a in consumed_indices or idx_b in consumed_indices or idx_c in consumed_indices:
                continue

            # Check repName pattern: LTR - int - LTR
            if not (
                not a['repName'].endswith('-int') and
                b['repName'].endswith('-int') and
                not c['repName'].endswith('-int')
            ):
                continue

            start = min(a['genoStart'], b['genoStart'], c['genoStart'])
            end = max(a['genoEnd'], b['genoEnd'], c['genoEnd'])
            span = end - start

            if span > max_span or span < min_span:
                continue
            
            valid_indices = [idx for idx in (idx_a, idx_b, idx_c) if idx in subfam_table.index]
            if valid_indices:
                
                subfam_table.loc[valid_indices, 'is_intact_ltr'] = True
                
                subfam_table.loc[valid_indices, 'ltr_pair_id'] = int(run_id)
            else:
                continue
            
            consumed_indices.update([idx_a, idx_b, idx_c])

    return subfam_table



def two_frags_handler(df, common_length, overlap_threshold):
    # Initialize columns if missing
    for col in ['is_complete', 'is_join', 'is_overlap', 'join_order']:
        if col not in df.columns:
            if col == 'join_order':
                df[col] = -1  # default invalid order
            else:
                df[col] = False
    two_ids = df[df['is_two']]['id'].unique()

    for block_id in two_ids:
        block_df = df[df['id'] == block_id].copy()
        strand = block_df['strand'].iloc[0]
        if strand == '-':
            block_df = block_df.sort_index(ascending=False)
        idxs = block_df.index
        first = block_df.iloc[0]
        second = block_df.iloc[1]
        
        if block_df['is_complete'].any():
            for idx in idxs:
                if not df.at[idx, 'is_complete']:
                    df.at[idx, 'is_singleton'] = True
        elif first['repEnd'] < second['repStart'] and \
             ((first['repLeft'] >= second['repEnd'] - second['repStart']) or second['repLeft'] == 0):
            # Simple join
            df.loc[[first.name, second.name], 'is_join'] = True
            df.loc[[first.name, second.name], 'join_order'] = [0, 1]
        elif (first['repEnd'] - second['repStart'] < overlap_threshold * common_length) and \
             (first['repStart'] < second['repStart']):
            # Overlap join
            df.loc[[first.name, second.name], ['is_join', 'is_overlap']] = True
            df.loc[[first.name, second.name], 'join_order'] = [0, 1]
    return df


def twoplus_frags_handler(df, common_length, overlap_threshold):

    twoplus_ids = df[df['is_twoplus']]['id'].unique()
    for col in ['is_complete', 'is_join', 'is_overlap', 'join_order']:
        if col not in df.columns:
            if col == 'join_order':
                df[col] = -1  # default invalid order
            else:
                df[col] = False
    joined = set()
    for block_id in twoplus_ids:

        block_df = df[df['id'] == block_id].copy()

        strand = block_df['strand'].iloc[0]
        block_df = block_df.sort_index(ascending=(strand != '-'))
        join_order_counter = 0
        for i in range(block_df.shape[0]-1):
            f1 = block_df.iloc[i]
            f2 = block_df.iloc[i+1]
            idxs = [f1.name, f2.name]

            if f1.name in joined or f2.name in joined:
                continue
    
            # Skip if any is already complete
            if df.loc[f1.name, 'is_complete']==True:
                continue

            if f1['repEnd'] < f2['repStart'] and \
               (f1['repLeft'] >= f2['repEnd'] - f2['repStart'] or f2['repLeft'] == 0):
                df.loc[idxs, 'is_join'] = True
                df.loc[idxs, 'join_order'] = [join_order_counter, join_order_counter + 1]
                joined.update(idxs)
                join_order_counter += 2

            elif (f1['repEnd'] - f2['repStart'] < overlap_threshold * common_length) and \
                 (f1['repStart'] < f2['repStart']):
                df.loc[idxs, ['is_join', 'is_overlap']] = True
                df.loc[idxs, 'join_order'] = [join_order_counter, join_order_counter + 1]
                joined.update(idxs)
                join_order_counter += 2

        remaining_idxs = set(block_df.index) - joined
        df.loc[list(remaining_idxs), 'is_singleton'] = True
    return df

def assign_internal_ids(df):
    df = df.copy()
    repname = df['repName'].unique()[0]
    te_class = df['repClass'].unique()[0]
    if 'LTR' in te_class and '-int' not in repname:    
        ltr_counts = df['ltr_pair_id'].value_counts()
        valid_ltr = ltr_counts[ltr_counts == 2].index
        df['valid_ltrpair'] = df['ltr_pair_id'].isin(valid_ltr)

        # Resolve ltr_pair_id per group (id, internal_index)
        group_keys = ['id', 'internal_index']
        resolved_ltr = (
        df[df['valid_ltrpair']]
        .groupby(group_keys)['ltr_pair_id']
        .first()
        .astype(str)  # convert to string early
        )
        df = df.join(resolved_ltr.rename('resolved_ltr_pair_id'), on=group_keys)
    else:
        df[['resolved_ltr_pair_id', 'valid_ltrpair']] = pd.NA
    # Fill missing with empty string for ease
    df['resolved_ltr_pair_id'] = df['resolved_ltr_pair_id'].fillna('')

    # Build internal_id vectorized
    rep = df['repName'].astype(str)
    block_id = df['id'].astype(str)
    internal_idx = df['internal_index'].astype(int).astype(str)
    ltr_id = df['resolved_ltr_pair_id'].astype(str)
    status = df['status'].astype(str)

    has_ltr = ltr_id != ''

    df['internal_id'] = ''
    df.loc[has_ltr, 'internal_id'] = rep[has_ltr] + '_' + block_id[has_ltr] + '_' + internal_idx[has_ltr] + '_' + ltr_id[has_ltr] + '_' + status[has_ltr]
    df.loc[~has_ltr, 'internal_id'] = rep[~has_ltr] + '_' + block_id[~has_ltr] + '_' + internal_idx[~has_ltr] + '__' + status[~has_ltr]

    df.drop(columns=['valid_ltrpair', 'resolved_ltr_pair_id'], inplace=True)
    return df


def assign_internal_index(df: pd.DataFrame) -> pd.DataFrame:
    block_id_list = df['id'].unique()
    for block_id in block_id_list:
        block_df = df[df['id'] == block_id]
        strand = block_df['strand'].iloc[0]

        sorted_block = block_df.sort_index(ascending=(strand != '-'))

        indices = []
        current_index = 0
        prev_joined = False

        for i, (idx, row) in enumerate(sorted_block.iterrows()):
            if i == 0:
                indices.append(current_index)
                prev_joined = row['is_join']
            else:
                if row['is_join'] and prev_joined:
                    indices.append(current_index)
                else:
                    current_index += 1
                    indices.append(current_index)
                prev_joined = row['is_join']

        df.loc[sorted_block.index, 'internal_index'] = indices

    return df

def assign_status(df: pd.DataFrame) -> pd.DataFrame:
    def get_status(row):
        if row.get('is_complete', False):
            return 'complete'
        elif row.get('is_join', False):
            return 'join'
        else:
            return 'singleton'
    df['status'] = df.apply(get_status, axis=1)
    return df
#%%
#detect repClass
def annotation_mending(subfamily, output_dir,overlap_threshold):
    print(f'process: {subfamily}')
    global repeatmasker_table
    output_filepath = f'{output_dir}/{subfamily}.internal_id.txt'
    if not os.path.exists(output_filepath):
        te_class =repeatmasker_table[repeatmasker_table['repName'] == subfamily]['repClass'].unique()[0]
        if 'LTR' in te_class and '-int' not in subfamily:
            subfam_table = ltr_handler(subfamily)
        else:
            cols_to_keep = ['id', 'repName', 'repClass', 'genoStart', 'genoEnd', 'repStart', 'repEnd', 'repLeft', 'strand']
            subfam_table = repeatmasker_table.loc[repeatmasker_table['repName'] == subfamily, cols_to_keep].copy()
        subfam_table['length'] = subfam_table['genoEnd'] - subfam_table['genoStart']
        subfam_table['repLength'] = subfam_table['repEnd'] - subfam_table['repStart'] + 1
        if (subfam_table['repLeft'] == 0).any():
            common_length = subfam_table[subfam_table.repLeft==0].repEnd.mode()[0]
        else: 
            common_length = (subfam_table.repLeft + subfam_table.repEnd).mode()[0]
        id_counts = subfam_table['id'].value_counts()
        subfam_table['is_singleton'] = subfam_table['id'].map(id_counts) == 1
        subfam_table['is_two'] = subfam_table['id'].map(id_counts) == 2
        subfam_table['is_twoplus'] = subfam_table['id'].map(id_counts) > 2
        subfam_table = tag_completeness(subfam_table)
        subfam_table = two_frags_handler(subfam_table, common_length, overlap_threshold)
        subfam_table = twoplus_frags_handler(subfam_table, common_length, overlap_threshold)
        subfam_table = assign_status(subfam_table)
        subfam_table = assign_internal_index(subfam_table)
        subfam_table = assign_internal_ids(subfam_table)
        internal_id_table = subfam_table[['internal_id']].rename_axis('rmsk_index').reset_index()
        internal_id_table.to_csv(f'{output_dir}/{subfamily}.internal_id.txt', sep='\t', index=False)
        return internal_id_table
    else:
        print(f'file already exists')   
    
# %%
def main(repeatmasker_filepath, subfamily_list, output_dir, overlap_threshold):
    global repeatmasker_table
    repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
    global ltr_by_chr_strand, ltr_only
    # Run this once at the start before calling ltr_handler repeatedly
    ltr_only = repeatmasker_table[
        repeatmasker_table['repClass'].str.contains('LTR', na=False)
    ].copy()

    ltr_only = ltr_only[ltr_only['genoName'].str.match(r'chr[\dXY]+$')]

    cols_to_keep = ['genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass']
    ltr_only = ltr_only[cols_to_keep]

    ltr_by_chr_strand = {
        (chrom, strand): df.sort_values('genoStart' if strand == '+' else 'genoEnd')
        for (chrom, strand), df in ltr_only.groupby(['genoName', 'strand'])
    }
    if subfamily_list is None:
        repname_counts = repeatmasker_table['repName'].value_counts().reset_index()
        repname_counts.columns = ['repName', 'count']
        subfamily_list=repname_counts['repName']
    for subfam in subfamily_list:
        annotation_mending(subfam, output_dir, overlap_threshold)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="annotation mending script")
    parser.add_argument("-r", "--repeatmasker",required=True,
                        help="Path to RepeatMasker output file (.tsv) with modified headers compatible with pandas.read_csv()")
    parser.add_argument("-o","--output", required=True,
                        help="Directory for internal ID tables")
    parser.add_argument("-s","--subfamily_list",nargs="+", default=None,
                        help="List of subfamily names to process (optional). If not given, the script will run through all subfamilies available on the RepeatMasker output table.",)
    parser.add_argument("-S", "--subfamily_file",  type=str,
                        help="Optional path to a text file with a list of subfamilies (one per line).")
    parser.add_argument("-t","--overlap_threshold", type=float, default=0.1,
                        help="Maximum overlap fraction allowed for merging fragments (default: 0.1)" )
    args = parser.parse_args()
    subfamily_list = None
    if args.subfamily_list:
        subfamily_list = args.subfamily_list
    elif args.subfamily_file:
        with open(args.subfamily_file) as f:
            subfamily_list = [line.strip() for line in f if line.strip()]
    main(
        repeatmasker_filepath=args.repeatmasker,
        subfamily_list=subfamily_list,
        output_dir=args.output,
        overlap_threshold=args.overlap_threshold,
    )
