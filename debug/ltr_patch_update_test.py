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
#%%

repeatmasker_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'

repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
#%%
def get_base_subfamily(subfamily):
    base = re.match(r'^[A-Za-z0-9]+', subfamily).group()
    # Remove last char only if it's a letter (a-z or A-Z)
    if base[-1].isalpha():
        base = base[:-1]
    return base

def ltr_handler(subfamily):
    global repeatmasker_table
    block_id_list = repeatmasker_table[repeatmasker_table['repName'] == subfamily]['id'].unique()
    base_subfamily = get_base_subfamily(subfamily)
    tagged_blocks = []

    for block_id in block_id_list:
        block_df = repeatmasker_table[repeatmasker_table['id'] == block_id].copy()
        block_df['tag'] = block_df['id'].copy()
        block_df['tag'] = block_df['tag'].astype(str)
        rep_names = list(block_df['repName'])

        # Match LTR variants (e.g. MER83, MER83A, MER83_, etc.)
        ltr_pattern = re.compile(rf"^{re.escape(base_subfamily)}([A-Z]|_)?$")
        ltr_positions = [i for i, name in enumerate(rep_names) if ltr_pattern.match(name)]

        # Match internal elements (e.g. MER83-int, MER83A-int, MER83_-int, etc.)
        int_pattern = re.compile(rf"^{re.escape(base_subfamily)}([A-Z]|_)?-int$")
        int_positions = [i for i, name in enumerate(rep_names) if int_pattern.match(name)]

        # Tag LTRs based on relative position to -int
        if int_positions:
            for i in ltr_positions:
                row_label = block_df.index[i]
                if any(pos > i for pos in int_positions):
                    block_df.at[row_label, 'tag'] = f'{block_id}_bINT'
                elif any(pos < i for pos in int_positions):
                    block_df.at[row_label, 'tag'] = f'{block_id}_aINT'
                else:
                    block_df.at[row_label, 'tag'] = f'{block_id}_nINT'
        else:
            # No -int found, tag all LTRs as nINT
            for i in ltr_positions:
                row_label = block_df.index[i]
                block_df.at[row_label, 'tag'] = f'{block_id}_nINT'

        # Filter only for the current subfamily (e.g., just MER83)
        block_df = block_df[block_df['repName'] == subfamily]
        tagged_blocks.append(block_df)

    subfam_table = pd.concat(tagged_blocks)
    subfam_table['tag'] = subfam_table['tag'].astype(str)  # Ensure string type for downstream operations
    return subfam_table

#%%
test_df=ltr_handler('THE1C')
#%%
#%%
LTR_list=repeatmasker_table[repeatmasker_table['repClass'].str.contains('LTR', na=False)]['repName'].unique()
# Step 1: Convert to set for fast lookup
LTR_set = set(LTR_list)
ltr_names = [name for name in LTR_list if not name.endswith('-int')]
ltr_only = repeatmasker_table[repeatmasker_table['repClass'].str.contains('LTR', na=False)].copy()
ltr_only = ltr_only[ltr_only['genoName'].str.match(r'chr[\dXY]+$')]
cols_to_keep = ['genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass']
ltr_only = ltr_only[cols_to_keep]
ltr_by_chr_strand = {
    (chrom, strand): df.sort_values('genoStart' if strand == '+' else 'genoEnd')
    for (chrom, strand), df in ltr_only.groupby(['genoName', 'strand'])
}
#%%
max_gap = 15000
max_span = 15000
min_span = 1000
intact_ltrs = []

subfam_table=repeatmasker_table[repeatmasker_table['repName']=='THE1C']
subfam_table = subfam_table[subfam_table['genoName'].str.match(r'chr[\dXY]+$')]
cols_to_keep = ['genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass']
subfam_trimmed = subfam_table[cols_to_keep]
consumed_indices = set()
for run_id in range(len(subfam_table)):
    
    row_idx = subfam_trimmed.iloc[run_id].name
    print(run_id,row_idx)
    current_row = ltr_only[ltr_only.index==int(row_idx)]

    # Skip if already part of an earlier trio
    if row_idx in consumed_indices:
        continue

    chrom = current_row['genoName'].values[0]
    strand = current_row['strand'].values[0]
    anchor_pos = current_row['genoStart'].values[0] if strand == '+' else current_row['genoEnd'].values[0]
    ltr_class = current_row['repClass'].values[0]

    try:
        candidates = ltr_by_chr_strand[(chrom, strand)]
    except KeyError:
        continue  # no LTRs on this chrom/strand

    # Get nearby entries within 10kb window
    if strand == '+':
        window_df = candidates[
            (candidates['repClass'] == ltr_class) &
            (candidates['genoStart'] >= anchor_pos) &
            (candidates['genoStart'] <= anchor_pos + max_gap)
        ].copy()
        window_df = window_df.sort_values('genoStart')
    else:
        window_df = candidates[
            (candidates['repClass'] == ltr_class) &
            (candidates['genoEnd'] <= anchor_pos) &
            (candidates['genoEnd'] >= anchor_pos - max_gap)
        ].copy()
        window_df = window_df.sort_values('genoEnd', ascending=False)
    if len(window_df) <3:
        continue

    # Try to find LTR -int LTR pattern
    for i in range(len(window_df) - 2):
        a, b, c = window_df.iloc[i], window_df.iloc[i+1], window_df.iloc[i+2]
        idx_a, idx_b, idx_c = a.name, b.name, c.name

        if idx_a in consumed_indices or idx_b in consumed_indices or idx_c in consumed_indices:
            continue

        # Check repName pattern: LTR - -int - LTR
        if not (
            not a['repName'].endswith('-int') and
            b['repName'].endswith('-int') and
            not c['repName'].endswith('-int')
        ):
            continue

        # Strand-aware span
        start = min(a['genoStart'], b['genoStart'], c['genoStart'])
        end   = max(a['genoEnd'],   b['genoEnd'],   c['genoEnd'])
        span  = end - start

        if span > max_span or span < min_span:
            continue

        # Passed all: record and mark as used
        trio = window_df.loc[[idx_a, idx_b, idx_c]]
        ltr_pair_id = f"{a['repName']}_{run_id}"
        trio['is_int'] = trio['repName'].str.endswith('-int')
        trio['ltr_pair_id'] = ltr_pair_id
        trio['rmsk_index'] = trio.index

        if a['strand'] == '+':
            trio.loc[idx_a, 'tag'] = f"{ltr_pair_id}_bINT"
            trio.loc[idx_b, 'tag'] = f"{ltr_pair_id}_INT"
            trio.loc[idx_c, 'tag'] = f"{ltr_pair_id}_aINT"
        else:
            # Reverse assignment for reverse strand
            trio.loc[idx_a, 'tag'] = f"{ltr_pair_id}_aINT"
            trio.loc[idx_b, 'tag'] = f"{ltr_pair_id}_INT"
            trio.loc[idx_c, 'tag'] = f"{ltr_pair_id}_bINT"

        intact_ltrs.append(trio)
        consumed_indices.update([idx_a, idx_b, idx_c])
        break  # move to next anchor LTR
# %%
intact_ltr_df = pd.concat(intact_ltrs)
tag_table=intact_ltr_df[['rmsk_index', 'tag']]
# %%
subfam_table_test=subfam_table.merge(tag_table,right_index=True, left_index=True, how='left').drop(columns=['rmsk_index'])
subfam_table_test['tag'] = subfam_table_test['tag'].fillna('nINT')
subfam_table_test['tag'] = subfam_table_test['id'].astype(str) + '_' + subfam_table_test['tag']
# %%
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
    else:
        print(f'file already exists')   
    return internal_id_table
#%%
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

#%%
def assign_internal_ids(df):
    df = df.copy()

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

#%%

#%%
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
#%%
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
#%%
def tag_completeness(df):
    def completeness_mask(group):
        return (group['repStart'] == 1) & (group['repLeft'] == 0)

    # Apply the completeness_mask per group by 'id' and assign back with transform
    df['is_complete'] = df.groupby('id', group_keys=False).apply(completeness_mask)

    return df

#%%


def ltr_handler(subfamily):
    global repeatmasker_table, ltr_by_chr_strand  # use the global precomputed dict

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

        current_row = ltr_only.loc[[row_idx]]  # you might want to access global ltr_only as well

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

# %%
ltr_handler('THE1C')
# %%
annotation_mending('THE1C', output_dir='./',overlap_threshold=0.1)
# %%
def apply_ltrint_correction(df, age_table_dir, mode):
    if df['ltrpair_id'].isnull().all():
        return df
    global ltr_only
    df_extend=df.merge(ltr_only, left_on='rmsk_index', right_index = True)
    grouped=df_extend.dropna(subset=['ltrpair_id']).groupby('ltrpair_id')

    corrected_ages = {}
    for group in grouped:
        if len(group) <2:
            continue
        chrom = group['genoName'].iloc[0]
        strand = group['strand'].iloc[0]
        start = group['genoStart'].min()
        end = group['genoEnd'].max()
        ltrpair_id = group['ltrpair_id'].iloc[0]
        family = group['repName'].iloc[0]
        te_class = group['repClass'].iloc[0]

        candidates = ltr_by_chr_strand.get((chrom, strand))
        if mode == 'ltr-int':
            target_name = family + '-int'
            internal_hit = candidates[
                (candidates['genoStart'] >= start) &
                (candidates['genoEnd'] <= end) &
                (candidates['repName']==target_name)
            ]
        elif mode == 'greedy':
            internal_hit = candidates[
                (candidates['genoStart'] >= start) &
                (candidates['genoEnd'] <= end) &
                (candidates['repClass'] == te_class)
            ]
        else:
            raise ValueError(f"Unsupported correction mode: {mode}")
        
        if len(internal_hit) != 1:
            raise ValueError(
                f"[INT MISMATCH] ltrpair_id '{ltrpair_id}' expected 1 internal, found {len(internal_hit)}"
            )
        int_entry = internal_hit.iloc[0]
        int_name = int_entry['repName']
        int_idx = int_entry.name  # rmsk_index

        age_path = f"{age_table_dir}/{int_name}.teatime.txt"
        try:
            int_ages = pd.read_csv(age_path, sep='\t', index_col=0)
            if int_idx not in int_ages.index:
                print(f"[WARNING] rmsk_index '{int_idx}' not found in {int_name}.teatime.txt. Falling back to LTR age only.")
                corrected_ages[ltrpair_id] = max(group['age'])
                continue
        except FileNotFoundError:
            print(f"[WARNING] Age file for internal '{int_name}' not found at {age_path}. Falling back to LTR age only.")
            corrected_ages[ltrpair_id] = max(group['te_age'])
            continue
        
        
        int_age = int_ages.loc[int_idx]['te_age']
        ltr_ages = group['te_age'].tolist()
        all_ages = ltr_ages + [int_age]
        corrected_ages[ltrpair_id] = max(all_ages)
    df['corrected_age'] = df['ltrpair_id'].map(corrected_ages).fillna(df['te_age'])
    return df
#%%
def unpack_internal_id(df):
    def parse_id(s):
        if '__' in s:
            # Pattern for IDs with __ separating status (and no explicit ltrpair_id)
            pattern = r'(?P<family>[^_]+)_(?P<ltr_id>\d+)_(?P<internal_index>\d+)__(?P<status>.+)'
        else:
            # Pattern for IDs with ltrpair_id but no __ separator
            pattern = r'(?P<family>[^_]+)_(?P<ltr_id>\d+)_(?P<internal_index>\d+)_(?P<ltrpair_id>\d+)_(?P<status>.+)'

        match = re.match(pattern, s)
        return match.groupdict() if match else {}


    unpacked = df['internal_id'].apply(parse_id).apply(pd.Series)
    return pd.concat([df, unpacked], axis=1)

def apply_strict_correction(df):
    if df['ltrpair_id'].isnull().all():
        return df
    paired = df.dropna(subset=['ltrpair_id']).copy()
    max_ages = paired.groupby('ltrpair_id')['age'].max()
    df['corrected_age'] = df['ltrpair_id'].map(max_ages)
    df['corrected_age'] = df['corrected_age'].fillna(df['age'])
    return df

def ltr_correction(subfamily, output_dir, internal_id_dir, age_table_dir, correction_mode):
    #load table
    output_filepath = f"{output_dir}/{subfamily}.ltr_fix.txt"
    if (os.path.isfile(output_filepath) == False):
        print('start',subfamily)

        internal_id_filepath = f'{internal_id_dir}/{subfamily}.internal_id.txt'
        internal_id=pd.read_csv(internal_id_filepath, sep='\t', index_col = 0)

        age_table_filepath = f'{age_table_dir}/{subfamily}.teatime.txt'
        age_table = pd.read_csv(age_table_filepath, sep='\t')
        internal_ids_age = internal_id.merge(age_table, on='internal_id')
        internal_ids_age = unpack_internal_id(internal_ids_age)

        if correction_mode == 'strict':
            updated_internal_ids_age = apply_strict_correction(internal_ids_age)
        else:
            updated_internal_ids_age = apply_ltrint_correction(internal_ids_age,correction_mode)
        
        ltr_patched_df = updated_internal_ids_age[['internal_id', 'te_age','corrected_age']]
        ltr_patched_df.to_csv(output_filepath, sep='\t', index=False)
        print('done',subfamily)
    else:
        print('already done', subfamily)

#%%
subfamily = 'THE1C'
internal_id_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline'
internal_id_filepath = f'{internal_id_dir}/{subfamily}.internal_id.txt'
internal_id=pd.read_csv(internal_id_filepath, sep='\t', index_col = 0)
age_table_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/pipeline'
age_table_filepath = f'{age_table_dir}/{subfamily}.teatime.txt'
age_table = pd.read_csv(age_table_filepath, sep='\t')
internal_ids_age = internal_id.merge(age_table, on='internal_id')
internal_ids_age = unpack_internal_id(internal_ids_age)
#%%
internal_ids_age = unpack_internal_id(internal_ids_age)
# %%
