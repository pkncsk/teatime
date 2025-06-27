#%% 
import pandas as pd
import itertools
from concurrent.futures import ProcessPoolExecutor
from ma_mapper import utility
import re
#%% INPUT PARAMETERS
output_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/annotation/' 
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility
global repeatmasker_table
repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
subfamily = 'LTR107_Mam'
#%%
def get_base_subfamily(subfamily):
    base = re.match(r'^[A-Za-z0-9]+', subfamily).group()
    # Remove last char only if it's a letter (a-z or A-Z)
    if base[-1].isalpha():
        base = base[:-1]
    return base


def ltr_handler(subfamily):
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


def tag_singletons(df):
    tag_counts = df['tag'].value_counts()
    df['tag'] = df.apply(
        lambda row: f"{row['tag']}_singleton" if tag_counts[row['tag']] == 1 else row['tag'],
        axis=1
    )
    return df

from numpy import diff, arctan, rad2deg
import numpy as np
def dynamic_nearcomplete_threshold(subfam_table, bins=100, angle_cutoff=-45, default_threshold=95):
    """
    Estimate near-complete cutoff (% of common length) using slope change in cumulative distribution.

    Parameters:
    - subfam_table: DataFrame with columns ['repStart', 'repEnd', 'repLeft']
    - bins: Number of histogram bins (default=100)
    - angle_cutoff: Angle in degrees to detect sharp drop (default=-45°)
    - default_threshold: Fallback if no sharp drop found (default=95%)

    Returns:
    - cutoff_pct: Estimated % of common length for near-complete threshold
    """

    # Step 1: Compute common length
    if subfam_table[subfam_table.repLeft == 0].shape[0] != 0:
        common_length = subfam_table[subfam_table.repLeft == 0]['repEnd'].mode()[0]
    else:
        common_length = (subfam_table['repEnd'] + subfam_table['repLeft']).mode()[0]

    # Step 2: Compute % of common length
    rep_length = subfam_table['repEnd'] - subfam_table['repStart'] + 1
    pct_length = (rep_length / common_length * 100).clip(0, 100)

    # Step 3: Histogram and cumulative sum
    counts, bin_edges = np.histogram(pct_length, bins=bins, range=(0, 100))
    accum = np.cumsum(counts)

    # Step 4: Normalize for slope analysis
    x_vals = bin_edges[:-1]
    y_vals = accum
    x_norm = (x_vals - x_vals.min()) / (x_vals.max() - x_vals.min())
    y_norm = (y_vals - y_vals.min()) / (y_vals.max() - y_vals.min())

    # Step 5: Calculate slope/angle
    slopes = diff(y_norm) / diff(x_norm)
    angles = rad2deg(arctan(slopes))

    # Step 6: Find sharp drop
    sharp_turn_idx = np.where(angles < angle_cutoff)[0]
    if len(sharp_turn_idx) > 0:
        cutoff_pct = x_vals[sharp_turn_idx[0]]
    else:
        cutoff_pct = default_threshold

    return cutoff_pct/100

def two_frags_handler(block_df, common_length):
    """
    Process a two-fragment block:
    - If one is _complete, tag the other as _close.
    - If they can be joined, tag as _open/_close.
    - If small overlap, still allow join.
    - Otherwise, tag both as _close.

    Returns:
        updated_df: The modified DataFrame
        category: A string indicating which list to append to
    """
    first = block_df.iloc[0]
    second = block_df.iloc[1]
    first_label = first.name
    second_label = second.name

    if block_df['tag'].str.contains('_complete').any():
        incomplete_rows = block_df[~block_df['tag'].str.contains('_complete')]
        block_df.loc[incomplete_rows.index, 'tag'] += '_close'
        return block_df, 'one_frag_with_complete'

    elif first['repEnd'] < second['repStart'] and (
        (first['repLeft'] >= second['repEnd'] - second['repStart']) or second['repLeft'] == 0
    ):
        block_df.at[first_label, 'tag'] += '_open'
        block_df.at[second_label, 'tag'] += '_close'
        return block_df, 'two_frags_simple_join'

    elif (first['repEnd'] - second['repStart'] < 0.1 * common_length) and (
        first['repStart'] < second['repStart']
    ):
        block_df.at[first_label, 'tag'] += '_overlap_open'
        block_df.at[second_label, 'tag'] += '_overlap_close'
        return block_df, 'two_frags_overlap_join'

    else:
        block_df.at[first_label, 'tag'] += '_close'
        block_df.at[second_label, 'tag'] += '_close'
        return block_df, 'two_frags_non_join'


#%% mocked table check
#%%
#detect repClass
te_class =repeatmasker_table[repeatmasker_table['repName'] == subfamily]['repClass'].unique()[0]
if 'LTR' in te_class and '-int' not in subfamily:
    subfam_table = ltr_handler(subfamily)
else:
    subfam_table = repeatmasker_table[repeatmasker_table['repName'] == subfamily].copy()
    subfam_table['tag'] = subfam_table['id'].copy()
    subfam_table['tag'] = subfam_table['id'].astype(str)
subfam_table['length'] = subfam_table['genoEnd'] - subfam_table['genoStart']
subfam_table['repLength'] = subfam_table['repEnd'] - subfam_table['repStart'] + 1
# Calculate common length
if subfam_table[subfam_table.repLeft==0].shape[0] != 0:
    common_length = subfam_table[subfam_table.repLeft==0].repEnd.mode()[0]
else: 
    common_length = (subfam_table.repLeft + subfam_table.repEnd).mode()[0]
nearcomplete_cutoff = dynamic_nearcomplete_threshold(subfam_table)
print(f"Estimated near-complete threshold: {nearcomplete_cutoff:.1f}% of common length")
#%%
#separate singletons from multis
subfam_table = tag_singletons(subfam_table)
#%%
singleton_table=subfam_table[subfam_table['tag'].str.contains('singleton', na=False)]
multis_table = subfam_table[~subfam_table['tag'].str.contains('singleton', na=False)]
#%%
all_complete_df=[]
partial_complete_df = []
block_id_list = multis_table['tag'].unique()
for block_id in block_id_list:
    block_df = multis_table[multis_table['tag']==block_id].copy()
    block_df['tag'] = block_df['tag'].astype(str)
    #update tag for complete and nearcomplete
    mask = (block_df['repStart'] == 1) & (block_df['repLeft'] == 0)
    block_df.loc[mask, 'tag'] = block_df.loc[mask, 'tag'] + '_complete'
    if block_df['tag'].str.contains('_complete').all():
        all_complete_df.append(block_df)
    else:
        partial_complete_df.append(block_df)
    #print(block_df[['id','genoStart','genoEnd','strand', 'repStart', 'repEnd', 'repLeft','tag']])

#separate complete block from calculation
if all_complete_df:
    all_complete_multis_table = pd.concat(all_complete_df)
else:
    all_complete_multis_table = pd.DataFrame()
#%%
one_frag_with_complete = []
two_frags_simple_join = []
more_than_two_frag = []
two_frags_overlap_join = []
two_frags_non_join = []
all_tags_pattern = '_open|_mid|_close|_singleton|_test|_complete' 
for block_df in partial_complete_df:
    print(block_df['id'].unique()[0])
    if block_df['strand'].unique()[0] == '-':
        block_df = block_df.sort_index(ascending=False)
    row_number = block_df.shape[0]
    if row_number >2:
        while not block_df['tag'].str.contains(all_tags_pattern).all():
            untagged = block_df[~block_df['tag'].str.contains(all_tags_pattern)].copy()
            # Limit untagged view to the top-most untagged segment (no mixing with already-tagged rows later in block)
            start_idx = untagged.index[0]
            following_rows = block_df.loc[start_idx:]

            # Cut off at first already-tagged row (if any)
            stop_idx = None
            for i, row in following_rows.iterrows():
                if re.search(all_tags_pattern, row['tag']):
                    stop_idx = i
                    break

            # Define working sub-block
            if stop_idx:
                working_block = block_df.loc[start_idx:stop_idx].iloc[:-1]  # exclude tagged row
            else:
                working_block = block_df.loc[start_idx:]
                # === CHAIN LOGIC ===
            if working_block.shape[0] == 1:
                print('singleton')
                block_df.at[working_block.index[0], 'tag'] += '_singleton'
                continue
            elif working_block.shape[0] == 2:
                print('duo')
                updated_df, _ = two_frags_handler(working_block, common_length)
                block_df.loc[updated_df.index, 'tag'] = updated_df['tag']
                continue
            else:
                chain = [0]
                overlap = False
                for i in range(1, len(working_block)):
                    prev = working_block.iloc[chain[-1]]
                    curr = working_block.iloc[i]
                    print(f'{prev.repStart}\t{prev.repEnd}::{curr.repStart}\t{curr.repEnd}')
                    
                    if (prev['repEnd'] < curr['repStart']) and \
                    ((prev['repLeft'] >= curr['repEnd'] - curr['repStart']) or curr['repLeft'] == 0):
                        chain.append(i)
                        print('chain')
                        if curr['repLeft'] == 0:
                            break
                    elif (prev['repEnd'] - curr['repStart'] < 0.1 * common_length) and (prev['repStart'] < curr['repStart']):
                        print('chain')
                        overlap=True
                        chain.append(i)
                        if curr['repLeft'] == 0:
                            break
                    else:
                        break

                if len(chain) >= 2:
                    # Valid chain: tag
                    col_idx = block_df.columns.get_loc('tag')
                    real_indices = working_block.iloc[chain].index
                    if overlap:
                        block_df.loc[working_block.iloc[chain].index, 'tag'] += '_overlap'
                    block_df.iloc[block_df.index.get_loc(real_indices[0]), col_idx] += '_open'
                    block_df.iloc[block_df.index.get_loc(real_indices[-1]), col_idx] += '_close'
                    for mid in real_indices[1:-1]:
                        block_df.iloc[block_df.index.get_loc(mid), col_idx] += '_mid'
                else:
                    # Fallback: tag all in working block with '_test'
                    #for idx in working_block.index:
                    block_df.at[working_block.index[0], 'tag'] += '_close'
        more_than_two_frag.append(block_df)
            
    else:
        #check if any row in the block_df, at 'tag' column contains '_complete', if yes, fill another row tag with row['tag']+'_close' 
        updated_df, category = two_frags_handler(block_df, common_length)
        if category == 'one_frag_with_complete':
            one_frag_with_complete.append(updated_df)
        elif category == 'two_frags_simple_join':
            two_frags_simple_join.append(updated_df)
        elif category == 'two_frags_overlap_join':
            two_frags_overlap_join.append(updated_df)
        else:
            two_frags_non_join.append(updated_df)
#%%
#merge all subtable
fragment_tables = [
    all_complete_multis_table,
    *more_than_two_frag,
    *two_frags_simple_join,
    *two_frags_overlap_join,
    *two_frags_non_join,
    *one_frag_with_complete,
    singleton_table,
]
merged_table = pd.concat(fragment_tables)
# %%
def extract_block_id_from_tag(tag: str) -> str:
    pattern = r"^(.*?)(?:_(complete|open|mid|close|singleton|test|overlap_open|overlap_mid|overlap_close))$"
    match = re.match(pattern, tag)
    if match:
        return match.group(1)
    return tag
#%%
import re
from collections import defaultdict
from typing import List

def convert_tags_to_internal_ids(tags: List[str], subfamily: str) -> List[str]:
    result = []
    individual_counter = defaultdict(int)
    join_group_id = 1
    join_stack = []
    current_suffix_type = None  # either 'join' or 'overlap_join'

    def flush_join_stack():
        nonlocal join_stack, join_group_id
        if join_stack:
            for t, b, suffix in join_stack:
                result.append(f"{subfamily}_{b}_{suffix}_{join_group_id}")
            join_group_id += 1
            join_stack = []

    for tag in tags:
        base = extract_block_id_from_tag(tag)

        # Normal join group
        if tag.endswith(('_open', '_mid', '_close')):
            suffix = 'join'

            if tag.endswith('_open'):
                # If a previous group was left open, flush it before starting a new one
                flush_join_stack()
                current_suffix_type = suffix
                join_stack.append((tag, base, suffix))
            elif tag.endswith('_mid'):
                if current_suffix_type == suffix:
                    join_stack.append((tag, base, suffix))
                else:
                    # mid without open — skip or handle as orphan
                    individual_counter[(base, 'mid')] += 1
                    count = individual_counter[(base, 'mid')]
                    result.append(f"{subfamily}_{base}_mid_{count}")
            elif tag.endswith('_close'):
                if current_suffix_type == suffix:
                    join_stack.append((tag, base, suffix))
                    flush_join_stack()
                    current_suffix_type = None
                else:
                    # Orphan close
                    individual_counter[(base, 'close')] += 1
                    count = individual_counter[(base, 'close')]
                    result.append(f"{subfamily}_{base}_singleton2_{count}")

        # Overlap join group
        elif tag.endswith(('_overlap_open', '_overlap_mid', '_overlap_close')):
            suffix = 'overlap_join'

            if tag.endswith('_overlap_open'):
                flush_join_stack()
                current_suffix_type = suffix
                join_stack.append((tag, base, suffix))
            elif tag.endswith('_overlap_mid'):
                if current_suffix_type == suffix:
                    join_stack.append((tag, base, suffix))
                else:
                    individual_counter[(base, 'overlap_mid')] += 1
                    count = individual_counter[(base, 'overlap_mid')]
                    result.append(f"{subfamily}_{base}_overlap_mid_{count}")
            elif tag.endswith('_overlap_close'):
                if current_suffix_type == suffix:
                    join_stack.append((tag, base, suffix))
                    flush_join_stack()
                    current_suffix_type = None
                else:
                    individual_counter[(base, 'overlap_close')] += 1
                    count = individual_counter[(base, 'overlap_close')]
                    result.append(f"{subfamily}_{base}_overlap_close_{count}")

        # Independent tags
        else:
            suffix = tag.split('_')[-1]
            individual_counter[(base, suffix)] += 1
            count = individual_counter[(base, suffix)]
            result.append(f"{subfamily}_{base}_{suffix}_{count}")

    # Flush any unfinished join
    flush_join_stack()

    return result
#%%
#make internal_id
# Group tags by their extracted block ID
merged_table['block_id'] = merged_table['tag'].map(extract_block_id_from_tag)
# Apply convert_tags_to_internal_ids to each block
def process_block(df_block):
    tags = df_block['tag'].tolist()
    internal_ids = convert_tags_to_internal_ids(tags, subfamily)
    return pd.Series(internal_ids, index=df_block.index)

merged_table['internal_id'] = (
    merged_table
    .groupby('block_id', group_keys=False, sort=False)
    [['tag']] # explicitly select relevant column
    .apply(process_block)
)

# Save result
internal_id_table = merged_table[['tag', 'internal_id']]
internal_id_table = internal_id_table.rename_axis('rmsk_index').reset_index()
#%%
internal_id_table.to_csv(f'{output_dir}/{subfamily}.internal_id.txt', sep='\t', index=False)
#%%
# extra
# %%
block_id_list=repeatmasker_table[repeatmasker_table['repName'] == subfamily]['id'].unique()
base_subfamily = re.match(r'^[A-Za-z0-9]+', subfamily).group()[:-1]
tagged_blocks = []
for block_id in block_id_list:
    print(block_id)
    block_df = repeatmasker_table[repeatmasker_table['id'] == block_id].copy()
    block_df['tag'] = None

    rep_names = list(block_df['repName'])

    # LTRs: same base family, not internal
    ltr_positions = [i for i, name in enumerate(rep_names)
                    if re.match(f"^{base_subfamily}[A-Z]$", name)]

    # Positions with matching -int
    
    int_positions = [i for i, name in enumerate(rep_names)
                    if re.match(f"^{base_subfamily}[A-Z]-int$", name)]
    if int_positions:
        print('pass1')
        for i in ltr_positions:
            row_label = block_df.index[i]
            if any(pos > i for pos in int_positions):
                block_df.at[row_label, 'tag'] = f'{block_id}_bINT'
            elif any(pos < i for pos in int_positions):
                block_df.at[row_label, 'tag'] = f'{block_id}_aINT'
            else:
                block_df.at[row_label, 'tag'] = f'{block_id}_nINT'
    else:
        row_label=block_df.index[0]
        block_df.loc[row_label:, 'tag'] = f'{block_id}_nINT'
    block_df = block_df[block_df['repName']==subfamily]
    tagged_blocks.append(block_df)

subfam_table = pd.concat(tagged_blocks)   
# %%
block_id_list=repeatmasker_table[repeatmasker_table['repName'] == subfamily]['id'].unique()
base_subfamily = get_base_subfamily(subfamily)
block_id = 3785581
block_df = repeatmasker_table[repeatmasker_table['id'] == block_id].copy()
#%%
block_df['tag'] = None

rep_names = list(block_df['repName'])

# LTRs: same base family, not internal
ltr_positions = [i for i, name in enumerate(rep_names)
                 if re.match(f"^{base_subfamily}[A-Z]_?$", name)]
#%%
# Positions with matching -int

int_positions = [i for i, name in enumerate(rep_names)
                if re.match(f"^{base_subfamily}[A-Z]-int$", name)]
if int_positions != []:
    for i in ltr_positions:
        row_label = block_df.index[i]
        print('pass1')
        if any(pos > i for pos in int_positions):
            block_df.at[row_label, 'tag'] = f'{block_id}_bINT'
        elif any(pos < i for pos in int_positions):
            block_df.at[row_label, 'tag'] = f'{block_id}_aINT'
        else:
            block_df.at[row_label, 'tag'] = f'{block_id}_nINT'
else:
    print('pass2')
    row_label=block_df.index[0]
    block_df.loc[:, 'tag'] = f'{block_id}_nINT'
block_df = block_df[block_df['repName']==subfamily]
tagged_blocks.append(block_df)
# %%
