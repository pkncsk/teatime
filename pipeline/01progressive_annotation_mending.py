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

def tag_singletons(df):
    tag_counts = df['tag'].value_counts()

    df['tag'] = df.apply(
        lambda row: f"{row['tag']}_singleton" if tag_counts[row['tag']] == 1 else row['tag'],
        axis=1
    )
    return df


def two_frags_handler(block_df, common_length, overlap_threshold):
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

    elif (first['repEnd'] - second['repStart'] < overlap_threshold * common_length) and (
        first['repStart'] < second['repStart']
    ):
        block_df.at[first_label, 'tag'] += '_overlap_open'
        block_df.at[second_label, 'tag'] += '_overlap_close'
        return block_df, 'two_frags_overlap_join'

    else:
        block_df.at[first_label, 'tag'] += '_close'
        block_df.at[second_label, 'tag'] += '_close'
        return block_df, 'two_frags_non_join'

def twoplus_frags_handler(block_df, common_length, overlap_threshold):
    all_tags_pattern = '_open|_mid|_close|_singleton|_test|_complete' 
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
            updated_df, _ = two_frags_handler(working_block, common_length, overlap_threshold)
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
                elif (prev['repEnd'] - curr['repStart'] < overlap_threshold * common_length) and (prev['repStart'] < curr['repStart']):
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
    return block_df

def extract_block_id_from_tag(tag: str) -> str:
    pattern = r"^(.*?)(?:_(complete|open|mid|close|singleton|test|overlap_open|overlap_mid|overlap_close))$"
    match = re.match(pattern, tag)
    if match:
        return match.group(1)
    return tag

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

def process_block(df_block, subfamily):
    tags = df_block['tag'].tolist()
    internal_ids = convert_tags_to_internal_ids(tags, subfamily)
    return pd.Series(internal_ids, index=df_block.index)
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
        #separate singletons from multis
        subfam_table = tag_singletons(subfam_table)
        singleton_table=subfam_table[subfam_table['tag'].str.contains('singleton', na=False)]
        multis_table = subfam_table[~subfam_table['tag'].str.contains('singleton', na=False)]

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
        one_frag_with_complete = []
        two_frags_simple_join = []
        more_than_two_frag = []
        two_frags_overlap_join = []
        two_frags_non_join = []
        
        for block_df in partial_complete_df:
            print(block_df['id'].unique()[0])
            if block_df['strand'].unique()[0] == '-':
                block_df = block_df.sort_index(ascending=False)
            row_number = block_df.shape[0]
            if row_number >2:
                update_df = twoplus_frags_handler(block_df, common_length, overlap_threshold)
                more_than_two_frag.append(update_df)
                    
            else:
                #check if any row in the block_df, at 'tag' column contains '_complete', if yes, fill another row tag with row['tag']+'_close' 
                updated_df, category = two_frags_handler(block_df, common_length, overlap_threshold)
                if category == 'one_frag_with_complete':
                    one_frag_with_complete.append(updated_df)
                elif category == 'two_frags_simple_join':
                    two_frags_simple_join.append(updated_df)
                elif category == 'two_frags_overlap_join':
                    two_frags_overlap_join.append(updated_df)
                else:
                    two_frags_non_join.append(updated_df)
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
        #make internal_id
        # Group tags by their extracted block ID
        merged_table['block_id'] = merged_table['tag'].map(extract_block_id_from_tag)
        # Apply convert_tags_to_internal_ids to each block
        merged_table['internal_id'] = (
            merged_table
            .groupby('block_id', group_keys=False, sort=False)
            [['tag']] # explicitly select relevant column
            .apply(lambda df_block: process_block(df_block, subfamily))
        )
        # Save result
        internal_id_table = merged_table[['tag', 'internal_id']]
        internal_id_table = internal_id_table.rename_axis('rmsk_index').reset_index()
        internal_id_table.to_csv(f'{output_dir}/{subfamily}.internal_id.txt', sep='\t', index=False)
    else:
        print(f'file already exists')
# %%
def main(repeatmasker_filepath, subfamily_list, output_dir, overlap_threshold):
    global repeatmasker_table
    repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
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
        repeatmasker_filepath=args.repeatmasker_filepath,
        subfamily_list=subfamily_list,
        output_dir=args.output_dir,
        overlap_threshold=args.overlap_threshold,
    )
