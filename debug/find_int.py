#%%
import pandas as pd
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
subfamily  = 'THE1C'
from ma_mapper import utility
repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
#%%
repeatmasker_table
# %%
subfamily_table=repeatmasker_table[repeatmasker_table['repName']==subfamily]
#%%
subfamily_block_id=subfamily_table['id'].unique()
#%%
multiple_te_blocks = []
for id in subfamily_block_id:
    block_id_table=repeatmasker_table[repeatmasker_table['id']==id]
    if block_id_table.shape[0] > 1:
        multiple_te_blocks.append(block_id_table)
# %%
repeatmasker_table[repeatmasker_table['id']==6363]
# %%
len(multiple_te_blocks)
# %%
multiple_te_blocks[2]
# %%
repeatmasker_table[repeatmasker_table.index.isin([6085,6086,6087,6088])]
# %%
# recover block with multiple TE first
# check pattern LTR, int, LTR
# probe +4 -4 rows, check LTR, int, LTR

# recover block with singular TE
# probe +4 -4 rows, check LTR, int, LTR
#%%
import pandas as pd

# Load RepeatMasker table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility

repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)

# Select LTR subfamily
subfamily = 'THE1C'
subfamily_table = repeatmasker_table[repeatmasker_table['repName'] == subfamily]

# Extract row indices
subfamily_indices = subfamily_table.index.tolist()

# Initialize categorization and tables
categorized_blocks = {
    'complete': [],
    'partial': [],
    'standalone': []
}

categorized_tables = {
    'complete': pd.DataFrame(),
    'partial': pd.DataFrame(),
    'standalone': pd.DataFrame()
}
#%%import pandas as pd

# Load RepeatMasker table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility

repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)

# Select LTR subfamily
subfamily = 'THE1C'
subfamily_table = repeatmasker_table[repeatmasker_table['repName'] == subfamily]

# Extract row indices
subfamily_indices = subfamily_table.index.tolist()

# Initialize categorization and tables
categorized_blocks = {
    'complete': [],
    'partial': [],
    'standalone': []
}

categorized_tables = {
    'complete': pd.DataFrame(),
    'partial': pd.DataFrame(),
    'standalone': pd.DataFrame()
}

# Check for patterns in a ±10 row window
for idx in subfamily_indices:
    # Define the ±10 row window
    start_idx = max(0, idx - 10)
    end_idx = min(len(repeatmasker_table) - 1, idx + 10)
    window = repeatmasker_table.loc[start_idx:end_idx]
    
    # Check strand consistency
    target_strand = repeatmasker_table.loc[idx, 'strand']
    window = window[window['strand'] == target_strand]

    # Get unique repNames in the window
    rep_names = list(window['repName'])

    # Locate occurrences of subfamily elements
    ltr_positions = [i for i, name in enumerate(rep_names) if name == subfamily]
    int_positions = [i for i, name in enumerate(rep_names) if name == f"{subfamily}-int"]

    # Check for complete pattern: 'xxx' - 'xxx-int' - 'xxx'
    if len(ltr_positions) >= 2 and len(int_positions) >= 1:
        first_ltr = ltr_positions[0]
        last_ltr = ltr_positions[-1]
        has_int_between = any(first_ltr < pos < last_ltr for pos in int_positions)

        if has_int_between:
            categorized_blocks['complete'].append(idx)
            extracted_table = window.iloc[first_ltr:last_ltr + 1]  # Include everything between first and last LTR
            categorized_tables['complete'] = pd.concat([categorized_tables['complete'], extracted_table])
            continue  # No need to check further, already categorized

    # Check for partial cases and extract all rows in between
    if len(ltr_positions) >= 1 and len(int_positions) >= 1:
        first = min(ltr_positions + int_positions)
        last = max(ltr_positions + int_positions)
        categorized_blocks['partial'].append(idx)
        extracted_table = window.iloc[first:last + 1]  # Include everything in between
        categorized_tables['partial'] = pd.concat([categorized_tables['partial'], extracted_table])
        continue  # No need to check standalone

    elif len(ltr_positions) > 1:  # Case: 'xxx' - 'xxx'
        first = ltr_positions[0]
        last = ltr_positions[-1]
        categorized_blocks['partial'].append(idx)
        extracted_table = window.iloc[first:last + 1]  # Include everything in between
        categorized_tables['partial'] = pd.concat([categorized_tables['partial'], extracted_table])
        continue  # No need to check standalone

    # If none of the above, classify as standalone
    categorized_blocks['standalone'].append(idx)
    categorized_tables['standalone'] = pd.concat([categorized_tables['standalone'], window[window['repName'] == subfamily]])

# Remove duplicates and reset index
for category in categorized_tables:
    categorized_tables[category] = categorized_tables[category].drop_duplicates().reset_index(drop=True)

# Print summary
for category, indices in categorized_blocks.items():
    print(f"{category}: {len(indices)} entries")

# Extracted tables
complete_table = categorized_tables['complete']
one_sided_table = categorized_tables['partial']
standalone_table = categorized_tables['standalone']

# Save if needed
# complete_table.to_csv("complete_ltr_blocks.csv", index=False)
# one_sided_table.to_csv("one_sided_ltr_blocks.csv", index=False)
# standalone_table.to_csv("standalone_ltr_blocks.csv", index=False)

#%%
for category in categorized_tables:
    categorized_tables[category] = categorized_tables[category].drop_duplicates().reset_index(drop=True)
#%%
for category, indices in categorized_blocks.items():
    print(f"{category}: {len(indices)} entries")

complete_table = categorized_tables['complete']
one_sided_table = categorized_tables['partial']
standalone_table = categorized_tables['standalone']

# complete_table.to_csv("complete_ltr_blocks.csv", index=False)
# one_sided_table.to_csv("one_sided_ltr_blocks.csv", index=False)
# standalone_table.to_csv("standalone_ltr_blocks.csv", index=False)
#%% debug
idx=194
start_idx = max(0, idx - 10)
end_idx = min(len(repeatmasker_table) - 1, idx + 10)
window = repeatmasker_table.loc[start_idx:end_idx]

# Check strand consistency
target_strand = repeatmasker_table.loc[idx, 'strand']
window = window[window['strand'] == target_strand]

# Get unique repNames in the window
rep_names = set(window['repName'])

# Identify `xxx` and `xxx-int`
has_ltr = subfamily in rep_names
has_int = f"{subfamily}-int" in rep_names
# %%
import pandas as pd

# Load RepeatMasker table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility

repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)

# Select LTR subfamily
subfamily = 'THE1C'
subfamily_table = repeatmasker_table[repeatmasker_table['repName'] == subfamily]

# Extract row indices
subfamily_indices = subfamily_table.index.tolist()

# Store blocks separately instead of merging
categorized_blocks = {
    'complete': [],
    'one-sided': [],
    'standalone': []
}

# Track used indices to prevent duplication
used_indices = set()

# Check for patterns in a ±10 row window
for idx in subfamily_indices:
    if idx in used_indices:  # Skip if already used
        continue

    # Define the ±10 row window
    start_idx = max(0, idx - 20)
    end_idx = min(len(repeatmasker_table) - 1, idx + 20)
    window = repeatmasker_table.loc[start_idx:end_idx]

    # Check strand consistency
    target_strand = repeatmasker_table.loc[idx, 'strand']
    window = window[window['strand'] == target_strand]

    # Get unique repNames in the window
    rep_names = list(window['repName'])

    # Locate occurrences of subfamily elements
    ltr_positions = [i for i, name in enumerate(rep_names) if name == subfamily]
    int_positions = [i for i, name in enumerate(rep_names) if name == f"{subfamily}-int"]

    # Check for complete pattern: 'xxx' - 'xxx-int' - 'xxx'
    if len(ltr_positions) >= 2 and len(int_positions) >= 1:
        first_ltr = ltr_positions[0]
        last_ltr = ltr_positions[-1]
        has_int_between = any(first_ltr < pos < last_ltr for pos in int_positions)

        if has_int_between:
            extracted_table = window.iloc[first_ltr:last_ltr + 1]  # Include everything between first and last LTR
            categorized_blocks['complete'].append(extracted_table)

            # Mark these indices as used
            used_indices.update(extracted_table.index)
            continue  # No need to check further, already categorized

    # Check for partial cases and extract all rows in between
    if len(ltr_positions) >= 1 and len(int_positions) >= 1:
        first = min(ltr_positions + int_positions)
        last = max(ltr_positions + int_positions)
        extracted_table = window.iloc[first:last + 1]  # Include everything in between
        categorized_blocks['one-sided'].append(extracted_table)

        # Mark these indices as used
        used_indices.update(extracted_table.index)
        continue  # No need to check standalone

    elif len(ltr_positions) > 1:  # Case: 'xxx' - 'xxx'
        first = ltr_positions[0]
        last = ltr_positions[-1]
        extracted_table = window.iloc[first:last + 1]  # Include everything in between
        categorized_blocks['one-sided'].append(extracted_table)

        # Mark these indices as used
        used_indices.update(extracted_table.index)
        continue  # No need to check standalone

    # If none of the above, classify as standalone
    extracted_table = window[window['repName'] == subfamily]

    # Make sure rows aren’t already used
    extracted_table = extracted_table[~extracted_table.index.isin(used_indices)]
    if not extracted_table.empty:
        categorized_blocks['standalone'].append(extracted_table)
        used_indices.update(extracted_table.index)

# Print summary
for category, blocks in categorized_blocks.items():
    print(f"{category}: {len(blocks)} blocks")

# Each block is a separate DataFrame in the list
complete_blocks = categorized_blocks['complete']
one_sided_blocks = categorized_blocks['one-sided']
standalone_blocks = categorized_blocks['standalone']

# %%
global_block_count = 1  # Global running counter

# Add block_id column to each block in categorized blocks
for category, blocks in categorized_blocks.items():
    for i, block in enumerate(blocks):
        block_id = f"{subfamily}_{category}_{global_block_count}"
        block['block_id'] = block_id  # Add block_id column
        global_block_count += 1  # Increment counter

# After adding the block IDs, you can concatenate the blocks as before
all_blocks = pd.concat([block for blocks in categorized_blocks.values() for block in blocks])

# Now `all_blocks` will contain a `block_id` column for each individual block

# %%
internal_id_table = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/{subfamily}.internal_id.txt'
age_table = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/{subfamily}.txt'
# %%
internal_id_df = pd.read_csv(internal_id_table, sep='\t')
age_df = pd.read_csv(age_table, sep='\t')
# %%
age_internal_id_df = pd.merge(age_df, internal_id_df, on='internal_id')
#%%
internal_id_table_int = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/{subfamily}-int.internal_id.txt'
age_table_int = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/{subfamily}-int.txt'
# %%
internal_id_df_int = pd.read_csv(internal_id_table_int, sep='\t')
age_df_int = pd.read_csv(age_table_int, sep='\t')
# %%
age_internal_id_df_int = pd.merge(age_df_int, internal_id_df_int, on='internal_id')
#%%
age_internal_id = pd.concat([age_internal_id_df, age_internal_id_df_int])
# %%
all_blocks_age = pd.merge(all_blocks, age_internal_id,left_index=True, right_on='rmsk_index', how='left')
# %%
