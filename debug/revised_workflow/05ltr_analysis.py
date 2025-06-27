#%%
import pandas as pd
#%%
from ma_mapper import utility
global repeatmasker_table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
repname_counts = repeatmasker_table['repName'].value_counts().reset_index()
repname_counts.columns = ['repName', 'count']
subfam_target=repname_counts[repname_counts['count']<1000]['repName'].unique()
subfamily = subfam_target
subfamily='THE1C'
the1c_patched = pd.read_csv(f'/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/revised_workflow/{subfamily}.ltr_fix.txt', sep='\t', index_col =0)
#%%
df=the1c_patched
#%%
import pandas as pd

# Assuming your DataFrame is `df`

grouped = df.groupby('block_id')

# Prepare dict to hold block_id and row count for each category
summary = {
    'Single entry': [],
    'Two entries': [],
    'Two entries + match': [],
    'Two entries + mismatch': [],
    'More than two entries': [],
    'More than two entries + match': [],
    'More than two entries + mismatch': [],
    'Total match': [],
    'Total mismatch': [],
}

for block_id, group in grouped:
    row_count = len(group)
    ages = group['te_age'].unique()
    age_match = len(ages) == 1

    # Categorize by number of rows
    if row_count == 1:
        summary['Single entry'].append((block_id, row_count))
    elif row_count == 2:
        summary['Two entries'].append((block_id, row_count))
        if age_match:
            summary['Two entries + match'].append((block_id, row_count))
        else:
            summary['Two entries + mismatch'].append((block_id, row_count))
    else:
        summary['More than two entries'].append((block_id, row_count))
        if age_match:
            summary['More than two entries + match'].append((block_id, row_count))
        else:
            summary['More than two entries + mismatch'].append((block_id, row_count))

    # Track total matches/mismatches regardless of row count
    if age_match:
        summary['Total match'].append((block_id, row_count))
    else:
        summary['Total mismatch'].append((block_id, row_count))

# Grand totals for blocks and rows
total_blocks = len(grouped)
total_rows = len(df)

# Build summary table rows
summary_rows = []
for category, entries in summary.items():
    block_count = len(entries)
    row_count = sum(rc for _, rc in entries)
    summary_rows.append([
        category,
        block_count,
        f"{100 * block_count / total_blocks:.2f}%",
        row_count,
        f"{100 * row_count / total_rows:.2f}%"
    ])

# Append grand total row
summary_rows.append([
    "Grand total",
    total_blocks,
    "100.00%",
    total_rows,
    "100.00%"
])

# Create DataFrame for better formatting
summary_df = pd.DataFrame(summary_rows, columns=[
    'Category', 'Block count', 'Block %', 'Row count', 'Row %'
])

print(summary_df)
#%%
import pandas as pd

# Assuming your DataFrame is `df`

grouped = df.groupby('block_id')

# Prepare dict to hold block_id and row count for each category
summary = {
    'Single entry': [],
    'Two entries': [],
    'Two entries + match': [],
    'Two entries + mismatch': [],
    'More than two entries': [],
    'More than two entries + match': [],
    'More than two entries + mismatch': [],
    'Total match': [],
    'Total mismatch': [],
}

for block_id, group in grouped:
    row_count = len(group)
    ages = group['corrected_te_age'].unique()
    age_match = len(ages) == 1

    # Categorize by number of rows
    if row_count == 1:
        summary['Single entry'].append((block_id, row_count))
    elif row_count == 2:
        summary['Two entries'].append((block_id, row_count))
        if age_match:
            summary['Two entries + match'].append((block_id, row_count))
        else:
            summary['Two entries + mismatch'].append((block_id, row_count))
    else:
        summary['More than two entries'].append((block_id, row_count))
        if age_match:
            summary['More than two entries + match'].append((block_id, row_count))
        else:
            summary['More than two entries + mismatch'].append((block_id, row_count))

    # Track total matches/mismatches regardless of row count
    if age_match:
        summary['Total match'].append((block_id, row_count))
    else:
        summary['Total mismatch'].append((block_id, row_count))

# Grand totals for blocks and rows
total_blocks = len(grouped)
total_rows = len(df)

# Build summary table rows
summary_rows = []
for category, entries in summary.items():
    block_count = len(entries)
    row_count = sum(rc for _, rc in entries)
    summary_rows.append([
        category,
        block_count,
        f"{100 * block_count / total_blocks:.2f}%",
        row_count,
        f"{100 * row_count / total_rows:.2f}%"
    ])

# Append grand total row
summary_rows.append([
    "Grand total",
    total_blocks,
    "100.00%",
    total_rows,
    "100.00%"
])

# Create DataFrame for better formatting
summary_df = pd.DataFrame(summary_rows, columns=[
    'Category', 'Block count', 'Block %', 'Row count', 'Row %'
])

print(summary_df)
# %%
block_with_complete=the1c_patched[the1c_patched['tag'].str.contains('complete')]['block_id'].unique()
# %%
the1c_patched[the1c_patched['block_id'].isin(block_with_complete)]
# %%
# Step 1: Group by block_id
df=the1c_patched
grouped = df.groupby('block_id')

# Step 2: Filter groups with exactly 2 rows
two_row_groups = grouped.filter(lambda x: len(x) == 2)

# Step 3: Further filter those with both aINT and bINT
desired_blocks = two_row_groups.groupby('block_id').filter(
    lambda x: set(x['int_type']) == {'aINT', 'bINT'}
)

# Optional: get unique block_ids
result_block_ids = desired_blocks['block_id'].unique()

# Display results
print(result_block_ids)

# %%
the1c_patched[the1c_patched['block_id'].isin(result_block_ids)]
#%%

# Load RepeatMasker table
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility

repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
#%%
repeatmasker_table[repeatmasker_table['id']==10428]
#%%
#chr1:7790171-7792441