# %%
import pandas as pd
#%%
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
from ma_mapper import utility
repeatmasker_table = utility.repeatmasker_prep(repeatmasker_filepath)
subfamily = 'THE1C'
subfamily_table = repeatmasker_table[repeatmasker_table['repName'] == subfamily]
subfamily_indices = subfamily_table.index.tolist()
#%%
block_te_counts = (
    repeatmasker_table
    .groupby(['id', 'repName', 'repClass'])
    .size()
    .reset_index(name='te_count')
)
#%%
#Define the summarization function
def summarize_subfamily(sub_df):
    total_blocks = len(sub_df)
    solo_blocks = (sub_df['te_count'] == 1).sum()
    multi_blocks = (sub_df['te_count'] > 1).sum()
    
    total_te = sub_df['te_count'].sum()
    solo_te = sub_df.loc[sub_df['te_count'] == 1, 'te_count'].sum()
    multi_te = sub_df.loc[sub_df['te_count'] > 1, 'te_count'].sum()
    
    return pd.Series({
        'total_blocks': total_blocks,
        'solo_blocks': solo_blocks,
        'multi_blocks': multi_blocks,
        'total_te': total_te,
        'solo_te': solo_te,
        'multi_te': multi_te,
        'solo_block_pct': round(100 * solo_blocks / total_blocks, 2) if total_blocks > 0 else 0,
        'multi_block_pct': round(100 * multi_blocks / total_blocks, 2) if total_blocks > 0 else 0,
        'solo_te_pct': round(100 * solo_te / total_te, 2) if total_te > 0 else 0,
        'multi_te_pct': round(100 * multi_te / total_te, 2) if total_te > 0 else 0,
    })

# Group by repName (subfamily) and summarize with include_groups=False
# Group and apply for repName summary
summary_by_subfamily = (
    block_te_counts
    .groupby('repName')
    .apply(summarize_subfamily, include_groups=False)
    .reset_index()
)

# Group and apply for repClass summary
summary_by_class = (
    block_te_counts
    .groupby('repClass')
    .apply(summarize_subfamily, include_groups=False)
    .reset_index()
)

# Optional: sort by most fragmented (repName or repClass)
summary_by_subfamily = summary_by_subfamily.sort_values(by='multi_te_pct', ascending=False)
summary_by_class = summary_by_class.sort_values(by='multi_te_pct', ascending=False)

# Preview
print("Summary by Subfamily:")
print(summary_by_subfamily.head())
print("\nSummary by RepClass:")
print(summary_by_class.head())
# %%
summary_by_subfamily.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/annotation_mending_revise/multientry_by_subfamily.txt', sep='\t')
# %%
summary_by_class.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/annotation_mending_revise/multientry_by_class.txt', sep='\t')
# %%
