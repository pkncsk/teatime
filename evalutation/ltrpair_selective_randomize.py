import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load all pair-level te_age table (must include: internal_id, block_id, element_family, te_age)
full_table = pd.read_csv(
    "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/final_LTR_sample_3000.teatime.txt",
    sep="\t"
)

# Extract block_id if not present
if 'block_id' not in full_table.columns:
    full_table['block_id'] = full_table['internal_id'].apply(lambda x: "_".join(x.split("_")[:2]))

# Extract element_family if not present
if 'element_family' not in full_table.columns:
    full_table['element_family'] = full_table['internal_id'].apply(lambda x: x.split("_")[0])

# Filter: only keep block_ids where both aINT and bINT have same te_age
equal_age_blocks = (
    full_table.groupby('block_id')
    .filter(lambda g: g['te_age'].nunique() == 1 and len(g) == 2)  # valid pair with equal age
)

# Extract unique one-row-per-pair with representative te_age
pair_ages = (
    equal_age_blocks
    .groupby(['block_id', 'element_family'])['te_age']
    .first()
    .reset_index()
)

# Plot: Histogram per subfamily with mean & median
subfamilies = pair_ages['element_family'].unique()

for sf in sorted(subfamilies):
    sub_df = pair_ages[pair_ages['element_family'] == sf]
    if len(sub_df) < 5:
        continue  # skip tiny families

    plt.figure(figsize=(6, 4))
    plt.hist(sub_df['te_age'], bins=30, edgecolor='black', color='skyblue')
    mean_val = sub_df['te_age'].mean()
    median_val = sub_df['te_age'].median()

    # Add mean/median lines
    plt.axvline(mean_val, color='red', linestyle='--', label=f"Mean: {mean_val:.1f}")
    plt.axvline(median_val, color='green', linestyle=':', label=f"Median: {median_val:.1f}")

    plt.title(f"{sf} — Age Distribution (n={len(sub_df)})")
    plt.xlabel("tAge (MYA)")
    plt.ylabel("Number of Pairs")
    plt.legend()
    plt.tight_layout()
    plt.show()
