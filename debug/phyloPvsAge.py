#%%
import pandas as pd
import os
from ma_mapper import utility
import numpy as np
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
#%% INPUT PARAMETERS
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
mean_phyloP_filepath  = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/phyloP.txt'
combined_table_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime/combined_age_div_lenient/'
#%% INITIATION
age_canon = utility.age_canon()[:-1]
mean_phyloP = pd.read_csv(mean_phyloP_filepath, sep='\t', index_col = 0).reset_index()
if combined_table_dir is None:
    combined_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
combined_te_age_filepath = f'{combined_table_dir}/all_subfamilies.itd.txt'
combined_te_age_df = pd.read_csv(combined_te_age_filepath, sep='\t')
#%%
#%% Extract subfamily name from internal_id
def extract_subfamily(internal_id):
    parts = internal_id.split('_')
    if len(parts) >= 3:
        return '_'.join(parts[:-2])
    return internal_id  # fallback if structure is unexpected

combined_te_age_df['subfamily'] = combined_te_age_df['internal_id'].apply(extract_subfamily)

#%% Group by subfamily and compute median te_age
median_age_per_subfamily = (
    combined_te_age_df
    .groupby('subfamily')['te_age']
    .median()
    .reset_index()
    .rename(columns={'te_age': 'median_te_age'})
)

# %%

#%% Map median age to nearest canonical age
median_age_per_subfamily['nearest_canon_age'] = median_age_per_subfamily['median_te_age'].apply(
    lambda x: find_nearest(age_canon, x)
)

#%% Output result
print(median_age_per_subfamily.head())
# %%
medage_meanphy=median_age_per_subfamily.merge(mean_phyloP, on='subfamily')
# %%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10,4))
flierprops = dict(marker='.', markersize=3, linestyle='none', markeredgecolor='black', alpha=0.6)
boxprops = dict(facecolor = 'white',alpha=0.5)
grouped_data = [medage_meanphy[medage_meanphy['nearest_canon_age'] == age]['mean_phyloP'] for age in age_canon]

for i, d in enumerate(grouped_data):
        y = d
        x = np.random.normal(age_canon[i], 0.4, len(y))  # Jitter added to the x-values
        ax.plot(x, y, 'bo', alpha=0.05, markersize=1)
ax.boxplot(grouped_data, positions=age_canon, widths=2, flierprops=flierprops,patch_artist=True, boxprops=boxprops)

subset=medage_meanphy
from scipy import stats
# First regression line (0 to 43.2)
res = stats.linregress(subset['nearest_canon_age'], subset['mean_phyloP'])
ax.plot(subset['nearest_canon_age'], res.intercept + res.slope * subset['nearest_canon_age'], color='red', label='Regression Line (0-96)', linewidth=1, alpha=0.5)

set_size=medage_meanphy[(medage_meanphy['nearest_canon_age'] > 0)].shape[0]
ax.text(0.99, 0.01, f'TE subfamilies n={set_size}', 
        transform=ax.transAxes, 
        fontsize=12, 
        verticalalignment='bottom', 
        horizontalalignment='right')
ax.axhline(y=0, color='grey', linewidth=1, alpha=0.5)
ax.axhline(y=3.342447251181431, color='blue',linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlim(0,110)
ax.set_xlabel('median TE age (MYA) estimated by TEA-TIME approach', fontsize=12)
ax.set_ylabel('mean phyloP', fontsize=12)
ax.set_title('Mean phyloP grouped by TEA-TIME', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.show()
# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import t

# Your data: assumes medage_meanphy with columns ['nearest_canon_age', 'mean_phyloP']
# and age_canon is a list/array of canonical ages

fig, ax = plt.subplots(figsize=(10, 4))
flierprops = dict(marker='.', markersize=3, linestyle='none', markeredgecolor='black', alpha=0.6)
boxprops = dict(facecolor='white', alpha=0.5)

# Group data by canonical age
grouped_data = [medage_meanphy[medage_meanphy['nearest_canon_age'] == age]['mean_phyloP'] for age in age_canon]


# Add boxplot
ax.boxplot(grouped_data, positions=age_canon, widths=2, flierprops=flierprops,
           patch_artist=True, boxprops=boxprops)
subset = medage_meanphy[(medage_meanphy['nearest_canon_age'] > 0) & (medage_meanphy['nearest_canon_age'] <= 96.0)]

# Subset for regression
#subset = medage_meanphy

# Linear regression
res = stats.linregress(subset['nearest_canon_age'], subset['mean_phyloP'])

# Create evenly spaced x values for prediction and CI
x_vals = np.linspace(subset['nearest_canon_age'].min(), subset['nearest_canon_age'].max(), 200)
y_pred = res.intercept + res.slope * x_vals

from scipy.stats import pearsonr, spearmanr

# === Compute stats ===
x = subset['nearest_canon_age']
y = subset['mean_phyloP']
n = len(x)

# Regression already done earlier with:
# res = stats.linregress(x, y)
slope = res.slope
intercept = res.intercept
r_value = res.rvalue
r2 = r_value ** 2
adj_r2 = 1 - (1 - r2) * (n - 1) / (n - 2)  # Adjusted R²

# Correlations
pearson_corr, pearson_p = pearsonr(x, y)
spearman_corr, spearman_p = spearmanr(x, y)
set_size = subset.shape[0]
# === Format stats box ===
textstr = (
    f"Linear regression (n = {n}):\n"
    f"y = {slope:.3f}x + {intercept:.3f}\n"
    f"R² = {r2:.3f}\n"
    f"Adjusted R² = {adj_r2:.3f}\n"
    f"TE subfamilies n={set_size}\n"
    f"Pearson r = {pearson_corr:.3f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_corr:.3f} (p = {spearman_p:.2e})"
)

# === Add text box to plot ===
props = dict(boxstyle='square, pad=0.5', facecolor='white', alpha=0.7)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment='left', bbox=props)


# Compute residuals and standard error
n = len(subset)
x = subset['nearest_canon_age']
y = subset['mean_phyloP']
y_fit = res.intercept + res.slope * x
residuals = y - y_fit
s_err = np.sqrt(np.sum(residuals**2) / (n - 2))

# Confidence interval
confidence = 0.95
t_val = t.ppf((1 + confidence) / 2., n - 2)
mean_x = np.mean(x)
SE = s_err * np.sqrt(1/n + (x_vals - mean_x)**2 / np.sum((x - mean_x)**2))
ci_upper = y_pred + t_val * SE
ci_lower = y_pred - t_val * SE

# Plot regression line
ax.plot(x_vals, y_pred, color='red', label='Regression Line (0–96)', linewidth=1.5)

# Plot confidence interval band
ax.fill_between(x_vals, ci_lower, ci_upper, color='red', alpha=0.2, label='95% CI')


# Styling
ax.axhline(y=0, color='grey', linewidth=1, alpha=0.5)
ax.axhline(y=3.342447251181431, color='blue', linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlim(0, 110)
ax.set_xlabel('median TE age (MYA) estimated by TEA-TIME approach', fontsize=12)
ax.set_ylabel('mean phyloP', fontsize=12)
ax.set_title('Mean phyloP grouped by TEA-TIME', fontsize=14)
plt.xticks(rotation=45, ha='right')
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/phylPvstaage_.png", dpi=300, bbox_inches="tight")

plt.show()

# %%
