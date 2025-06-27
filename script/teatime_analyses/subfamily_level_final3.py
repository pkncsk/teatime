#%% initialtion
import pandas as pd
import numpy as np
from ma_mapper import utility
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
#%% INPUT PARAMETERS
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
combined_table_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/combined_age_div_lenient'
#%% table prep
age_canon = utility.age_canon() #list of TEATIME age 
age_ref_table_template = utility.age_reference_table() #table for visualization component
age_ref_table = pd.DataFrame(data=age_ref_table_template).iloc[0:-1]
repeatmasker_table  = utility.repeatmasker_prep(repeatmasker_filepath) # repeatmasker
repeatmasker_table['rmsk_index'] = repeatmasker_table.index
if combined_table_dir is None:
    combined_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
combined_te_age_filepath = f'{combined_table_dir}/all_subfamilies.itd.txt'
combined_te_age_df = pd.read_csv(combined_te_age_filepath, sep='\t') #precomputed TEATIME age
repeatmasker_update=repeatmasker_table.merge(combined_te_age_df, how = 'left', on='rmsk_index') #combining two tables for information needed
age_df = repeatmasker_update[~repeatmasker_update['te_age'].isna()].copy() #filter out entries without TEATIME
##convert kimura distances into estimate age
age_df['div_percent'] = age_df['te_div']
age_df['div_fraction'] = age_df['div_percent']/100
age_df['div_age'] = age_df['div_fraction']/(2*2e-9)
age_df['div_age_mya'] = age_df['div_age']/1e+6
age_df['length'] =  age_df.genoEnd - age_df.genoStart
##binning kimuradistance 
bins = [-0.1]
genomesize = 3049315783
for i in range(0,62):
    bins.append(i+0.9)
age_df['binned'] = pd.cut(age_df['div_percent'], bins=bins, labels=list(range(0,62)))
#binning derived age by 1
bins = [-1]
for i in range(0,158):
    bins.append(i+0.9)
age_df['mya_binned'] = pd.cut(age_df['div_age_mya'], bins=bins, labels=list(range(0,158)))
#binning derived age by 10
"""
bins = [-1]
for i in range(0,160,10):
    bins.append(i+9)
age_df['mya_10binned'] = pd.cut(age_df['div_age_mya'], bins=bins, labels=list(range(0,160,10)))
"""
#%%
bins = [-1]
for i in range(0,160):
    bins.append(i+0.9)
age_df['mya_10binned'] = pd.cut(age_df['div_age_mya'], bins=bins, labels=list(range(0,160)))
#add genome coverage
age_df['percent_coverage'] = age_df.length/genomesize*100
# Create a custom mapping for specific subclasses
comparable_cat = ['SINE/MIR','SINE/tRNA-Deu','SINE/tRNA-RTE','SINE/tRNA','SINE/Alu','SINE/5S-Deu-L2','Retroposon/SVA','LINE/Penelope','LINE/Dong-R4','LINE/Jockey','LINE/L2','LINE/CR1','LINE/RTE-X','LINE/RTE-BovB','LINE/L1','LINE/L1-Tx1', 'LTR/ERVK','LTR/ERV1','LTR','LTR/ERVL','LTR/ERVL-MaLR','LTR/Gypsy','RC/Helitron','DNA/TcMar','DNA/TcMar-Mariner','DNA/TcMar-Pogo','DNA/TcMar-Tc1','DNA/TcMar-Tc2','DNA/TcMar-Tigger','DNA/PiggyBac','DNA/MULE-MuDR','DNA/Merlin','DNA','DNA/Kolobok','DNA/hAT','DNA/hAT-Ac','DNA/hAT-Blackjack','DNA/hAT-Charlie','DNA/hAT-Tag1','DNA/hAT-Tip100','DNA/PIF-Harbinger','Unknown']
comparable_cat_color = ['#D7B4F8','#CE9BF7','#C481F5','#B966F4','#B358F3','#A637F1','#FF4D4D','#ACD8E5','#99B3D7','#8FA1CF','#625CB1','#483AA2','#38299A','#38299A','#00008B','#00008B','#90ED90','#73CD70','#65BD61','#57AE51','#57AE51','#489E42','#FF00FF','#FF936C','#FF936C','#FF936C','#FF936C','#FF936C','#FF936C','#FF865E','#FF7850','#FF7149','#FF6A42','#FF5A34','#FF512D','#FF512D','#FF512D','#FF512D','#FF512D','#FF512D','#FF4825','#999999']
comparable_cat.reverse()
comparable_cat_color.reverse()
custom_groups = {
    'L1HS': 'L1PA',
    'L1PA2': 'L1PA',
    'L1PA3': 'L1PA',
    'L1PA4': 'L1PA',
    'L1PA5': 'L1PA',
    'L1PA6': 'L1PA',
    'L1PA7': 'L1PA',
    'L1PA8': 'L1PA',
    'L1PA10': 'L1PA',
    'L1PA11': 'L1PA',
    'L1PA12': 'L1PA',
    'L1PA13': 'L1PA',
    'L1PA14': 'L1PA',
    'L1PA15': 'L1PA',
    'L1PA16': 'L1PA',
    'L1PA17': 'L1PA',
    'L1PA8A': 'L1PA',
    'L1PA15-16': 'L1PA',
    # Add more mappings as needed
}
#%%

# Apply the custom group to the repClass_ column
age_df['custom_group'] = age_df['repName'].map(custom_groups).fillna(age_df['repName'])
#%%
#summarize the table by age/and 10bin
age_count_bysubfam = age_df.groupby(['te_age', 'mya_10binned','custom_group']).count().reset_index()[['te_age','mya_10binned','custom_group','genoName']].rename(columns={'genoName':'count'})
subfamily_of_interest = 'L1PA'
subfamily_of_interest_filename = subfamily_of_interest.replace('/', '_')
#age_count_subfam=age_count_bysubfam[age_count_bysubfam['custom_group']==subfamily_of_interest]
age_count_subfam=age_count_bysubfam[age_count_bysubfam['custom_group']==subfamily_of_interest]
global_xlim = 73.8
global_ylim = 70
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_subfam["te_age"].unique())  # Original x-labels
y_labels = sorted(age_count_subfam["mya_10binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_subfam.iterrows():
    x_idx = x_labels.index(row["te_age"])
    y_idx = y_labels.index(row["mya_10binned"])
    heatmap_array[y_idx, x_idx] = row["count"]
# Mask zero or negative values
masked_array = np.ma.masked_less_equal(heatmap_array, 0)

# Set color for masked values
cmap = plt.cm.viridis.copy()
cmap.set_bad(color='lightgray')  # Replace with desired hash-like color

# Use LogNorm
norm = mcolors.LogNorm(vmin=1e-4, vmax=np.nanmax(heatmap_array))

# Define log-normalization for colors (ignoring zero values)
norm = mcolors.LogNorm(vmin=max(1, np.nanmin(heatmap_array)), vmax=np.nanmax(heatmap_array))

# 🔹 Step 5: Plot using `pcolormesh()` with custom bin edges
fig, ax = plt.subplots(figsize=(8, 6))
cax = ax.pcolormesh(x_edges, y_edges, masked_array, cmap=cmap, norm=norm, shading="flat")

for x in x_edges:
    ax.axvline(x, color="white", linewidth=0.5, alpha=0.5)
for y in y_edges:
    ax.axhline(y, color="white", linewidth=0.5, alpha=0.5)

# Add colorbar
cbar = fig.colorbar(cax)
cbar.set_label("Count")

# 🔹 Step 6: Set x-ticks at bin edges
ax.set_xticks(x_edges)
ax.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in y_tick_positions])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120
ax.set_ylim(0, global_ylim)
ax.set_xlim(0, global_xlim)
ax.set_xlabel("TEATIME age (MYA)")
ax.set_ylabel("estimated age from kimuradistance (MYA)")
ax.set_title(f"{subfamily_of_interest} TE counts grouped by TEATIME and Kimura distance derived age", fontdict={'fontsize':10})

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig(f"./fig/{subfamily_of_interest_filename}_te_heatmap_zoom.png", dpi=300, bbox_inches="tight")
plt.show()

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_subfam["te_age"].unique())  # Original x-labels
y_labels = sorted(age_count_subfam["mya_10binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_subfam.iterrows():
    x_idx = x_labels.index(row["te_age"])
    y_idx = y_labels.index(row["mya_10binned"])
    heatmap_array[y_idx, x_idx] = row["count"]

# Calculate per-column mean and std ignoring NaNs
col_means = np.nanmean(heatmap_array, axis=0)
col_stds = np.nanstd(heatmap_array, axis=0)
col_stds[col_stds == 0] = 1  # avoid division by zero

# Z-center columns
heatmap_z = (heatmap_array - col_means) / col_stds

vmax = np.nanmax(np.abs(heatmap_z))
norm = mcolors.Normalize(vmin=-4, vmax=4)

# Mask NaNs and actual zeros
masked_z = np.ma.masked_where((heatmap_array == 0) | np.isnan(heatmap_array), heatmap_z)

# Set color for masked values
cmap = plt.cm.seismic.copy()
cmap.set_bad(color='lightgray') 

# 🔹 Step 5: Plot using `pcolormesh()` with custom bin edges
fig, ax = plt.subplots(figsize=(8, 6))
cax = ax.pcolormesh(x_edges, y_edges, masked_z, cmap=cmap, norm=norm, shading="flat")

for x in x_edges:
    ax.axvline(x, color="white", linewidth=0.5, alpha=0.5)
for y in y_edges:
    ax.axhline(y, color="white", linewidth=0.5, alpha=0.5)

# Add colorbar
cbar = fig.colorbar(cax)
cbar.set_label("Z-score")

# 🔹 Step 6: Set x-ticks at bin edges
ax.set_xticks(x_edges)
ax.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in y_tick_positions])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120
ax.set_ylim(0, global_ylim)
ax.set_xlim(0, global_xlim)
ax.set_xlabel("TEATIME age (MYA)")
ax.set_ylabel("estimated age from kimuradistance (MYA)")
ax.set_title(f"Per column Z-scores of {subfamily_of_interest} TE count grouped by TEATIME and Kimura distance derived age", fontdict={'fontsize':10})

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig(f"./fig/{subfamily_of_interest_filename}_te_heatmap_z_zoom.png", dpi=300, bbox_inches="tight")
plt.show()
#%%
#calculate fraction
total_count = age_count_subfam["count"].sum()
age_count_subfam["fraction"] = age_count_subfam["count"] / total_count

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_subfam["te_age"].unique())  # Original x-labels
y_labels = sorted(age_count_subfam["mya_10binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_subfam.iterrows():
    x_idx = x_labels.index(row["te_age"])
    y_idx = y_labels.index(row["mya_10binned"])
    heatmap_array[y_idx, x_idx] = row["fraction"]
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

masked_fraction = np.ma.masked_where((heatmap_array == 0) | np.isnan(heatmap_array), heatmap_array)

# Step 7: Set colormap and color for masked values
cmap = plt.cm.viridis.copy()
cmap.set_bad('lightgray')  # or 'white', etc.

# Step 8: Normalize (no log)
norm = mcolors.Normalize(vmin=0, vmax=np.nanmax(heatmap_array))  # or fixed vmax like 0.4
# 🔹 Step 5: Plot using `pcolormesh()` with custom bin edges
fig, ax = plt.subplots(figsize=(8, 6))
cax = ax.pcolormesh(x_edges, y_edges, masked_fraction, cmap=cmap, norm=norm, shading="flat")

for x in x_edges:
    ax.axvline(x, color="white", linewidth=0.5, alpha=0.5)
for y in y_edges:
    ax.axhline(y, color="white", linewidth=0.5, alpha=0.5)

# Add colorbar
cbar = fig.colorbar(cax)
cbar.set_label("Fraction of grand total")

# 🔹 Step 6: Set x-ticks at bin edges
ax.set_xticks(x_edges)
ax.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in y_tick_positions])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120
ax.set_ylim(0, global_ylim)
ax.set_xlim(0, global_xlim)
ax.set_xlabel("TEATIME age (MYA)")
ax.set_ylabel("Estimated age from kimuradistance (MYA)")
ax.set_title(f"Fraction of {subfamily_of_interest} TE grouped by TEATIME and Kimura distance derived age", fontdict={'fontsize':10})

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig(f"./fig/{subfamily_of_interest_filename}_te_heatmap_f_zoom.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_subfam["te_age"].unique())  # Original x-labels
y_labels = sorted(age_count_subfam["mya_10binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_subfam.iterrows():
    x_idx = x_labels.index(row["te_age"])
    y_idx = y_labels.index(row["mya_10binned"])
    heatmap_array[y_idx, x_idx] = row["count"]

# Normalize each column by its sum (fraction of column total)
column_sums = np.nansum(heatmap_array, axis=0)
column_sums[column_sums == 0] = np.nan  # Avoid division by zero
heatmap_frac = heatmap_array / column_sums  # Broadcasting

# Normalize color range between 0 and 1 (fractions)
norm = mcolors.Normalize(vmin=0, vmax=1)

norm = mcolors.Normalize(vmin=0, vmax=np.nanmax(heatmap_frac))  # or fixed vmax like 0.4

# Mask NaNs and actual zeros
masked_z = np.ma.masked_where((heatmap_frac == 0) | np.isnan(heatmap_frac), heatmap_frac)

# Set color for masked values
cmap = plt.cm.plasma.copy()
cmap.set_bad(color='lightgray') 

# 🔹 Step 5: Plot using `pcolormesh()` with custom bin edges
fig, ax = plt.subplots(figsize=(8, 6))
cax = ax.pcolormesh(x_edges, y_edges, masked_z, cmap=cmap, norm=norm, shading="flat")

for x in x_edges:
    ax.axvline(x, color="white", linewidth=0.5, alpha=0.5)
for y in y_edges:
    ax.axhline(y, color="white", linewidth=0.5, alpha=0.5)

# Add colorbar
cbar = fig.colorbar(cax)
cbar.set_label("Fraction of column total")

# 🔹 Step 6: Set x-ticks at bin edges
ax.set_xticks(x_edges)
ax.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)

# Set y-ticks normally
ax.set_yticks(y_labels)
y_tick_positions = np.arange(0, len(y_labels), 10)
ax.set_yticks(y_tick_positions)
ax.set_yticklabels([y_labels[i] for i in y_tick_positions])
#ax.set_yticklabels(y_labels)
# Truncate y-axis to exclude outlier above 120
ax.set_ylim(0, global_ylim)
ax.set_xlim(0, global_xlim)
ax.set_xlabel("TEATIME age (MYA)")
ax.set_ylabel("estimated age from kimuradistance (MYA)")
ax.set_title(f"Per column fractions of {subfamily_of_interest} TE count grouped by TEATIME and Kimura distance derived age", fontdict={'fontsize':10})

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig(f"./fig/{subfamily_of_interest_filename}_te_heatmap_z2_zoom.png", dpi=300, bbox_inches="tight")
plt.show()
#%%