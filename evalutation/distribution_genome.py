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
    'SINE/Alu': 'SINE/Alu',
    'SINE/MIR': 'SINE/MIR',
    'Retroposon/SVA': 'Retroposon/SVA',
    'LINE/L1': 'LINE/L1',
    'LINE/L2': 'LINE/L2',
    'LINE/CR1': 'LINE/CR1',
    'LTR/ERV1': 'LTR/ERV1',
    'LTR/ERVL': 'LTR/ERVL',
    'LTR/ERVL-MaLR': 'LTR/ERVL',
    'RC/Helitron':'RC/Helitron',
    'DNA/TcMar':'DNA/TcMar',
    'DNA/TcMar-Mariner':'DNA/TcMar',
    'DNA/TcMar-Pogo':'DNA/TcMar',
    'DNA/TcMar-Tc1':'DNA/TcMar',
    'DNA/TcMar-Tc2':'DNA/TcMar',
    'DNA/TcMar-Tigger':'DNA/TcMar',
    'DNA/hAT': 'DNA/hAT',
    'DNA/hAT-Ac': 'DNA/hAT',
    'DNA/hAT-Blackjack': 'DNA/hAT',
    'DNA/hAT-Charlie': 'DNA/hAT',
    'DNA/hAT-Tag1': 'DNA/hAT',
    'DNA/hAT-Tip100': 'DNA/hAT',
    'Unknown': 'Unknown'
    # Add more mappings as needed
}
# Group other subclasses into their respective major class followed by "/others"
for subclass in comparable_cat:
    if subclass not in custom_groups:
        major_class = subclass.split('/')[0]
        custom_groups[subclass] = f'{major_class}/others'

# Apply the custom group to the repClass_ column
age_df['custom_group'] = age_df['repClass'].map(custom_groups).fillna(age_df['repClass'])
#%%
#summarize the table by age/and 10bin
age_count_overall = age_df.groupby(['te_age', 'mya_10binned']).count().reset_index()[['te_age','mya_10binned','genoName']].rename(columns={'genoName':'count'})
#%%
x_raw = age_df['te_age'].astype(float)
y_raw = age_df['mya_10binned'].astype(int)
#%%
def percentile_bounds(series, lower_pct=0.01, upper_pct=0.99):
    lower = series.quantile(lower_pct)
    upper = series.quantile(upper_pct)
    return lower, upper
x_min, x_max = percentile_bounds(x_raw, 0.01, 0.99)
y_min, y_max = percentile_bounds(y_raw, 0.01, 0.99)
print("x:",x_min,x_max,"y:",y_min,y_max)
# %%
x_min, x_max = percentile_bounds(x_raw, 0.10, 0.90)
y_min, y_max = percentile_bounds(y_raw, 0.10, 0.90)
print("x:",x_min,x_max,"y:",y_min,y_max)
#%%
global_ylim=110
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_overall["te_age"].unique())  # Original x-labels
y_labels = sorted(age_count_overall["mya_10binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_overall.iterrows():
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

ax.set_xlabel("tAge (MYA)", fontdict={'fontsize':10})
ax.set_ylabel("kAge (MYA)",fontdict={'fontsize':10})
ax.set_title("TE counts grouped by TEATIME and Kimura distance derived age")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig("overall_te_heatmap.png", dpi=300, bbox_inches="tight")
plt.show()

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_overall["te_age"].unique())  # Original x-labels
y_labels = sorted(age_count_overall["mya_10binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

# 🔹 Step 4: Fill the matrix with counts
for _, row in age_count_overall.iterrows():
    x_idx = x_labels.index(row["te_age"])
    y_idx = y_labels.index(row["mya_10binned"])
    heatmap_array[y_idx, x_idx] = row["count"]

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

ax.set_xlabel("TEATIME age (MYA)")
ax.set_ylabel("estimated age from kimuradistance (MYA)")
ax.set_title("Per column Z-scores of TE count grouped by TEATIME and Kimura distance derived age")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig("overall_te_heatmap_z.png", dpi=300, bbox_inches="tight")
plt.show()
#%%

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get unique x (te_age) and y (mya_10binned) values
x_labels = sorted(age_count_overall["te_age"].unique())  # Original x-labels
y_labels = sorted(age_count_overall["mya_10binned"].unique())

# 🔹 Step 1: Define bin edges for x-axis based on real values
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])  # Add extra edge for the last cell

# 🔹 Step 2: Define bin edges for y-axis (uniform in this case)
y_edges = np.array(y_labels + [y_labels[-1] + (y_labels[-1] - y_labels[-2])])

# 🔹 Step 3: Create 2D matrix filled with NaN
heatmap_array = np.full((len(y_labels), len(x_labels)), np.nan)

for _, row in age_count_overall.iterrows():
    x_idx = x_labels.index(row["te_age"])
    y_idx = y_labels.index(row["mya_10binned"])
    heatmap_array[y_idx, x_idx] = row["count"]
# Calculate per-column mean and std ignoring NaNs
column_sums = np.nansum(heatmap_array, axis=0)
column_sums[column_sums == 0] = np.nan  # Avoid division by zero
heatmap_frac = heatmap_array / column_sums  # Broadcasting

# Normalize color range between 0 and 1 (fractions)
norm = mcolors.Normalize(vmin=0, vmax=0.05)

#norm = mcolors.Normalize(vmin=0, vmax=np.nanmax(heatmap_frac))  # or fixed vmax like 0.4

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

ax.set_xlabel("tAge (MYA)", fontdict={'fontsize':12})
ax.set_ylabel("kAge (MYA)", fontdict={'fontsize':12})
ax.set_title("Per column Fraction of TE grouped by TEATIME and Kimura distance derived age")

# 🔹 Step 7: Add gridlines to match custom bin sizes
ax.grid(True, which="both", color="white", linestyle="-", linewidth=0.5, alpha=0.5)
fig.patch.set_facecolor("white")
plt.savefig("overall_te_heatmap_f.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
# Flatten the matrix and get indices of top 10 values
flat_indices = np.argpartition(heatmap_frac.flatten(), -10)[-10:]
top10_flat = flat_indices[np.argsort(heatmap_array.flatten()[flat_indices])[::-1]]

# Convert flat indices back to 2D coordinates
top10_coords = [np.unravel_index(idx, heatmap_frac.shape) for idx in top10_flat]

# Show value + coordinates
top10_values = [(coord, heatmap_frac[coord]) for coord in top10_coords]

# Print nicely
for (i, j), val in top10_values:
    print(f"Row {i} (y: {y_labels[i]}), Col {j} (x: {x_labels[j]}): Count = {val}")
# %%
col_dict={
    'Unknown': '#999999',
    'DNA/others': '#FF6A42',
    'DNA/TcMar':'#FF512D',
    'DNA/hAT': '#FF936C',
    'RC/Helitron': '#FF00FF',
    'LTR/others': '#90ED90',
    'LTR/ERVL': '#57AE51',
    'LTR/ERV1': '#73CD70',
    'LINE/others': '#ACD8E5',
    'LINE/CR1':'#99B3D7',
    'LINE/L1': '#00008B',
    'LINE/L2': '#625CB1',
    'Retroposon/SVA': '#FF4D4D',
    'SINE/others': '#C481F5',
    'SINE/Alu': '#B358F3',
    'SINE/MIR': '#D7B4F8'
}
#%%
# Count entries per te_age and custom group
age_count_by_class = (
    age_df.groupby(['te_age', 'custom_group']).size().unstack(fill_value=0).sort_index()
)

# Make sure all groups are in the color dict
age_count_by_class = age_count_by_class.reindex(columns=sorted(col_dict.keys()), fill_value=0)

# Prepare x (age) and bottom for stacking
x = age_count_by_class.index
bottom = np.zeros(len(x))
# %%
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(12,12))
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
main = fig.add_subplot(grid[0:50,20:70])
legend = fig.add_subplot(grid[0:30,72:75])

x_labels = sorted(age_count_by_class.index)  # should match your heatmap x_labels
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])
bar_widths = np.diff(x_edges)  # Width of each bin

# Initialize bottom properly as an array of zeros with same length as bins
bottom = np.zeros_like(x_labels, dtype=float)

for label in age_count_by_class.columns:
    values = age_count_by_class.loc[x_labels, label].values  # force alignment
    main.bar(x_edges[:-1], values, width=bar_widths, bottom=bottom,
           label=label, color=col_dict[label], align='edge')
    bottom += values

# Matching heatmap tick logic
main.set_xticks(x_edges)
main.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)
main.set_xlim(x_edges[0], x_edges[-1])

# X-axis formatting
main.set_xlabel("tAge (MYA)", fontdict={'fontsize':11})
main.set_ylabel("counts", fontdict={'fontsize':11})
main.set_title("Stacked TE count by TEATIME age and subclass")
#main.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')

from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib

labels = list(col_dict.keys())
colors = list(col_dict.values())

cmap = ListedColormap(colors)
bounds = range(len(colors) + 1)
norm = BoundaryNorm(bounds, cmap.N)

data_legend = fig.colorbar(
    matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
    cax=legend,
    orientation='vertical'
)
data_legend.set_ticks([(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds) - 1)])
data_legend.set_ticklabels(labels, fontsize='small')
data_legend.ax.set_title('TE Class', fontsize='small')

plt.tight_layout()
plt.savefig("stacked_te_age_by_group.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

fig = plt.figure(figsize=(10,10))
grid = fig.add_gridspec(nrows=100, ncols=100, hspace=0)
main = fig.add_subplot(grid[0:50, 20:70])
legend = fig.add_subplot(grid[0:30, 72:75])

x_labels = sorted(age_count_by_class.index)  # should match your heatmap x_labels
x_edges = np.array(x_labels + [x_labels[-1] + (x_labels[-1] - x_labels[-2])])
bar_widths = np.diff(x_edges)  # Width of each bin

bottom = np.zeros_like(x_labels, dtype=float)

for label in age_count_by_class.columns:
    values = age_count_by_class.loc[x_labels, label].values  # force alignment
    main.bar(x_edges[:-1], values, width=bar_widths, bottom=bottom,
             label=label, color=col_dict[label], align='edge')
    bottom += values

# Set x-axis ticks and labels matching heatmap style
main.set_xticks(x_edges)
main.set_xticklabels([str(x) for x in x_labels] + ["105+"], rotation=45, ha="right", fontsize=8)
main.set_xlim(x_edges[0], x_edges[-1])
main.set_xlabel("tAge (MYA)")
main.set_ylabel("Count")
main.set_title("TE count by TEATIME age estimation")

# ======== Manual legend with colored rectangles ========
legend.axis('off')  # Hide axes

n_labels = len(col_dict)
height_per_rect = 1.0  # arbitrary unit height per color block
legend.set_ylim(0, n_labels * height_per_rect)
legend.set_xlim(0, 2)  # extra width for text spacing

for i, (label, color) in enumerate(reversed(list(col_dict.items()))):
    y = n_labels * height_per_rect - (i + 1) * height_per_rect  # reversed order from bottom to top
    rect = Rectangle((0, y), 1, height_per_rect * 0.9, facecolor=color,
                     edgecolor='black', linewidth=0.8)
    legend.add_patch(rect)
    legend.text(1.1, y + height_per_rect*0.45, label, va='center', ha='left', fontsize=8)

plt.tight_layout()
plt.savefig("stacked_te_age_by_group_manual_legend_v2.png", dpi=300, bbox_inches="tight")
plt.show()


# %%
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# Example DataFrame with 'te_age' discrete and 'kimura_age' continuous
# age_df = pd.DataFrame({'te_age': [...], 'kimura_age': [...]})

# Aggregate y by x
agg = age_df.groupby('te_age').agg(
    mean_div=('div_age_mya', 'mean'),
    median_div=('div_age_mya', 'median'),
    count=('te_age', 'count')
).reset_index()

X = agg['te_age'].values.reshape(-1, 1)
y = agg['mean_div'].values
weights = agg['count'].values
#%%
# Fit weighted linear regression
model = LinearRegression()
model.fit(X, y, sample_weight=weights)

y_pred = model.predict(X)

# Calculate weighted R²
def weighted_r2(y_true, y_pred, weights):
    ss_res = np.sum(weights * (y_true - y_pred) ** 2)
    ss_tot = np.sum(weights * (y_true - np.average(y_true, weights=weights)) ** 2)
    return 1 - ss_res / ss_tot

r2 = weighted_r2(y, y_pred, weights)
print(f"Weighted R²: {r2:.3f}")
print(f"Regression slope: {model.coef_[0]:.3f}")
print(f"Regression intercept: {model.intercept_:.3f}")

# %%
from scipy.stats import pearsonr, spearmanr

# Pearson (linear) correlation on aggregated mean data
pearson_corr, pearson_p = pearsonr(agg['te_age'], agg['mean_div'])

# Spearman (rank-based) correlation, useful if non-linear but monotonic
spearman_corr, spearman_p = spearmanr(agg['te_age'], agg['mean_div'])

print(f"Pearson r: {pearson_corr:.3f} (p={pearson_p:.2e})")
print(f"Spearman ρ: {spearman_corr:.3f} (p={spearman_p:.2e})")
# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(6,4))
plt.scatter(X, y, s=weights/10000, alpha=0.7, label="Mean Kimura Age (weighted)")
plt.plot(X, y_pred, color='red', linewidth=2, label="Weighted Regression")
plt.xlabel("TEA-TIME age (MYA)")
plt.ylabel("Mean Kimura-derived age (MYA)")
plt.title("Kimura vs TEA-TIME Age Estimates")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.3)
plt.tight_layout()
plt.show()
# %% non mean
from scipy.stats import pearsonr, spearmanr

x_raw = age_df['te_age'].astype(float)           # TEA-TIME
y_raw = age_df['div_age_mya'].astype(float)      # Kimura

# Pearson correlation (linear)
pearson_r, pearson_p = pearsonr(x_raw, y_raw)

# Spearman correlation (rank-based)
spearman_rho, spearman_p = spearmanr(x_raw, y_raw)

print(f"Pearson r: {pearson_r:.3f} (p={pearson_p:.2e})")
print(f"Spearman ρ: {spearman_rho:.3f} (p={spearman_p:.2e})")


# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression

x_raw = age_df['te_age'].astype(float)
y_raw = age_df['div_age_mya'].astype(float)

# Correlations
pearson_r, pearson_p = pearsonr(x_raw, y_raw)
spearman_rho, spearman_p = spearmanr(x_raw, y_raw)

# Linear regression fit
X = x_raw.values.reshape(-1, 1)
model = LinearRegression()
model.fit(X, y_raw)
y_pred = model.predict(X)

# R² and adjusted R²
r2 = model.score(X, y_raw)
n = len(y_raw)
p = 1  # number of predictors
adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p - 1)

# Prepare regression equation text
slope = model.coef_[0]
intercept = model.intercept_
equation_text = f"y = {slope:.3f}x + {intercept:.3f}"

# Plot
plt.figure(figsize=(6, 6))
plt.scatter(x_raw, y_raw, s=5, alpha=0.3, label='Data points')
plt.plot(x_raw, y_pred, color='red', linewidth=2, label='Linear fit')

# Annotation text
textstr = (
    f"Regression equation:\n{equation_text}\n"
    f"R² = {r2:.3f}\n"
    f"Adjusted R² = {adj_r2:.3f}\n\n"
    f"Pearson r = {pearson_r:.3f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_rho:.3f} (p = {spearman_p:.2e})"
)
plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
               verticalalignment='top', bbox=dict(boxstyle="square,pad=0.5", facecolor='white', alpha=0.7))

plt.ylim(0, global_ylim)
plt.xlim(0, None)
plt.xlabel('tAge (MYA)', fontdict={'fontsize': 12})
plt.ylabel('kAge (MYA)', fontdict={'fontsize': 12})
plt.title('Raw Data and Linear Regression Fit')
plt.legend()
plt.tight_layout()
plt.savefig("overall_te_fit.png", dpi=300, bbox_inches="tight")
plt.show()



# %%

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression

# Aggregate (mean y per x), also collect count as weight
agg = age_df.groupby('te_age').agg(
    mean_kAge=('div_age_mya', 'mean'),
    count=('div_age_mya', 'count')
).reset_index()

X_mean = agg['te_age'].values.reshape(-1, 1)
y_mean = agg['mean_kAge'].values
weights = agg['count'].values

# Weighted linear regression
model = LinearRegression()
model.fit(X_mean, y_mean, sample_weight=weights)
y_pred = model.predict(X_mean)

# Weighted R²
def weighted_r2(y_true, y_pred, weights):
    ss_res = np.sum(weights * (y_true - y_pred) ** 2)
    ss_tot = np.sum(weights * (y_true - np.average(y_true, weights=weights)) ** 2)
    return 1 - ss_res / ss_tot

r2 = weighted_r2(y_mean, y_pred, weights)

# Adjusted R² using effective sample size
n_eff = (np.sum(weights))**2 / np.sum(weights**2)
p = 1  # number of predictors
adj_r2 = 1 - (1 - r2) * (n_eff - 1) / (n_eff - p - 1)

slope = model.coef_[0]
intercept = model.intercept_

# Correlations (Pearson/Spearman)
pearson_corr, pearson_p = pearsonr(agg['te_age'], agg['mean_kAge'])
spearman_corr, spearman_p = spearmanr(agg['te_age'], agg['mean_kAge'])

# Plot
plt.figure(figsize=(6, 6))
plt.plot(X_mean, y_pred, color='red', linewidth=2, label='Weighted fit on means')
plt.scatter(X_mean, y_mean, color='black', s=10, label='Mean kAge per tAge')

# Annotation text
textstr = (
    f"Weighted regression:\ny = {slope:.3f}x + {intercept:.3f}\n"
    f"Weighted R² = {r2:.3f}\n"
    f"Adjusted R² = {adj_r2:.3f}\n\n"
    f"Pearson r = {pearson_corr:.3f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_corr:.3f} (p = {spearman_p:.2e})"
)

plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
               verticalalignment='top', bbox=dict(boxstyle="square, pad=0.5", facecolor='white', alpha=0.7))

plt.ylim(0, global_ylim)
plt.xlim(0, None)
plt.xlabel('tAge (MYA)', fontsize=12)
plt.ylabel('kAge (MYA)', fontsize=12)
plt.title('Aggregated Regression Fit and Raw Data Correlation')
plt.legend()
plt.tight_layout()
plt.savefig("overall_te_fitw.png", dpi=300, bbox_inches="tight")
plt.show()

## %%

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression

# Group by te_age and get median (unweighted)
agg_median = age_df.groupby('te_age').agg(
    median_kAge=('div_age_mya', 'median')
).reset_index()

x_median = agg_median['te_age'].astype(float)
y_median = agg_median['median_kAge'].astype(float)

# Correlations (note: this is between medians, not raw)
pearson_r, pearson_p = pearsonr(x_median, y_median)
spearman_rho, spearman_p = spearmanr(x_median, y_median)

# Linear regression (unweighted)
X = x_median.values.reshape(-1, 1)
model = LinearRegression()
model.fit(X, y_median)
y_pred = model.predict(X)

# R² and adjusted R²
r2 = model.score(X, y_median)
n = len(y_median)
p = 1
adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p - 1)

# Regression line equation
slope = model.coef_[0]
intercept = model.intercept_
equation_text = f"y = {slope:.3f}x + {intercept:.3f}"

# Plot
plt.figure(figsize=(6, 6))
plt.scatter(x_median, y_median, s=10, color='black', label='Median kAge per tAge')
plt.plot(x_median, y_pred, color='red', linewidth=2, label='Unweighted fit on medians')

# Annotation
textstr = (
    f"Regression equation:\n{equation_text}\n"
    f"R² = {r2:.3f}\n"
    f"Adjusted R² = {adj_r2:.3f}\n\n"
    f"Pearson r = {pearson_r:.3f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_rho:.3f} (p = {spearman_p:.2e})"
)
plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
               verticalalignment='top', bbox=dict(boxstyle="square,pad=0.5", facecolor='white', alpha=0.7))

plt.ylim(0, global_ylim)
plt.xlim(0, None)
plt.xlabel('tAge (MYA)', fontsize=12)
plt.ylabel('kAge (MYA)', fontsize=12)
plt.title('Unweighted Linear Regression on Median kAge')
plt.legend()
plt.tight_layout()
plt.savefig("overall_te_fit_median.png", dpi=300, bbox_inches="tight")
plt.show()

# %%
