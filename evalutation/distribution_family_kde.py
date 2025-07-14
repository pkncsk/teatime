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
custom_group = {
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
#age_df['repName'] = age_df['repName'].map(custom_group).fillna(age_df['repName'])
#%%
#summarize the table by age/and 10bin
#age_count_bysubfam = age_df.groupby(['te_age', 'mya_10binned','custom_group']).count().reset_index()[['te_age','mya_10binned','custom_group','genoName']].rename(columns={'genoName':'count'})
#subfamily_of_interest = 'L1PA'
#subfamily_of_interest_filename = subfamily_of_interest.replace('/', '_')
#age_count_subfam=age_count_bysubfam[age_count_bysubfam['custom_group']==subfamily_of_interest]
#age_count_subfam=age_count_bysubfam[age_count_bysubfam['custom_group']==subfamily_of_interest]
#%%
#%%
subfamily_of_interest = 'L1PA12'
global_xlim = 73.8
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from matplotlib.lines import Line2D

# Step 1: Filter and prepare
df_subfam = age_df[age_df['repName'] == subfamily_of_interest].copy()
df_subfam['te_age'] = df_subfam['te_age'].astype(float)
df_subfam['div_age_mya'] = df_subfam['div_age_mya'].astype(float)

te_age_bins = sorted(df_subfam['te_age'].unique())
te_age_bins = [t for t in te_age_bins if 0 <= t <= 70]

# Step 2: Precompute all KDEs and get global max density
kde_data = []
global_max_density = 0
x_grid = np.linspace(0, global_xlim, 500)

for t in te_age_bins:
    subset = df_subfam[df_subfam['te_age'] == t]['div_age_mya']
    if len(subset) < 5:
        continue

    kde = gaussian_kde(subset, bw_method=0.3)
    y_vals = kde(x_grid)
    kde_data.append((t, y_vals, subset))
    global_max_density = max(global_max_density, y_vals.max())

# Step 3: Plot
fig, ax = plt.subplots(figsize=(6, 6))
n_bins = len(kde_data)
row_spacing = 1.0  # fixed spacing per row (raw y units), just offsets
plot_height = n_bins * row_spacing

for i, (t, y_vals, subset) in enumerate(kde_data):
    y_offset = i * row_spacing
    scaled_y_vals = y_vals / global_max_density  # scale all by *same* max

    ax.fill_between(x_grid, y_offset, y_offset + scaled_y_vals, alpha=0.6, color='grey')

    median = subset.median()
    mean = subset.mean()
    ax.vlines(median, y_offset, y_offset + (scaled_y_vals.max()), color='red', linestyle='-', linewidth=0.5)
    ax.vlines(mean, y_offset, y_offset + (scaled_y_vals.max()), color='black', linestyle='--', linewidth=0.5)

    ax.text(global_xlim + 1, y_offset + 0.05, f"{t} MYA", va='bottom', fontsize=10)

# Final formatting
ax.set_xlim(0, global_xlim)
ax.set_ylim(0, plot_height)
ax.set_yticks([])
ax.set_xlabel("kAge (MYA)")
ax.set_title(f"KDE of kAge per tAge: {subfamily_of_interest}", fontsize=12)

custom_lines = [
    Line2D([0], [0], color='red', lw=0.5, label='Median'),
    Line2D([0], [0], color='black', linestyle='--', lw=0.5, label='Mean'),
]
ax.legend(handles=custom_lines, loc='upper right', fontsize=8)

plt.tight_layout()
plt.savefig(f'{subfamily_of_interest}_kde.png', dpi=300, bbox_inches="tight")
plt.show()
# %%
subfamily_of_interest = 'L1HS'
global_xlim = 73.8
global_ylim = 70
# Copy and jitter te_age
df_subfam = age_df[age_df['repName'] == subfamily_of_interest].copy()
df_plot = df_subfam.copy()
jitter_strength = 1
df_plot['te_age_jitter'] = df_plot['te_age'] + np.random.uniform(-jitter_strength, jitter_strength, size=len(df_plot))

# Set up plot
fig, ax = plt.subplots(figsize=(6, 6))

# Scatterplot with jittered x
ax.scatter(
    df_plot['te_age_jitter'],
    df_plot['div_age_mya'],
    s=10,
    alpha=0.3,
    color='darkblue',
    edgecolor='none'
)

# Control axes
ax.set_xlim(0, global_xlim)
ax.set_ylim(0, global_ylim)

# Labeling
ax.set_xlabel("tAge (MYA)")
ax.set_ylabel("kAge (MYA)")
ax.set_title(f"Jittered scatter: {subfamily_of_interest}")
plt.savefig(f'{subfamily_of_interest}_jitter.png', dpi=300, bbox_inches="tight")
plt.show()

# %%
subfamily_of_interest = 'L1PA12'
global_xlim = 73.8
global_ylim = 70
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
df_subfam = age_df[age_df['repName'] == subfamily_of_interest].copy()
# Copy and round
df_plot = df_subfam.copy()
x_col = 'te_age'
y_col = 'div_age_mya'

round_factor = 0.1  # Adjust resolution
df_plot['x_round'] = (df_plot[x_col] / round_factor).round() * round_factor
df_plot['y_round'] = (df_plot[y_col] / round_factor).round() * round_factor

# Count duplicates
collapsed = df_plot.groupby(['x_round', 'y_round']).size().reset_index(name='count')

# Set up axis
fig, ax = plt.subplots(figsize=(6, 6))

# Scatter with size by count
sc = ax.scatter(
    collapsed['x_round'],
    collapsed['y_round'],
    s=collapsed['count'] * 5,  # Scale size (tweak this)
    alpha=0.6,
    color='darkblue',
    edgecolor='black',
    linewidth=0.2
)
from matplotlib.lines import Line2D

# Representative count sizes
size_legend = [1, 5, 10, 20, 50, 100]
handles = [
    Line2D([], [], marker='o', linestyle='None', markersize=np.sqrt(c * 5), 
           markerfacecolor='darkblue', alpha=0.6, label=f'{c}', 
           markeredgecolor='black', markeredgewidth=0.3)
    for c in size_legend
]

# Add legend
ax.legend(handles=handles, title='Count per point', loc='upper right', frameon=False)
# Axes control
ax.set_xlim(0, global_xlim)
ax.set_ylim(0, global_ylim)
ax.set_xlabel("tAge (MYA)")
ax.set_ylabel("kAge (MYA)")
ax.set_title(f"Aggregated scatter (dot size = count): {subfamily_of_interest}")
plt.savefig(f'{subfamily_of_interest}_aggregated_scatter.png', dpi=300, bbox_inches="tight")
plt.show()
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib import colormaps
# Input parameters
subfamily_of_interest = 'L1HS'
global_xlim = 73.8
global_ylim = 70
bin_width = 1.0

# Prepare data
df_subfam = age_df[age_df['repName'] == subfamily_of_interest].copy()
bins = np.arange(0, global_ylim + bin_width, bin_width)
bin_centers = (bins[:-1] + bins[1:]) / 2
te_age_bins = sorted(df_subfam['te_age'].unique())
te_age_bins = [t for t in te_age_bins if 0 <= t <= global_xlim]

# Precompute smoothed histograms
smoothed_lines = []
max_y = 0
for t in te_age_bins:
    subset = df_subfam[df_subfam['te_age'] == t]['div_age_mya'].dropna()
    if len(subset) < 5:
        smoothed_lines.append(None)
        continue
    counts, _ = np.histogram(subset, bins=bins)
    smooth = savgol_filter(counts, 7, 2) if len(counts) >= 7 else counts
    smoothed_lines.append(smooth)
    max_y = max(max_y, max(smooth))

global_te_ages = sorted(age_df['te_age'].dropna().unique())
global_te_ages = [t for t in global_te_ages if 0 <= t <= 73.8]
base_cmap = colormaps.get_cmap('tab20').resampled(len(global_te_ages))
te_age_color = {t: base_cmap(i) for i, t in enumerate(global_te_ages)}
# Set up color map for distinct colors


# Plot setup
fig, ax = plt.subplots(figsize=(6, 6))
row_width = 1.0
x_spacing = 0.1

legend_elements = []

for i, (t, smoothed) in enumerate(zip(te_age_bins, smoothed_lines)):
    if smoothed is None:
        continue

    x_base = i * (row_width + x_spacing)
    scaled = (smoothed / max_y) * row_width
    color = te_age_color[t]

    # Optional: light box background
    rect = Rectangle((x_base, 0), row_width, global_ylim, color='lightgray', alpha=0.1)
    ax.add_patch(rect)

    # Fill under curve
    ax.fill_betweenx(bin_centers, x_base, x_base + scaled, color=color, alpha=0.85)

    # Outline
    ax.plot(x_base + scaled, bin_centers, lw=1.2, color=color)

    # Label at the bottom


    # (optional) collect for legend
    legend_elements.append(Line2D([0], [0], color=color, lw=4, label=f"{t} MYA"))

# X-tick positions = centers of boxes
xtick_positions = [
    i * (row_width + x_spacing)
    for i in range(len(te_age_bins))
]
ax.set_xticks(xtick_positions)
ax.set_xticklabels([str(t) for t in te_age_bins], fontsize=9)
# Axis labels and formatting
ax.set_ylim(0, global_ylim)
ax.set_xlim(0, len(te_age_bins) * (row_width + x_spacing))
ax.set_ylabel("kAge (MYA)", fontdict={'fontsize':10})
ax.set_xlabel("tAge (MYA)", fontdict={'fontsize':10})
#ax.set_xticks([])  # custom labels only
ax.set_title(f"Smoothed histogram of {subfamily_of_interest}", fontsize=12)

# Adjust layout to fit bottom labels
plt.tight_layout(rect=[0.05, 0.1, 1, 1])

# Save and show
plt.savefig(f"{subfamily_of_interest}_smoothed_histogram_transposed.png", dpi=300, bbox_inches="tight")
plt.show()
# %%
