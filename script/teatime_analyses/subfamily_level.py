#%%
import pandas as pd
import numpy as np
from ma_mapper import utility
#%% INPUT PARAMETERS
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
combined_table_dir = None
internal_id_dir = None
age_table_dir = None
te_alignment_dir = None
#%%
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
def consensus(matrix, cutoff=0.5):
    num_to_nucleotide = {0: "-", 1: "A", 2: "C", 3: "T", 4: "G"}
    consensus = []
    transposed_matrix = matrix.T
    
    for column in transposed_matrix:
        counts = np.bincount(column.astype(int), minlength=5)
        total = np.sum(counts)
        
        # Identify the most frequent nucleotide and its frequency
        max_count = np.max(counts)
        most_frequent_nucleotide = np.argmax(counts)
        
        # Determine if it meets the cutoff ratio
        if max_count / total >= cutoff:
            consensus.append(num_to_nucleotide[most_frequent_nucleotide])
        else:
            consensus.append("N")  # Use 'N' to indicate no clear consensus
    
    return ''.join(consensus)

def degenerate_consensus(matrix, cutoff=0.3):

    iupac_code = {
        frozenset(['A']): "A", frozenset(['C']): "C",
        frozenset(['G']): "G", frozenset(['T']): "T",
        frozenset(['A', 'G']): "R", frozenset(['C', 'T']): "Y",
        frozenset(['G', 'C']): "S", frozenset(['A', 'T']): "W",
        frozenset(['G', 'T']): "K", frozenset(['A', 'C']): "M",
        frozenset(['A', 'C', 'G']): "V", frozenset(['A', 'C', 'T']): "H",
        frozenset(['A', 'G', 'T']): "D", frozenset(['C', 'G', 'T']): "B",
        frozenset(['A', 'C', 'G', 'T']): "N",
        frozenset(['-']): "-",  # Handle positions with gaps only
        frozenset(['-', 'A']): "A", frozenset(['-', 'C']): "C",
        frozenset(['-', 'G']): "G", frozenset(['-', 'T']): "T",
        frozenset(['-', 'A', 'G']): "R", frozenset(['-', 'C', 'T']): "Y",
        frozenset(['-', 'G', 'C']): "S", frozenset(['-', 'A', 'T']): "W",
        frozenset(['-', 'G', 'T']): "K", frozenset(['-', 'A', 'C']): "M",
        frozenset(['-', 'A', 'C', 'G']): "V", frozenset(['-', 'A', 'C', 'T']): "H",
        frozenset(['-', 'A', 'G', 'T']): "D", frozenset(['-', 'C', 'G', 'T']): "B",
        frozenset(['-', 'A', 'C', 'G', 'T']): "N"
    }
    
    num_to_nucleotide = {0: "-", 1: "A", 2: "C", 3: "T", 4: "G"}
    consensus = []
    transposed_matrix = matrix.T
    
    for column in transposed_matrix:
        counts = np.bincount(column.astype(int), minlength=5)
        total = np.sum(counts)
        
        # Get nucleotides (including gaps) that meet the cutoff ratio
        nucleotides = [num_to_nucleotide[i] for i in range(5) if counts[i] / total >= cutoff]
        
        # Convert the set of nucleotides to a frozenset and use IUPAC code
        consensus.append(iupac_code[frozenset(nucleotides)] if nucleotides else "N")
    
    return ''.join(consensus)

import math

def count_transitions_transversions(seq1, seq2):
    transitions = 0
    transversions = 0

    transitions_pairs = {'AG', 'GA', 'CT', 'TC'}
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}

    for a, b in zip(seq1, seq2):
        # Skip if either sequence has 'N' or one sequence has a gap ('-')
        if a == 'N' or b == 'N' or a == '-' or b == '-':
            continue
        if a == b:
            continue
        pair = a + b
        if pair in transitions_pairs:
            transitions += 1
        elif (a in purines and b in pyrimidines) or (a in pyrimidines and b in purines):
            transversions += 1

    return transitions, transversions

def div_calc(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length.")

    transitions, transversions = count_transitions_transversions(seq1, seq2)
    valid_positions = sum(1 for a, b in zip(seq1, seq2) if a != 'N' and b != 'N' and a != '-' and b != '-')

    if valid_positions == 0:
        return float('inf')  # Infinite distance if no valid positions are left

    p = transitions / valid_positions
    q = transversions / valid_positions

    # Avoid domain errors in log and square root calculations
    if (1 - 2 * p - q) <= 0 or (1 - 2 * q) <= 0:
        return float('inf')  # Infinite distance if the formula goes out of bounds

    try:
        dist = -0.5 * math.log((1 - 2 * p - q) * math.sqrt(1 - 2 * q))
    except ValueError:
        dist = float('inf')  # Infinite distance in case of math domain errors

    return dist

def num_to_seq(sequence):
    num_to_nucleotide = {0: "-", 1: "A", 2: "C", 3: "T", 4: "G"}
    return ''.join(num_to_nucleotide.get(int(num), 'N') for num in sequence)

#%% INITIATION
age_canon = utility.age_canon()
age_ref_table_template = utility.age_reference_table()
age_ref_table = pd.DataFrame(data=age_ref_table_template).iloc[0:-1]
repeatmasker_table  = utility.repeatmasker_prep(repeatmasker_filepath)
repeatmasker_table['rmsk_index'] = repeatmasker_table.index
if combined_table_dir is None:
    combined_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
combined_te_age_filepath = f'{combined_table_dir}/all_subfamilies.itd.txt'
combined_te_age_df = pd.read_csv(combined_te_age_filepath, sep='\t')
#%%
repeatmasker_update=repeatmasker_table.merge(combined_te_age_df, how = 'left', on='rmsk_index')
#repeatmasker_update_path = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/repeatmasker_update/age_hg38.txt'
#repeatmasker_update = pd.read_csv(repeatmasker_update_path, sep='\t', low_memory=False,index_col=0)
#%%
test = repeatmasker_update[repeatmasker_update['repName']=='L1HS']
#%%
age_df = repeatmasker_update
#%%
age_df['div_percent'] = age_df['te_div']
age_df['div_fraction'] = age_df['div_percent']/100
age_df['div_age'] = age_df['div_fraction']/(2*2e-9)
age_df['div_age_mya'] = age_df['div_age']/1e+6
age_df['length'] =  age_df.genoEnd - age_df.genoStart
bins = [-0.1]
genomesize = 3049315783
for i in range(0,62):
    bins.append(i+0.9)
age_df['binned'] = pd.cut(age_df['div_percent'], bins=bins, labels=list(range(0,62)))
bins = [-1]
for i in range(0,158):
    bins.append(i+0.9)
age_df['mya_binned'] = pd.cut(age_df['div_age_mya'], bins=bins, labels=list(range(0,158)))
age_df['percent_coverage'] = age_df.length/genomesize*100
comparable_cat = ['SINE/MIR','SINE/tRNA-Deu','SINE/tRNA-RTE','SINE/tRNA','SINE/Alu','SINE/5S-Deu-L2','Retroposon/SVA','LINE/Penelope','LINE/Dong-R4','LINE/Jockey','LINE/L2','LINE/CR1','LINE/RTE-X','LINE/RTE-BovB','LINE/L1','LINE/L1-Tx1', 'LTR/ERVK','LTR/ERV1','LTR','LTR/ERVL','LTR/ERVL-MaLR','LTR/Gypsy','RC/Helitron','DNA/TcMar','DNA/TcMar-Mariner','DNA/TcMar-Pogo','DNA/TcMar-Tc1','DNA/TcMar-Tc2','DNA/TcMar-Tigger','DNA/PiggyBac','DNA/MULE-MuDR','DNA/Merlin','DNA','DNA/Kolobok','DNA/hAT','DNA/hAT-Ac','DNA/hAT-Blackjack','DNA/hAT-Charlie','DNA/hAT-Tag1','DNA/hAT-Tip100','DNA/PIF-Harbinger','Unknown']
comparable_cat_color = ['#D7B4F8','#CE9BF7','#C481F5','#B966F4','#B358F3','#A637F1','#FF4D4D','#ACD8E5','#99B3D7','#8FA1CF','#625CB1','#483AA2','#38299A','#38299A','#00008B','#00008B','#90ED90','#73CD70','#65BD61','#57AE51','#57AE51','#489E42','#FF00FF','#FF936C','#FF936C','#FF936C','#FF936C','#FF936C','#FF936C','#FF865E','#FF7850','#FF7149','#FF6A42','#FF5A34','#FF512D','#FF512D','#FF512D','#FF512D','#FF512D','#FF512D','#FF4825','#999999']
comparable_cat.reverse()
comparable_cat_color.reverse()
# %%
# Create a custom mapping for specific subclasses
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
# %% select subfamily
age_df_sub=age_df[age_df['repName']=='L1PA8']
#AGE SPREAD
age_df_sub=age_df_sub[~age_df_sub['te_age'].isna()]
age_set=utility.age_canon()
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(2,20))
grid = fig.add_gridspec(nrows = 300, ncols = 100, hspace=0)
bar_axes = []
# Calculate the maximum y-limit across all histograms
max_count = 0
for age in age_set[:7]:
    counts, _ = np.histogram(age_df_sub[age_df_sub['te_age'] == age]['div_age_mya'], bins=bins)
    if counts.max() > max_count:
        max_count = counts.max()

for idx, age in enumerate(age_set[:7]):
    bar_axes.append(fig.add_subplot(grid[0+20*idx:20+20*idx,0:100]))
    bar_axes[idx].hist(age_df_sub[age_df_sub['te_age']==age]['div_age_mya'],bins= bins)
    age_anno=age_ref_table[age_ref_table['age']==age]['representative'].values[0]
    bar_axes[idx].text(0.95, 0.85, f'Age: {age} MYA\n{age_anno}', 
            transform=bar_axes[idx].transAxes, 
            fontsize=8, 
            verticalalignment='top', 
            horizontalalignment='right')
    if idx < len(age_set[:7])-1:
        bar_axes[idx].set_xticks([])
    yticks = bar_axes[idx].get_yticks()
    # Remove the 0 from y-ticks
    yticks = yticks[yticks != 0]
    # Set the new y-ticks without 0
    bar_axes[idx].set_yticks(yticks)
    bar_axes[idx].set_xlim(0, 110)
#%%
#import statsmodels.api as sm
#from scipy import stats
age_df_sub=age_df[age_df['repName']=='L1HS']
#%%
age_df_sub=age_df_sub[~age_df_sub['te_age'].isna()]
# Calculate Q1 (25th percentile) and Q3 (75th percentile) for column 'A'
Q1 = age_df_sub['te_age'].quantile(0.25)
Q3 = age_df_sub['te_age'].quantile(0.75)
IQR = Q3 - Q1
# Define bounds for outliers
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR
# Filter out rows where 'A' is outside of the bounds
age_df_sub = age_df_sub[age_df_sub['te_age'] <= upper_bound]
#age_df_sub.fillna(0, inplace=True)
age_set=age_df_sub['te_age'].sort_values().unique()
age_count=age_df_sub.groupby(['div_age_mya', 'te_age']).count().reset_index()[['div_age_mya','te_age','genoName']].rename(columns={'genoName':'count'})

#X = sm.add_constant(age_count['te_age'])
#model = sm.WLS(age_count['binned'], X, weights=age_count['count'])
#results = model.fit()
#%%
import matplotlib.pyplot as plt

plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,10))
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
plot_axes=fig.add_subplot(grid[0:100,0:100])
plot_axes.scatter(age_count['te_age'], age_count['div_age_mya'], s=age_count['count'], alpha=0.5)
#plot_axes.plot(X, results.predict(X), color='red', linewidth=1)

plot_axes.set_ylabel('div')
plot_axes.set_xlabel('te_age (MYA)')
#%%
import pandas as pd
from ma_mapper import mapper
import numpy as np
subfamily= 'THE1C'
subfamily_filename = subfamily.replace('/','_')
if te_alignment_dir is None:
    te_alignment_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
alignment_file = f'{te_alignment_dir}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
if age_table_dir is None:
    age_table_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
teatime_filepath = f'{age_table_dir}/{subfamily_filename}.teatime.txt'
age_tbl_df = pd.read_csv(teatime_filepath, sep='\t')

if internal_id_dir is None:
    internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
internal_id_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
internal_id_df = pd.read_csv(internal_id_filepath, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')

te_age_internal_id=internal_id_sort.merge(age_tbl_df, on='internal_id', how='left')
age_default_id = pd.DataFrame()
te_age_internal_id=te_age_internal_id.merge(age_df['binned'], how='left',left_on='rmsk_index', right_index=True)
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id[['te_age','binned']] = te_age_internal_id[['te_age','binned']]
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
#%%
metadata_age.fillna(0, inplace=True)
age_set=metadata_age['te_age'].sort_values().unique()
#%%
consensus_list=[]
degenerate_consensus_list=[]
for age in age_set:
    metadata_sub=metadata_age[metadata_age['te_age']==age]
    sub_index = metadata_sub.index
    aligment_sub = alignment_filtered[sub_index]
    consensus_list.append(consensus(aligment_sub))
    degenerate_consensus_list.append(degenerate_consensus(aligment_sub))
#%%
for idx, age in enumerate(age_set):
    main_seq = degenerate_consensus_list[idx]
    for idx2, age2 in enumerate(age_set):
        test_seq = degenerate_consensus_list[idx2]
        div=div_calc(main_seq, test_seq)
        print(f'TEage\t{age}\tvs\tTEage\t{age2}:\t{div*100}')
#%%
overall_consensus=degenerate_consensus(alignment_filtered)
div_list = []
for row in alignment_filtered:
    div_list.append(div_calc(overall_consensus,num_to_seq(row))*100)
#%%
metadata_age['div_test'] = div_list
#num_to_seq(alignment_filtered[9729])
#div_calc(overall_consensus,alignment_filtered[9729])*100
#%%
bins = [-0.1]
for i in range(0,51):
    bins.append(i+0.9)
metadata_age['binned_test'] = pd.cut(metadata_age['div_test'], bins=bins, labels=list(range(0,51)))
#%%
metadata_age.fillna(0, inplace=True)
age_set=metadata_age['te_age'].sort_values().unique()
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(5,20))
grid = fig.add_gridspec(nrows = 300, ncols = 100, hspace=0)
bar_axes = []
# Calculate the maximum y-limit across all histograms
max_count = 0
for age in age_set:
    counts, _ = np.histogram(metadata_age[metadata_age['te_age'] == age]['binned'], bins=bins)
    if counts.max() > max_count:
        max_count = counts.max()

for idx, age in enumerate(age_set):
    bar_axes.append(fig.add_subplot(grid[0+20*idx:20+20*idx,0:100]))
    bar_axes[idx].hist(metadata_age[metadata_age['te_age']==age]['binned_test'],bins= bins)
    age_anno=age_ref_table[age_ref_table['age']==age]['representative'].values[0]
    bar_axes[idx].text(0.95, 0.85, f'Age: {age} MYA\n{age_anno}', 
            transform=bar_axes[idx].transAxes, 
            fontsize=10, 
            verticalalignment='top', 
            horizontalalignment='right')
    bar_axes[idx].set_ylim(0, max_count)
#%%
