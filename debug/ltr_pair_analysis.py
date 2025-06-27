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
categorized_blocks = {
    'complete': [],
    'partial': [],
    'standalone': []
}
used_indices = set()
for idx in subfamily_indices:
    if idx in used_indices:  # Skip if already used
        continue

    start_idx = max(0, idx - 20)
    end_idx = min(len(repeatmasker_table) - 1, idx + 20)
    window = repeatmasker_table.loc[start_idx:end_idx]

    target_strand = repeatmasker_table.loc[idx, 'strand']
    window = window[window['strand'] == target_strand]

    rep_names = list(window['repName'])

    ltr_positions = [i for i, name in enumerate(rep_names) if name == subfamily]
    int_positions = [i for i, name in enumerate(rep_names) if name == f"{subfamily}-int"]

    if len(ltr_positions) >= 2 and len(int_positions) >= 1:
        first_ltr = ltr_positions[0]
        last_ltr = ltr_positions[-1]
        has_int_between = any(first_ltr < pos < last_ltr for pos in int_positions)

        if has_int_between:
            extracted_table = window.iloc[first_ltr:last_ltr + 1] 
            categorized_blocks['complete'].append(extracted_table)
            used_indices.update(extracted_table.index)
            continue  

    if len(ltr_positions) >= 1 and len(int_positions) >= 1:
        first = min(ltr_positions + int_positions)
        last = max(ltr_positions + int_positions)
        extracted_table = window.iloc[first:last + 1]  
        categorized_blocks['partial'].append(extracted_table)

        used_indices.update(extracted_table.index)
        continue  

    elif len(ltr_positions) > 1: 
        first = ltr_positions[0]
        last = ltr_positions[-1]
        extracted_table = window.iloc[first:last + 1] 
        categorized_blocks['partial'].append(extracted_table)

        used_indices.update(extracted_table.index)
        continue 

    extracted_table = window[window['repName'] == subfamily]

    extracted_table = extracted_table[~extracted_table.index.isin(used_indices)]
    if not extracted_table.empty:
        categorized_blocks['standalone'].append(extracted_table)
        used_indices.update(extracted_table.index)


for category, blocks in categorized_blocks.items():
    print(f"{category}: {len(blocks)} blocks")
#%%
complete_blocks = categorized_blocks['complete']
one_sided_blocks = categorized_blocks['partial']
standalone_blocks = categorized_blocks['standalone']

global_block_count = 1 

for category, blocks in categorized_blocks.items():
    for i, block in enumerate(blocks):
        block_id = f"{subfamily}_{category}_{global_block_count}"
        block['block_id'] = block_id  # Add block_id column
        global_block_count += 1  # Increment counter
all_blocks = pd.concat([block for blocks in categorized_blocks.values() for block in blocks])
#%%
internal_id_table = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/{subfamily}.internal_id.txt'
age_table = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_table, sep='\t')
age_df = pd.read_csv(age_table, sep='\t')
age_internal_id_df = pd.merge(age_df, internal_id_df, on='internal_id')
internal_id_table_int = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/{subfamily}-int.internal_id.txt'
age_table_int = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/{subfamily}-int.txt'
internal_id_df_int = pd.read_csv(internal_id_table_int, sep='\t')
age_df_int = pd.read_csv(age_table_int, sep='\t')
age_internal_id_df_int = pd.merge(age_df_int, internal_id_df_int, on='internal_id')
age_internal_id = pd.concat([age_internal_id_df, age_internal_id_df_int])
all_blocks_age = pd.merge(all_blocks, age_internal_id,left_index=True, right_on='rmsk_index', how='left')
#%%
from ma_mapper import sequence_alignment
from Bio import SeqIO
from itertools import combinations
source_fasta = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa'
records = SeqIO.to_dict(SeqIO.parse(open(source_fasta), 'fasta'))
#%%
alignment_results = []
for block_id in all_blocks_age['block_id'].unique():
    print(block_id)
    block_df = all_blocks_age[all_blocks_age['block_id'] == block_id].copy()

    valid_internal_ids = block_df['internal_id'].dropna()
    valid_internal_ids = valid_internal_ids[valid_internal_ids.apply(lambda x: isinstance(x, str))].unique()

    merged_sequences = {}
    div_values = {}
    age_values = {}
    #print(valid_internal_ids)
    for internal_id in valid_internal_ids:
        sub_df = block_df[block_df['internal_id'] == internal_id]

        merged_seq = "".join(
            sequence_alignment.extract_sequence(row['genoName'], row['genoStart'], row['genoEnd'], row['strand'], records)
            for _, row in sub_df.iterrows()
        )
        merged_sequences[internal_id] = merged_seq
        age_values[internal_id] = sub_df['te_age'].values[0] 

    ltr_sequences = {k: v for k, v in merged_sequences.items() if "int" not in k.lower()}
    int_age = [round(v,2) for k, v in age_values.items() if "int" in k.lower()]
    
    for (id1, seq1), (id2, seq2) in combinations(ltr_sequences.items(), 2):
        age1 = age_values[id1]
        age2 = age_values[id2]
        alignment_result={
            "Block ID": block_id,
            "LTR1 ID": id1,
            "LTR2 ID": id2,
            "LTR1 age": age1,
            "LTR2 age": age2,
            "int Age": int_age,
        }

        alignment_results.append(alignment_result)
        """
        #TODO: fragment mergeing algo isn't perfect, some case it joins fragments at the opposite side of the -int, we need to come up with a better system, length aware or int aware merging 
        if block_id == 'THE1C_complete_37':
            break
        """
alignment_df = pd.DataFrame(alignment_results)
#%%
complete_block_df=alignment_df[alignment_df['Block ID'].str.contains('complete')]
block_num=len(complete_block_df['Block ID'].unique())
filtered_df = complete_block_df[complete_block_df['Block ID'].duplicated(keep=False) == False]
twoLTR_num=len(filtered_df['Block ID'].unique())
equal_teatime_df=filtered_df[filtered_df['LTR1 age']==filtered_df['LTR2 age']]
equalTEATIME_num = len(equal_teatime_df['Block ID'].unique())
print(f'subfamily: {subfamily}')
print(f'all complete LTR blocks:{block_num}')
print(f'blocks with two LTRs:{twoLTR_num}')
perc_hit=equalTEATIME_num/twoLTR_num
print(f'Two LTR with equal teatime age:{equalTEATIME_num}')
print(f'simple case with equal TEATIME/LTR block with two LTRs: {perc_hit}')
#%%
unequal_teatime_df=filtered_df[filtered_df['LTR1 age']!=filtered_df['LTR2 age']]
# %%
import pycurl
import certifi
from io import BytesIO
import os

def send_reqeust_to_gbrowser(config_request):
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.WRITEDATA, buffer)
    c.setopt(c.CAINFO, certifi.where())
    c.setopt(pycurl.POST, 1)
    c.setopt(pycurl.POSTFIELDS, config_request)
    c.setopt(pycurl.URL, 'https://genome.ucsc.edu/cgi-bin/hgTracks')
    c.setopt(pycurl.VERBOSE, 0)
    c.perform()
    c.close()

def image_reqest(image_request_command, output_dir):
    cmd = ['curl', '\''+image_request_command+'\'','>', output_dir]
    print(' '.join(cmd))
    os.system(' '.join(cmd))
hgsid = '2532104992_Y0xzL3iaoaQzHwhKKiL4n43D0urW'
import numpy as np
#%%
for _, row in unequal_teatime_df.iterrows():
    #find rmsk index
    #testrow=unequal_teatime_df.iloc[1]
    block_id = row['Block ID']
    ltr1_internal_id=row['LTR1 ID']
    ltr1_rmsk_index=internal_id_df[internal_id_df['internal_id']==ltr1_internal_id]['rmsk_index'].values[0]
    ltr1_row=subfamily_table[subfamily_table.index == ltr1_rmsk_index]
    ltr1_chrom = ltr1_row['genoName'].values[0]
    ltr1_start = ltr1_row['genoStart'].values[0]
    ltr1_end   = ltr1_row['genoEnd'].values[0]
    ltr1_strand = ltr1_row['strand'].values[0]
    ltr2_internal_id=row['LTR2 ID']
    ltr2_rmsk_index=internal_id_df[internal_id_df['internal_id']==ltr2_internal_id]['rmsk_index'].values[0]
    ltr2_row=subfamily_table[subfamily_table.index == ltr2_rmsk_index]
    ltr2_chrom = ltr2_row['genoName'].values[0]
    ltr2_start = ltr2_row['genoStart'].values[0]
    ltr2_end   = ltr2_row['genoEnd'].values[0]
    ltr2_strand = ltr2_row['strand'].values[0]
    #bs filter
    if (ltr1_chrom == ltr1_chrom) and (ltr1_strand == ltr2_strand):
        frame_start = np.min([ltr1_start, ltr2_start])
        frame_end = np.max([ltr1_end,ltr2_end])
        ltr1_coords = f'{ltr1_chrom}:{ltr1_start}-{ltr1_end}'
        print(f'ltr1 coords: {ltr1_coords}')
        ltr2_coords = f'{ltr2_chrom}:{ltr2_start}-{ltr2_end}'
        print(f'ltr1 coords: {ltr2_coords}')
        probe_coords = f'{ltr1_chrom}:{frame_start}-{frame_end}'
        extend_coords = f'{ltr1_chrom}:{frame_start-500}-{frame_end+500}'
        print(f'coords: {probe_coords}')
    else: 
        print('unexpected behavior, different chromosome and strand')
    #NOTE: for detailed customization (not needed until michael requested it)
    """
    config_request = f'hgsid={hgsid}&cons241wayViewalign.showCfg=off&cons241wayViewphyloP.showCfg=off&g=cons241way&hgsid={hgsid}&g=cons241way&displaySubtracks=all'
    send_reqeust_to_gbrowser(config_request)
    """
    image_request_command = f'https://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg38&position={extend_coords}&hgsid={hgsid}&highlight=hg38.{ltr1_coords}%23FFC107%7Chg38.{ltr2_coords}%23FFC107'
    output_filepath = f'/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/unequal_ltr_pair_img_extend/{block_id}_{probe_coords}.png'
    image_reqest(image_request_command, output_filepath)
# %%
#find missing block
all_blocks_age[all_blocks_age['block_id']=='THE1C_complete_37']
# %%
