#%%
import pandas as pd
import re
import argparse
import os
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table

def strict_correction(block_df):
    strict_rows = block_df[block_df['int_type'].isin(['aINT', 'bINT'])]
    if not strict_rows.empty:
        return strict_rows['te_age'].max()
    else:
        return block_df['te_age'].max()  # fallback to any available TE age
    
def preload_internal_te_ages(subfamily_internal, internal_id_dir, age_table_dir):
    """
    Load and merge internal ID and age table for a given subfamily (e.g., THE1C-int),
    then return a dictionary mapping block_id to max TE age.
    """
    try:
        int_id_path = f"{internal_id_dir}/{subfamily_internal}.internal_id.txt"
        int_age_path = f"{age_table_dir}/{subfamily_internal}.teatime.txt"

        int_id_table = pd.read_csv(int_id_path, sep="\t")
        int_id_table['block_id'] = int_id_table['internal_id'].str.extract(r'^(.*?_\d+)_')[0]

        int_age_table = pd.read_csv(int_age_path, sep="\t")

        int_full = int_id_table.merge(int_age_table, on='internal_id')

        return int_full.groupby('block_id')['te_age'].max().to_dict()

    except FileNotFoundError:
        print(f"[Warning] Cannot preload internal data for {subfamily_internal}")
        return None   

def get_base_subfamily(subfamily):
    base = re.match(r'^[A-Za-z0-9]+', subfamily).group()
    # Remove last char only if it's a letter (a-z or A-Z)
    if base[-1].isalpha():
        base = base[:-1]
    return base

def get_related_internal_te_ages(block_id, subfamily, repeatmasker_table,
                                 internal_id_dir, age_table_dir):
    """
    Scan a block for any related -int elements (same root), and return all matching te_age values.
    """
    base_subfamily = get_base_subfamily(subfamily)

    block_rmsk = repeatmasker_table[repeatmasker_table['id'] == block_id]
    repnames = block_rmsk['repName'].unique()

    # Match internal variants like THE1*, THE1C-int, THE1A-int, etc.
    int_pattern = re.compile(rf"^{re.escape(base_subfamily)}([A-Z]|_)?-int$")
    matching_ints = [name for name in repnames if int_pattern.match(name)]

    collected_ages = []

    for int_name in matching_ints:
        try:
            id_path = f"{internal_id_dir}/{int_name}.internal_id.txt"
            age_path = f"{age_table_dir}/{int_name}.teatime.txt"

            int_ids = pd.read_csv(id_path, sep='\t')
            int_ids['block_id'] = int_ids['internal_id'].str.extract(r'^(.*?_\d+)_')[0]

            int_ages = pd.read_csv(age_path, sep='\t')
            merged = int_ids.merge(int_ages, on='internal_id')

            block_ages = merged[merged['block_id'] == block_id]['te_age'].tolist()
            collected_ages.extend(block_ages)
        except FileNotFoundError:
            continue  # Skip missing files silently

    return collected_ages
#%%
def ltr_correction(subfamily, output_dir, internal_id_dir, age_table_dir, correction_mode):
    #load table
    output_filepath = f"{output_dir}/{subfamily}.ltr_fix.txt"
    if (os.path.isfile(output_filepath) == False):
        print('start',subfamily)

        internal_id_filepath = f'{internal_id_dir}/{subfamily}.internal_id.txt'
        internal_id=pd.read_csv(internal_id_filepath, sep='\t', index_col = 0)
        age_table_filepath = f'{age_table_dir}/{subfamily}.teatime.txt'
        age_table = pd.read_csv(age_table_filepath, sep='\t')
        internal_ids_age = internal_id.merge(age_table, on='internal_id')
        # Extract LTR status from ID strings
        internal_ids_age['int_type'] = internal_ids_age['internal_id'].str.extract(r'_(aINT|bINT|nINT|singleton)')
        # extract block_id
        internal_ids_age['block_id'] = internal_ids_age['internal_id'].str.extract(r'^(.*?_\d+)_')[0]
        # Identify candidate LTR blocks
        int_counts = internal_ids_age.pivot_table(index='block_id', columns='int_type', aggfunc='size', fill_value=0)
        #eligible_blocks = int_counts[(int_counts.get('aINT', 0) ==1) & (int_counts.get('bINT', 0) ==1)].index
        # Drop duplicate internal_id entries — one row per LTR piece
        dedup = internal_ids_age.drop_duplicates(subset='internal_id')

        # Now do the groupby eligibility on deduplicated data
        valid_blocks = (
            dedup.groupby(['block_id', 'int_type'])['internal_id']
            .nunique()
            .unstack(fill_value=0)
        )

        # Now filter eligible blocks
        eligible_blocks = valid_blocks[
            (valid_blocks.get('aINT', 0) == 1) & 
            (valid_blocks.get('bINT', 0) == 1)
        ].index

        # Prepare corrected age storage

        int_te_age_lookup = {}
        if correction_mode == 'ltr-int':
            expected_internal = f"{subfamily}-int"
            int_te_age_lookup = preload_internal_te_ages(expected_internal, internal_id_dir, age_table_dir)

        corrected_ages = {}
        for block_id in eligible_blocks:
            block_df = internal_ids_age[internal_ids_age['block_id'] == block_id]
            if correction_mode == 'strict' or int_te_age_lookup is None:
                # Use age from current subfamily only
                corrected_ages[block_id]=strict_correction(block_df)
            
            elif correction_mode == 'ltr-int':
                # Only accept exact internal match for this subfamily (e.g., THE1C-int)

                block_rmsk = repeatmasker_table[repeatmasker_table['id'] == block_id]
                internal_rmsk = block_rmsk[block_rmsk['repName'] == expected_internal]
                
                if not internal_rmsk.empty:
                    # Gather internal aINT/bINT ages from current subfamily
                    self_ages = block_df[block_df['int_type'].isin(['aINT', 'bINT'])]['te_age']
                    int_ages = int_te_age_lookup.get(block_id)
                    # Combine and take max
                    combined_ages = pd.concat([self_ages, int_ages])
                    corrected_ages[block_id] = combined_ages.max()

                else:
                    print('internal part with the same name not found: fallback to "strict" method')
                    corrected_ages[block_id] = strict_correction(block_df)

            elif correction_mode == 'greedy':
                self_ages = block_df[block_df['int_type'].isin(['aINT', 'bINT'])]['te_age']

                int_ages = get_related_internal_te_ages(
                    block_id, subfamily, repeatmasker_table,
                    internal_id_dir, age_table_dir
                )

                combined_ages = pd.concat([self_ages, pd.Series(int_ages)]) if int_ages else self_ages
                corrected_ages[block_id] = combined_ages.max()
        
        # 1. Extract block_id from internal_id
        age_table['block_id'] = age_table['internal_id'].str.extract(r'^(.*?_\d+)_')[0]

        # 2. Add corrected_te_age column
        age_table['corrected_te_age'] = age_table.apply(
            lambda row: corrected_ages.get(row['block_id'], row['te_age']),
            axis=1
        )
        age_table.to_csv(output_filepath, sep='\t', index=False)
        print('done',subfamily)
    else:
        print('already done', subfamily)
#%%
def main(internal_id_dir,
         age_table_dir,
         subfamily_list,
         repeatmasker_filepath,
         corrected_age_table_dir,
         correction_mode
         ):
    global repeatmasker_table
    repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
    if subfamily_list is None:
        repname_counts = repeatmasker_table['repName'].value_counts().reset_index()
        repname_counts.columns = ['repName', 'count']
        subfamily_list=repname_counts[repname_counts['count']<1000]['repName'].unique()
    for subfamily in subfamily_list:
        ltr_correction(subfamily,corrected_age_table_dir,internal_id_dir, age_table_dir,correction_mode)
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct TE ages for LTR elements using various correction modes." )
    parser.add_argument('-d', '--internal_id_dir', required=True,
                        help="Directory containing internal ID tables")
    parser.add_argument("-a", "--age_table_dir", required=True,
                        help="Directory containing TE age tables")
    parser.add_argument('-r', '--repeatmasker', required=True,
                        help="Path to RepeatMasker output file (.tsv) with modified headers compatible with pandas.read_csv()")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for corrected TE age tables")
    parser.add_argument('-s', '--subfamily_list', nargs='+', default=None,
                        help="List of subfamily names to process (optional). If not given, the script will run through all subfamilies available on the RepeatMasker output table")
    parser.add_argument("-S", "--subfamily_file",  type=str,
                        help="Optional path to a text file with a list of subfamilies (one per line).")
    parser.add_argument("-m","--correction_mode", choices=['strict', 'ltr-int', 'greedy'], default='strict',
                        help="Correction strategy: 'strict' (default): LTR only , 'ltr-int': find -int that match subfamily name, or 'greedy': find -int from subfamilies in the same family group.")
    
    args = parser.parse_args()

    subfamily_list = None
    if args.subfamily_list:
        subfamily_list = args.subfamily_list
    elif args.subfamily_file:
        with open(args.subfamily_file) as f:
            subfamily_list = [line.strip() for line in f if line.strip()]
    main(
        internal_id_dir=args.internal_id_dir,
        age_table_dir=args.age_table_dir,
        subfamily_list=subfamily_list,
        repeatmasker_filepath=args.repeatmasker,
        corrected_age_table_dir=args.output,
        correction_mode=args.correction_mode
    )
