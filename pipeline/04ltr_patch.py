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

def apply_ltrint_correction(df, internal_id_dir,age_table_dir, mode):
    if df['ltrpair_id'].isnull().all():
        return df
    global ltr_only
    df_extend=df.merge(ltr_only, left_on='rmsk_index', right_index = True)
    grouped=df_extend.dropna(subset=['ltrpair_id']).groupby('ltrpair_id')

    corrected_ages = {}
    for run_id, group in grouped:
        if len(group) <2:
            continue
        chrom = group['genoName'].iloc[0]
        strand = group['strand'].iloc[0]
        start = group['genoStart'].min()
        end = group['genoEnd'].max()
        ltrpair_id = group['ltrpair_id'].iloc[0]
        family = group['repName'].iloc[0]
        te_class = group['repClass'].iloc[0]

        candidates = ltr_by_chr_strand.get((chrom, strand))
        if mode == 'ltr-int':
            target_name = family + '-int'
            internal_hit = candidates[
                (candidates['genoStart'] >= start) &
                (candidates['genoEnd'] <= end) &
                (candidates['repName']==target_name)
            ]
        elif mode == 'greedy':
            internal_hit = candidates[
                (candidates['genoStart'] >= start) &
                (candidates['genoEnd'] <= end) &
                (candidates['repClass'] == te_class)
            ]
        else:
            raise ValueError(f"Unsupported correction mode: {mode}")
        
        if len(internal_hit) != 1:
            print(f"[INT MISMATCH] ltrpair_id '{ltrpair_id}' expected 1 internal, found {len(internal_hit)}. falling back to LTR age only")
            corrected_ages[ltrpair_id] = max(group['te_age'])
            continue
        int_entry = internal_hit.iloc[0]
        int_name = int_entry['repName']
        int_idx = int_entry.name  # rmsk_index
        int_id_path = f"{internal_id_dir}/{int_name}.internal_id.txt"
        age_path = f"{age_table_dir}/{int_name}.teatime.txt"
        
        try:
            int_id=pd.read_csv(int_id_path, sep='\t')
            
            int_ages = pd.read_csv(age_path, sep='\t')
            int_id_age = int_id.merge(int_ages, on='internal_id')

            if int_idx not in int_id_age['rmsk_index'].values:
                print(f"[WARNING] rmsk_index '{int_idx}' not found in {int_name}.teatime.txt. Falling back to LTR age only.")
                corrected_ages[ltrpair_id] = max(group['te_age'])
                continue
        except FileNotFoundError:
            print(f"[WARNING] Age file for internal '{int_name}' not found at {age_path}. Falling back to LTR age only.")
            corrected_ages[ltrpair_id] = max(group['te_age'])
            continue
        
        
        int_age = int_id_age[int_id_age['rmsk_index']==int_idx]['te_age'].values[0]
        ltr_ages = group['te_age'].tolist()
        all_ages = ltr_ages + [int_age]
    
        corrected_ages[ltrpair_id] = max(all_ages)
    df['corrected_age'] = df['ltrpair_id'].map(corrected_ages).fillna(df['te_age'])
    return df

def unpack_internal_id(df):

    def parse_id(s):
        if '__' in s:
            # Pattern for IDs with __ separating status (and no explicit ltrpair_id)
            pattern = r'(?P<family>[^_]+)_(?P<ltr_id>\d+)_(?P<internal_index>\d+)__(?P<status>.+)'
        else:
            # Pattern for IDs with ltrpair_id but no __ separator
            pattern = r'(?P<family>[^_]+)_(?P<ltr_id>\d+)_(?P<internal_index>\d+)_(?P<ltrpair_id>\d+)_(?P<status>.+)'

        match = re.match(pattern, s)
        return match.groupdict() if match else {}

    unpacked = df['internal_id'].apply(parse_id).apply(pd.Series)
    return pd.concat([df, unpacked], axis=1)

def apply_strict_correction(df):
    if df['ltrpair_id'].isnull().all():
        return df
    paired = df.dropna(subset=['ltrpair_id']).copy()
    print(paired)
    max_ages = paired.groupby('ltrpair_id')['te_age'].max()
    print(max_ages)
    df['corrected_age'] = df['ltrpair_id'].map(max_ages)
    df['corrected_age'] = df['corrected_age'].fillna(df['te_age'])
    print(df)
    return df
#%%
def ltr_correction(subfamily, output_dir, internal_id_dir, age_table_dir, correction_mode):
    #load table
    output_filepath = f"{output_dir}/{subfamily}.ltr_fix.txt"
    if (os.path.isfile(output_filepath) == False):
        print('start',subfamily)

        internal_id_filepath = f'{internal_id_dir}/{subfamily}.internal_id.txt'
        internal_id=pd.read_csv(internal_id_filepath, sep='\t')
        age_table_filepath = f'{age_table_dir}/{subfamily}.teatime.txt'
        age_table = pd.read_csv(age_table_filepath, sep='\t')
        internal_ids_age = internal_id.merge(age_table, on='internal_id')
        internal_ids_age = unpack_internal_id(internal_ids_age)

        if correction_mode == 'strict':
            updated_internal_ids_age = apply_strict_correction(internal_ids_age)
        else:
            updated_internal_ids_age = apply_ltrint_correction(internal_ids_age, internal_id_dir,age_table_dir,correction_mode)
        
        ltr_patched_df = updated_internal_ids_age[['internal_id', 'te_age','corrected_age']]
        ltr_patched_df.to_csv(output_filepath, sep='\t', index=False)
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
    global ltr_by_chr_strand, ltr_only
    # Run this once at the start before calling ltr_handler repeatedly
    ltr_only = repeatmasker_table[
        repeatmasker_table['repClass'].str.contains('LTR', na=False)
    ].copy()

    ltr_only = ltr_only[ltr_only['genoName'].str.match(r'chr[\dXY]+$')]

    cols_to_keep = ['genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass']
    ltr_only = ltr_only[cols_to_keep]

    ltr_by_chr_strand = {
        (chrom, strand): df.sort_values('genoStart' if strand == '+' else 'genoEnd')
        for (chrom, strand), df in ltr_only.groupby(['genoName', 'strand'])
    }
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
