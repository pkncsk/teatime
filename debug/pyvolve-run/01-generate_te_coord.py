#%%
import os
import pandas as pd
import numpy as np
work_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-run/simTHE1C'
# Simulated genome size
genome_length = 10000

# Parameters for TEs
num_tes = 1
te_min_length = 2269
te_max_length = 2270

# Generate TE positions and lengths
# Step 1: Randomize the seed using a "master seed" (or from system entropy)
np.random.seed(575)
seed_list = np.random.randint(0, 2**32, size=100, dtype=np.uint32) 

for idx, seed in enumerate(seed_list):
    print(idx)
    np.random.seed(seed)
    starts = np.random.randint(0, genome_length - te_max_length, size=num_tes)
    lengths = np.random.randint(te_min_length, te_max_length + 1, size=num_tes)
    ends = np.clip(starts + lengths, None, genome_length)

    # Generate internal_id
    subfamily = "simTHE1C"
    state = "COMPLETE"
    internal_ids = [f"{subfamily}_{state}_{i}" for i in range(num_tes)]

    # Build reference table (no rmsk_index here)
    ref_df = pd.DataFrame({
        'internal_id': internal_ids,
        'genoName': ['chr1'] * num_tes,
        'genoStart': starts,
        'genoEnd': ends,
        'repName': [subfamily] * num_tes,
        'strand': np.random.choice(['+', '-'], size=num_tes)
    })

    # Sort and reset index
    ref_df = ref_df.sort_values('genoStart').reset_index(drop=True)

    # Save ref table
    
    os.makedirs(f"{work_dir}/{idx}", exist_ok=True)
    ref_path = f'{work_dir}/{idx}/simTHE1C_{idx}.ref.txt'
    ref_df.to_csv(ref_path, sep='\t', header=True, index=True)

    # Build internal_id table: index becomes rmsk_index
    internal_df = pd.DataFrame({
        'rmsk_index': ref_df.index,
        'internal_id': ref_df['internal_id']
    })

    # Save internal_id mapping
    internal_path = ref_path.replace('ref', 'internal_id')
    internal_df.to_csv(internal_path, sep='\t', header=True, index=False)

# %%
