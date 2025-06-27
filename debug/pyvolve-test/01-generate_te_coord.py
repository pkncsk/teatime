#%%
import pandas as pd
import numpy as np

# Simulated genome size
genome_length = 10000

# Parameters for TEs
num_tes = 1
te_min_length = 2269
te_max_length = 2270

# Generate TE positions and lengths
np.random.seed(42)
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
    'genoName': ['simchr'] * num_tes,
    'genoStart': starts,
    'genoEnd': ends,
    'repName': [subfamily] * num_tes,
    'strand': np.random.choice(['+', '-'], size=num_tes)
})

# Sort and reset index
ref_df = ref_df.sort_values('genoStart').reset_index(drop=True)

# Save ref table
ref_path = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/pyvolve-test/simTHE1C/simulated_te_ref.txt'
ref_df.to_csv(ref_path, sep='\t', header=True, index=True)

# Build internal_id table: index becomes rmsk_index
internal_df = pd.DataFrame({
    'rmsk_index': ref_df.index,
    'internal_id': ref_df['internal_id']
})

# Save internal_id mapping
internal_path = ref_path.replace('ref.txt', 'internal_id_table.txt')
internal_df.to_csv(internal_path, sep='\t', header=True, index=False)

# %%
