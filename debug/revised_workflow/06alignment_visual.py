#load package
from numpy import repeat
import pandas as pd
from Bio.AlignIO import MafIO
from concurrent.futures import ProcessPoolExecutor
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#%%
#%% INPUT PARAMETERS
target_species = 'hg38'
divergence_table_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_c.txt'
subfamily_list = ['final_LTR_sample_3000']
subfamily='final_LTR_sample_3000'
repeatmasker_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
maf_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/'
internal_id_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/evalutation/'
e_table_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/e_value_simple/'
age_table_dir = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime447/teatime_simple/'
#%% function load
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table

def divergence_table_prep(divergence_table_filepath):
    species = pd.read_table(divergence_table_filepath, sep='\t')
    species['meta_name']=species['track_name']

    return species

def get_spliced_mod(self, starts, ends, strand=1):
    # Dictionary for IUPAC ambiguity codes for 2-base combinations
    iupac_code = {
        frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S', frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }
    from typing import List, Tuple
    def trim_overlapping_fragments(starts: List[int], ends: List[int]) -> Tuple[List[int], List[int]]:
        """
        Trim overlapping TE fragments by adjusting the start positions to avoid base duplication.

        This is used to **prevent per-base duplication** when processing overlapping TE fragments,
        especially in multi-species alignments. It ensures that each base in the reference genome
        is only used once, even if it is part of multiple fragment calls.

        This function is necessary for **multi-species sequence alignments** where 
            overlapping TE fragments (e.g. from RepeatMasker) can introduce duplicated sequence 
            content if not properly merged. Without this step, fragments with physical overlap 
            in one species can be duplicated in the alignment, causing:

                - False inflation of sequence length
                - Artificially high TE conservation across species
                - Misleading phylogenetic signals
                - Incorrect ancestral reconstruction

            In contrast, if working in **single-genome**, retaining overlaps 
            may preserve internal deletions or tandem repeat events and should be handled differently.

        Parameters:
        -----------
        starts : List[int]
            List of start positions.
        ends : List[int]
            List of end positions (must be same length as `starts`).

        Returns:
        --------
        Tuple[List[int], List[int]]
            Adjusted (trimmed) start and end positions with no overlapping bases.
        """
        if not starts or not ends or len(starts) != len(ends):
            raise ValueError("Starts and ends must be non-empty and of equal length.")

        trimmed_starts = [starts[0]]
        trimmed_ends = [ends[0]]

        for i in range(1, len(starts)):
            prev_end = trimmed_ends[-1]
            curr_start = starts[i]
            curr_end = ends[i]

            # Trim current start to be >= previous end
            adjusted_start = max(curr_start, prev_end)
            if adjusted_start >= curr_end:
                # Skip completely overlapped region
                continue
            trimmed_starts.append(adjusted_start)
            trimmed_ends.append(curr_end)

        return trimmed_starts, trimmed_ends

    def convert_to_iupac(sequence):
        unique_bases = frozenset(sequence)
        if len(unique_bases) == 1:
            return sequence[0].upper()  
        return iupac_code.get(unique_bases, 'N')  # Default to 'N' for any unhandled cases
    
    def process_sequence_localized(sequence):
        sequence = sequence.upper()
        filtered_sequence = [base for base in sequence if base != '-']

        #base_counts = Counter(sequence)
        #most_common_bases = base_counts.most_common()
        #max_count = most_common_bases[0][1]
        #consensus_bases = [base for base, count in most_common_bases if count == max_count]
        new_base = convert_to_iupac(filtered_sequence)
        return new_base

    if strand not in (1, -1): 
        raise ValueError("Strand must be 1 or -1, got %s" % strand)
    
    starts, ends = trim_overlapping_fragments(starts, ends)
    fetched = list(self.search(starts, ends))
    #return fetched
    expected_letters = sum(end - start for start, end in zip(starts, ends))
    if len(fetched) == 0:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})

    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    
    exon_positions = []
    for start, end in zip(starts, ends):
        
        exon_positions.extend(range(start, end))
    exon_positions = sorted(set(exon_positions))  # Ensure unique & sorted
    
    total_rec_length = 0
    ref_first_strand = None
    split_by_position = defaultdict(lambda: defaultdict(str))
    for multiseq in fetched:
        if len(multiseq) == 1 and multiseq[0].id == self._target_seqname:
            # Fast-path for human-only alignment block
            rec = multiseq[0]
            rec_start = rec.annotations["start"]
            gapped_seq = rec.seq
            if ref_first_strand is None:
                ref_first_strand = rec.annotations["strand"]
            if "-" not in gapped_seq:
                #print('debug-humanonly: superfast route')
                # Fully ungapped, directly assign bases by offset
                for exonstart, exonend in zip(starts, ends):
                    # Clamp the range to only what's within the alignment block
                    effective_start = max(exonstart, rec_start)
                    effective_end = min(exonend, rec_start + len(gapped_seq))
           
                    if effective_start >= effective_end:
                        # Exon lies completely outside this block — skip
                        continue

                    start_idx = effective_start - rec_start
                    end_idx = effective_end - rec_start
                    for offset, base in zip(range(effective_start, effective_end), gapped_seq[start_idx:end_idx]):
                        split_by_position[rec.id][offset] = base
            else:
                #print('debug-humanonly: fast route')
                # Gapped version: build position map efficiently
                valid_indexes = []
                positions = []
                realpos = rec_start
                for i, base in enumerate(gapped_seq):
                    if base != "-":
                        valid_indexes.append((i, realpos))
                        positions.append(realpos)
                        realpos += 1

                for exonstart, exonend in zip(starts, ends):
                    i_start = bisect_left(positions, exonstart)
                    i_end = bisect_right(positions, exonend - 1)
                 
                    for i in range(i_start, i_end):
                        aln_idx, real_pos = valid_indexes[i]
                    
                        split_by_position[rec.id][real_pos] += gapped_seq[aln_idx]

            continue  # Skip rest of loop for this block
        # Find the human reference sequence
        ref = next((r for r in multiseq if r.id == self._target_seqname), None)
        if ref is None:
            raise ValueError(
                "Did not find %s in alignment bundle" % (self._target_seqname,)
            )
        
        try:
            if ref_first_strand is None:
                ref_first_strand = ref.annotations["strand"]
                if ref_first_strand not in (1, -1):
                    raise ValueError("Strand must be 1 or -1")
            elif ref_first_strand != ref.annotations["strand"]:
                raise ValueError(
                    "Encountered strand='%s' on target seqname, expected '%s'"
                    % (strand, ref_first_strand)
                )
        except KeyError:
            raise ValueError(
                "No strand information for target seqname (%s)" % self._target_seqname
            ) from None

        ref_start = ref.annotations["start"]
        ref_seq = ref.seq

        # Gapped path: build ref_pos_map
        #print("debug-humanref: fast route")
        ref_pos_map = []
        realpos = ref_start
        for i, base in enumerate(ref_seq):
            if base != "-":
                ref_pos_map.append((i, realpos))
                realpos += 1
        
        block_positions  = [pos for _, pos in ref_pos_map]
     
        # Now extract aligned fragments for all sequences using the same alignment coordinates
        for exonstart, exonend in zip(starts, ends):
            i_start = bisect_left(block_positions, exonstart)
            i_end = bisect_right(block_positions, exonend - 1)
   
            if i_start == i_end:
                continue  # This exon not represented in the alignment block

            aln_indexes = [ref_pos_map[i][0] for i in range(i_start, i_end)]
            ref_coords = block_positions[i_start:i_end]  # genomic coords covered in this exon by this block
      
            for rec in multiseq:
                rec_id = rec.id
                gapped_seq = rec.seq
                frag = "".join(gapped_seq[i] for i in aln_indexes)

                for pos, base in zip(ref_coords, frag):
                    split_by_position[rec_id][pos] += base
          
            # Count how many times each seqrec.id appears (detect duplicates)
            id_counts = Counter(seqrec.id for seqrec in multiseq)
            duplicate_ids = {seq_id for seq_id, count in id_counts.items() if count > 1}

            # Process duplicates only
            for seq_id in duplicate_ids:
                for pos, seq_str in split_by_position[seq_id].items():
                    # If sequence length > 1, means concatenation happened (duplication)
                    if len(seq_str) > 1:
                        new_seq = process_sequence_localized(seq_str)
                        split_by_position[seq_id][pos] = new_seq        
    realpos_to_len = {
        pos: len(gapped_fragment)
        for pos, gapped_fragment in split_by_position[self._target_seqname].items()
        if len(gapped_fragment) > 1
    }

    seqid_list = []
    seq_list = []
    
    for seqid in all_seqnames:
        seq_split = split_by_position[seqid]
        seq_splice = []
        filler_char = "N" if seqid == self._target_seqname else "-"
        append = seq_splice.append

        for exonstart, exonend in zip(starts, ends):
            for real_pos in range(exonstart, exonend):
                if real_pos in seq_split:
                    append(seq_split[real_pos])
                elif real_pos in realpos_to_len:
                    append(filler_char * realpos_to_len[real_pos])
                else:
                    append(filler_char)

        seqid_list.append(seqid)
        seq_list.append(Seq("".join(seq_splice))) 
    target_index = seqid_list.index(self._target_seqname)
    if len(seq_list[target_index].replace("-", "")) != expected_letters:
        raise ValueError(
            "coordinates %s , %s Returning %s letters for target seqname (%s), expected %s"
            % (
                starts,
                ends,
                len(seq_list[target_index].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(seq_list[target_index])
    for seqid, seq in zip(seqid_list, seq_list):
        if len(seq) != ref_subseq_len:
            raise ValueError(
                "Returning length %s for %s, expected %s"
                % (len(seq), seqid, ref_subseq_len)
            )
  
    # Create a DataFrame
    df = pd.DataFrame({
        'seqid': seqid_list,
        'seq': [seq.reverse_complement() if strand != ref_first_strand else seq for seq in seq_list]
    })
    return df

MafIO.MafIndex.get_spliced = get_spliced_mod
#%% set up information table
subfamily_filename = subfamily.replace('/','_') 
if internal_id_dir is None:
        internal_id_dir = '/'.join(str.split(repeatmasker_filepath, sep ='/')[:-1])
input_filepath = f'{internal_id_dir}/{subfamily_filename}.internal_id.txt'
internal_id_tbl = pd.read_csv(input_filepath, sep='\t')
repeatmasker_table = repeatmasker_prep(repeatmasker_filepath)
divergence_table=divergence_table_prep(divergence_table_filepath)
age_table_filepath = f'{age_table_dir}/{subfamily_filename}.teatime.txt'
age_tbl=pd.read_csv(age_table_filepath, sep='\t')
#%% merge information tables
idt_tbl=internal_id_tbl.merge(age_tbl, on='internal_id')
idtc_table=repeatmasker_table.merge(idt_tbl, left_index = True, right_on='rmsk_index')
#%%
def extract_block_id(internal_id):
    return "_".join(internal_id.split("_")[:2])

idtc_table['block_id'] = idtc_table['internal_id'].apply(extract_block_id)
# Identify block_ids with ANY NaN in te_age
bad_blocks = idtc_table[idtc_table['te_age'].isna()]['block_id'].unique()
filtered_age_df = idtc_table[~idtc_table['block_id'].isin(bad_blocks)].copy()
#%%
age_diff_df = (
    filtered_age_df.groupby('block_id')['te_age']
    .agg(['min', 'max'])
    .assign(age_diff=lambda x: x['max'] - x['min'])
    .reset_index()
)
same_age_pairs = age_diff_df[age_diff_df['age_diff'] == 0].copy()
ltr0_df = same_age_pairs[same_age_pairs['min']==0]
#%%
filtered_te_df=idtc_table
extender = 5000
#internal_id = 'THE1C_4280464_nINT_singleton_1'
te_list=filtered_te_df['internal_id'].unique()

for internal_id in te_list:
    print(internal_id)
    outliers_subset= filtered_te_df[filtered_te_df['internal_id']==internal_id]
    strand = outliers_subset['strand'].unique()[0]
    if strand=='-':
        strand = -1
    else:
        strand = 1
    start_list=outliers_subset['genoStart'].to_list()
    end_list=outliers_subset['genoEnd'].to_list()
    chrom = outliers_subset['genoName'].unique()[0]
    maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
    mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
    target_species = 'hg38'
    target_chrom = f'{target_species}.{chrom}'
    if extender:
        start_flanked=[min(start_list)-extender] + start_list + [max(end_list)]
        end_flanked = [min(start_list)] + end_list + [max(end_list)+extender]
    else:
        start_flanked = start_list
        end_flanked = end_list
    te_age = outliers_subset['te_age'].unique()[0]
    # extract maf
    maf_min = min(start_flanked)
    maf_max = max(end_flanked)
    index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
    spliced_maf_full =index_maf.get_spliced([maf_min],[maf_max],strand)
    # reorder maf extraction according to phylotree
    from Bio import Phylo
    phylotree_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/phylogenetic_tree/hg38.447way.nh.txt'
    # Read the tree file
    tree = Phylo.read(phylotree_filepath, "newick")
    leaf_order = [leaf.name for leaf in tree.get_terminals()]
    # Ensure index is aligned and reorder to match the tree
    spliced_maf_full['track'] = spliced_maf_full['seqid'].str.split('.').str[0]
    spliced_maf_full=spliced_maf_full.merge(divergence_table, left_on='track', right_on='track_name')
    # Create a categorical type based on tree order
    spliced_maf_full['tree_name'] = pd.Categorical(spliced_maf_full['tree_name'], categories=leaf_order, ordered=True)
    # Sort df by species and (optionally) chromosome name
    df_sorted = spliced_maf_full.sort_values(['tree_name', 'seqid'])  # or replace seq_id with 'chrom' if you have that separately
    # make annotation
    from itertools import chain

    # Flatten all regions (assume same chromosome)
    maf_ranges = list(zip(start_flanked, end_flanked))
    # Get min/max for union

    ref_seq = df_sorted[df_sorted['meta_name']==target_species]['seq'].values[0] # the aligned reference sequence

    ref_map = {}
    aln_pos = 0

    genomic_pos = maf_min
    while genomic_pos < maf_max:
        # Skip gaps in the alignment when building the map
        if ref_seq[aln_pos] != '-':
            ref_map[genomic_pos] = aln_pos
            genomic_pos += 1
        aln_pos += 1

    # Filter relevant repeats overlapping the alignment
    repeats_in_region = repeatmasker_table[
        (repeatmasker_table['genoName'] == chrom) &
        (repeatmasker_table['genoStart'] > maf_min) &
        (repeatmasker_table['genoEnd'] < maf_max)
    ].copy()
    #
    import numpy as np
    for idx, row in repeats_in_region.iterrows():
        repeats_in_region.at[idx, 'aln_start'] = ref_map.get(row['genoStart'], np.nan)
        repeats_in_region.at[idx, 'aln_end'] = ref_map.get(row['genoEnd'], np.nan)

    #
    import re

    # Assume `subfamily` is your target TE name, e.g., 'THE1C'
    # Step 1: Extract the alphabetic prefix before the final letter/digit using regex
    match = re.match(r'([A-Za-z]+)', subfamily)
    if match:
        te_prefix = match.group(1)
    else:
        raise ValueError(f"Could not extract prefix from subfamily name: {subfamily}")

    # Step 2: Filter repeats that start with this prefix
    te_family_repeats = repeats_in_region[
        repeats_in_region['repName'].str.startswith(te_prefix)
    ].copy()
    # Add annotation column: 'target' for exact match, 'related' for prefix match
    for idx, row in te_family_repeats.iterrows():
        if row['repName'] == subfamily:
            te_family_repeats.at[idx, 'annotation'] = 'target'
        else:
            te_family_repeats.at[idx, 'annotation'] = 'related'

    #plot
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    import matplotlib.patches as patches
    import matplotlib.gridspec as gridspec
    import numpy as np

    raster_mode = 'annotate'
    # Define colormap based on mode
    if raster_mode == 'detail':
        colors = ['green', 'yellow', 'red', 'blue', 'white', 'gray']  # ACGT- N
    elif raster_mode == 'annotate':
        colors = ['black', 'black', 'black', 'black', 'white', 'gray']  # ACGT = black, - = white, N = gray
    else:
        raise ValueError("Invalid raster_mode: choose 'detail' or 'annotate'")
    # Encode sequences to integers for imshow
    char_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4, 'N': 5}
    df_sorted['seq_array'] = df_sorted['seq'].apply(lambda s: [char_map.get(c.upper(), 5) for c in s])
    array = np.array(df_sorted['seq_array'].tolist())

    # Colormap for sequence letters
    colors = ['black', 'black', 'black', 'black', 'white', 'gray']
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries=np.arange(len(colors) + 1) - 0.5, ncolors=len(colors))


    n_tracks = len(df_sorted)
    annotation_rows = 1  # One row per annotation level
    total_rows = annotation_rows + n_tracks

    # Create GridSpec layout
    fig = plt.figure(figsize=(6, 0.2 * total_rows))
    grid = fig.add_gridspec(nrows=total_rows+3, ncols=1000, hspace=0)
    # Reserve top 40 rows (0–39) for annotations
    ax0 = fig.add_subplot(grid[0:1, 250:])

    # Remaining rows for MSA (40–500)
    ax1 = fig.add_subplot(grid[1:, 250:])

    # Annotation axis (for TE bars)

    ax0.set_xlim(0, array.shape[1])
    ax0.set_ylim(0, annotation_rows)
    ax0.axis('off')
    for idx, row in te_family_repeats.iterrows():
        if row['annotation'] == 'target':
            ax0.axvspan(row['aln_start'], row['aln_end'], color='red', alpha=0.2)
    # Genome bar (white with black border)
    ax0.add_patch(
        patches.Rectangle(
            (maf_min, 0),                    # x, y
            maf_max - maf_min,              # width
            1,                              # height
            edgecolor='black',
            facecolor='white'
        )
    )

    # Overlay TEs
    for _, row in te_family_repeats.iterrows():
        color = 'red' if row['annotation'] == 'target' else 'blue'
        ax0.add_patch(
            patches.Rectangle(
                (row['genoStart'], 0),                        # still in genome coords
                row['genoEnd'] - row['genoStart'],            # width
                1,
                edgecolor='none',
                facecolor=color,
                alpha=0.8
            )
        )

    # Set limits and style
    ax0.set_xlim(maf_min, maf_max)
    ax0.set_ylim(0, 1.2)



    region_strs = [f"{s}-{e}" for s, e in zip(start_list, end_list)]
    region_part = ":".join(region_strs)

    # Final title
    title = f"{internal_id}:age:{te_age}:::{chrom}:{region_part}{strand}"
    # Set title
    ax0.set_title(title, fontsize=10)# Adjust y as needed for spacing

    # Alignment axis (for MSA heatmap)
    ax1.imshow(array, cmap=cmap, norm=norm, aspect='auto', interpolation='none')
    ax1.set_yticks(np.arange(n_tracks))
    ax1.set_yticklabels(df_sorted['seqid'], fontsize = 8)
    ax1.set_xlabel('Alignment Position')
    ax1.set_xlim(0, array.shape[1])
    for idx, row in te_family_repeats.iterrows():
        if row['annotation'] == 'target':
            ax1.axvspan(row['aln_start'], row['aln_end'], color='red', alpha=0.2)

    # Add red vertical span to BOTH annotation and alignment axes

    
    #plt.tight_layout()
    #plt.show()
    image_filepath = f'/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/debug/revised_workflow/ltr0/outlier_age{te_age}_{chrom}_{maf_min}-{maf_max}.jpeg'
    fig.savefig(image_filepath,format='jpg',bbox_inches = 'tight')