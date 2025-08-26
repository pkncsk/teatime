def block_overlaps(query_start, query_end, block_start, block_length):
    block_end = block_start + block_length
    return not (block_end < query_start or block_start > query_end)

def extract_maf_region(maf_input, maf_output, target_id, region_start, region_end, verbose=False):
    block_idv = []
    block_state = False

    with open(maf_input, 'r') as file_input, open(maf_output, 'w') as file_output:
        for line in file_input:
            if line.startswith('##'):
                file_output.write(line)
                continue

            linesplit = line.strip().split()

            if not linesplit:
                if block_idv:
                    if verbose:
                        print(f"Writing block with {len(block_idv)} lines")
                    file_output.write('\n'.join(block_idv) + '\n\n')
                    block_idv = []
                block_state = False
                continue

            if linesplit[0] == 'a':
                block_state = True

            elif linesplit[0] == 's':
                if linesplit[1] == target_id:
                    start = int(linesplit[2])
                    size = int(linesplit[3])
                    if start > region_end:
                        if verbose:
                            print(f"Start {start} > region_end {region_end}, exiting loop early.")
                        return  # Exit cleanly here
                    if block_overlaps(region_start, region_end, start, size):
                        block_state = True
                    else:
                        block_state = False
                        block_idv = []

            if block_state:
                block_idv.append(line.strip())

if __name__ == '__main__':
    maf_input = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chr1.maf"
    maf_output = "/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/dev/teatime/test/chr1_region_cut.maf"
    target_id = "hg38.chr1"
    region_start = 5146456
    region_end = 7798441
    extract_maf_region(maf_input, maf_output, target_id, region_start, region_end, verbose=True)
