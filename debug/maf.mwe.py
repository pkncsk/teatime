import sys
from typing import Optional

def block_overlaps(qstart: int, qend: int, bstart: int, bsize: int) -> bool:
    return not (bstart + bsize < qstart or bstart > qend)

def filter_maf_pristine(
    maf_in: str, maf_out: str, ref_prefix: str, chrom: str,
    qstart: int, qend: int, verbose: bool = False
):
    """
    Extract raw MAF blocks overlapping [qstart, qend] on reference species,
    writing them verbatim to maf_out.
    """
    with open(maf_in) as fin, open(maf_out, 'w') as fout:
        # copy header lines (##...)
        while True:
            pos = fin.tell()
            line = fin.readline()
            if not line or not line.startswith('##'):
                fin.seek(pos)
                break
            fout.write(line)

        block_num = 0
        block_lines = []
        for raw in fin:
            if raw.strip() == "":
                # end of one block
                if block_lines:
                    block_num += 1
                    header = block_lines[0].strip().split()
                    if header and header[0] == 'a':
                        # find reference 's' line
                        s_lines = [l for l in block_lines if l.startswith('s ')]
                        kept = False
                        for s in s_lines:
                            parts = s.strip().split()
                            sid = parts[1]
                            if sid.startswith(f"{ref_prefix}.{chrom}"):
                                try:
                                    bstart = int(parts[2])
                                    bsize = int(parts[3])
                                except ValueError:
                                    continue
                                if block_overlaps(qstart, qend, bstart, bsize):
                                    fout.writelines(block_lines)
                                    fout.write('\n')
                                    kept = True
                                    break
                        if verbose:
                            print(f"Block {block_num}: {'KEPT' if kept else 'SKIPPED'}")
                    else:
                        if verbose:
                            print(f"Block {block_num}: no 'a' line, SKIPPED")
                    block_lines = []
                continue
            block_lines.append(raw)

        # last block if not trailing newline
        if block_lines:
            block_num += 1
            header = block_lines[0].strip().split()
            if header and header[0] == 'a':
                s_lines = [l for l in block_lines if l.startswith('s ')]
                for s in s_lines:
                    parts = s.strip().split()
                    sid = parts[1]
                    if sid.startswith(f"{ref_prefix}.{chrom}"):
                        try:
                            bstart = int(parts[2])
                            bsize = int(parts[3])
                        except ValueError:
                            continue
                        if block_overlaps(qstart, qend, bstart, bsize):
                            fout.writelines(block_lines)
                            fout.write('\n')
                            if verbose:
                                print(f"Block {block_num}: KEPT")
                            break
            else:
                if verbose:
                    print(f"Block {block_num}: no 'a' line, SKIPPED")

if __name__ == "__main__":
    maf_input = sys.argv[1]
    maf_output = sys.argv[2]
    ref_species = sys.argv[3]        # e.g., "hg38"
    chr_name = sys.argv[4]           # e.g., "chr1"
    start = int(sys.argv[5])
    end = int(sys.argv[6])
    verbose = ("--verbose" in sys.argv)

    filter_maf_pristine(maf_input, maf_output, ref_species, chr_name, start, end, verbose)
    print("Done.")
