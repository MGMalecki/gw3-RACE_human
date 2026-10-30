fasta_file = "fasta/genome/GRCh38.primary_assembly.genome.fa.split/GRCh38.primary_assembly.genome.part_chr1.fa"
bed_file = "AT_mismatch_20nt_genome_1.bed"

def find_loose_homopolymers(seq, target_base="A", min_len=20, max_mismatch=1):
    stretches = []
    i = 0
    while i < len(seq):
        start = i
        count = 0
        mismatches = 0
        while i < len(seq):
            if seq[i] == target_base:
                count += 1
            elif seq[i] in "ACGT" and mismatches < max_mismatch:
                mismatches += 1
                count += 1
            else:
                break
            i += 1
        if count >= min_len:
            stretches.append((start, start + count))
        else:
            i = start + 1  # Slide window by 1
    return stretches

with open(fasta_file, "r") as fasta, open(bed_file, "w") as bed_out:
    seq_id = None
    sequence = []

    for line in fasta:
        line = line.strip()
        if line.startswith(">"):
            if seq_id and sequence:
                seq_str = "".join(sequence).upper()
                for base in ["A", "T"]:
                    matches = find_loose_homopolymers(seq_str, target_base=base, min_len=20, max_mismatch=1)
                    for start, end in matches:
                        bed_out.write(f"{seq_id}\t{start}\t{end}\n")
            seq_id = line[1:].split()[0]
            sequence = []
        else:
            sequence.append(line)

    # Final sequence
    if seq_id and sequence:
        seq_str = "".join(sequence).upper()
        for base in ["A", "T"]:
            matches = find_loose_homopolymers(seq_str, target_base=base, min_len=20, max_mismatch=1)
            for start, end in matches:
                bed_out.write(f"{seq_id}\t{start}\t{end}\n")

print(f"BED file saved: {bed_file}")
