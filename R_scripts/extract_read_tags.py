import argparse

# Reads the bowtie2 END-TO-END alignment of R2 against the UNMASKED genome and writes one
# row per aligned read. bowtie2 runs with --no-unal, so presence in this file IS the
# "aligned end-to-end" signal and absence means the read did not align that way.
#
# This alignment never contributes a tail measurement -- the tail comes exclusively from
# the soft clip of the STAR alignment to the AT-masked genome. Its purpose is the opposite:
# an end-to-end alignment means the WHOLE read, putative tail included, matched the genome,
# so a non-zero tail_len on such a read is templated sequence rather than a 3' addition.
#
# Columns: 1 read_id   2 unique|multi   3 NM (edit distance, "NA" if absent)
# NM is new as of 2026-09-30; run_tail_analysis.R widens defensively, so two-column files
# written earlier still read back (NM becomes NA) and reproduce the old result exactly.


def parse_sam(sam_path, output_path):
    with open(sam_path, "r") as sam_file, open(output_path, "w") as out_file:
        for line in sam_file:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            read_id = fields[0]
            tags = fields[11:]

            # Inspect Bowtie2 optional tags for multi-mapping XS:i: flag
            is_multimapped = any(field.startswith("XS:i:") for field in tags)
            tag = "multi" if is_multimapped else "unique"

            # Edit distance of the end-to-end alignment. A read that matches the genome
            # perfectly over its full length (NM 0) is the strongest possible evidence
            # that an apparent tail on it is genomic.
            nm = "NA"
            for field in tags:
                if field.startswith("NM:i:"):
                    nm = field[5:]
                    break

            out_file.write(f"{read_id}\t{tag}\t{nm}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract mapped read names and classify as unique or multi-mappers.")
    parser.add_argument("-i", "--input", required=True, help="Input SAM file")
    parser.add_argument("-o", "--output", required=True, help="Output file to write read names and mapping status")

    args = parser.parse_args()
    parse_sam(args.input, args.output)
