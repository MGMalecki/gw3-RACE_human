# Changelog

Development history of this pipeline. It was built in the `pipeline_Z` repository and split out
into `gw3-RACE_human` on 2026-10-30; entries dated before then were written under the old name.
Displaced here from `README.md` when that became a user manual.
Nothing operational lives in this file; it is a record of how the current design was arrived at.

## 2026-10-21 — environment narrowed to the pipeline

`environment_pipeline.yml` (env `pipeline_z`, 13 packages) replaces the previous general-purpose
`seq` environment as the repo's declared dependency set. The workflow calls snakemake, fastp, seqtk,
bowtie2 and STAR, plus R with dplyr, tidyr, stringr and stringi; bedtools and seqkit are needed only
for the one-off reference build. The retired `seq` environment additionally carried samtools, fastqc,
fastq-pair, pysam, rseqc, bedops, modkit, jbrowse2, r-readr, r-yaml and four UCSC tools, none of
which are invoked anywhere in `Snakefile`, `scripts/` or `R_scripts/`.

Exact builds the pipeline was validated against, recorded here because the lock file they came from
(`environment.lock.txt`) is retired with the `seq` environment:

```
bedtools-2.30.0-h7d7f7ad_2      r-base-4.2.0-h1ae530e_0
bowtie2-2.4.1-py38he513fc3_0    r-dplyr-1.0.9-r42h884c59f_0
fastp-0.23.2-hb7a2d85_2         r-stringi-1.7.6-r42h884c59f_0
python-3.8.12-hb7a2778_2        r-stringr-1.4.0-r42h6115d3f_0
seqkit-2.10.0-h9ee0642_0        r-tidyr-1.2.0-r42h884c59f_0
seqtk-1.3-h5bf99c6_3            snakemake-7.32.4-hdfd78af_1
star-2.7.10b-h9ee0642_0
```

`environment_pipeline.yml` states these as minimums rather than pins, with one deliberate
divergence: it requires `python >=3.9` where the validated run used 3.8.12. The only Python in the
workflow is `R_scripts/extract_read_tags.py`, which imports `argparse` alone and is unaffected;
3.8 is end-of-life and current snakemake builds require 3.9 or newer.

## 2026-10-21 — documentation restructure

`README.md` rewritten as a manual rather than a narrative. New `docs/OUTPUT_COLUMNS.md` documents all
46 `_final_small.rds` columns and the 5 extra `_final_big.rds` columns. AT-masking recipe recovered
and documented: `scripts/mark_AT_20_genome.py` is the original script, and the masking character was
verified as `N` by decoding chr21 out of the STAR index — 3,294 short N runs clustered at 20–49 nt
(median 23), separable from 30 assembly-gap runs ≥1 kb.

Noted at the same time: neither the masked nor the unmasked genome FASTA is still on disk. The STAR
index is the only surviving copy of the masked sequence.

## 2026-10-20 — flag columns and schema stamping ("Phase A")

Six columns added to `_final_small.rds`, taking it from 40 to 46: `tail_route`, `r2_genomic_mapped`,
`softclip_3p_len`, `grep_route_reason`, `flag_grep_misrouted`, `flag_genomic_tail`. Routing and every
threshold left untouched, so the retained columns are bit-identical to the previous schema — verified
across all twelve libraries. `schema_version = 2026-10-20-flag1`.

The motivating bug: `loc_status` answers *"did R2 land where R1's gene says it should?"* and was being
read as *"does this read have an alignment at all?"*. Reads that aligned cleanly but landed elsewhere
were handed to the sequence-only peeler, which calls a genomic A-stretch a tail. The two flag columns
measure the size of that population rather than silently changing the routing.

Rebuild was idempotent (one library byte-identical across two runs), peak RSS 16.0 GiB, 58 min for
twelve libraries. Manifest with md5s at `logs/rebuild/manifest.tsv`; QC tables under
`logs/rebuild/qc/`.

## 2026-09-30 — STAR acceptance thresholds lowered to 0.50

`outFilterScoreMinOverLread` and `outFilterMatchNminOverLread` moved from STAR's default 0.66 to
0.50. These compare the aligned portion to the *full* read length, so on an 81 nt R2 they cap the
largest reportable anchored tail: 0.66 → 28 nt, 0.50 → 40 nt, 0.30 → 56 nt.

0.50 was adopted because it had zero collateral — none of the 304,592 already-anchored reads changed
their `tail_len` or class — while admitting 1,469 reads with tails ≥30 nt, 70.37% of them confirmed.
0.30 was rejected: 93.78% of its additional rescue is multimapping.

## 2026-09-25 — multimapper handling

`outFilterMultimapNmax` / `outSAMmultNmax` set to `20 / 1`. Previously `1 / 1`, under which any read
at more than one locus was discarded as unmapped (`uT:A:3`), lost its CIGAR, and fell through to the
sequence-only branch. At `20 / 1` the read keeps its best alignment and a real CIGAR, with `NH`
recording how many loci it hit.

## 2026-09-24 → 09-30 — twelve-library production run

All twelve libraries processed on the settings above, each output carrying its original FASTQ stem.
`scripts/run_all_samples.sh` added as a sequential, resumable driver, with `logs/run_status.tsv` as
its state.

## The fastp quality filter

Removed entirely rather than left configurable. `-q`/`-u` had been passed alongside `-Q` — where
fastp ignores them — while `enable_quality_filter` defaulted to `True`, so anyone running without the
config override silently got the filter on.

The reason it must stay off: Illumina quality collapses immediately after a homopolymer. Mean Phred
of the 20 nt following a tail falls from Q36.4 (no tail) to Q16.9 (57–80 nt tail) through cluster
dephasing, while the homopolymer itself holds Q34–36. Measured with `-q 30 -u 10` on 94,198
adapter-bearing reads: 34.6% of all reads dropped, and 99.4–99.9% of reads carrying a tail ≥20 nt.
Output verified byte-identical before and after removal.

## Rerun-trigger fixes

`R_scripts/extract_read_tags.py` and the two R helpers were added as declared rule inputs so
Snakemake hashes them. Before that, the `code` rerun-trigger covered only the shell command string,
which does not change when a script's contents do — an edit to a helper left `snakemake -n` reporting
"Nothing to be done" while results reflected stale code.
