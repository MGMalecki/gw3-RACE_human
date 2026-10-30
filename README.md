# gw3-RACE_human

A Snakemake pipeline that turns paired-end 3'-end RACE sequencing (3' RACE-seq) into per-read
measurements of non-templated 3' tails. For every molecule it reports whether the transcript carries
a tail, how long that tail is, what it is made of — poly(A), uridylation, or mixed — and which gene
it came from. Input is two FASTQ files per library; output is two tidy R tables.

**Contents**

- [What it answers, and what it does not](#what-it-answers-and-what-it-does-not)
- [Input requirements](#input-requirements)
- [Installation and references](#installation-and-references)
- [Running it on your data](#running-it-on-your-data)
- [How the pipeline works](#how-the-pipeline-works)
- [Output](#output)
- [QC](#qc)
- [Caveats and limitations](#caveats-and-limitations)
- [Worked example: the reference dataset](#worked-example-the-reference-dataset)
- [Repository layout](#repository-layout)

---

## What it answers, and what it does not

Per gene, across a library, it supports:

- the **tailed fraction** — what share of molecules carry a non-templated 3' extension;
- the **tail length distribution**;
- the **tail composition classes** — pure poly(A), oligo-U, A+U mixtures, and so on.

It is **not** an isoform quantifier, not a general expression pipeline, and not a polyadenylation-site
caller. Gene assignment is deliberately collapsed to the gene level
([why](docs/OUTPUT_COLUMNS.md#1-read-identity-and-provenance)), and 3'-end positions are reported
relative to gene-level annotation bounds, not transcript ends.

## Input requirements

### The mate contract

This is the assumption that silently ruins results if it is wrong:

| Mate | Role |
|---|---|
| **R1** | gene-identifying. Aligned to the transcriptome; its only job is to say which gene the molecule came from |
| **R2** | 3'-end. Carries the non-templated tail. Every tail measurement comes from this mate |

A library with the two swapped produces output that looks entirely normal and means nothing. Check
this before anything else.

### Read structure

- R2 must contain the protocol's adapter signature — 15 nt, then `GTCAG`. Reads without it are
  discarded early.
- A **15 nt UMI** is expected at the start of R2 and is moved into the read header.
- The first 5 nt after the UMI are trimmed.
- R2 must be **≥80 nt** after trimming; shorter reads are dropped.

All four are in the `fastp_filtering_R2` rule of the `Snakefile`. A different UMI length or adapter
means editing that rule — they are not currently config keys.

### References

| Reference | Used for |
|---|---|
| Primary genome FASTA | the bowtie2 unmasked genome index |
| **AT-masked** primary genome FASTA | the STAR index — this is what makes tails measurable |
| GENCODE annotation GTF | splice junctions in the STAR index |
| Transcriptome FASTA | the bowtie2 index R1 is aligned to |
| Gene annotation table (`.rds`) | gene coordinates, strand and biotype for the R stage |

The reference set was built from GRCh38 primary assembly with GENCODE v47 (junctions) and a v50 gene
table. Substituting other versions is fine; see [Caveats](#caveats-and-limitations) for what the
mismatch costs.

## Installation and references

```bash
git clone https://github.com/MGMalecki/gw3-RACE_human.git && cd gw3-RACE_human
conda env create -f environment_pipeline.yml     # creates env "pipeline_z"
conda activate pipeline_z
```

That file is the pipeline's complete dependency set — 13 packages. The workflow calls **snakemake**,
**fastp**, **seqtk**, **bowtie2** and **STAR**, plus R with **dplyr**, **tidyr**, **stringr** and
**stringi**; `bedtools` and `seqkit` are there for the one-off reference build and can be commented
out if you already have the indexes. Nothing else is used anywhere in `Snakefile`, `scripts/` or
`R_scripts/`.

The validated version set is fastp 0.23.2, seqtk 1.3, bowtie2 2.4.1, STAR 2.7.10b, bedtools 2.30.0,
seqkit 2.10.0, snakemake 7.32.4 — recorded with exact build strings in
[`docs/CHANGELOG.md`](docs/CHANGELOG.md). `environment_pipeline.yml` states these as minimums rather
than pins; **the one place to be careful is STAR**, because an index can only be read by a compatible
STAR version — build your own and you are self-consistent, inherit someone else's and you should
match their version (check `versionGenome` in the index's `genomeParameters.txt`).

**Disk and memory, before you start:** the STAR index alone is ~23.6 GB on disk and ~24 GB resident
when loaded, which is why libraries are processed one at a time. Intermediate `R_input/` TSVs run
0.2–3 GB per library. Budget ~100 GB plus your FASTQs.

### Building the AT-masked genome

R2's tail is measured as the 3' soft clip of its alignment. Against an unmasked genome a poly(A) or
poly(U) tail aligns straight into any nearby genomic A- or T-stretch: the aligner extends through the
tail rather than clipping it, and the measurement disappears. Masking long A and T runs out of the
reference removes those landing sites, so the tail has nowhere to go but the soft clip. **Every
anchored tail call rests on this.**

`scripts/mark_AT_20_genome.py` finds the runs and writes them as BED. The rule it implements is a
*loose* homopolymer with one mismatch allowed:

- A-runs and T-runs are scanned separately, over the uppercased sequence, so soft-masked repeat
  regions are included.
- A run may contain at most **one** non-target base; the second ends it, so `ATATAT…` never
  accumulates.
- `min_len = 20` is the length of the whole run **including** that permitted mismatch. 19 target
  bases plus one other base therefore qualifies. Do not "fix" this to a strict `A{20,}` — it would
  change the reference.
- `N` ends a run, so assembly gaps are never bridged.

One mismatch rather than an exact match is deliberate: genomic A-stretches are mostly imperfect, and
a tail can still align into an interrupted stretch, so an exact rule would leave the commonest
landing sites unmasked.

```bash
# 1. Locate the A/T runs. Paths are hardcoded at the top of the script; it was run per
#    chromosome over a `seqkit split` directory, then the BEDs concatenated. The scan is a
#    pure-Python per-base loop -- minutes per large chromosome, so split rather than running
#    it over the whole assembly at once.
python scripts/mark_AT_20_genome.py

# 2. Apply them. Masking must be to N: soft (lowercase) masking is uppercased by STAR at
#    index build and would silently do nothing.
bedtools maskfasta -fi <GRCh38.primary_assembly.genome.fa> \
                   -bed AT_mismatch_20nt_genome.bed \
                   -fo GRCh38.primary_masked_AT_new.fa
```

The BED is unsorted and unmerged and A-intervals may overlap T-intervals; `bedtools maskfasta`
accepts that. Both `seqkit` and `bedtools` are in `environment_pipeline.yml`.

### Building the indexes

```bash
# STAR, masked genome
STAR --runMode genomeGenerate \
     --genomeDir fasta/index_masked_genome_star/ \
     --genomeFastaFiles GRCh38.primary_masked_AT_new.fa \
     --sjdbGTFfile <gencode.vNN.primary_assembly.basic.annotation.gtf> \
     --sjdbOverhang 79 --genomeSAindexNbases 14 --runThreadN 16

# bowtie2, unmasked genome (the templated-tail check)
bowtie2-build --threads 16 <GRCh38_primary_assembly.fa> \
  fasta/index_bowtie_genome/GRCh38_primary_genome

# bowtie2, transcriptome (R1 gene assignment)
bowtie2-build --threads 16 <full_transcriptome.fa> \
  fasta/index_bowtie_full_transcriptome/full_trans_
```

`--sjdbOverhang 79` matches the 80 nt read-length floor; change it if you change that.

## Running it on your data

### 1. Put your FASTQs where the pipeline can see them

**R1 and R2 go in two different directories.** That is the convention the repo uses and what
`.gitignore` already covers:

```
fastq_read1/           <- R1 files only, the gene-identifying mate
fastq_read2/           <- R2 files only, the 3'-end mate that carries the tail
```

So one library contributes one file to each:

```
fastq_read1/MYSAMPLE_A_R1_001.fastq.gz
fastq_read2/MYSAMPLE_A_R2_001.fastq.gz
```

Gzipped is fine. The two directories are a convention, not a requirement — the paths in
`config/samples.tsv` are read as-is, so any readable location works — but keeping R1 and R2 apart
makes it much harder to swap the mates by accident, which matters because
[swapped mates are undetectable from the output](#the-mate-contract).

### 2. List your libraries — `config/samples.tsv`

Tab-separated, one row per library, with this header. A two-library example:

```
output_name	pipeline_name	genotype	replicate	sample_index	fastq_r1	fastq_r2
wt_1	MYSAMPLE_A	WT	1	S01	fastq_read1/MYSAMPLE_A_R1_001.fastq.gz	fastq_read2/MYSAMPLE_A_R2_001.fastq.gz
ko_1	MYSAMPLE_B	KO	1	S02	fastq_read1/MYSAMPLE_B_R1_001.fastq.gz	fastq_read2/MYSAMPLE_B_R2_001.fastq.gz
```

(Real tabs between fields, not spaces. The shipped file has twelve rows you can use as a model.)

| Column | Required | What it does |
|---|---|---|
| `output_name` | label only | short handle for your own use; also accepted by `ONLY=` |
| `pipeline_name` | **yes** | **the stem** — see below |
| `genotype` | label only | grouping label, carried into the QC report |
| `replicate` | label only | — |
| `sample_index` | label only | — |
| `fastq_r1` | **yes** | path to that library's R1 file |
| `fastq_r2` | **yes** | path to that library's R2 file |

### What "the stem" is, and where you set it

**The stem is whatever you put in the `pipeline_name` column.** It is the one identifier you choose,
it is set in that column and nowhere else, and everything downstream is named from it. For
`pipeline_name = MYSAMPLE_A`:

| | |
|---|---|
| its generated config | `config/generated/MYSAMPLE_A.yaml` |
| its outputs | `R_output/MYSAMPLE_A_final_small.rds`, `..._final_big.rds`, `..._deduplication.csv` |
| its logs | `logs/run_MYSAMPLE_A_<timestamp>.log`, a row keyed `MYSAMPLE_A` in `logs/run_status.tsv` |
| what you pass to run one library | `ONLY=MYSAMPLE_A` |
| what you pass to re-finalise it without re-aligning | `Rscript scripts/rebuild/rebuild_one.R MYSAMPLE_A` |

It does **not** have to match your FASTQ filenames — the paths are given explicitly in the two
`fastq_*` columns. Matching them is just convenient. Anything filesystem-safe works; avoid spaces.

Everywhere `<stem>` appears below, substitute your own `pipeline_name`.

### 3. Point at your references — `config/config_template.yaml`

`config/` holds exactly **two files you edit** plus one directory you don't:

| Path | Yours to edit? |
|---|---|
| `config/samples.tsv` | **yes** — your libraries |
| `config/config_template.yaml` | **yes** — reference paths and STAR parameters |
| `config/generated/<stem>.yaml` | **no** — written by the driver, one per library, as a record of the exact config each run used |

(`config/generated/` currently holds 17 files: the twelve libraries plus five parameter-trial
variants. A top-level `config.yaml` also exists and is legacy — see below.)

The driver copies the template once per library, substituting `__SAMPLE__`, `__R1__` and `__R2__`
from the `samples.tsv` row. You edit everything else. The `reference:` block is the part you must
change:

```yaml
reference:
  transcriptome_bowtie2_index:  fasta/index_bowtie_full_transcriptome/full_trans_
  genome_bowtie2_index:         fasta/index_bowtie_genome/GRCh38_primary_genome
  genome_star_masked_index_dir: fasta/index_masked_genome_star/
  gene_annotation_rds:          annotation/gene_annotations_v50.rds
```

These are index *prefixes* for bowtie2 and a *directory* for STAR, matching what
[Building the indexes](#building-the-indexes) produces. Leave the STAR parameters below that block
alone unless you have read [stage 3](#3-r2--at-masked-genome-star) — they change what is measurable,
not just how fast it runs.

### 4. Run

**All libraries**, sequentially — two concurrent STAR jobs would want ~48 GB, so the driver never
runs them in parallel:

```bash
conda activate pipeline_z
CORES=16 bash scripts/run_all_samples.sh
```

**One library** — same driver, filtered by stem (or by `output_name`). This writes that library's
config for you before running it:

```bash
ONLY=MYSAMPLE_A bash scripts/run_all_samples.sh
```

**Dry run**, to see the plan without executing it. This needs the library's generated config to
exist already, so either run `ONLY=` once or copy the template by hand:

```bash
snakemake --configfile config/generated/MYSAMPLE_A.yaml -n
```

Progress goes to `logs/run_status.tsv`, one row per library.

**Re-running is safe and resumes.** Any library already recorded with exit 0 is skipped, so an
interruption costs at most one library of compute. A failure is recorded, keeps its intermediates for
diagnosis, and the loop carries on to the next library.

> **To force a full re-run, move `logs/run_status.tsv` aside first** — otherwise the driver reports
> "already completed, skipping" for every library and exits in seconds having done nothing.

### 5. Nothing — the tables are finished

`snakemake` writes `_final_small.rds` already typed, with fixed factor levels and the
`schema_version` stamp. There is no second command.

If you ever need to re-apply the schema to one library *without* re-running the alignments —
after editing the class vocabulary, say — `scripts/rebuild/rebuild_one.R` does exactly that from
the existing `R_input/` tables:

```bash
Rscript scripts/rebuild/rebuild_one.R MYSAMPLE_A                  # in place
Rscript scripts/rebuild/rebuild_one.R MYSAMPLE_A --scratch tmp/   # or to a scratch dir
```

It calls the same two functions the workflow calls, in the same order, so the two paths cannot
drift. It does not rewrite `_final_big.rds`.

### Calling snakemake directly

The driver is a convenience wrapper, not a requirement. If you invoke `snakemake` yourself, two
flags matter:

```bash
snakemake --configfile config/generated/<stem>.yaml --cores 16 \
  --default-resources "tmpdir='$PWD/tmp'"
```

- **`--default-resources "tmpdir='$PWD/tmp'"`** — without it Snakemake points `resources.tmpdir` at
  `/tmp`, and the `sort` in `extract_read1_names` can spill several GB there. The driver sets this;
  a bare `snakemake` call does not.
- **`--notemp`** keeps the `temp()` intermediates (per-stage SAM and FASTQ, ~12 GB per library) when
  you are debugging a stage. Omit it otherwise.

> **A config is mandatory and there is no default.** The workflow reads `sample`, `R1` and `R2` from
> `--configfile` only; with none supplied it stops immediately and tells you what to pass. (A
> hardcoded `configfile: "config.yaml"` used to sit at the top of the `Snakefile`, which meant the
> workflow refused to start unless that file existed and silently merged it under whatever you
> passed. Both behaviours are gone.) `run_pipeline.sh`, the top-level `config.yaml` and `RUN_NOTES.md`
> are leftovers from an early single-sample run and are not used.

## How the pipeline works

Nine rules per library. Each stage below says *why*, because several choices look like mistakes until
you know the reason.

### 1. R1 → gene assignment (bowtie2, transcriptome)

R1 is aligned to the transcriptome and the transcript ID is collapsed to its gene. No MAPQ or
multimapper filter is applied, deliberately: a read from a shared exon multimaps across isoforms of
its own gene, which is ambiguity that does not affect the gene call. Filtering it would
systematically deplete multi-isoform genes.

Align against the **full** transcriptome, not protein-coding only. The Snakefile's fallback default
is the protein-coding index; `config_template.yaml` overrides it. This matters: with rRNA present on
the R1 side, rRNA fragments get a gene assignment and their R2 mates land in the rDNA repeat arrays,
which is why R2 multimapper handling has such a large effect on the results.

### 2. R2 → adapter filter, UMI extraction, trimming (grep + fastp)

Reads without the adapter signature are dropped, the UMI moves to the header, 5 nt are trimmed, and
reads under 80 nt go.

> **The quality filter is disabled on purpose (`-Q`), and it is not a config key.** Illumina quality
> collapses immediately after a homopolymer: on the reference set the mean Phred of the 20 nt
> following a tail falls from Q36.4 (no tail) to Q16.9 (57–80 nt tail) through cluster dephasing,
> while the homopolymer itself stays at Q34–36. Measured with `-q 30 -u 10` on 94,198 adapter-bearing
> reads: **34.6% of all reads dropped, and 99.4–99.9% of reads carrying a tail ≥20 nt.** A quality
> filter removes precisely the molecules the experiment is about.

### 3. R2 → AT-masked genome (STAR)

The alignment that produces the measurement. Its 3' soft clip **is** the tail. See
[Building the AT-masked genome](#building-the-at-masked-genome) for why masking is what makes this
work.

Two parameter groups are worth understanding before changing them:

- **`outFilterScoreMinOverLread` / `MatchNminOverLread`** compare the aligned portion to the *full
  read length*, so they cap the largest measurable tail. On an 81 nt read: 0.66 (STAR's default) caps
  anchored tails at 28 nt; 0.50 caps them at 40. The reference set uses 0.50, which admitted 1,469
  reads with tails ≥30 nt and changed the tail length or class of **none** of the 304,592 reads that
  were already anchored. Going to 0.30 would raise the cap to 56 nt, but 93.78% of the extra reads are
  multimappers.
- **`outFilterMultimapNmax` / `outSAMmultNmax`** at `20 / 1` keep the best alignment of a multimapper
  so the read retains a real CIGAR, with `NH` recording how many loci it hit. At `1 / 1` any
  multimapper is discarded as unmapped and falls through to the sequence-only branch instead.

### 4. R2 → unmasked genome, end-to-end (bowtie2)

Contributes no tail and no coordinate. It is a single-purpose check: if the whole read *including*
the putative tail matched the genome, the tail may be templated. Read the resulting ratio the right
way round — details in [the column reference](docs/OUTPUT_COLUMNS.md#8-independent-templated-tail-check-bowtie2-unmasked-genome-end-to-end).

### 5. UMI deduplication

Rows are deduplicated on gene + R1 transcript position + UMI, so **one output row is one molecule**.
`R_output/<stem>_deduplication.csv` records the before and after counts; their ratio is the first
number to look at in QC.

### 6. Tail calling — two branches

| Branch | How the tail is obtained | Evidence |
|---|---|---|
| `cigar` | 3' soft clip of a confirmed, nearby alignment | strong; anchored and length-capped |
| `grep` | homopolymer prefix matching on the sequence alone | weak; needed for coverage, not length-capped |

A read goes to `grep` when it has no usable alignment — unmapped, on the wrong chromosome, ≥50 kb
from its gene, or the gene is unannotated. `grep_route_reason` records which.

`tail_route` keeps the two separable so you can pool them or not. Pooling is the normal thing to do.

The one case where the branches are not equivalent is the `oU` (oligo-U) class, which the grep
branch essentially cannot produce — 550,807 cigar-route reads against 14 grep-route across the
twelve reference libraries. Branch mix is therefore worth a glance in QC before pooling `oU`
specifically: if it tracks the grouping you are comparing, that class will too. Within
protein-coding genes the ratio is comparatively stable, so this is mostly a concern for
whole-transcriptome or biotype-split comparisons.

### 7. Routing flags

`loc_status` answers *"did R2 land where R1's gene says it should?"* — a gene-assignment verdict. It
must never be read as *"does this read have an alignment at all?"*, which is `r2_genomic_mapped`.
Conflating them hands cleanly-aligned reads to the sequence-only peeler, which then calls a genomic
A-stretch a tail.

Two columns measure that: `flag_grep_misrouted` (routed to grep despite having an alignment) and
`flag_genomic_tail` (the alignment positively contradicts the tail — zero 3' soft clip). **The
pipeline flags and does not drop.** You choose. Filters are in
[the column reference](docs/OUTPUT_COLUMNS.md#6-routing-and-the-two-branches).

### 8. Classification

`type` and `class` are assigned by `R_scripts/assign_type_and_class.R`. Class vocabulary and the full
derivation are in [`docs/TAIL_CATEGORIES.md`](docs/TAIL_CATEGORIES.md).

Both length floors — grep-branch and cigar-branch — are applied **downstream, not here.** Tails are
present in the output down to 1 nt, so changing your mind about a minimum length costs a filter, not
a re-run. [What is baked in and what is not](docs/OUTPUT_COLUMNS.md#thresholds-what-is-baked-in-and-what-is-not).

## Output

Per library, in `R_output/`:

| File | What |
|---|---|
| `<stem>_final_small.rds` | 46 columns, one row per molecule — the analysis table |
| `<stem>_final_big.rds` | 51 columns — the same plus R2's sequence, the tail string and R1's alignment detail |
| `<stem>_deduplication.csv` | read counts before and after UMI deduplication |

Every column is documented in **[`docs/OUTPUT_COLUMNS.md`](docs/OUTPUT_COLUMNS.md)**.

```r
d <- readRDS("R_output/<stem>_final_small.rds")

# tailed fraction per gene, anchored branch only, contradicted tails excluded
library(dplyr)
d %>%
  filter(tail_route == "cigar", !flag_genomic_tail) %>%
  group_by(geneR1) %>%
  summarise(molecules = n(),
            tailed    = sum(tail_len > 0),
            frac      = tailed / molecules) %>%
  filter(molecules >= 50) %>%
  arrange(desc(frac))
```

> **The schema is applied by the workflow.** The R stage calls `retype_small_table()` and
> `write_small_atomic()` from `R_scripts/finalise_small_table.R`: fixed factor levels, integer typing,
> the `schema_version` attribute, and an atomic write that reads the file back and checks rows,
> columns and stamp before renaming it into place. Every run also appends a provenance row — rows,
> columns, schema, size, md5 — to `logs/pipeline_manifest.tsv`.
>
> Fixed levels are what let you `rbind` libraries without factor levels silently reordering, and
> without `chrR2 == seqnames` failing with *"level sets of factors are different"*. If you are
> handed tables from an unknown vintage, check before analysing:
>
> ```r
> stopifnot(identical(attr(d, "schema_version"), "2026-10-20-flag1"))
> ```
>
> The rule also self-checks before writing: `tail_route` must partition every row with no `NA`, and
> `softclip_3p_len` must equal `tail_len` on every cigar-route row. Either failing aborts the rule
> rather than producing a table that looks fine.
>
> `_final_big.rds` is **not** stamped — it is the pre-finalise superset, written before the typing.

## QC

`reports/run_report.qmd` renders a per-run report (requires Quarto, which is **not** in
`environment_pipeline.yml` — the pipeline runs without it):

```bash
quarto render reports/run_report.qmd
```

Four numbers to look at first:

1. **Dedup ratio** — from the `_deduplication.csv` files.
2. **R2 unique vs multi share** — `r2_status`.
3. **Branch mix** — the share of `tail_route == "grep"`.
4. **Flag rates** — `flag_grep_misrouted` and `flag_genomic_tail`.

Branch mix is a QC number rather than an analysis one. Worth knowing where it sits (54.4–81.2% grep
across the reference set, varying with mappability and therefore with genotype) so that a shift
between libraries is recognised as a library property rather than read as biology. Keeping the QC
tables split by branch is what makes that visible; the analysis tables can be pooled.

## Caveats and limitations

**The `oU` class is cigar-branch only.** Pooling branches is fine in general, but `oU` counts carry
the library's branch mix with them. See [§6](#6-tail-calling--two-branches).

**Gene-level TTS bounds.** The annotation table is one row per gene, so `start`/`end` are the union
extent across all isoforms and the annotated 3' end is the most distal transcript end of any biotype.
Measured against per-transcript coordinates, the median multi-isoform gene's *typical* isoform end
sits 648 nt proximal to that bound; 52.5% of multi-isoform genes are more than 500 nt proximal, and
the effect is ~22× worse for lncRNA genes (median 2,962 nt) than protein-coding ones (132 nt). Any
distance-to-TTS measure inherits this, biased one way, strongest in the richest-annotated genes.

**Annotation version mismatch.** The reference STAR index was built with GENCODE v47 junctions while
the R stage joins a v50 gene table — a three-release gap. The join fails for 3.17–6.24% of rows
(mean 4.05%). Note `loc_status == "no_annotation"` is narrower (0.47–0.87%) because it is only
assigned to reads that mapped; an unmapped read lacking annotation is recorded as `r2_unmapped`.
Build both from one release to avoid this.

**Short tails are the artefact-prone bin.** 2–3 nt tails are the least reliable. On NovaSeq, the
two-channel chemistry's dark-G state inflates `AC` over `AG` in poly-T contexts, a composition
artefact specific to short tails. Treat the 2–3 nt bin separately or exclude it.

**`templated_suspect` does not work on the grep branch.** A genuine poly(A) tail and a genomic
poly-A stretch are the same sequence. Informative on the anchored branch, ~85% false-positive rate on
the unanchored one.

**No isoform resolution.** By design — see
[the column reference](docs/OUTPUT_COLUMNS.md#1-read-identity-and-provenance).

**`_final_big.rds` carries no `schema_version`.** Only the small table is typed and stamped; see [Output](#output).

**Protocol-specific assumptions.** The UMI length, adapter signature, trim length and read-length
floor are hardcoded in the fastp rule, not configurable. The pipeline has only been exercised on
libraries from one protocol; the mate contract in particular is an assumption it cannot check.

## Worked example: the reference dataset

Numbers quoted throughout come from twelve HEK-type knockout libraries (TUT4, TUT7, TUT4/7, MOV10,
DIS3L2, wild-type; two replicates each), sequenced paired-end 2×81 nt on NovaSeq.

All twelve are summarised in `logs/rebuild/qc/summary.tsv`.

| | |
|---|---|
| Molecules per library after dedup | 2,428,569 – 12,970,298 |
| Runtime per library | 12–52 min, tracking depth |
| Full twelve-library run | ~4 h sequential |
| Peak RSS, R stage | 16.0 GiB worst case |
| Branch mix | 54.4–81.2% grep |
| R2 uniquely mapped | 19.1–34.1% |
| R2 multimapped | 4.2–27.7% |
| `flag_grep_misrouted` | 4.5–7.7% of all molecules |
| `flag_genomic_tail` | 1.7–3.2% of all molecules, i.e. 37.5–44.8% of the misrouted set |
| Max anchored tail length | 40 nt in all twelve (the 0.50 threshold cap) |

Use these as a sizing guide only. This is the one section where dataset-specific material belongs.

## Repository layout

| Path | What |
|---|---|
| `Snakefile` | the workflow: nine rules plus `all` |
| `config/samples.tsv` | one row per library — the single source of truth |
| `config/config_template.yaml` | reference paths and STAR parameters |
| `config/generated/<stem>.yaml` | the exact config used per library; provenance |
| `scripts/run_all_samples.sh` | the driver: all libraries, sequential, resumable |
| `scripts/mark_AT_20_genome.py` | AT-run finder for building the masked reference |
| `R_scripts/run_tail_analysis.R` | Snakemake `script:` target for the R stage |
| `R_scripts/run_AT_tail_analysis_df.R` | tail measurement, routing and flags |
| `R_scripts/assign_type_and_class.R` | `type` and `class` assignment |
| `R_scripts/finalise_small_table.R` | typing, fixed factor levels, `schema_version` stamp, atomic write — called by the workflow |
| `scripts/rebuild/rebuild_one.R` | re-applies the schema to one library from existing `R_input/`, without re-aligning |
| `reports/run_report.qmd` | QC report (needs Quarto) |
| `docs/OUTPUT_COLUMNS.md` | **every output column explained** |
| `docs/TAIL_CATEGORIES.md` | class vocabulary and derivation, cited to file:line |
| `environment_pipeline.yml` | **the pipeline's dependencies** — 13 packages, env `pipeline_z` |
| `logs/run_status.tsv` | per-library start/finish/exit; drives the driver's resume logic |

Analysis downstream of these tables — thresholds, pooling, biotype comparisons, reports — lives in a
separate repository. This one ends at the RDS files.
