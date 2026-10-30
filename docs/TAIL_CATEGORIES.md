# Tail categories in `R_output/<stem>_final_*.rds`

**Written 2026-09-30.** Every statement below is transcribed from the code, with file and line
references. Counts are from `T7_14_S326_MM3` (the 0.50 trial: `star_multimap_nmax: 20`,
`star_score_min_over_lread: 0.50`), 3,453,168 rows, unless another run is named.

---

## 1. The three alignments, and which one produces a tail

| # | Alignment | Index | Mode | Output | Role |
|---|---|---|---|---|---|
| 1 | **R1 → transcriptome** | `fasta/index_bowtie_full_transcriptome/full_trans_` — the **full** transcriptome, 670,670 sequences, all biotypes (653 rRNA/Mt_rRNA), **not** the 112,218-sequence `prot_coding_` index | bowtie2, end-to-end | `R_input/<s>_R1_aligned.tsv` | Identifies the **gene** (`geneR1`) and supplies the dedup key (`mapR1`, the transcript coordinate) |
| 2 | **R2 → AT-masked genome** | `fasta/index_masked_genome_star/` | **STAR**, soft-clipping | `R_input/<s>_R2_aligned_to_masked_genome.tsv` | **The only source of a tail measurement.** The soft clip *is* the tail |
| 3 | **R2 → unmasked genome** | `fasta/index_bowtie_genome/GRCh38_primary_genome` | bowtie2, **end-to-end** (default mode, no `--local`) | `R_input/<s>_R2_genome_mapped.tsv` | **Flag only.** Contributes no tail and no coordinate |

### Yes — the CIGAR tail comes only from the masked genome

`Snakefile:199-204` writes the R2 table from the **STAR masked-genome SAM** only:

```awk
$1 !~ /^@/ { nh="NA"; if (match($0, /NH:i:[0-9]+/))   nh=substr($0, RSTART+5, RLENGTH-5);
             ut="NA"; if (match($0, /uT:A:[0-9]/))    ut=substr($0, RSTART+5, 1);
             as="NA"; if (match($0, /AS:i:-?[0-9]+/)) as=substr($0, RSTART+5, RLENGTH-5);
             nm="NA"; if (match($0, /nM:i:[0-9]+/))   nm=substr($0, RSTART+5, RLENGTH-5);
             print $1, $2, $3, $4, $6, $10, nh, ut, as, nm }
```

Field `$6` (CIGAR) becomes `V5`, which `run_tail_analysis.R:32-35` renames to **`CIGAR_R2`**. The
soft clip is parsed from it at `run_AT_tail_analysis_df.R:83-84`:

```r
softclip_str[idx_rev_map] <- str_extract(cigar[idx_rev_map], "\\d+(?=S$)")   # flag 16
softclip_str[idx_fwd_map] <- str_extract(cigar[idx_fwd_map], "^\\d+(?=S)")   # flag 0
```

That is the only place a soft-clip length enters the pipeline. **No code path derives a tail from
alignment 3.**

### Yes — the unmasked alignment is purely a false-positive flag

`R_scripts/extract_read_tags.py` writes **three columns** — read id, `unique`/`multi`, and `NM`
(the edit distance, added 2026-09-30; it was two columns before, and the reader widens defensively
so older files still load with `NM = NA`):

```
A00805:339:HGMF5DRX5:1:1101:24343:1000:GAAATCGGAGAACCT	unique	0
A00805:339:HGMF5DRX5:1:1101:32045:1000:TAGGGAGAAATCGCA	multi	2
```

No coordinate and no CIGAR — so it still cannot contribute a tail. bowtie2 runs with `--no-unal`
(`Snakefile:145`), so **presence in the file means "aligned end-to-end to the unmasked genome"** and
absence means it did not.

`al_ete_bt` is built from it, unchanged, in `run_tail_analysis.R`:

```r
R2_2_unique <- R2_2 %>% filter(V2 == "unique")
grep_aligned <- dedup_dat %>%
    filter(tail_len < 40, type == "grep", read_name %in% R2_2_unique$V1) %>%
    pull(read_name)
dedup_dat$al_ete_bt <- ifelse(dedup_dat$read_name %in% grep_aligned, "yes", "no")
```

**It is a flag, not a filter** — no row is dropped or altered. The logic is sound: an end-to-end
alignment means the *whole* read, putative tail included, matched the genome, so a non-zero
`tail_len` on such a read is templated sequence rather than a 3′ addition.

**Two blind spots follow from the predicate** — both addressed by `bt_nm_over_tail` /
`templated_suspect` in §7b, which leave `al_ete_bt` itself untouched:

- It fires only for `type == "grep"`, so **no `cigar` read is ever flagged** — measured:
  `al_ete_bt == "yes"` is 66,466 reads, **all** of them `grep`, **zero** `cigar`.
- It fires only below `tail_len < 40`, so long unanchored tails get no check — even though a long
  A-run peeled into A-rich genomic sequence is exactly the failure mode the A-regex permits (§5).

### What is masked, and the one thing that is not recoverable

The STAR index was built from **`genome/GRCh38.primary_masked_AT_new.fa`**
(`fasta/index_masked_genome_star/genomeParameters.txt`). The name says **AT-masked**, which is the
mechanism the `cigar` branch depends on: if genomically encoded A-rich sequence were alignable, it
would be absorbed into the `M` block and the tail silently under-called.

> **Open item.** *How* the AT masking was produced — tool, window, threshold — is **not recoverable
> from this repository.** The source FASTA is not on disk and the index records only its filename.
> `README.md:58-59` already carries this as a TODO. Every `cigar` tail's interpretation rests on
> this, so it is worth writing down from memory or lab notes while it is still known.

> **Version mismatch, noted not fixed.** The STAR index splice junctions come from
> **GENCODE v47** (`sjdbGTFfile genome/gtf/gencode.v47.primary_assembly.basic.annotation.gtf`),
> while the gene coordinates used for the location check come from
> **`annotation/gene_annotations_v50.rds`**. A three-release gap between the junctions R2 is
> aligned against and the gene bounds `loc_status` is measured from.

---

## 2. `loc_status` — where R2 landed, and it decides the branch

Assigned at `run_AT_tail_analysis_df.R:57-71`. Categories are mutually exclusive and cover every
row. `stop_pos` is the gene's 3′ end — `end` on `+` strand, `start` on `−` (lines 36-37).

| `loc_status` | Condition (line) | Count | Branch taken |
|---|---|---|---|
| `confirmed` | same chromosome as `geneR1` **and** `abs(stop_pos − posR2) < 50000` (61, 68) | 909,121 | **anchored** |
| `multi_locus` | mapped **and** `NH > 1` — assigned **last**, overrides any verdict above (58, 71) | 439,069 | **anchored** |
| `r2_unmapped` | `CIGAR_R2 == "*"` or `chrR2 == "*"` or `NA` (57, 64) | 1,837,747 | fallback |
| `far_from_gene` | right chromosome, **≥ 50 kb** from the gene's 3′ end (67) | 147,382 | fallback |
| `diff_chrom` | mapped and annotated, but `chrR2 != seqnames` (66) | 89,804 | fallback |
| `no_annotation` | `geneR1` absent from the annotation table (59, 65) | 30,045 | fallback |

**The 50 kb window is the `far_from_gene` cutoff** (line 61). `NH` is `0`/`1`/`>1` from STAR, or
`NA` for tables written before 2026-09-25.

### The branch gate, and a trap worth knowing

```r
mapped_case <- loc_status %in% c("confirmed", "multi_locus")      # line 75
```

**`diff_chrom` and `far_from_gene` reads have a perfectly good STAR CIGAR, and the pipeline
deliberately refuses to use it.** They are routed to the unanchored scan. So:

> **anchored branch = `confirmed` ∪ `multi_locus`, and nothing else.**

Measured, `loc_status` × `type` — a clean partition, no overlap:

|  | `cigar` | `grep` | `no_tail` |
|---|---|---|---|
| `confirmed` | 283,299 | 0 | 625,822 |
| `multi_locus` | 41,972 | 0 | 397,097 |
| `diff_chrom` | 0 | 89,804 | 0 |
| `far_from_gene` | 0 | 147,382 | 0 |
| `no_annotation` | 0 | 30,045 | 0 |
| `r2_unmapped` | 0 | 1,837,747 | 0 |

`multi_locus` is kept **out** of `confirmed` on purpose: the soft-clip *offset* is stable across
paralogous copies, so the tail call survives, but which locus STAR reported is arbitrary among
equal-scoring ties (lines 51-56).

---

## 3. `type` — how the tail was called

From `tail_source` (`run_AT_tail_analysis_df.R:97-100`) via `assign_type_and_class.R:19-24`.
**This is a method label, not a biological one.**

| `type` | `tail_source` | Meaning | Count | Median `tail_len` |
|---|---|---|---|---|
| **`cigar`** | `cigar` | Anchored: soft clip > 0 on a `confirmed`/`multi_locus` alignment | 325,271 (9.4%) | **1** |
| **`no_tail`** | `full_match` | Anchored: alignment present, **no soft clip** → no tail | 1,022,919 (29.6%) | 0 |
| **`grep`** | `fallback` | Unanchored: whole oriented read scanned from position 1 | 2,104,978 (61.0%) | **57** |
| `dubious` | `NA` | Empty/invalid sequence | **0** | — |

The `cigar` median of 1 against the `grep` median of 57 is **not** biology — see §6.

**Orientation** (lines 29-32), which is why peeling from position 1 is legitimate:

```r
is_rev <- !is.na(flag) & flag == 16
seq_oriented <- ifelse(is_rev, stri_reverse(seq),
                                chartr("ATCGatcg","TAGCtagc", seq))
```

Reverse-flag reads are **reversed**; forward-flag reads are **complemented**. Either way the 3′ end
of the original RNA ends up at index 1.

---

## 4. `class` — what the tail is made of

The peel runs **T → G → C → A**, each from the current start
(`run_AT_tail_analysis_df.R:123-147`):

```r
t_len: ^T+      g_len: ^G+      c_len: ^C+      a_match: ^A+(?:[^A]A{5,})*
```

`A_count = str_count(a_match, "A")`; `T_count`, `G_count`, `C_count` are the run lengths.
Classified at `assign_type_and_class.R:25-43`:

| `class` | Condition | Reading | Count | % |
|---|---|---|---|---|
| `none` | `tail_len == 0` or all four counts 0 | No detected 3′ addition | 1,022,930 | 29.62 |
| `A` | A only | Canonical poly(A) | 2,207,553 | 63.93 |
| `AU` | A **and** T, no G/C | **Uridylated poly(A)** — the TUT4/7 product | 61,683 | 1.79 |
| `oU` | T only, no A/G/C | **Uridylation of an A-less 3′ end** | 45,647 | 1.32 |
| `mixed` | `(T_count + G_count + C_count) > 1` after the above fail | ≥2 non-A species | 45,374 | 1.31 |
| `AC` | A **and** C, no T/G | Cytidylated poly(A) | 25,937 | 0.75 |
| `AG` | A **and** G, no T/C | Guanylated poly(A) | 20,347 | 0.59 |
| `C` | C only | **At the error floor — do not interpret** (§7) | 13,796 | 0.40 |
| `G` | G only | **At the error floor — do not interpret** (§7) | 9,901 | 0.29 |
| `check` | fallback | Should be unreachable | **0** | — |

### `oU` is structurally reachable only from the anchored branch

Measured: **45,647 `oU` reads, every one of them `cigar`. Zero from `grep`.** This is not a
sampling accident, it is forced by the code:

- For `cigar`, `work_seq` is **truncated to the soft clip** (line 107). If the whole clip is T's,
  the peel yields `T_count > 0` and nothing else → `oU`.
- For `grep`, `work_seq` is the **whole oriented read** (line 115). After the T-run the peel
  continues into genomic sequence, which essentially always contributes a G, C or A — so
  `A_count == 0 & G_count == 0 & C_count == 0` is unreachable in practice.

**Consequence: uridylation of a deadenylated end can only be seen on reads that aligned to the
masked genome.** Anything that moves reads from `grep` into `cigar` raises measured `oU`.

> **This makes `oU` as a share of all rows a confounded statistic — do not use it.** Across the
> twelve samples, relaxing the threshold raised `oU` as % of the table by 49–79%, but on the
> **anchored branch** the rise was much smaller, and in the **high-confidence set** (anchored +
> `confirmed` + `NH == 1`) the U≥2 rate moved by only ±0.3 percentage points. Almost all of the
> apparent increase was the anchored branch getting bigger, not more uridylation being detected.
> Always normalise on `type == "cigar"` at minimum, and preferably on the strict set.

---

## 5. `tail_len` is a **span**, not a base count

`A_count + T_count + G_count + C_count ≤ tail_len`, and the two differ for a **different reason in
each branch**. This trips people up, so it is worth stating plainly:

| Branch | `tail_len` is | Line | Rows where `tail_len ≠ Σcounts` |
|---|---|---|---|
| `cigar` | the **soft-clip length**, whatever it contains | 109 | 22,174 (**6.8%** of `cigar`); excess median 2, max 39 |
| `grep` | `a_len + t_len + g_len + c_len`, where `a_len = nchar(a_match)` — **the whole A-regex span** | 150-151 | 351,530 (**16.7%** of `grep`) |
| `no_tail` | 0 | 103 | 0 (100% agree) |

The `grep` gap is the A-regex: `^A+(?:[^A]A{5,})*` matches across interruptions, so
`a_len` counts the interrupting non-A bases while `A_count` does not.

The `cigar` gap means part of the clipped segment was not accounted for by any of the four peels.

**So `tail_len` answers "how many bases are non-templated", and `A_count` answers "how many of them
are A". Do not use one as a proxy for the other.**

---

## 6. The measurement ceiling — the single most important caveat

`outFilterScoreMinOverLread` / `outFilterMatchNminOverLread` compare the aligned portion to the
**full read length**. A soft-clipped tail adds nothing to the numerator but stays in the
denominator, so for an 81 nt read (fastp `-l 80`) the largest reportable anchored tail is

```
ceiling = read_len × (1 − threshold)
```

**Measured, exactly as predicted, on the same sample and FASTQ:**

| threshold | predicted ceiling | **measured max `tail_len` (`cigar`)** |
|---|---|---|
| 0.66 (STAR default) | 27.5 | **28** |
| 0.50 | 40.5 | **40** |
| 0.30 | 56.7 | **56** |

The distribution does not taper — it stops dead. At 0.66 there were 146 reads at 28 nt and **zero**
at 29+.

**Where the unanchored tails actually sit** — `grep`-branch `tail_len` at 0.66, n = 2,127,432,
**median 57**:

| `tail_len` | n | % | Reachable by the anchored branch? |
|---|---|---|---|
| 0–4 | 238,022 | 11.19 | already |
| 5–27 | 164,603 | 7.74 | already at 0.66 |
| 28–40 | 301,840 | 14.19 | needs 0.50 |
| 41–56 | 357,093 | 16.79 | needs 0.30 |
| **57–80** | **451,789** | **21.24** | **no threshold reaches these** |
| **81+** | **614,085** | **28.87** | **right-censored at read length** |

**Half of all unanchored tails (50.1%) are ≥ 57 nt and structurally unreachable at any threshold
setting, and 28.87% are jammed against the 81 nt read-length wall.** Tail length in this assay is
limited by read length far more than by any aligner parameter. Means and quantiles above ~60 are
not tail lengths.

---

## 6a. WHY long tails cannot be anchored — a double penalty, measured

Read length is only half the story. The other half is **sequencing quality after the
homopolymer**, and it is the larger effect.

Raw R2 layout is `[15 nt UMI][GTCAG][81 nt insert]`, and the tail sits at the insert's **5′
start** (it is sequenced as poly-T and complemented downstream). So the genomic anchor is
whatever follows the homopolymer. Measured over 2,000,000 reads of `T7_14_S326`:

| leading T-run (= tail) | reads | mean Q **of the homopolymer** | mean Q **of the 20 nt anchor after it** | anchor nt remaining |
|---|---|---|---|---|
| 0 | 49.97% | — | **36.4** | 81.0 |
| 1–4 | 17.30% | 36.2 | **36.2** | 79.7 |
| 5–9 | 0.86% | 34.0 | **33.7** | 74.5 |
| 10–19 | 1.44% | 34.6 | **30.5** | 66.2 |
| 20–27 | 2.01% | 35.4 | **24.7** | 57.3 |
| 28–40 | 4.88% | 35.2 | **20.8** | 46.8 |
| 41–56 | 6.10% | 34.9 | **18.9** | 32.4 |
| 57–80 | 7.28% | 34.4 | **16.9** | 13.5 |
| 81 (all T) | 10.16% | 35.2 | — | 0.0 |

**The homopolymer itself stays at Q34–36 throughout. Only what follows it collapses, Q36.4 →
Q16.9** — from ~0.023% to ~2% expected error per base, roughly **80×**. That is the signature of
Illumina cluster dephasing: while every molecule is on T the run calls cleanly, but once
synchrony is lost the heterogeneous bases afterwards receive mixed signal.

So a long-tailed read is penalised twice: **fewer anchor bases** (81 → 13.5) **and those bases far
less reliable**.

### This is visible inside the tables, and it explains the filter behaviour

Anchored (`cigar`) reads, `T7_14_S326` at 0.50 — mismatch count rises **12×** with tail length:

| `tail_len` | n | mean `nM` | % with `nM ≥ 2` | mean `score_frac` |
|---|---|---|---|---|
| 1–4 | 233,081 | 0.226 | 3.96 | 0.950 |
| 5–9 | 22,522 | 0.692 | 15.19 | 0.864 |
| 10–19 | 47,320 | 1.461 | 41.09 | 0.747 |
| 20–27 | 18,010 | 2.709 | **78.36** | 0.625 |
| 28–40 | 4,340 | 2.679 | **76.59** | 0.537 |

**Both STAR filters therefore penalise long tails, not just the score filter.**
`outFilterMismatchNoverReadLmax 0.05` allows ~4 mismatches in 81 nt, and by `tail_len` 20–40
around 77% of reads already carry ≥2 — they are approaching that wall as well.

It also predicts the observed cost of relaxing the threshold. STAR's mismatch rate per base was
**0.28% at 0.66** and **1.17% at 0.30**; the anchors admitted at 0.30 sit at Q18.9–20.8, i.e.
0.8–1.3% expected error. Predicted and observed agree. Relaxing the threshold does not merely
admit shorter anchors — it admits **dephased** ones.

> **The consequence for interpretation: measured tail length and location confidence are in
> direct conflict, and no aligner parameter resolves it.** The information was not recorded by
> the sequencer. Longer reads, or a protocol that does not read through the homopolymer before
> the anchor, is the only real fix.

**One practical note, and it is not optional.** fastp runs with **`-Q`** — quality filtering
**off** — hardcoded in the Snakefile since 2026-09-30, so these low-quality anchors are retained
rather than discarded. It has to be off: a quality filter removes precisely the reads this assay
is about. Measured with fastp's own `-q 30 -u 10` on 94,198 adapter-bearing reads of
`T7_14_S326`, against the leading T-run (= the tail):

| leading T-run | would be **discarded** |
|---|---|
| 0 | 12.2% |
| 1–4 | 12.0% |
| 5–9 | 45.1% |
| 10–19 | 87.4% |
| **20–27** | **99.4%** |
| **28–40** | **99.9%** |
| **41–56** | **99.9%** |
| **57–80** | **96.7%** |

34.6% of all reads, and **99.4–99.9% of everything with a tail ≥ 20 nt**. Before 2026-09-30 the
`enable_quality_filter` key defaulted to `True`, so a run without the config override would have
silently done this; the key is gone and `-Q` is no longer tunable.

Quality therefore has to be handled at analysis time instead, which is what `nM`, `score_frac`
and `match_frac` are for.
That is defensible — dropping them would discard precisely the long-tailed reads of interest —
but it means quality has to be handled at analysis time, which is what `nM`, `score_frac` and
`match_frac` are for.

---

## 7. Modification identity is far more robust than length

### Identity survives anchoring; length does not

From the 0.30 trial, 373,352 reads have **both** an unanchored and an anchored call, so the
unanchored call can be scored against anchored ground truth:

| `class` (unanchored) | n | % surviving anchoring | `confirmed` only |
|---|---|---|---|
| `A` | 359,078 | **100.0** | 100.0 |
| `C` | 376 | 100.0 | 100.0 |
| `G` | 186 | 100.0 | 100.0 |
| `AG` | 2,246 | 98.4 | 98.9 |
| `AU` | 5,762 | 96.3 | 98.7 |
| `AC` | 3,837 | 94.6 | 95.6 |
| **`mixed`** | 1,867 | **77.1** | 81.7 |

Tail **length** over the same reads was over-called by a median of 18 nt.

**Mechanism:** the peel reads T → G → C → A *from position 1*, so the non-A modification is read
**before** the A-run can over-extend. Over-run corrupts the A count, not the identity of the
terminal non-A bases. Every misclassification loses the A component (`AU` → `oU` 3.7%,
`AC` → `C` 5.4%, `mixed` → `oU` 22.5%).

> **"Is there a U?" is reliable. "Is the U on a poly(A) tail?" is less so. `mixed` is the one
> genuinely unstable class — re-derive it from the counts.**

### Mono-G and mono-C are a free internal error floor

If a terminal non-A call were pure sequencing error, U, G and C would appear at similar rates. They
do not. Anchored reads at 0.66, n = 678,555:

| base | any length | exactly 1 nt | **≥ 2 nt** | ≥ 3 nt |
|---|---|---|---|---|
| **U** | 70,655 (10.41%) | 51,297 (7.56%) | **19,358 (2.85%)** | 7,841 (1.16%) |
| G | 7,117 (1.05%) | 6,350 (0.94%) | 767 (0.11%) | 176 (0.03%) |
| C | 9,391 (1.38%) | 8,408 (1.24%) | 983 (0.14%) | 288 (0.04%) |

Signal-to-floor (floor = mean of G and C): any length **8.6×**, mono **7.0×**, **≥2 nt 22.1×**.
Requiring ≥2 nt roughly triples specificity. Strictest set — anchored + `confirmed` + `NH == 1`:
U ≥2 nt **38.4×** above the floor, residual contamination ≈ 2.6%.

---

## 7a. `r2_status` — why R2 is, or is not, placed

`loc_status` collapses every kind of alignment failure into `r2_unmapped`. STAR records the actual
reason on each unmapped record as `uT:A:<n>`, captured as column 8 since 2026-09-30 and exposed as
`r2_status`. Measured on `T7_14_S326` at 0.50, and it reconciles exactly to the 3,453,168 rows:

| `r2_status` | `uT` | n | % | meaning |
|---|---|---|---|---|
| `unique` | — | 1,176,359 | 34.07 | placed at one locus (`NH == 1`) |
| `multi` | — | 439,072 | 12.71 | placed at several (`NH > 1`) |
| **`dropped_short`** | **1** | **1,183,135** | **34.26** | **an alignment was found, then rejected for short aligned length** |
| `no_seed` | 0 | 654,448 | 18.95 | no acceptable seed: no alignment at all |
| `dropped_multi` | 3 | 154 | 0.00 | more loci than `star_multimap_nmax` allows |
| `dropped_mismatch` | 2 | 0 | — | too many mismatches |

It is **orthogonal to `loc_status` and `type`** — it reports the aligner's verdict, not the
location check. The mapping is clean: `multi` → all `multi_locus`; `unique` → `confirmed` +
`diff_chrom` + `far_from_gene` + `no_annotation`; all three `dropped_*`/`no_seed` → `r2_unmapped`.

**So `r2_status` splits the previously opaque 1,837,747-read `r2_unmapped` bucket into
`dropped_short` (64.4%) and `no_seed` (35.6%).** These are very different populations:

| `r2_status` (within `grep`) | n | median `tail_len` | % ≥ 57 | % ≥ 81 |
|---|---|---|---|---|
| `dropped_short` | 1,183,135 | **48** | 33.57 | 2.11 |
| `no_seed` | 654,448 | **81** | **99.80** | **89.39** |

> **`no_seed` reads have a median `tail_len` of 81 — the entire read — and 89.4% are at the
> 81 nt wall.** These are reads with no genomic anchor whatsoever (10.16% of raw reads are
> *entirely* poly-T). Their `tail_len` is a **censored lower bound, not a measurement**, and they
> carry no location information at all. They are the bulk of the `81+` bin in §6 and should be
> excluded from any tail-length distribution.
>
> `dropped_short` is the genuinely informative bucket: a real anchor, a real tail of ~48 nt
> median, rejected for aligned length. **This is the category that is causally entangled with the
> measurement** — a read is short-aligned largely *because* much of it is non-templated.

### Reconstructing a stricter threshold without re-aligning

`AS` and `nM` are carried through, so the fractions STAR itself thresholds on are recomputed as
`score_frac = AS / read_len` and `match_frac = (ΣM − nM) / read_len`, and exposed as
`passes_030`, `passes_040`, `passes_050`, `passes_066`.

**Validated against ground truth:** `passes_066` computed from the 0.50 run, compared with whether
STAR actually placed each read in the independent 0.66 run — **99.833% agreement, with zero false
positives** (5,768 false negatives in 1,554,383). Within the 0.50 run, `passes_050` is TRUE for
99.36% of placed reads and **0.000%** of `dropped_short`.

**Caveat:** the filters also gate which loci survive, so under a looser setting a read can return
at a higher `NH` or at a different locus among equal-scoring ties. Reliable for *"would this read
have been reported at all"*; weaker for *"where"*.

---

## 7b. `templated_suspect` and `bt_nm_over_tail` — using the unmasked alignment properly

Mapping end-to-end is **not** on its own evidence that a tail is templated: bowtie2 will force a
read through by scoring every tail base as a mismatch. Measured on the anchored branch,
`bt_ete_NM / tail_len` has **median exactly 1.00** — the typical end-to-end alignment mismatches
through precisely the tail, which is the signature of a **genuine** non-templated addition.

A tail is suspect only when that ratio is **low**, meaning most of the tail actually matched the
genome. `bt_nm_over_tail` is emitted as a continuous ratio so the cut is an analysis choice;
`templated_suspect` is the convenience boolean at ≤ 0.5.

| `tail_len` (anchored) | n | % `templated_suspect` | median `bt_nm_over_tail` |
|---|---|---|---|
| 1–4 | 233,081 | 2.64 | 1.00 |
| 5–9 | 22,522 | 5.54 | 0.83 |
| 10–19 | 47,320 | 8.72 | 0.82 |
| 20–27 | 18,010 | 11.51 | 0.37 |
| 28–40 | 4,340 | 9.17 | 0.25 |

Reassuringly, the most damning combination — a long anchored tail that *also* aligns end-to-end at
`NM == 0` — is essentially empty: **1 read out of 22,350** with `tail_len ≥ 20`. The anchored tail
calls are not templated artefacts.

> **The confound, and why this is a flag and never a filter.** A genuine poly(A) tail and a
> genomic poly-A stretch are *the same sequence*. An all-A read aligns end-to-end at `NM 0` to any
> long genomic A-run whether its tail is real or the product of internal priming. On the
> unanchored branch, where the peel can consume the whole read, the flag therefore cannot
> discriminate and fires on ~85% of end-to-end mappers (753,756 `grep` reads). **Informative on
> the anchored branch; suggestive at best on the unanchored one.**

---

## 8. Recommended filters

```r
d <- readRDS("R_output/<stem>_final_small.rds")

# --- anchored vs not ---------------------------------------------------------
anchored <- subset(d, type %in% c("cigar", "no_tail"))   # == confirmed | multi_locus
conf     <- subset(d, loc_status == "confirmed")         # anchored AND location verified
ambig    <- subset(d, !is.na(NH) & NH > 1)               # placement arbitrary
severe   <- subset(d, !is.na(NH) & NH > 5)               # heavily repetitive

# --- high-confidence uridylation (~97% specific) -----------------------------
uri <- subset(d, type == "cigar" & loc_status == "confirmed" &
                 !is.na(NH) & NH == 1 & T_count >= 2)

# --- mono-U: usable, but subtract the floor ----------------------------------
uri_mono   <- subset(d, type == "cigar" & T_count == 1)
floor_rate <- mean(c(mean(d$G_count >= 2), mean(d$C_count >= 2)))   # free negative control

# --- templated-sequence suspicion (the unmasked end-to-end check) ------------
#  al_ete_bt == "yes" only ever marks grep reads with tail_len < 40 (kept unchanged).
#  Prefer the ratio: low bt_nm_over_tail means the "tail" matched the genome.
table(d$al_ete_bt, d$type)
susp <- subset(d, type == "cigar" & templated_suspect)         # ~4% of anchored reads

# --- why is R2 unplaced, and would a stricter threshold have kept it? --------
table(d$r2_status, d$type)
as_if_066 <- subset(d, passes_066)        # reproduces the 0.66 result, no re-alignment

# --- EXCLUDE the uninformative bucket before any tail-length summary ---------
#  no_seed reads have no anchor at all and median tail_len 81 -- a censored bound.
len_ok <- subset(d, r2_status != "no_seed")

# --- why is R2 on another chromosome? ----------------------------------------
d$chrom_mismatch <- !is.na(d$seqnames) & d$chrR2 != d$seqnames
table(d$chrom_mismatch, d$loc_status)
```

| Call | Verdict |
|---|---|
| U ≥2 nt, `cigar`, `confirmed`, `NH == 1` | **Trust** — ~2.6% residual |
| U ≥2 nt, `cigar`, any location | Trust; ~4.5% floor |
| mono-U, `cigar` | Usable; subtract the G/C floor |
| mono-G / mono-C | **Do not interpret** — at the error floor |
| `mixed` | **Re-derive** from counts; ~23% reclassifies |
| any `class` where `type == "grep"` | Identity OK; **`A_count` and `tail_len` unreliable** |
| `oU` | Only ever from `cigar` (§4) — absence in `grep` is structural, not biological |
| `tail_len` ≥ 57 | Unreachable by the anchored branch at any threshold (§6) |
| `tail_len` ≥ 81 | **Right-censored** — a read-length artefact, not a length |

---

## 9. Column reference

`run_tail_analysis.R:75` drops five columns from the small table:

```r
cols_to_remove <- c("seq_R2", "transR1", "CIGAR_R1", "tail", "flagR1")
```

| Column | small | big | Source |
|---|---|---|---|
| `read_name`, `geneR1`, `mapR1` | ✓ | ✓ | R1 alignment; `mapR1` is the dedup coordinate |
| `flagR2`, `chrR2`, `posR2`, `CIGAR_R2` | ✓ | ✓ | STAR masked-genome alignment |
| `UMI` | ✓ | ✓ | 15 nt, from the read header |
| `NH` | ✓ | ✓ | STAR: loci placed. `0`/`1`/`>1`; `NA` if pre-2026-09-25 |
| `seqnames`, `start`, `end`, `strand` | ✓ | ✓ | Gene bounds, GENCODE v50 |
| `A_count`, `T_count`, `G_count`, `C_count` | ✓ | ✓ | Peel counts |
| `tail_len` | ✓ | ✓ | **Span**, see §5 |
| `high_AT_mix` | ✓ | ✓ | AT content of R2 > 93% |
| `tail_source`, `type`, `class`, `loc_status` | ✓ | ✓ | §2–§4 |
| `al_ete_bt` | ✓ | ✓ | Unmasked end-to-end flag, §1 — unchanged, narrow predicate |
| `r2_status` | ✓ | ✓ | Why R2 is/is not placed, from `uT`. §7a |
| `uT`, `AS`, `nM` | ✓ | ✓ | Raw STAR tags behind `r2_status` and the `passes_*` flags |
| `read_len`, `score_frac`, `match_frac` | ✓ | ✓ | The fractions STAR thresholds on. §7a |
| `passes_030/040/050/066` | ✓ | ✓ | Would this read pass a stricter threshold. §7a |
| `bt_ete_mapped`, `bt_ete_unique`, `bt_ete_NM` | ✓ | ✓ | Unmasked end-to-end result, per read. §7b |
| `bt_nm_over_tail`, `templated_suspect` | ✓ | ✓ | Templated-tail suspicion. §7b |
| **`tail`** | — | ✓ | The tail **sequence**. Needed to re-derive a class |
| **`seq_R2`** | — | ✓ | Full R2 read |
| `transR1`, `CIGAR_R1`, `flagR1` | — | ✓ | R1 detail |

**To re-derive `class` from scratch you need the big table** (`tail`), not the small one.
Base qualities are not carried in either.

---

## 10. Known limits

| Limit | Evidence |
|---|---|
| **Anchor quality collapses after the homopolymer, Q36.4 → Q16.9 (~80× error rate)** | §6a — the dominant reason long tails cannot be located, and not fixable by any parameter |
| **`no_seed` reads (18.95%) have median `tail_len` 81 and no anchor at all** | §7a — a censored bound, not a measurement. Exclude from length summaries |
| `tail_len` right-censored at read length (81 nt) | 28.87% of `grep` reads sit in the `81+` bin |
| Anchored ceiling = `read_len × (1 − threshold)` | Measured 28 / 40 / 56 at 0.66 / 0.50 / 0.30 |
| 50.1% of unanchored tails (≥57 nt) are unreachable at any threshold | §6 |
| A 1–4 nt tail can be absorbed as mismatches rather than soft-clipped | `outFilterMismatchNoverReadLmax 0.05` ≈ 4 mismatches per 81 nt. Short tails stay ambiguous at **every** threshold; long tails cannot be absorbed this way |
| `^A+(?:[^A]A{5,})*` tolerates interruptions followed by ≥5 A, so it can run into A-rich genomic sequence | 61.0% of rows are `grep` |
| `oU` is unreachable from `grep` | §4 — structural, and it biases any `grep`-vs-`cigar` comparison of uridylation |
| `no_annotation` is silent | 30,045 rows (0.87%), no warning emitted |
| `posR2` is raw SAM `POS`, so its meaning flips with strand | `PIPELINE_EVALUATION.md` finding 5b |
| `al_ete_bt` never marks a `cigar` read | §1 — measured 0 |
| AT-masking procedure unrecorded | §1 |
| STAR junctions GENCODE v47 vs gene bounds v50 | §1 |
