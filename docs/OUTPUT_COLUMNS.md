# Output columns

Reference for the two per-read tables the pipeline writes per library:

| File | Rows | Columns | Contents |
|---|---|---|---|
| `R_output/<stem>_final_small.rds` | one per deduplicated read | **46** | the analysis table |
| `R_output/<stem>_final_big.rds` | identical rows | **51** | the same 46 plus five raw columns (below) |

Both are plain `data.frame`s; load with `readRDS()`. `_final_small.rds` carries
`attr(df, "schema_version")` when it has been through `R_scripts/finalise_small_table.R`
(see [Typing and the schema stamp](#typing-and-the-schema-stamp)).

One row is one *molecule*, not one sequenced read: UMI deduplication happens before these
tables are written, so counting rows counts molecules.

---

## 1. Read identity and provenance

| Column | Type | Meaning |
|---|---|---|
| `read_name` | character | Illumina read header with the 15 nt UMI suffix stripped. Joins R1 and R2 |
| `UMI` | character | The 15 nt UMI read off the header |
| `geneR1` | factor | **Gene assignment**, Ensembl gene ID, version suffix stripped. Comes from R1's transcriptome alignment, collapsed transcript → gene |
| `mapR1` | integer | R1's alignment position in **transcript** coordinates of whichever transcript bowtie2 reported. Used only as part of the deduplication key; it is not a genomic coordinate |

**Isoforms are not resolved.** R1 reads from a shared exon are compatible with every isoform of
their gene; bowtie2 reports one arbitrarily among equal-scoring ties and the pipeline keeps only the
gene. There is no isoform-level readout, by design — 3'-RACE R1 carries no information that could
break the tie.

## 2. Gene annotation (joined from the annotation RDS on `geneR1`)

| Column | Type | Meaning |
|---|---|---|
| `seqnames` | factor | Gene's chromosome |
| `start`, `end` | integer | Gene's outermost bounds, i.e. the union extent across **all** its annotated transcripts |
| `strand` | factor | `+` / `-` |

`NA` in all four means the gene ID was absent from the annotation table. On the reference set that
is 3–6% of rows.

**`start`/`end` are gene-level, not transcript-level.** For a `+` gene the annotated 3' end is
`end`, which is the most distal transcript end of any isoform or biotype. A read terminating at a
genuine proximal polyadenylation site therefore scores as far upstream. If you compute distance to
the TTS, this is the dominant source of error, and it is worse for lncRNA genes than coding ones.

## 3. R2 alignment to the AT-masked genome (STAR)

R2 is the mate that carries the tail. These are read straight from the SAM.

| Column | Type | Meaning |
|---|---|---|
| `flagR2` | integer | SAM flag. `16` = reverse strand; this is what orients the tail |
| `chrR2` | factor | R2's chromosome, or `*` when unmapped |
| `posR2` | integer | First aligned base, 1-based. `0` when unmapped |
| `CIGAR_R2` | character | R2's CIGAR, or `*`. The 3' soft clip here **is** the tail measurement |
| `NH` | integer | Number of loci STAR placed the read at. `>1` = multimapper |
| `uT` | factor | STAR's reason code for an *unmapped* record, so `NA` whenever the read aligned. `0` no alignment, `1` too short, `3` multimapping beyond the limit |
| `AS` | integer | Alignment score |
| `nM` | integer | Mismatch count |

## 4. Recomputable acceptance thresholds

STAR's acceptance tests compare the aligned portion to the **full read length**, so they cap the
largest reportable tail. These columns let you re-apply any stricter threshold without re-aligning.

| Column | Type | Meaning |
|---|---|---|
| `read_len` | integer | R2 length after fastp. Fixed by the protocol (81 nt on the reference set) |
| `score_frac` | numeric | `AS / read_len` |
| `match_frac` | numeric | matched bases / `read_len` |
| `passes_030`, `passes_040`, `passes_050`, `passes_066` | logical | Would this read survive `outFilterScoreMinOverLread` / `MatchNminOverLread` at 0.30 / 0.40 / 0.50 / 0.66 |

The run threshold is in `config/config_template.yaml`. At 0.50 on an 81 nt read the largest
anchored tail that can be reported is 40 nt; at STAR's default 0.66 it is 28. **`passes_*` can only
be used to tighten, never to loosen** — reads rejected by the threshold actually used were never
written with a CIGAR in the first place.

## 5. Tail measurement

| Column | Type | Meaning |
|---|---|---|
| `tail_len` | integer | Tail length in nt. `0` = no tail |
| `A_count`, `T_count`, `G_count`, `C_count` | integer | Base composition of the tail. `T_count` is uridylation — the table is in DNA alphabet |
| `high_AT_mix` | logical | R2 is >93% A+T over its whole length. A low-complexity warning, not a tail property |
| `tail_source` | factor | How the tail was obtained: `cigar` (3' soft clip of a real alignment), `full_match` (aligned, no soft clip, so no tail), `fallback` (sequence-only scan) |

## 6. Routing and the two branches

This is the most important section for interpreting anything.

| Column | Type | Meaning |
|---|---|---|
| `tail_route` | factor | `cigar` or `grep` — which branch produced `tail_len`. Partitions every row, never `NA` |
| `r2_genomic_mapped` | logical | Does this read have an R2 alignment at all |
| `softclip_3p_len` | integer | 3' soft clip measured for **every** row with a real CIGAR, independent of routing. `NA` when there is no CIGAR; `0` is meaningful and means the aligner clipped nothing |
| `grep_route_reason` | factor | Why a read went to the grep branch: `diff_chrom`, `far_from_gene`, `no_annotation`, `r2_unmapped`, `no_seed`, `dropped_short`, `dropped_multi`, `dropped_mismatch`, `dropped_mate`. `NA` for cigar-route reads |
| `flag_grep_misrouted` | logical | Sent to the grep branch **despite** having an alignment |
| `flag_genomic_tail` | logical | The above **and** `softclip_3p_len == 0` — the alignment positively contradicts the tail |

`flag_genomic_tail` is a strict subset of `flag_grep_misrouted`. **The pipeline flags and never
drops**; which rows to exclude is your decision. The one-line filters:

```r
d <- readRDS("R_output/<stem>_final_small.rds")
clean   <- subset(d, !flag_genomic_tail)              # drop contradicted tails only
stricter <- subset(d, !flag_grep_misrouted)           # drop everything misrouted
anchored <- subset(d, tail_route == "cigar")          # strongest evidence only
```

### What differs between the branches

The `cigar` branch anchors the tail to an alignment and is length-capped by the read geometry (§4).
The `grep` branch reads the tail off the sequence for reads with no usable alignment: weaker
evidence, and not length-capped, so it can consume the whole read. `tail_route` keeps them
separable; pooling is the normal choice.

One class behaves differently. **`oU` is effectively cigar-branch only** — 22,608 cigar-route reads
against 1 grep-route read in a single reference library, 550,807 against 14 across twelve. `oU`
counts therefore carry the library's branch mix with them, so check branch mix is not confounded
with your grouping before comparing that class. Branch mix is more stable within protein-coding
genes than across all biotypes.

## 7. Status classifications

| Column | Type | Levels |
|---|---|---|
| `loc_status` | factor | `confirmed`, `multi_locus`, `diff_chrom`, `far_from_gene`, `no_annotation`, `r2_unmapped` |
| `r2_status` | factor | `unique`, `multi`, `dropped_short`, `dropped_multi`, `dropped_mismatch`, `no_seed`, `dropped_mate` |
| `type` | factor | `cigar`, `grep`, `no_tail`, `dubious` |
| `class` | factor | `A`, `AC`, `AG`, `AU`, `C`, `G`, `mixed`, `none`, `oU`, `check` |

> **`loc_status` answers "did R2 land where R1's gene says it should?" — a gene-assignment verdict.
> It is not "does this read have an alignment?"** For alignment existence use `r2_genomic_mapped`.
> Confusing the two routes cleanly-aligned reads to the sequence-only peeler, which then calls a
> genomic A-stretch a tail. That is exactly what `flag_grep_misrouted` measures.

`far_from_gene` means the right chromosome but ≥50 kb from the gene's annotated 3' end.

## 8. Independent templated-tail check (bowtie2, unmasked genome, end-to-end)

A second alignment of R2 against the **unmasked** genome in end-to-end mode. It contributes no tail
and no coordinate. Its only job: if the whole read including the putative tail matched the genome,
the "tail" may be templated sequence.

| Column | Type | Meaning |
|---|---|---|
| `bt_ete_mapped` | logical | R2 aligned end-to-end to the unmasked genome |
| `bt_ete_unique` | logical | …and uniquely |
| `bt_ete_NM` | integer | Mismatch count of that alignment. `NA` when unmapped |
| `bt_nm_over_tail` | numeric | `bt_ete_NM / tail_len`. `NA` when unmapped or `tail_len == 0` |
| `templated_suspect` | logical | `bt_nm_over_tail <= 0.5` |
| `al_ete_bt` | factor | `yes`/`no`. Legacy version of the same idea, fires only for grep reads with `tail_len < 40` |

**Read the ratio the right way round.** A *high* `bt_nm_over_tail` is reassuring: median 1.00 on the
anchored branch means the end-to-end alignment mismatches through precisely the tail, the signature
of genuine non-templated addition. A *low* ratio is the suspicious case.

**`templated_suspect` is informative on the anchored branch and near-useless on the unanchored
one.** A real poly(A) tail and a genomic poly-A stretch are the same sequence, so an all-A read
aligns end-to-end at `NM` 0 to any long genomic A-run regardless of origin. On the grep branch it
fires on ~85% of end-to-end mappers. Prefer `al_ete_bt`'s successor columns and treat the flag as
suggestive, never as a filter.

## 9. In `_final_big.rds` only

| Column | Type | Meaning |
|---|---|---|
| `seq_R2` | character | Full R2 sequence after fastp |
| `tail` | character | The tail as a string, or `"full match"` |
| `transR1` | character | The transcript ID bowtie2 reported for R1 — an arbitrary pick among isoform ties, see §1 |
| `flagR1` | integer | R1's SAM flag |
| `CIGAR_R1` | character | R1's CIGAR |

Use the big table when you need the sequence itself — re-deriving a tail under a different rule, or
auditing a classification. Use the small table for everything else; it is roughly a third smaller.

## Typing and the schema stamp

`R_scripts/finalise_small_table.R` converts the character columns to factors with **fixed** level
sets, coerces the counts to integer, and stamps `attr(df, "schema_version")`. Fixed levels are what
make libraries bindable: without them, `rbind` across libraries silently reorders factor levels and
`chrR2 == seqnames` fails with *"level sets of factors are different"*.

> **The workflow applies this.** The R stage calls `retype_small_table()` then
> `write_small_atomic()`, so `_final_small.rds` arrives typed and stamped, with a provenance row in
> `logs/pipeline_manifest.tsv`. `scripts/rebuild/rebuild_one.R` re-applies the same two functions to
> one library from the existing `R_input/` tables when you need to re-finalise without re-aligning.
>
> `_final_big.rds` is **not** stamped — it is written before the typing, as the pre-finalise superset.
> If you are handed tables of unknown vintage:
>
> ```r
> stopifnot(identical(attr(d, "schema_version"), "2026-10-20-flag1"))
> ```

## Thresholds: what is baked in and what is not

| Threshold | Where | Changing it needs |
|---|---|---|
| STAR score/match fraction | `config/config_template.yaml` | a re-run — it decides whether a CIGAR exists at all |
| AT-masking rule (≥20 nt, 1 mismatch) | `scripts/mark_AT_20_genome.py` | a reference rebuild |
| UMI length, read-length floor | `Snakefile`, fastp rule | a re-run |
| **grep-branch minimum tail length** | **downstream only** | **nothing — filter the table** |
| **CIGAR-branch minimum tail length** | **downstream only** | **nothing — filter the table** |
| class vocabulary | `R_scripts/assign_type_and_class.R` | re-running the R stage only |

Short tails are present in the output down to 1 nt. Both length floors are applied by whatever
consumes these tables, which is why you can change your mind about them for free. The 2–3 nt bin is
the most artefact-prone; see the caveats in the main README before trusting it.
