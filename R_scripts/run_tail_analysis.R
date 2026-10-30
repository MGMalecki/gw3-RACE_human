library(stringr)
library(stringi)
library(tidyr)
library(dplyr)

source("R_scripts/run_AT_tail_analysis_df.R")
source("R_scripts/assign_type_and_class.R")
source("R_scripts/finalise_small_table.R")

# Read input alignment tables
R1 <- read.delim(snakemake@input[["tabR1"]], header = FALSE, stringsAsFactors = FALSE)
R2 <- read.delim(snakemake@input[["tabR2"]], header = FALSE, stringsAsFactors = FALSE)

# NH:i: -- how many loci STAR placed this read at -- is captured as a 7th column only
# by Snakefiles that ask the awk for it. R_input tables written before 2026-09-25 have
# six columns; NH is NA for those, and every NH-dependent branch downstream is written
# to be a no-op when it is NA, so re-running this script on an old table reproduces the
# old result exactly.
R2$NH <- if (ncol(R2) >= 7) suppressWarnings(as.integer(R2$V7)) else NA_integer_

# uT:A: / AS:i: / nM:i: are captured as columns 8-10 (Snakefile rule
# extract_columns_from_R1_R2). They are what make two things possible:
#   uT     -- why an unmapped read is unmapped. STAR writes this only for unmapped
#             records, so it is NA for every read that aligned; r2_status below reads
#             it together with NH to give one status per read.
#   AS, nM -- the acceptance thresholds, reconstructable per read. STAR compares the
#             alignment score and the matched-base count to the READ length, so both
#             fractions can be recomputed in R and any stricter threshold applied after
#             the fact instead of by re-aligning.
# Tables written before 2026-09-30 have seven columns or fewer; all three become NA and
# every branch keyed on them downstream is a no-op, so re-running this script on an old
# table reproduces the old result exactly.
R2$uT <- if (ncol(R2) >= 8)  as.character(R2$V8)                    else NA_character_
R2$AS <- if (ncol(R2) >= 9)  suppressWarnings(as.integer(R2$V9))    else NA_integer_
R2$nM <- if (ncol(R2) >= 10) suppressWarnings(as.integer(R2$V10))   else NA_integer_
R2$uT[!is.na(R2$uT) & R2$uT == "NA"] <- NA_character_

# Extract UMIs and clean read names for joining
R2$UMI <- stringr::str_extract(R2$V1, "(?<=:)[^:]{15}$")
R1$V1 <- sub(":[^:]{15}$", "", R1$V1)
R2$V1 <- sub(":[^:]{15}$", "", R2$V1)

# Join mates on read headers
merged_dat <- dplyr::inner_join(
    dplyr::select(R1, V1, V2, V3, V4, V5, V7),
    dplyr::select(R2, V1, V2, V3, V4, V5, V6, UMI, NH, uT, AS, nM),
    by = "V1"
)

colnames(merged_dat)[1:11] <- c(
    "read_name", "geneR1", "transR1", "flagR1", "mapR1", "CIGAR_R1",
    "flagR2", "chrR2", "posR2", "CIGAR_R2", "seq_R2"
)

# Strip transcript version suffixes before deduplication
merged_dat$geneR1 <- sub("\\..*", "", merged_dat$geneR1)

# Deduplicate based on Gene ID, Transcription Coordinate, and UMI
dedup_dat <- merged_dat %>%
    distinct(geneR1, mapR1, UMI, .keep_all = TRUE)

dat <- data.frame(
    metric = c("before_dedup", "after_dedup"),
    count = c(nrow(merged_dat), nrow(dedup_dat))
)

# Merge reference annotations (path now comes from config.yaml via the Snakefile)
annotation <- readRDS(snakemake@input[["annotation"]])
dedup_dat <- dplyr::left_join(dedup_dat, annotation[, 1:5], by = c("geneR1" = "gene_id"))

# Analyze tail lengths and assign structural classes
dedup_dat <- run_AT_tail_analysis_df(dedup_dat)
dedup_dat <- assign_type_and_class(dedup_dat)

# --- The unmasked end-to-end check -------------------------------------------------
# bowtie2 ran against the UNMASKED genome in end-to-end mode with --no-unal, so this
# table holds one row per read that aligned over its FULL length. That is the opposite
# of a tail measurement and it is the only independent evidence the pipeline has that a
# called tail is really non-templated: if the whole read including the putative tail
# matched the genome, the "tail" is genomic sequence.
R2_2 <- read.delim(snakemake@input[["tabR2_2"]], header = FALSE, stringsAsFactors = FALSE)
R2_2$V1 <- sub(":[^:]{15}$", "", R2_2$V1)
R2_2$NM <- if (ncol(R2_2) >= 3) suppressWarnings(as.integer(R2_2$V3)) else NA_integer_
R2_2_unique <- R2_2 %>% filter(V2 == "unique")

# al_ete_bt -- UNCHANGED, so nothing downstream that already uses it breaks. Note the
# narrow predicate: it fires only for grep reads below 40 nt, so no cigar read is ever
# marked and long unanchored tails get no check at all.
grep_aligned <- dedup_dat %>%
    filter(tail_len < 40, type == "grep", read_name %in% R2_2_unique$V1) %>%
    pull(read_name)

dedup_dat$al_ete_bt <- ifelse(dedup_dat$read_name %in% grep_aligned, "yes", "no")

# The same evidence, as a plain per-read annotation covering EVERY read and every branch,
# with no length cutoff baked in. Filter on these at analysis time instead of inheriting
# the hard-coded 40 above.
bt_nm <- setNames(R2_2$NM, R2_2$V1)
dedup_dat$bt_ete_mapped <- dedup_dat$read_name %in% R2_2$V1
dedup_dat$bt_ete_unique <- dedup_dat$read_name %in% R2_2_unique$V1
dedup_dat$bt_ete_NM     <- unname(bt_nm[dedup_dat$read_name])

# Mapping end-to-end is NOT on its own evidence that a tail is templated: bowtie2 will
# force a read through by scoring every tail base as a mismatch. Measured on the anchored
# branch, bt_ete_NM / tail_len has median exactly 1.00 -- i.e. the typical end-to-end
# alignment mismatches through precisely the tail, which is the signature of a GENUINE
# non-templated addition. A tail is suspect only when that ratio is LOW, meaning most of
# the tail actually matched the genome.
#
# Emitted as a continuous ratio so the cut is an analysis choice, not baked in here.
dedup_dat$bt_nm_over_tail <- ifelse(
    dedup_dat$bt_ete_mapped & dedup_dat$tail_len > 0 & !is.na(dedup_dat$bt_ete_NM),
    dedup_dat$bt_ete_NM / dedup_dat$tail_len, NA_real_)

# The convenience boolean at half the tail length. Among anchored reads with tails >= 10 nt
# this flags ~28%; among 1-4 nt tails almost none, which is the right behaviour -- a 1 nt
# tail cannot be distinguished from a single mismatch either way.
#
# IMPORTANT CONFOUND, and the reason this is a flag and never a filter: a genuine poly(A)
# tail and a genomic poly-A stretch are the same sequence. An all-A read will align
# end-to-end at NM 0 to any long genomic A-run whether its tail is real or the product of
# internal priming. For the unanchored (grep) branch, where the peel can consume the whole
# read, this flag therefore cannot discriminate and fires on ~85% of end-to-end mappers.
# It is informative on the anchored branch and suggestive at best on the unanchored one.
dedup_dat$templated_suspect <- !is.na(dedup_dat$bt_nm_over_tail) &
                               dedup_dat$bt_nm_over_tail <= 0.5

# Write output files
saveRDS(dedup_dat, snakemake@output[["big_table"]])

# Columns kept in the big table only. read_name, chrR2/posR2/flagR2/CIGAR_R2 are
# deliberately retained in the small table: posR2 is the first aligned base of R2
# and the flag/CIGAR are what make it interpretable (strand, and the soft-clip
# offset to the read's real start) when referring it to an annotated TSS.
cols_to_remove <- c("seq_R2", "transR1", "CIGAR_R1", "tail", "flagR1")
small_table <- dedup_dat %>% select(-any_of(cols_to_remove))
message("small_table: dropped -> ",
        paste(setdiff(colnames(dedup_dat), colnames(small_table)), collapse = ", "))
message("small_table: kept    -> ", paste(colnames(small_table), collapse = ", "))

# --- Finalise: typing, fixed factor levels, schema_version stamp -------------------
# retype_small_table() and write_small_atomic() come from finalise_small_table.R and
# are the SAME functions scripts/rebuild/rebuild_one.R calls, in the same order. That
# is deliberate: the pipeline path and the standalone re-finalise path cannot drift,
# because there is only one implementation of the schema.
#
# Without this the table has the right 46 columns with the right contents but is
# untyped and unstamped, and binding libraries together then silently reorders factor
# levels -- `chrR2 == seqnames` fails with "level sets of factors are different".
small_table <- retype_small_table(small_table)

# --- Self-checks: cheap, and they fail the rule rather than writing a bad table -----
# Both verified to hold exactly on all twelve reference libraries.
#   1. tail_route partitions every row: no NA, nothing outside {cigar, grep}.
#   2. On the anchored branch the tail IS the 3' soft clip, so the two columns that
#      record it must agree. A mismatch means the routing and the measurement have
#      come apart, which is the bug the flag columns exist to detect.
rt <- as.character(small_table$tail_route)
stopifnot(!anyNA(rt), all(rt %in% c("cigar", "grep")))
cig <- which(rt == "cigar")
stopifnot(!anyNA(small_table$softclip_3p_len[cig]),
          all(small_table$softclip_3p_len[cig] == small_table$tail_len[cig]))
message("self-checks passed: tail_route partitions ", length(rt), " rows; ",
        "softclip_3p_len == tail_len on all ", length(cig), " cigar-route rows")

# write_small_atomic() writes to a .tmp, reads it back and checks nrow/ncol/stamp
# before renaming into place, then appends a provenance row (rows, columns, schema,
# size, md5) to the manifest. `stem` defaults to the basename of the output path.
write_small_atomic(small_table, snakemake@output[["small_table"]],
                   manifest = "logs/pipeline_manifest.tsv")

write.csv(dat, snakemake@output[["deduplication"]], row.names = FALSE)