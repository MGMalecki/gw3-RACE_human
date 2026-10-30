#!/usr/bin/env Rscript
# Rebuild ONE library's _final_small.rds in place, with the flag group and the retyped
# schema, asserting identity against the existing file BEFORE it is replaced.
#
#   Rscript scripts/rebuild/rebuild_one.R <stem> [--scratch <dir>]
#
# --scratch writes to <dir>/<stem>_final_small.rds and never touches R_output/.
#
# Reproduces the Snakemake rule `tail_analysis` exactly (same inputs, same two helpers,
# same read.delim parsing) minus the big table and the dedup csv, which this pass does
# not touch. _final_big.rds is NOT written.

t0 <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: rebuild_one.R <stem> [--scratch <dir>] [--nrows <n>]")
stem <- args[1]
scratch <- if ("--scratch" %in% args) args[which(args == "--scratch") + 1] else NA_character_
nrows   <- if ("--nrows"   %in% args) as.integer(args[which(args == "--nrows") + 1]) else -1L

# --- RAM discipline (implementation_plan.md s2.3) -----------------------------------
Sys.setenv(OMP_NUM_THREADS = "2", MKL_NUM_THREADS = "2", OPENBLAS_NUM_THREADS = "2")
options(mc.cores = 1)
suppressPackageStartupMessages({
  library(stringr); library(stringi); library(tidyr); library(dplyr)
})
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(2)

source("R_scripts/run_AT_tail_analysis_df.R")
source("R_scripts/assign_type_and_class.R")
source("R_scripts/finalise_small_table.R")

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
msg("stem:", stem, "| schema:", SCHEMA_VERSION, "| nrows:", nrows)

tabR1       <- sprintf("R_input/%s_R1_aligned.tsv", stem)
tabR2       <- sprintf("R_input/%s_R2_aligned_to_masked_genome.tsv", stem)
tabR2_2     <- sprintf("R_input/%s_R2_genome_mapped.tsv", stem)
annot_path  <- "annotation/gene_annotations_v50.rds"   # config/generated/*.yaml
for (f in c(tabR1, tabR2, tabR2_2, annot_path)) if (!file.exists(f)) stop("missing input: ", f)

out_path <- if (is.na(scratch)) {
  sprintf("R_output/%s_final_small.rds", stem)
} else {
  dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
  file.path(scratch, sprintf("%s_final_small.rds", stem))
}
log_path <- if (is.na(scratch)) {
  sprintf("logs/rebuild/%s_compare.txt", stem)
} else {
  file.path(scratch, sprintf("%s_compare.txt", stem))
}
dir.create("logs/rebuild", showWarnings = FALSE, recursive = TRUE)

# ---------------- everything below mirrors R_scripts/run_tail_analysis.R -------------
msg("reading R1")
R1 <- read.delim(tabR1, header = FALSE, stringsAsFactors = FALSE, nrows = nrows)
msg("reading R2")
R2 <- read.delim(tabR2, header = FALSE, stringsAsFactors = FALSE, nrows = nrows)

R2$NH <- if (ncol(R2) >= 7)  suppressWarnings(as.integer(R2$V7)) else NA_integer_
R2$uT <- if (ncol(R2) >= 8)  as.character(R2$V8)                 else NA_character_
R2$AS <- if (ncol(R2) >= 9)  suppressWarnings(as.integer(R2$V9)) else NA_integer_
R2$nM <- if (ncol(R2) >= 10) suppressWarnings(as.integer(R2$V10)) else NA_integer_
R2$uT[!is.na(R2$uT) & R2$uT == "NA"] <- NA_character_

R2$UMI <- stringr::str_extract(R2$V1, "(?<=:)[^:]{15}$")
R1$V1 <- sub(":[^:]{15}$", "", R1$V1)
R2$V1 <- sub(":[^:]{15}$", "", R2$V1)

msg("joining mates")
merged_dat <- dplyr::inner_join(
  dplyr::select(R1, V1, V2, V3, V4, V5, V7),
  dplyr::select(R2, V1, V2, V3, V4, V5, V6, UMI, NH, uT, AS, nM),
  by = "V1")
rm(R1, R2); gc(full = TRUE)

colnames(merged_dat)[1:11] <- c(
  "read_name", "geneR1", "transR1", "flagR1", "mapR1", "CIGAR_R1",
  "flagR2", "chrR2", "posR2", "CIGAR_R2", "seq_R2")
merged_dat$geneR1 <- sub("\\..*", "", merged_dat$geneR1)

dedup_dat <- merged_dat %>% distinct(geneR1, mapR1, UMI, .keep_all = TRUE)
msg("dedup:", nrow(merged_dat), "->", nrow(dedup_dat))
rm(merged_dat); gc(full = TRUE)

annotation <- readRDS(annot_path)
dedup_dat <- dplyr::left_join(dedup_dat, annotation[, 1:5], by = c("geneR1" = "gene_id"))

msg("tail analysis")
dedup_dat <- run_AT_tail_analysis_df(dedup_dat)
dedup_dat <- assign_type_and_class(dedup_dat)

msg("bowtie2 end-to-end evidence")
R2_2 <- read.delim(tabR2_2, header = FALSE, stringsAsFactors = FALSE)
R2_2$V1 <- sub(":[^:]{15}$", "", R2_2$V1)
R2_2$NM <- if (ncol(R2_2) >= 3) suppressWarnings(as.integer(R2_2$V3)) else NA_integer_
R2_2_unique <- R2_2 %>% filter(V2 == "unique")

grep_aligned <- dedup_dat %>%
  filter(tail_len < 40, type == "grep", read_name %in% R2_2_unique$V1) %>%
  pull(read_name)
dedup_dat$al_ete_bt <- ifelse(dedup_dat$read_name %in% grep_aligned, "yes", "no")

bt_nm <- setNames(R2_2$NM, R2_2$V1)
dedup_dat$bt_ete_mapped <- dedup_dat$read_name %in% R2_2$V1
dedup_dat$bt_ete_unique <- dedup_dat$read_name %in% R2_2_unique$V1
dedup_dat$bt_ete_NM     <- unname(bt_nm[dedup_dat$read_name])
rm(R2_2, R2_2_unique, bt_nm, grep_aligned); gc(full = TRUE)

dedup_dat$bt_nm_over_tail <- ifelse(
  dedup_dat$bt_ete_mapped & dedup_dat$tail_len > 0 & !is.na(dedup_dat$bt_ete_NM),
  dedup_dat$bt_ete_NM / dedup_dat$tail_len, NA_real_)
dedup_dat$templated_suspect <- !is.na(dedup_dat$bt_nm_over_tail) &
                               dedup_dat$bt_nm_over_tail <= 0.5

cols_to_remove <- c("seq_R2", "transR1", "CIGAR_R1", "tail", "flagR1")
small_table <- dedup_dat %>% select(-any_of(cols_to_remove))
rm(dedup_dat); gc(full = TRUE)
# ---------------- end of the mirrored section ----------------------------------------

msg("retyping")
small_table <- retype_small_table(small_table)

cat("\n=== ", stem, " @ ", format(Sys.time()), " ===\n", sep = "", file = log_path, append = TRUE)
msg("comparing against the file about to be replaced")
ok <- if (nrows > 0) {
  msg("  (subset run: full-file comparison skipped by design)"); TRUE
} else {
  compare_to_existing(small_table, sprintf("R_output/%s_final_small.rds", stem), log_path)
}
if (!ok) stop("COMPARISON FAILED for ", stem, " -- nothing written, loop must abort")

msg("flag summary")
fs <- small_table %>%
  summarise(n = n(),
            route_cigar = sum(tail_route == "cigar", na.rm = TRUE),
            route_grep  = sum(tail_route == "grep",  na.rm = TRUE),
            route_na    = sum(is.na(tail_route)),
            misrouted   = sum(flag_grep_misrouted),
            genomic     = sum(flag_genomic_tail))
print(as.data.frame(fs))
cat(sprintf("flags\tn=%d\troute_cigar=%d\troute_grep=%d\troute_na=%d\tmisrouted=%d\tgenomic=%d\n",
            fs$n, fs$route_cigar, fs$route_grep, fs$route_na, fs$misrouted, fs$genomic),
    file = log_path, append = TRUE)
stopifnot(fs$route_cigar + fs$route_grep + fs$route_na == fs$n)

msg("writing", out_path)
write_small_atomic(small_table, out_path,
                   manifest = if (is.na(scratch)) "logs/rebuild/manifest.tsv" else NULL,
                   stem = stem)
msg("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min",
    "| peak R RSS", round(sum(gc()[, "max used"] * c(56, 8)) / 1048576, 1), "MB")
