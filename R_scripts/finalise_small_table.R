# Typing, schema stamping and the safe in-place writer for <stem>_final_small.rds.
#
# Sourced by both the Snakemake path (R_scripts/run_tail_analysis.R) and the standalone
# rebuild driver (scripts/rebuild/rebuild_one.R), so the two cannot drift apart.

SCHEMA_VERSION <- "2026-10-20-flag1"

# Fixed level sets. Every library gets the same levels in the same order for these
# columns whether or not it uses them all, so rbind() and bind_rows() both behave.
# geneR1 / chrR2 / seqnames are library-specific by nature and are left data-driven.
SCHEMA_LEVELS <- list(
  strand            = c("-", "+"),
  tail_source       = c("cigar", "full_match", "fallback"),
  loc_status        = c("confirmed", "multi_locus", "diff_chrom", "far_from_gene",
                        "no_annotation", "r2_unmapped"),
  r2_status         = c("unique", "multi", "dropped_short", "dropped_multi",
                        "dropped_mismatch", "no_seed", "dropped_mate"),
  type              = c("cigar", "grep", "no_tail", "dubious"),
  class             = c("A", "AC", "AG", "AU", "C", "G", "mixed", "none", "oU", "check"),
  al_ete_bt         = c("no", "yes"),
  uT                = c("0", "1", "2", "3", "4"),
  tail_route        = c("cigar", "grep"),
  grep_route_reason = c("diff_chrom", "far_from_gene", "no_annotation",
                        "no_seed", "dropped_short", "dropped_multi",
                        "dropped_mismatch", "dropped_mate", "r2_unmapped")
)

# Columns factored with whatever levels the library actually contains.
SCHEMA_FACTOR_DATA_DRIVEN <- c("geneR1", "chrR2", "seqnames")

SCHEMA_INTEGER <- c("tail_len", "A_count", "T_count", "G_count", "C_count",
                    "read_len", "NH", "start", "end", "posR2", "mapR1", "flagR2",
                    "AS", "nM", "bt_ete_NM", "softclip_3p_len")

retype_small_table <- function(df) {
  for (nm in names(SCHEMA_LEVELS)) {
    if (!nm %in% names(df)) next
    v <- as.character(df[[nm]])
    extra <- setdiff(unique(v[!is.na(v)]), SCHEMA_LEVELS[[nm]])
    if (length(extra)) {
      stop(sprintf("retype_small_table: column '%s' has unexpected value(s): %s",
                   nm, paste(extra, collapse = ", ")))
    }
    df[[nm]] <- factor(v, levels = SCHEMA_LEVELS[[nm]])
  }
  for (nm in intersect(SCHEMA_FACTOR_DATA_DRIVEN, names(df))) {
    df[[nm]] <- factor(as.character(df[[nm]]))
  }
  for (nm in intersect(SCHEMA_INTEGER, names(df))) {
    v <- df[[nm]]
    if (is.integer(v)) next
    iv <- suppressWarnings(as.integer(v))
    # Lossless or refuse: never silently truncate.
    bad <- which(!is.na(v) & (is.na(iv) | iv != v))
    if (length(bad)) {
      stop(sprintf("retype_small_table: integer downcast of '%s' is lossy at %d row(s)",
                   nm, length(bad)))
    }
    df[[nm]] <- iv
  }
  attr(df, "schema_version") <- SCHEMA_VERSION
  df
}

# Every column that existed before this change. Values must be bit-identical, modulo
# the factor/integer encoding. This is the test that licenses the in-place overwrite.
compare_to_existing <- function(new_df, old_path, log_path = NULL) {
  say <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n", sep = "")
    if (!is.null(log_path)) cat(msg, "\n", sep = "", file = log_path, append = TRUE)
    invisible(NULL)
  }
  if (!file.exists(old_path)) {
    say("COMPARE SKIPPED: no existing file at ", old_path)
    return(TRUE)
  }
  old <- readRDS(old_path)
  ok <- TRUE
  say("compare: ", basename(old_path))
  say("  old nrow=", nrow(old), " ncol=", ncol(old),
      " | new nrow=", nrow(new_df), " ncol=", ncol(new_df))
  if (nrow(old) != nrow(new_df)) {
    say("  FAIL: row count changed"); return(FALSE)
  }
  new_cols <- setdiff(names(new_df), names(old))
  gone     <- setdiff(names(old), names(new_df))
  say("  added columns (", length(new_cols), "): ", paste(new_cols, collapse = ", "))
  if (length(gone)) { say("  FAIL: columns disappeared: ", paste(gone, collapse = ", ")); ok <- FALSE }

  for (nm in names(old)) {
    if (!nm %in% names(new_df)) next
    a <- old[[nm]]; b <- new_df[[nm]]
    if (is.factor(b)) b <- as.character(b)
    if (is.factor(a)) a <- as.character(a)
    if (is.integer(b) && is.numeric(a) && !is.integer(a)) a <- as.integer(a)
    same <- isTRUE(all.equal(a, b, check.attributes = FALSE))
    if (!same) {
      n_diff <- sum(!(is.na(a) & is.na(b)) & (is.na(a) | is.na(b) | a != b))
      say("  FAIL ", nm, ": ", n_diff, " differing row(s)")
      ok <- FALSE
    }
  }
  say("  VERDICT: ", if (ok) "IDENTICAL on every retained column" else "MISMATCH")
  rm(old); gc(full = TRUE)
  ok
}

# Write to <path>.tmp then rename. An interrupted run can never leave a truncated final.
write_small_atomic <- function(df, path, manifest = "logs/rebuild/manifest.tsv",
                               stem = NULL) {
  tmp <- paste0(path, ".tmp")
  saveRDS(df, tmp, compress = "gzip")
  chk <- readRDS(tmp)
  if (nrow(chk) != nrow(df) || ncol(chk) != ncol(df) ||
      !identical(attr(chk, "schema_version"), SCHEMA_VERSION)) {
    stop("write_small_atomic: readback of ", tmp, " does not match what was written")
  }
  rm(chk); gc(full = TRUE)
  if (!file.rename(tmp, path)) stop("write_small_atomic: rename failed for ", path)

  if (!is.null(manifest)) {
    dir.create(dirname(manifest), showWarnings = FALSE, recursive = TRUE)
    if (!file.exists(manifest)) {
      cat("stem\tnrow\tncol\tschema_version\tmtime\tsize_mb\tmd5\n", file = manifest)
    }
    cat(sprintf("%s\t%d\t%d\t%s\t%s\t%.1f\t%s\n",
                if (is.null(stem)) sub("_final_small\\.rds$", "", basename(path)) else stem,
                nrow(df), ncol(df), SCHEMA_VERSION,
                format(file.mtime(path), "%Y-%m-%dT%H:%M:%S"),
                file.size(path) / 1048576,
                unname(tools::md5sum(path))),
        file = manifest, append = TRUE)
  }
  invisible(path)
}
