run_AT_tail_analysis_df <- function(df) {
  library(stringr)
  library(stringi)

  n <- nrow(df)
  seq      <- df$seq_R2
  flag     <- df$flagR2
  cigar    <- df$CIGAR_R2
  chrR2    <- df$chrR2
  seqnames <- df$seqnames
  strand   <- df$strand
  start    <- df$start
  end      <- df$end
  posR2    <- df$posR2
  # NH:i: from STAR -- the number of loci this read was placed at. Absent (all NA) for
  # tables built before 2026-09-25 and always 1 when STAR ran with
  # --outFilterMultimapNmax 1, so every branch keyed on it below is inert in both cases.
  nh       <- if ("NH" %in% names(df)) df$NH else rep(NA_integer_, nrow(df))
  # uT:A: reason code, alignment score and mismatch count from STAR. All NA for tables
  # written before 2026-09-30; every branch below is written to be inert in that case.
  uT       <- if ("uT" %in% names(df)) as.character(df$uT) else rep(NA_character_, nrow(df))
  as_score <- if ("AS" %in% names(df)) as.integer(df$AS)   else rep(NA_integer_, nrow(df))
  nmis     <- if ("nM" %in% names(df)) as.integer(df$nM)   else rep(NA_integer_, nrow(df))

  total_length <- nchar(seq)
  valid <- !is.na(seq) & !is.na(total_length) & total_length > 0

  # --- AT percentage ---
  at_content <- (str_count(seq, "[AT]") / total_length) * 100
  high_AT_mix <- at_content > 93
  high_AT_mix[is.na(high_AT_mix)] <- FALSE

  # --- Orient sequence so tail sits at index 1 ---
  is_rev <- !is.na(flag) & flag == 16
  seq_oriented <- ifelse(is_rev,
                          stri_reverse(seq),
                          chartr("ATCGatcg", "TAGCtagc", seq))
  seq2 <- seq_oriented

  # --- Decide mapped vs fallback branch, per row ---
  stop_pos <- ifelse(!is.na(strand) & strand == "+", end,
               ifelse(!is.na(strand) & strand == "-", start, NA))

  # mapped_case is a single boolean, so four quite different failures used to be
  # indistinguishable downstream -- an R2 that never mapped, one that mapped to a
  # different chromosome than the R1 gene, one that mapped far from it, and a gene
  # missing from the annotation all ended up as type == "grep". loc_status records
  # which of them it was. Categories are mutually exclusive and cover every row.
  #   confirmed     = R2 mapped, annotated, same chromosome as the gene, within 50 kb
  #   r2_unmapped   = no R2 alignment (CIGAR "*" / RNAME "*"; STAR keeps these via
  #                   --outSAMunmapped Within)
  #   no_annotation = gene ID absent from the annotation table, so there is nothing
  #                   to compare against
  #   diff_chrom    = R2 mapped, but to a different chromosome than the R1 gene
  #   far_from_gene = right chromosome, but >= 50 kb from the gene's 3' end
  #   multi_locus   = R2 aligned, but to more than one locus (NH > 1). Which locus STAR
  #                   reported is arbitrary among equal-scoring ties, so the location
  #                   check is meaningless for these and is not attempted -- but the
  #                   soft-clip boundary is not, because paralogous copies of the same
  #                   sequence clip at the same offset. They therefore keep a CIGAR-based
  #                   tail call while being excluded from "confirmed".
  r2_unmapped   <- is.na(cigar) | cigar == "*" | is.na(chrR2) | chrR2 == "*"
  multi_locus   <- !r2_unmapped & !is.na(nh) & nh > 1
  no_annotation <- !r2_unmapped & (is.na(seqnames) | is.na(stop_pos))
  same_chrom    <- !r2_unmapped & !no_annotation & chrR2 == seqnames
  near_gene     <- same_chrom & !is.na(posR2) & abs(stop_pos - posR2) < 50000

  loc_status <- rep("unclassified", n)   # defensive: should never survive
  loc_status[r2_unmapped]                                <- "r2_unmapped"
  loc_status[no_annotation]                              <- "no_annotation"
  loc_status[!r2_unmapped & !no_annotation & !same_chrom] <- "diff_chrom"
  loc_status[same_chrom & !near_gene]                    <- "far_from_gene"
  loc_status[near_gene]                                  <- "confirmed"
  # Assigned last so it overrides whichever location verdict the arbitrary reported
  # locus happened to produce. Inert when NH is NA or 1, so old inputs are unaffected.
  loc_status[multi_locus]                                <- "multi_locus"

  # Both branches have a usable CIGAR. "confirmed" also has a verified location;
  # "multi_locus" does not, and that is what loc_status records.
  mapped_case <- loc_status %in% c("confirmed", "multi_locus")
  mapped_case[is.na(mapped_case)] <- FALSE
  mapped_case <- mapped_case & valid

  # --- Extract soft-clip length for mapped rows (direction depends on strand flag) ---
  softclip_str <- rep(NA_character_, n)
  idx_rev_map  <- mapped_case & is_rev
  idx_fwd_map  <- mapped_case & !is_rev
  softclip_str[idx_rev_map] <- str_extract(cigar[idx_rev_map], "\\d+(?=S$)")
  softclip_str[idx_fwd_map] <- str_extract(cigar[idx_fwd_map], "^\\d+(?=S)")
  softclip_len <- suppressWarnings(as.integer(softclip_str))

  has_softclip  <- mapped_case & !is.na(softclip_len) & softclip_len > 0
  full_match    <- mapped_case & !has_softclip           # cigar mapped, no usable softclip
  fallback_case <- valid & !mapped_case                  # unmapped / not on same chr / too far

  # --- Soft clip for EVERY row that has a real CIGAR, not only mapped_case rows -------
  # softclip_len above is deliberately left gated on mapped_case, because it drives the
  # tail call and must keep producing byte-identical results. This is a second, ungated
  # measurement used only for the diagnostic flags below.
  #
  # Why it is needed: loc_status answers "did R2 land where R1's gene says it should?",
  # a GENE-ASSIGNMENT verdict. mapped_case uses it as if it answered "does this read have
  # an alignment at all?", an ALIGNMENT-EXISTENCE verdict. A read whose R2 aligned
  # uniquely and end-to-end but landed on another chromosome, >= 50 kb away, or in
  # unannotated space is therefore handed to the sequence-only peeler even though its
  # CIGAR is right there. These columns measure that, and nothing else: the routing is
  # NOT changed here (see implementation_plan.md D1/§2.1).
  #
  # Same strand convention as above -- flag 16 means the mRNA 3' end is the RIGHT-hand
  # end of the alignment, so the clip is the trailing \d+S; otherwise it is the leading
  # one. NA only where there is no CIGAR to read; a real CIGAR with no clip on the 3'
  # side is 0, which is the whole point of the flag.
  real_cigar <- !is.na(cigar) & cigar != "*" & !is.na(chrR2) & chrR2 != "*"
  softclip_3p <- rep(NA_integer_, n)
  softclip_3p[real_cigar &  is_rev] <- suppressWarnings(as.integer(
    str_extract(cigar[real_cigar &  is_rev], "\\d+(?=S$)")))
  softclip_3p[real_cigar & !is_rev] <- suppressWarnings(as.integer(
    str_extract(cigar[real_cigar & !is_rev], "^\\d+(?=S)")))
  softclip_3p[real_cigar & is.na(softclip_3p)] <- 0L

  # --- Explicit branch label: this is the ground truth for downstream type assignment ---
  # "cigar"      = soft-clip parsed from a confirmed, nearby genomic alignment (high confidence)
  # "full_match" = mapped + confirmed location, but no soft-clip -> effectively no tail
  # "fallback"   = sequence-only scan; either unmapped, wrong chromosome, or too far from
  #                the annotated 3' end -- genomic location NOT confirmed for the tail call
  # NA           = invalid/empty input sequence
  tail_source <- rep(NA_character_, n)
  tail_source[has_softclip]  <- "cigar"
  tail_source[full_match]    <- "full_match"
  tail_source[fallback_case] <- "fallback"

  work_seq <- rep(NA_character_, n)
  tail_len <- rep(0L, n)
  tail_out <- rep(NA_character_, n)

  # mapped rows with a real softclip: tail = first softclip_len chars of oriented seq
  work_seq[has_softclip] <- str_sub(seq2[has_softclip], 1, softclip_len[has_softclip])
  tail_out[has_softclip] <- work_seq[has_softclip]
  tail_len[has_softclip] <- softclip_len[has_softclip]

  # mapped rows with no usable softclip
  tail_out[full_match] <- "full match"

  # fallback rows: peel straight from the whole oriented sequence
  work_seq[fallback_case] <- seq2[fallback_case]

  # --- Vectorized T -> G -> C -> A peeling wherever work_seq is set ---
  peel_idx <- !is.na(work_seq)
  cur <- work_seq

  t_len <- rep(0L, n); g_len <- rep(0L, n); c_len <- rep(0L, n); a_count <- rep(0L, n)

  starts_T <- peel_idx & str_sub(cur, 1, 1) == "T"
  if (any(starts_T)) {
    t_match <- str_extract(cur[starts_T], "^T+")
    t_len[starts_T] <- nchar(t_match)
    cur[starts_T] <- str_sub(cur[starts_T], t_len[starts_T] + 1)
  }

  starts_G <- peel_idx & str_sub(cur, 1, 1) == "G"
  if (any(starts_G)) {
    g_match <- str_extract(cur[starts_G], "^G+")
    g_len[starts_G] <- nchar(g_match)
    cur[starts_G] <- str_sub(cur[starts_G], g_len[starts_G] + 1)
  }

  starts_C <- peel_idx & str_sub(cur, 1, 1) == "C"
  if (any(starts_C)) {
    c_match <- str_extract(cur[starts_C], "^C+")
    c_len[starts_C] <- nchar(c_match)
    cur[starts_C] <- str_sub(cur[starts_C], c_len[starts_C] + 1)
  }

  a_match <- rep(NA_character_, n)
  a_match[peel_idx] <- str_extract(cur[peel_idx], "^A+(?:[^A]A{5,})*")
  a_len <- ifelse(!is.na(a_match), nchar(a_match), 0L)
  a_count[peel_idx] <- ifelse(!is.na(a_match[peel_idx]), str_count(a_match[peel_idx], "A"), 0L)

  # tail_len / tail for fallback rows (mirrors original: sum of peeled pieces)
  tail_len[fallback_case] <- a_len[fallback_case] + t_len[fallback_case] +
                              g_len[fallback_case] + c_len[fallback_case]
  tail_out[fallback_case] <- str_sub(seq2[fallback_case], 1, tail_len[fallback_case])

  # --- Reset invalid (NA/empty seq) rows to match original early-return behavior ---
  a_count[!valid] <- 0L; t_len[!valid] <- 0L; g_len[!valid] <- 0L; c_len[!valid] <- 0L
  tail_len[!valid] <- 0L; tail_out[!valid] <- NA_character_; high_AT_mix[!valid] <- FALSE
  tail_source[!valid] <- NA_character_

  df$A_count     <- a_count
  df$T_count     <- t_len
  df$G_count     <- g_len
  df$C_len       <- c_len
  df$tail_len    <- tail_len
  df$high_AT_mix <- high_AT_mix
  df$tail        <- tail_out
  df$tail_source <- tail_source
  # Kept for every row, including invalid-sequence ones: it describes where R2
  # landed relative to the gene, which is independent of whether the sequence
  # was usable for tail calling.
  df$loc_status  <- loc_status
  # Always present so the schema is the same whether or not the input carried NH.
  # NA = the R_input table predates NH capture; 1 = uniquely placed; >1 = multi_locus.
  df$NH          <- nh

  # --- r2_status: WHY this read is, or is not, placed on the masked genome -----------
  # loc_status collapses every kind of failure into "r2_unmapped". STAR actually records
  # the reason on each unmapped record as uT:A:<n>, and it matters which one it is:
  #
  #   dropped_short  (uT:A:1)  An alignment WAS found and then rejected because the
  #                            aligned portion was too short a fraction of the read.
  #                            This bucket is entangled with the measurement itself --
  #                            a read is short-aligned largely BECAUSE much of it is
  #                            non-templated tail -- so these are reads whose tail is
  #                            real but whose location was unrecoverable at the
  #                            threshold in force. Raising outFilterScoreMinOverLread
  #                            moves reads out of here and into "unique"/"multi".
  #   dropped_multi  (uT:A:3)  Aligned fine, to more places than outFilterMultimapNmax
  #                            allows. Nothing to do with the tail.
  #   no_seed        (uT:A:0)  No acceptable seed or window: no alignment at all.
  #   dropped_mismatch (uT:A:2) Too many mismatches.
  #   dropped_mate   (uT:A:4)  Unmapped mate of a pair (not applicable here).
  #
  # Orthogonal to loc_status and to type: it describes the aligner's verdict, not the
  # location check or the tail-calling branch. NA only when the input predates uT capture.
  df$r2_status <- dplyr::case_when(
    !is.na(nh) & nh >  1 ~ "multi",
    !is.na(nh) & nh == 1 ~ "unique",
    uT == "1"            ~ "dropped_short",
    uT == "3"            ~ "dropped_multi",
    uT == "2"            ~ "dropped_mismatch",
    uT == "0"            ~ "no_seed",
    uT == "4"            ~ "dropped_mate",
    TRUE                 ~ NA_character_
  )

  # --- The genomic-tail flag group ---------------------------------------------------
  # Six diagnostic columns. Nothing is dropped, nothing is reclassified, no routing
  # changes: every one of type / class / tail_len / A_count / T_count / G_count /
  # C_len is produced by the code above, untouched.
  #
  # tail_route collapses the three-way tail_source onto the two real peelers. Note that
  # "full_match" (type == "no_tail") is a CIGAR-route read that HAS no tail, not a third
  # branch, so it is route "cigar". Any pooled denominator has to decide explicitly
  # whether those reads are in it.
  df$tail_route <- ifelse(is.na(tail_source), NA_character_,
                   ifelse(tail_source == "fallback", "grep", "cigar"))

  # Does this read have an alignment at all -- the question loc_status does NOT answer.
  df$r2_genomic_mapped <- !is.na(df$r2_status) & df$r2_status %in% c("unique", "multi")

  df$softclip_3p_len <- softclip_3p

  # Why this read was sent to the sequence-only peeler. loc_status collapses every
  # alignment failure into "r2_unmapped"; where that is the reason, r2_status says which
  # failure it was, so the two are composed here rather than losing the distinction.
  gr <- !is.na(df$tail_route) & df$tail_route == "grep"
  reason <- rep(NA_character_, n)
  reason[gr] <- loc_status[gr]
  unmapped_gr <- gr & loc_status == "r2_unmapped"
  reason[unmapped_gr] <- ifelse(is.na(df$r2_status[unmapped_gr]),
                                "r2_unmapped", df$r2_status[unmapped_gr])
  df$grep_route_reason <- reason

  # TWO flags, because the two rules are not the same set and do not make the same claim.
  # Measured on T7_16_S324: 154,781 reads are grep-routed despite having an alignment,
  # but only 66,122 of them (42.7%) have a zero-length 3' soft clip. The other 88,659
  # have a clip of 1-40 nt -- their alignment AGREES there is non-templated 3' material
  # and merely disagrees about how much. See docs/SCHEMA_final_small.md s1.
  #
  #   flag_grep_misrouted  the routing bug itself. Identical to the grep_genomic
  #                        workaround in Data_analysis_michal.R::data_filter(), so that
  #                        workaround can be checked as an exact no-op and retired.
  #   flag_genomic_tail    the narrower, better-founded claim: the alignment positively
  #                        CONTRADICTS the sequence-derived tail -- nothing was clipped,
  #                        yet grep reported a tail.
  #
  # Neither is conditioned on tail_len. They mean one thing each; length thresholds stay
  # in TAIL_FILTERS where all the other thresholds live.
  df$flag_grep_misrouted <- gr & df$r2_genomic_mapped
  df$flag_genomic_tail   <- df$flag_grep_misrouted &
                            !is.na(df$softclip_3p_len) & df$softclip_3p_len == 0L

  # --- Reconstruct any stricter acceptance threshold, without re-aligning ------------
  # STAR compares both the alignment score and the number of matched bases to the READ
  # length. Both numerators are recoverable per read: AS directly, and matched bases as
  # (sum of M ops) - nM. So the fractions STAR itself thresholds on can be recomputed,
  # and a stricter setting applied after the fact.
  #
  # CAVEAT, and it is a real one: this is not exactly equivalent to having run STAR at
  # the stricter value. The filters also gate which loci survive, so under a looser
  # setting a read can come back with a higher NH or at a different locus among
  # equal-scoring ties. Treat these flags as reliable for "would this read have been
  # reported at all" and weaker for "where".
  m_ops <- stringr::str_extract_all(ifelse(is.na(cigar), "", cigar), "\\d+(?=M)")
  cig_M <- vapply(m_ops, function(v) if (length(v)) sum(as.integer(v)) else 0L, integer(1))
  matched <- cig_M - ifelse(is.na(nmis), 0L, nmis)

  df$read_len   <- total_length
  df$score_frac <- ifelse(total_length > 0, as_score / total_length, NA_real_)
  df$match_frac <- ifelse(total_length > 0, matched  / total_length, NA_real_)
  for (thr in c(0.30, 0.40, 0.50, 0.66)) {
    df[[sprintf("passes_%03d", round(thr * 100))]] <-
      !is.na(df$score_frac) & !is.na(df$match_frac) &
      df$score_frac >= thr & df$match_frac >= thr
  }

  return(df)
}