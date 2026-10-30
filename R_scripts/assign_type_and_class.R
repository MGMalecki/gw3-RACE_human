assign_type_and_class <- function(df) {
  # Normalize C_count column name
  if (!"C_count" %in% colnames(df) && "C_len" %in% colnames(df)) {
    df <- df %>% rename(C_count = C_len)
  }

  if (!"tail_source" %in% colnames(df)) {
    stop("assign_type_and_class: expected a 'tail_source' column from ",
         "run_AT_tail_analysis_df() (cigar / full_match / fallback / NA). ",
         "Update run_AT_tail_analysis_df.R or check the input df.")
  }

  df <- df %>%
    mutate(
      # type now comes directly from the branch run_AT_tail_analysis_df actually
      # took for this read, instead of being re-derived from tail_len/CIGAR_R2
      # (which could mislabel a mapped-but-wrong-chromosome/too-far read as
      # "cigar" even though its tail was really called from the fallback scan).
      type = case_when(
        tail_source == "cigar"      ~ "cigar",     # soft-clip, confirmed nearby genomic location
        tail_source == "full_match" ~ "no_tail",    # confirmed location, no soft-clip -> no tail
        tail_source == "fallback"   ~ "grep",       # sequence-only scan, location NOT confirmed
        TRUE                        ~ "dubious"     # NA tail_source (invalid/empty input seq)
      ),
      class = case_when(
        tail_len == 0 | (A_count == 0 & T_count == 0 & G_count == 0 & C_count == 0) ~ "none",

        # Poly(A) and modified poly(A) classes
        A_count > 0  & T_count == 0 & G_count == 0 & C_count == 0 ~ "A",
        A_count > 0  & T_count > 0  & G_count == 0 & C_count == 0 ~ "AU",
        A_count > 0  & G_count > 0  & T_count == 0 & C_count == 0 ~ "AG",
        A_count > 0  & C_count > 0  & T_count == 0 & G_count == 0 ~ "AC",

        # Non-polyadenylated classes
        A_count == 0 & T_count > 0  & G_count == 0 & C_count == 0 ~ "oU",
        A_count == 0 & G_count > 0  & T_count == 0 & C_count == 0 ~ "G",
        A_count == 0 & C_count > 0  & T_count == 0 & G_count == 0 ~ "C",

        # Complex mixed non-A additions
        (T_count + G_count + C_count) > 1 ~ "mixed",

        TRUE ~ "check"
      )
    )

  return(df)
}