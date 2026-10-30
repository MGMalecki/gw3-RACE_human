# Config comes from --configfile only; there is no default config file.
#
# A hardcoded `configfile: "config.yaml"` used to sit here, which made the workflow
# refuse to start unless that file existed, and silently merged it underneath whatever
# --configfile supplied. Per-library configs are generated from
# config/config_template.yaml by scripts/run_all_samples.sh, so a repo-level default is
# both unnecessary and a trap.
_required = ("sample", "R1", "R2")
_missing = [k for k in _required if k not in config]
if _missing:
    raise WorkflowError(
        "Missing config key(s): " + ", ".join(_missing) + ".\n"
        "Pass a per-library config, e.g.\n"
        "  snakemake --configfile config/generated/<stem>.yaml --cores 16\n"
        "or let the driver generate one:  ONLY=<stem> bash scripts/run_all_samples.sh\n"
        "Config files are built from config/config_template.yaml; see README.md."
    )

SAMPLE = config["sample"]
R1_PATH = config["R1"]
R2_PATH = config["R2"]

# --- Reference paths (overridable via config.yaml -> reference:) ---
REF = config.get("reference", {})
TRANSCRIPTOME_BOWTIE2_INDEX = REF.get(
    "transcriptome_bowtie2_index", "fasta/index_bowtie_transcriptome/prot_coding_"
)
GENOME_BOWTIE2_INDEX = REF.get(
    "genome_bowtie2_index", "fasta/index_bowtie_genome/GRCh38_primary_genome"
)
GENOME_STAR_MASKED_INDEX_DIR = REF.get(
    "genome_star_masked_index_dir", "fasta/index_masked_genome_star/"
)
GENE_ANNOTATION_RDS = REF.get(
    "gene_annotation_rds", "annotation/gene_annotations.rds"
)

# --- STAR multimapper handling ---
# Defaults reproduce the 2026-09-24 behaviour exactly. outFilterMultimapNmax 1 makes
# STAR discard any read mapping to more than one locus; it is written to the SAM as
# unmapped (uT:A:3), lands in R_input with CIGAR "*", and is tail-called by the
# unanchored fallback scan. Raising the filter and setting outSAMmultNmax 1 keeps
# only the best alignment, so the read count per SAM line is unchanged but the read
# gains a real CIGAR. NH:i: records how many loci it actually hit.
STAR_MULTIMAP_NMAX = config.get("star_multimap_nmax", 1)
STAR_SAM_MULT_NMAX = config.get("star_sam_mult_nmax", 1)

# --- STAR acceptance thresholds ---
# Defaults reproduce the 2026-09-24 command exactly.
#   score/match 0.66 is STAR's own default and is LENGTH-DEPENDENT in the variable
#   this assay measures: R2 is 81 nt after fastp, so a read whose non-templated tail
#   exceeds ~27 nt cannot cover 66% of its length and is rejected as uT:A:1.
#   mismatch 0.05 gives 4 mismatches per 81 nt; the observed nM distribution is
#   truncated exactly there (see implementation_plan.md §14).
STAR_SCORE_MIN_OVER_LREAD    = config.get("star_score_min_over_lread", 0.66)
STAR_MATCH_NMIN_OVER_LREAD   = config.get("star_match_nmin_over_lread", 0.66)
STAR_MISMATCH_NOVER_READLMAX = config.get("star_mismatch_nover_readlmax", 0.05)

rule all:
    input:
        f"R_output/{SAMPLE}_final_big.rds",
        f"R_output/{SAMPLE}_final_small.rds",
        f"R_output/{SAMPLE}_deduplication.csv"

rule align_read1:
    input:
        fastq = R1_PATH
    output:
        sam = temp("sam/{sample}_R1_aligned.sam")
    params:
        index = TRANSCRIPTOME_BOWTIE2_INDEX
    log:
        "logs/align/{sample}_bowtie2_transcriptome.log"
    threads: 8
    shell:
        """
        bowtie2 -p {threads} -x {params.index} -U {input.fastq} -S {output.sam} --no-unal 2> {log}
        """

rule extract_read1_names:
    input:
        "sam/{sample}_R1_aligned.sam"
    output:
        names = temp("sam/{sample}_R1_names.txt")
    shell:
        """
        awk '{{if ($1 !~ /^@/) print $1}}' {input} | sort | uniq > {output.names}
        """

rule filter_read2_mates:
    input:
        fastqR2 = R2_PATH,
        fastqR1 = R1_PATH,
        names = "sam/{sample}_R1_names.txt"
    output:
        read2_aligned = temp("fastq_read2/{sample}_R2_R1aligned.fq"),
        read1_aligned = temp("fastq_read1/{sample}_R1_R1aligned.fq")
    shell:
        """
        seqtk subseq {input.fastqR2} {input.names} > {output.read2_aligned}
        seqtk subseq {input.fastqR1} {input.names} > {output.read1_aligned}
        """

rule filter_read2_with_adapter_sequence:
    input:
        "fastq_read2/{sample}_R2_R1aligned.fq"
    output:
        temp("fastq_read2/{sample}_R2_adapter.fq")
    shell:
        """
        (grep -E -A 2 -B 1 --no-group-separator "^[ATCGN]{{15}}GTCAG" {input} || true) > {output}
        """

rule fastp_filtering_R2:
    input:
        read2 = "fastq_read2/{sample}_R2_adapter.fq"
    output:
        read2_fp = temp("fastq_read2/{sample}_R2_fp.fq"),
        report_html = "fastq_read2/{sample}_fp_report.html",
        report_json = "fastq_read2/{sample}_fp_report.json"
    threads: 4
    shell:
        # -Q disables quality filtering, and that is NOT tunable here on purpose.
        #
        # This assay measures non-templated 3' tails, and Illumina quality collapses
        # immediately after a homopolymer: the mean Phred of the 20 nt following the tail
        # falls from Q36.4 (no tail) to Q16.9 (57-80 nt tail) through cluster dephasing,
        # while the homopolymer itself stays at Q34-36. A quality filter therefore removes
        # exactly the reads the experiment is about. Measured with fastp's own -q 30 -u 10
        # on 94,198 adapter-bearing reads of T7_14_S326: 34.6% of all reads dropped, and
        # among reads with a tail >= 20 nt, 99.4-99.9%.
        #
        # -q/-u used to be passed alongside -Q, where fastp ignores them, and
        # enable_quality_filter defaulted to True -- so anyone running without the config
        # override silently got the filter ON. Removing all three closes that trap.
        # Verified byte-identical output before and after (same md5, same read count).
        """
        fastp \
            -w {threads} \
            -i {input.read2} \
            -o {output.read2_fp} \
            -Q \
            --umi \
            --umi_loc=read1 \
            --umi_len=15 \
            --trim_front1=5 \
            -l 80 \
            -h {output.report_html} \
            -j {output.report_json}
        """

rule align_R2_bowtie_genome_end_to_end:
    input:
        reads = "fastq_read2/{sample}_R2_fp.fq",
        # Declared so Snakemake hashes it. The `code` rerun-trigger covers only the shell
        # command string, which does not change when the script's contents do -- the same
        # trap as PIPELINE_EVALUATION.md finding 6.
        tagger = "R_scripts/extract_read_tags.py"
    output:
        sam = temp("sam/{sample}_R2_aligned_to_genome.sam"),
        tsv = "R_input/{sample}_R2_genome_mapped.tsv"
    params:
        index = GENOME_BOWTIE2_INDEX
    log:
        "logs/align/{sample}_bowtie2_genome.log"
    threads: 16
    shell:
        """
        bowtie2 --no-unal --threads {threads} -x {params.index} -U {input.reads} -S {output.sam} 2> {log}
        python {input.tagger} -i {output.sam} -o {output.tsv}
        """

rule align_R2_to_masked_genome:
    input:
        "fastq_read2/{sample}_R2_fp.fq"
    output:
        sam = temp("R2_aligned_to_masked_genome/{sample}_Aligned.out.sam")
    params:
        genome_dir = GENOME_STAR_MASKED_INDEX_DIR,
        multimap_nmax = STAR_MULTIMAP_NMAX,
        sam_mult_nmax = STAR_SAM_MULT_NMAX,
        score_min     = STAR_SCORE_MIN_OVER_LREAD,
        match_nmin    = STAR_MATCH_NMIN_OVER_LREAD,
        mismatch_nov  = STAR_MISMATCH_NOVER_READLMAX
    threads: 12
    shell:
        """
        STAR \
          --genomeDir {params.genome_dir} \
          --readFilesIn {input} \
          --runThreadN {threads} \
          --outFileNamePrefix R2_aligned_to_masked_genome/{wildcards.sample}_ \
          --outSAMtype SAM \
          --outFilterMismatchNoverReadLmax {params.mismatch_nov} \
          --outFilterScoreMinOverLread {params.score_min} \
          --outFilterMatchNminOverLread {params.match_nmin} \
          --alignIntronMax 50000 \
          --outFilterMultimapNmax {params.multimap_nmax} \
          --outSAMmultNmax {params.sam_mult_nmax} \
          --outSAMprimaryFlag OneBestScore \
          --outSAMattributes NH HI AS nM \
          --outSAMunmapped Within
        """

# R2 columns written here, in order:
#   1 QNAME  2 FLAG  3 RNAME  4 POS  5 CIGAR  6 SEQ  7 NH  8 uT  9 AS  10 nM
# 7-10 are optional SAM tags and must be pulled by regex, not by column index: they sit
# past field 11 at variable positions. uT:A: is STAR's reason code for an unmapped record
# (--outSAMunmapped Within), and is the only way to tell a read that aligned but was
# rejected for short alignment length (uT:A:1) from one with no alignment at all
# (uT:A:0) or one suppressed as a multimapper (uT:A:3). AS and nM make the acceptance
# thresholds reconstructable per read without re-aligning -- see run_AT_tail_analysis_df.R.
rule extract_columns_from_R1_R2:
    input:
        read1_sam = "sam/{sample}_R1_aligned.sam",
        read2_sam = "R2_aligned_to_masked_genome/{sample}_Aligned.out.sam"
    output:
        read1_cols = "R_input/{sample}_R1_aligned.tsv",
        read2_cols = "R_input/{sample}_R2_aligned_to_masked_genome.tsv"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} !/^@/ {{split($3, a, "|"); print $1, a[2], a[1], $2, $4, $5, $6}}' {input.read1_sam} > {output.read1_cols}
        awk 'BEGIN {{OFS="\\t"}} $1 !~ /^@/ {{ \
             nh="NA"; if (match($0, /NH:i:[0-9]+/))    nh=substr($0, RSTART+5, RLENGTH-5); \
             ut="NA"; if (match($0, /uT:A:[0-9]/))     ut=substr($0, RSTART+5, 1); \
             as="NA"; if (match($0, /AS:i:-?[0-9]+/))  as=substr($0, RSTART+5, RLENGTH-5); \
             nm="NA"; if (match($0, /nM:i:[0-9]+/))    nm=substr($0, RSTART+5, RLENGTH-5); \
             print $1, $2, $3, $4, $6, $10, nh, ut, as, nm }}' {input.read2_sam} > {output.read2_cols}
        """

rule tail_analysis:
    input:
        tabR1 = "R_input/{sample}_R1_aligned.tsv",
        tabR2 = "R_input/{sample}_R2_aligned_to_masked_genome.tsv",
        tabR2_2 = "R_input/{sample}_R2_genome_mapped.tsv",
        annotation = GENE_ANNOTATION_RDS,
        # Declared as inputs so Snakemake hashes them. The `code` rerun-trigger only
        # covers the rule's own `script:`, so before this was added an edit to either
        # helper left `snakemake -n` reporting "Nothing to be done" and results could
        # silently reflect stale code (PIPELINE_EVALUATION.md finding 6).
        helper_tail  = "R_scripts/run_AT_tail_analysis_df.R",
        helper_class = "R_scripts/assign_type_and_class.R",
        # The finaliser: typing, fixed factor levels and the schema_version stamp.
        # Declared here for the same reason as the two above -- an edit to it must
        # invalidate the output, and the `code` rerun-trigger does not see sourced files.
        helper_final = "R_scripts/finalise_small_table.R"
    output:
        big_table = "R_output/{sample}_final_big.rds",
        small_table = "R_output/{sample}_final_small.rds",
        deduplication = "R_output/{sample}_deduplication.csv"
    script:
        "R_scripts/run_tail_analysis.R"