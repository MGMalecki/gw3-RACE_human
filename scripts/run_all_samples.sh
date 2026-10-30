#!/usr/bin/env bash
# Run every sample in config/samples.tsv through the pipeline, sequentially.
#
# Naming: each sample runs under its ORIGINAL FASTQ stem (T7_14_S326), so every
# file the pipeline writes - intermediates, alignment logs, final outputs - is
# already named after the sample it came from. Renaming the finals to friendly
# labels is deliberately deferred until after the run (see implementation_plan.md
# Step 9); the output_name column here is carried as documentation only.
#
# Intermediates: deleted by Snakemake itself via the Snakefile's temp() markers.
# This script issues no rm. The 2026-09-21 run kept them only because
# run_pipeline.sh passed --notemp; that flag is deliberately absent here.
#
# Sequential by design: the STAR index is ~24 GB resident, so two concurrent
# STAR jobs would want ~48 GB plus sort buffers. --cores still parallelises
# within a sample.
#
# Env:
#   ONLY=<stem|label>   run just that one sample
#   CORES=<n>           snakemake --cores (default 16)
#   CONDA_ENV=<name>    conda env to activate (default seq)
#   CONDA_SH=<path>     path to etc/profile.d/conda.sh, if autodetection fails
#   REPO=<path>         repo root (default: the parent of this script)

set -uo pipefail

# Repo root is derived from this script's own location, so a clone works
# anywhere. REPO=<path> overrides it if you call the script from elsewhere.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
cd "$REPO" || exit 1

# The env's activate.d/activate-binutils_linux-64.sh references unset vars
# (ADDR2LINE and friends), which is fatal under `set -u`. Drop -u across the
# activation only, then restore it.
set +u
# conda.sh location differs per install (miniconda3/anaconda3/mambaforge, system
# or user). Honour CONDA_SH if set, else try the usual places, else fail with a
# message naming the variable rather than a bare "conda: command not found".
if [ -n "${CONDA_SH:-}" ]; then
  source "$CONDA_SH"
else
  for c in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" \
           "$HOME/miniforge3" /opt/conda; do
    if [ -f "$c/etc/profile.d/conda.sh" ]; then source "$c/etc/profile.d/conda.sh"; break; fi
  done
fi
command -v conda >/dev/null || { echo "conda not found: set CONDA_SH=/path/to/etc/profile.d/conda.sh" >&2; exit 1; }
conda activate "${CONDA_ENV:-seq}"
set -u

# Containment: every scratch file stays inside the working folder. Snakemake
# otherwise defaults resources.tmpdir to /tmp, and the sort in
# extract_read1_names can spill several GB there.
mkdir -p "$REPO/tmp" "$REPO/logs/align" "$REPO/config/generated"
export TMPDIR="$REPO/tmp"

CORES=${CORES:-16}
STATUS=logs/run_status.tsv
[ -f "$STATUS" ] || printf 'pipeline_name\tplanned_label\tstarted\tfinished\texit\tnotes\n' > "$STATUS"

# STAR writes these as plain files, not declared outputs, so Snakemake never
# removes them. Copy (never move) into logs/align/ so all alignment stats for a
# sample sit together under one name.
collect_align_logs() {
    local pn="$1"
    local suf f
    for suf in Log.final.out SJ.out.tab Log.out; do
        f="R2_aligned_to_masked_genome/${pn}_${suf}"
        [ -f "$f" ] && cp -p "$f" logs/align/ && echo "  log: logs/align/${pn}_${suf}"
    done
    for f in "logs/align/${pn}_bowtie2_transcriptome.log" \
             "logs/align/${pn}_bowtie2_genome.log"; do
        [ -s "$f" ] && echo "  log: $f" || echo "  WARNING: missing or empty $f"
    done
}

# A zero exit code is not proof: Snakemake can exit 0 having skipped a rule.
outputs_ok() {
    local pn="$1" ok=0 f sz
    for f in "R_output/${pn}_final_big.rds" "R_output/${pn}_final_small.rds" \
             "R_output/${pn}_deduplication.csv"; do
        if [ -s "$f" ]; then
            echo "  present: $f ($(du -h "$f" | cut -f1))"
        else
            echo "  MISSING/EMPTY: $f"; ok=1
        fi
    done
    sz=$(stat -c %s "R_output/${pn}_final_big.rds" 2>/dev/null || echo 0)
    if [ "$sz" -le 10000000 ]; then
        echo "  final_big.rds implausibly small (${sz} bytes)"; ok=1
    fi
    return $ok
}

echo "=== driver started: $(date -Is)  host: $(hostname)  cores: $CORES ==="
echo "=== commit: $(git rev-parse --short HEAD) ==="
echo "=== conda env: ${CONDA_DEFAULT_ENV} ==="
df -h "$REPO" | tail -1

while IFS=$'\t' read -r out pn geno rep idx r1 r2; do
    [ "$out" = "output_name" ] && continue
    [ -z "${pn:-}" ] && continue
    if [ -n "${ONLY:-}" ] && [ "$out" != "$ONLY" ] && [ "$pn" != "$ONLY" ]; then
        continue
    fi

    # Resume: skip anything already finished with exit 0.
    if awk -F'\t' -v p="$pn" 'NR>1 && $1==p && $5==0 {found=1} END{exit !found}' "$STATUS"; then
        echo "[$pn] already completed, skipping"
        continue
    fi

    echo ""
    echo "=== [$pn]  genotype: $geno  planned label: $out  $(date -Is) ==="

    sed -e "s|__SAMPLE__|${pn}|" -e "s|__R1__|${r1}|" -e "s|__R2__|${r2}|" \
        config/config_template.yaml > "config/generated/${pn}.yaml"

    runlog="logs/run_${pn}_$(date +%Y%m%d_%H%M%S).log"
    t0=$(date -Is)
    # No --notemp: Snakemake deletes its own temp() intermediates as it goes.
    snakemake \
        --configfile "config/generated/${pn}.yaml" \
        --cores "$CORES" \
        --printshellcmds \
        --rerun-incomplete \
        --default-resources "tmpdir='$REPO/tmp'" \
        > "$runlog" 2>&1
    rc=$?
    t1=$(date -Is)
    echo "  snakemake exit=$rc  log=$runlog"

    collect_align_logs "$pn"

    if [ $rc -eq 0 ] && outputs_ok "$pn"; then
        note=ok
    else
        note="FAILED rc=$rc - see $runlog"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$pn" "$out" "$t0" "$t1" "$rc" "$note" >> "$STATUS"
    echo "  status: $note"
    df -h "$REPO" | tail -1
done < config/samples.tsv

echo ""
echo "=== driver finished: $(date -Is) ==="
column -t "$STATUS"
