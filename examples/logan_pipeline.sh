#!/usr/bin/env bash
# hoki on Logan/serratus DIAMOND output.
#
# Pipeline: DIAMOND blastx TSV  ->  ingest (fused spill dir)  ->  merge-shard (per-shard .lhg)
#           ->  merge-index (.lhgx over all shards)  ->  query (freq / saav).
#
# A "shard" = one DIAMOND file. A file may concatenate many accessions; `ingest -a auto`
# keys them apart by the qseqid prefix, so one file -> one shard .lhg.
set -euo pipefail

HOKI=${HOKI:-hoki}
BUCKET=${BUCKET:-s3://serratus-rayan/beetles/logan_jun9_26_run/diamond}
WORK=${WORK:-./hoki_run}
THREADS=${THREADS:-24}
FLUSH=${FLUSH:-24G}        # peak-RSS target per stage (default unset = 70% of cgroup/SLURM mem)
JOBS=${JOBS:-8}           # shards built concurrently
ACCS=${1:?usage: $0 accs.txt   # one shard/accession id per line}

mkdir -p "$WORK/spill" "$WORK/lhg"

# ── Phase 0 (upstream, for reference): how to produce the TSV hoki expects ──
# 14 columns, in this order:
#   qseqid qstart qend qlen qstrand sseqid sstart send slen pident evalue cigar qseq_translated full_qseq
#   qseqid = Logan unitig, "ACC_N_ka:f:..."  (accession = prefix before the first '_')
#   sseqid = OMA HOG id,    "...|N0.HOG...."
#
#   diamond blastx -d hogs.dmnd -q reads.fa --outfmt 6 \
#     qseqid qstart qend qlen qstrand sseqid sstart send slen pident evalue cigar \
#     qseq_translated full_qseq

# ── Phase 1+2: per shard, stream from S3 -> fused ingest -> merge-shard -> shard .lhg ──
build_shard() {
  local id=$1
  local sp="$WORK/spill/$id"
  aws s3 cp "$BUCKET/$id/$id.diamond.jun9_26.txt" - \
    | "$HOKI" ingest -a auto -t "$THREADS" --flush "$FLUSH" - "$sp"
  "$HOKI" merge-shard -t "$THREADS" --flush "$FLUSH" \
    "$WORK/lhg/$id.lhg" "$WORK/lhg/$id.lhgi" "$sp"
  rm -rf "$sp"              # spill is consumed by merge-shard; drop to reclaim disk
}
export -f build_shard
export HOKI BUCKET WORK THREADS FLUSH

parallel -j "$JOBS" build_shard :::: "$ACCS"
# SLURM equivalent: one array task per line of $ACCS, each calling build_shard "$ID".

# ── Phase 3: one cross-shard index over every per-shard .lhg (no data merge) ──
"$HOKI" merge-index "$WORK/global.lhgx" "$WORK/lhg/"

# ── Query across all shards via the .lhgx ──
#   $HOKI freq "$WORK/global.lhgx" N0.HOG0001527
#   $HOKI saav "$WORK/global.lhgx" N0.HOG0001527 100 --min-pident 90
echo "done: $WORK/global.lhgx  ($(ls "$WORK"/lhg/*.lhg | wc -l) shards)"
