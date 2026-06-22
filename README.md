# hoki

Position-centric inverted index over diamond blastx codon observations, keyed by OMA HOG.

**Deps**: C++17, zstd ≥1.4, cmake ≥3.16.

```
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)
```

---

## Pipeline

### Small-scale (all `.lhb` files in memory)

```
hoki convert -a ACC in.tsv out.lhb     # one per accession, fully parallel
hoki merge out.lhg out.lhgi *.lhb      # single inversion pass
```

### Large-scale (S3 / NFS — `ingest` path, no intermediate .lhb)

One job per input file; each writes its own isolated partition dir.
`merge-shard` then accepts all dirs at once and merges their acc registries automatically.

```bash
# Phase 1: one ingest job per TSV (parallelise with GNU parallel, Slurm, AWS Batch, …)
parallel -j 32 'hoki ingest -a auto {} parts/{/.}/' ::: inputs/*.tsv

# S3 example (stdin supported):
aws s3 cp s3://serratus-rayan/beetles/logan_jun9_26_run/diamond/DRR000001/DRR000001.diamond.jun9_26.txt - \
    | hoki ingest -a auto - parts/DRR000001/

# Phase 2: per-batch merge — pass the partition dir; opendir internally, no ARG_MAX
hoki merge-shard -t 192 out.lhg out.lhgi parts/

# Phase 3: global merge across all batch .lhg files — pass the dir, same convention
hoki merge -t 192 global.lhg global.lhgi batch_outputs/
```

`--hog-range START END` restricts HOGs processed (for cluster jobs splitting the HOG list across nodes).

### Two-phase path (via .lhb — retained for compatibility)

```bash
# Phase 1: convert
for f in batches/*.tsv.gz; do zcat $f | hoki convert -a auto - lhbs/${f%.tsv.gz}/; done

# Phase 2: partition (pass file list via stdin to avoid ARG_MAX)
find lhbs/ -name '*.lhb' | hoki partition -t 192 part/ -

# Phase 3: merge-shard
hoki merge-shard -t 192 out.lhg out.lhgi part/
```

---

## `hoki convert`

```
hoki convert -a ACC    [-z LVL=3] [-p MINPID=0] [-e MAXEV=1.0] [-v] in.tsv out.lhb
hoki convert -a auto   [-z LVL=3] [-p MINPID=0] [-e MAXEV=1.0] [-v] in.tsv out_dir/
```

Input: 14-col diamond blastx TSV (`qseqid qstart qend qlen qstrand sseqid sstart
send slen pident evalue cigar qseq_translated full_qseq`).

`sseqid` must be an OMA HOG ID (`bpgv2|N0.HOG...` or `panbarley.bpgv2|N0.HOG...`).
`qseqid` identifies the source assembly unitig (Logan format: `ACC_N_ka:f:...`).
`full_qseq` is consumed for ACGT validation and discarded.

**`-a auto`** — multi-accession mode: extracts the accession from the `qseqid` prefix
(everything before the first `_`). `out_dir/` is created automatically; one `ACC.lhb`
file is written per detected accession. Useful when input TSVs are concatenated across
many accessions (e.g. Rayan/Logan batches of thousands of accessions per file).

Skips: non-ACGT codon triplets; round-trip failures `codon_to_aa(pack(nt3)) != observed_aa`.

Writes one ShardBlock per HOG per accession into `.lhb`, zstd-compressed at level `LVL`.

## `hoki partition`

```
hoki partition [-t N] [-z LVL=3] out_dir/ input1.lhb [input2.lhb ...]
```

Reads all `.lhb` inputs and partitions observations by HOG across `N` worker threads.
Each thread writes a single sequential file `out_dir/tN.lhp` (no per-HOG files).
On completion writes `out_dir/partition.idx` (LHPI binary index) and
`out_dir/acc.registry` (LHGA accession registry).

Replaces the old approach that wrote one `.lhp` file per HOG — eliminates the
file-per-HOG metadata storm (catastrophic on NFS/EFS at scale).

## `hoki ingest`

```
hoki ingest -a ACC|auto [-z LVL=3] [-p MINPID=0] [-e MAXEV=1.0] [-v] in.tsv out_dir/
```

TSV → partition dir in a single pass — no intermediate `.lhb`. Equivalent to
`convert | partition` but without the intermediate file and with no ARG\_MAX risk.

`-a auto` extracts the accession from the `qseqid` prefix (before the first `_`), so one
TSV containing many accessions produces a correctly partitioned output dir.

Writes `out_dir/t0.lhp`, `out_dir/partition.idx`, and `out_dir/acc.registry`.
Each parallel job writes to its own isolated output dir; `merge-shard` accepts
multiple dirs and reconciles accession registries automatically.

## `hoki merge-shard`

```
hoki merge-shard [-t N=nproc] [--hog-range START END] \
                 out.lhg out.lhgi part_dir/ [part_dir2/ ...|-]
```

Accepts one or more partition dirs (from `ingest` or `partition`).
Three calling conventions, all avoiding shell ARG\_MAX:

```bash
# 1. Parent dir (recommended) — merge-shard calls opendir internally, no shell expansion
hoki merge-shard -t 192 out.lhg out.lhgi parts/

# 2. Stdin — one dir path per line
cat manifest.txt | hoki merge-shard -t 192 out.lhg out.lhgi -

# 3. Explicit list — fine for small N
hoki merge-shard -t 4 out.lhg out.lhgi parts/A/ parts/B/
```

Loads all `acc.registry` files, builds a merged sorted global accession list,
remaps per-dir local acc indices to global indices on the fly while reading extents.

Opens all `tN.lhp` files across all dirs (flat fd list), then inverts in parallel
(`N` threads, default = hardware concurrency). HOG payloads written in completion order;
final sort pass orders index entries by HOG ID before writing LHGI.

`--hog-range START END` restricts processing to HOGs in `[START, END]` (for parallel
cluster jobs splitting the HOG list across nodes).

## `hoki merge`

```
hoki merge [-t N] [-z LVL=3] [-zo LVL=3] [--buckets N=1]
           [--acc-registry file.acc] [--hog-range START END]
           [--profile] [--hot-threshold N=100]
           out.lhg out.lhgi input1 [input2 ...]
```

Inputs may be `.lhb` or `.lhg`, mixed. If a single directory is given, all `*.lhg`
files in it are collected via `opendir` (no shell expansion, no ARG\_MAX):

```bash
hoki merge -t 192 global.lhg global.lhgi parts/
```

Scan: all inputs scanned in parallel across `-t` threads (worker-stealing). Each input
produces a sorted list of `(hog_id, block_ref)` entries; `.lhb` files are buffered in
RAM (≤ 64 MB each).

Merge: K-way heap merge of the per-file sorted lists into a single HOG-ordered sequence.

Invert + write: for each HOG group — decompress blocks, accumulate observations,
counting-sort by position, serialize column-oriented inverted block, compress, write
`LHHE` entry. Groups are dispatched largest-first to balance threads. A dedicated
writer thread serializes output in HOG-ID order.

`--buckets N` writes N sharded `.lhg`/`.lhgi` pairs.
`--hot-threshold N` sets the minimum block count for parallel decode of a single HOG.

## `hoki saav`

```
hoki saav lhg lhgi HOG_ID POS [AA] [--min-pident N]
```

TSV stdout: `acc_id \t unitig_id \t hog_pos \t obs_aa \t codon \t pident`

`POS` is 0-based `sstart`. `AA` is a single-letter filter.

`unitig_id` is `acc_id_N` where `N` is the numeric field immediately after the first
`_` in `qseqid` (e.g. `DRR000001_42_ka:f:...` → `unitig_id = DRR000001_42`).

## `hoki freq`

```
hoki freq lhg lhgi HOG_ID [--min-pident N]
```

TSV stdout: `hog_pos \t codon \t obs_aa \t n_accessions`

Counts distinct `acc_idx` per codon (multiple unitigs from the same accession at the
same position count once).

---

## Wire formats

### `.lhb` — batch file (LHB_VERSION 1, SHARD_BLOCK_VERSION 4)

```
File header:
  [4]  "LHGB"
  [1]  version = 1
  [1]  flags   = 0
  [2]  pad
  varint(len) + acc_id

Repeated until EOF:
  [4]  "LHBT"  block tag
  varint(len) + hog_id
  ShardBlockHeader (28 bytes):
    [4]  "LHSB"
    [1]  version = 4
    [1]  flags   = 0
    [2]  pad
    [4]  compressed_sz  uint32 LE
    [4]  raw_sz         uint32 LE
    [4]  n_records
    [4]  min_sstart
    [4]  max_send
  zstd payload (compressed_sz bytes) →
    varint n_contigs
    n_contigs × (varint(len) + bytes)        // qseqid strings (unitig IDs)
    VarNT block:
      6 header varints:
        n_records / contig_col_bytes / sstart_col_bytes /
        span_col_bytes / bmp_bytes / n_codons
      contig_idx  [contig_col_bytes]          // varint per record, index into contig list
      sstart      [sstart_col_bytes]          // varint delta per record; first absolute
      span        [span_col_bytes]            // varint per record; send − sstart + 1
      qframe      [n_records × 1 byte]        // query frame
      pident      [n_records × 2 bytes]       // fp16-encoded percent identity
      evalue      [n_records × 2 bytes]       // fp16-encoded e-value
      bitmap      [bmp_bytes total]           // per record: ceil(span/8) bytes, bit i = obs at sstart+i
      codons      [n_codons × 1 byte]         // codon_idx 0-63; one per set bitmap bit, in order

[4]  "LHBE"  EOF sentinel
```

Codon encoding: A=0 C=1 G=2 T=3; `codon_idx = (nt0 << 4) | (nt1 << 2) | nt2`.

### `tN.lhp` — per-thread partition file (LHP_VERSION 2, no file header)

One file per worker thread. Records are appended in arrival order (no HOG sorting
within the file; the index handles lookup).

```
Repeated until EOF:
  [2]  hog_id_len  uint16 LE
  [hog_id_len bytes]  hog_id string
  PartitionEntry (16 bytes, packed):
    [4]  acc_idx        uint32 LE   // global accession index
    [4]  compressed_sz  uint32 LE
    [4]  raw_sz         uint32 LE
    [4]  n_records      uint32 LE
  [compressed_sz bytes]  zstd-compressed ShardBlock payload
```

### `partition.idx` — sorted HOG→extent index (LHPI_VERSION 1)

```
[4]  "LHPI"
[1]  version = 1
[3]  pad
[4]  n_threads  uint32 LE
[4]  n_hogs     uint32 LE
per HOG (sorted lex by hog_id):
  varint(len) + hog_id
  varint n_extents
  n_extents × PartitionIndexExtent (16 bytes):
    [4]  thread_idx    uint32 LE
    [8]  entry_offset  uint64 LE   // byte offset of PartitionEntry within tN.lhp
    [4]  entry_len     uint32 LE   // sizeof(PartitionEntry) + compressed_sz
```

### `acc.registry` — standalone accession registry

```
[4]  "LHGA"
[4]  n_accs  uint32 LE
per accession (sorted; position = acc_idx):
  varint(len) + acc_id
```

### `.lhg` — global inverted index (LHG_VERSION 8)

```
File header (16 bytes):
  [4]  "LHGG"
  [1]  version = 8
  [1]  flags   = 0
  [2]  pad
  [8]  index_offset   uint64 LE  // byte offset of LHGI section from file start

Per-HOG entry (sorted lex by HOG ID):
  [4]  "LHHE"
  [4]  stored_sz  uint32 LE
         bit 31 = 1  → raw (uncompressed); payload = stored_sz & 0x7FFFFFFF bytes
         bit 30 = 1  → zstd frame
  payload (stored_sz & 0x3FFFFFFF bytes)

Decompressed HOG payload:
  [local acc dict]
    varint n_local
    n_local × (varint delta_gacc, uint8 pident_u8)  // sorted global acc_idx; min pident

  [coverage dict]
    varint hog_length                // protein length in AA (0 = unknown)
    n_local × varint covered_aa      // AA positions covered per local acc

  [position headers]
    varint n_positions
    n_positions × (varint pos_delta, varint n_obs)   // pos_delta from previous; first absolute

  [acc column — all positions concatenated]
    total_obs × varint delta_local_acc_idx   // delta from previous within each position; first absolute

  [cnum column — all positions concatenated]
    total_obs × varint cnum          // numeric field after first '_' in qseqid (unitig number)

  [codon column — per position]
    for each position:
      uint8 consensus_codon_idx      // 0-63; 0xFF = no consensus (all explicit)
      if consensus != 0xFF:
        varint n_var
        n_var × varint delta_obs_ordinal   // ordinal of variant obs within this position
        n_var × uint8 codon_idx
      else:
        n_obs × uint8 codon_idx

Index section (at index_offset; also written standalone as .lhgi):
  [4]  "LHGI"
  [4]  n_hogs  uint32 LE
  per HOG (sorted):
    varint(len) + hog_id
    [8]  data_offset   uint64 LE
    [8]  data_length   uint64 LE
    [4]  n_accessions  uint32 LE
  [4]  "LHGA"
  [4]  n_accs  uint32 LE
  per accession (sorted; position = acc_idx):
    varint(len) + acc_id
```

### `.lhgi` — standalone index

Same bytes as the LHGI section embedded at `index_offset` in `.lhg`. Load once;
range-GET `.lhg[data_offset : data_offset+data_length]` per HOG.

---

## Varint

Unsigned LEB128: 7 bits/byte, LSB first, high bit = continuation.
