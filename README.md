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

### Large-scale (S3 / NFS — fused `ingest` path, no intermediate .lhb)

`ingest` decodes the TSV **once** and writes a self-contained *spill dir*: per-worker
compressed spill buckets plus the HOG and accession registries. `merge-shard` detects the
spill dir and builds the `.lhg` straight from the buckets — it does **not** re-decode or
re-scatter (the old partition path's pass 1). The input is decoded a single time across the
whole pipeline.

```bash
# Phase 1: one ingest job per shard's TSV (GNU parallel, Slurm, AWS Batch, …).
# Each writes an isolated spill dir.
parallel -j 32 'hoki ingest -a auto -t 24 {} spill/{/.}/' ::: inputs/*.tsv

# S3 example (stdin supported):
aws s3 cp s3://serratus-rayan/beetles/logan_jun9_26_run/diamond/DRR000001/DRR000001.diamond.jun9_26.txt - \
    | hoki ingest -a auto - spill/DRR000001/

# Phase 2: build each shard's .lhg from its spill dir (no decode/scatter pass)
hoki merge-shard -t 24 shard.lhg shard.lhgi spill/DRR000001/

# Phase 3: combine shards. Prefer one cross-shard index over N per-shard .lhg:
hoki merge-index global.lhgx shards/      # dir of per-shard .lhg (or list them)
```

Peak RSS for both `ingest` and `merge-shard` is bounded by `--flush` (default = 70% of the
cgroup/SLURM memory limit) regardless of shard size. `--hog-range START END` (partition
path only) restricts HOGs processed for cluster jobs splitting the HOG list across nodes.

The legacy `convert | partition | merge-shard` path (below) still exists and `merge-shard`
still accepts old partition dirs; only `ingest`'s own output changed to the fused spill dir.

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
hoki ingest -a ACC|auto [-t N] [-z LVL=3] [-p MINPID=0] [-e MAXEV=1.0] [-v] \
            [--flush N[KMGT]] [--buckets B] in.tsv out_dir/
```

TSV → fused **spill dir** in a single decode pass — no intermediate `.lhb`, no partition.
`-t N` runs N parallel workers. `-a auto` extracts the accession from the `qseqid` prefix
(before the first `_`), so one TSV containing many accessions produces a correctly keyed dir.

Each worker emits its observations as spill records bucketed by `hash(hog_name) % B` and
zstd-compresses them on the (otherwise decode-idle) worker cores. The dir contains:

- `tN.bucket.b` — worker `N`'s compressed spill for bucket `b` (`B` buckets per worker)
- `tN.hog.registry` — worker `N`'s local-HOG-id → name table
- `acc.registry` — global accession registry (LHGA)
- `spill.meta` — bucket count `B` and worker count `N`

`B` is auto-sized from the input size and the flush budget so each bucket fits `merge-shard`'s
per-bucket working set; override with `--buckets`. The spill is the only intermediate and is
consumed by `merge-shard`.

**Memory management** — peak RSS is bounded by `--flush` regardless of input size. The
threshold is auto-detected from the cgroup/SLURM memory limit and set to 70% of it; override
with `--flush`:

```bash
hoki ingest -a auto -t 24 in.tsv out/             # auto: 70% of cgroup/SLURM limit
hoki ingest -a auto -t 24 --flush 24G in.tsv out/ # explicit 24 GiB
```

## `hoki merge-shard`

```
hoki merge-shard [-t N=nproc] [--flush N[KMGT]] [--hog-range START END] \
                 out.lhg out.lhgi dir/ [dir2/ ...|-]
```

Builds a shard `.lhg` from `ingest` (fused spill dir) or `partition`/`convert` (partition dir).
The path is chosen by the input dir:

- **Fused spill dir** (has `spill.meta`): reads the per-worker compressed buckets, unions the
  worker HOG registries into global ids, and builds each HOG's inverted block directly — **no
  decode/scatter pass**, since `ingest` already produced the spill. This is the fast path
  (~15–28% faster and lighter on memory than the partition path: the uncompressed-spill
  round-trip the partition path does internally is gone).
- **Partition dir(s)** (have `partition.idx`): the legacy path — decode each `tN.lhp` frame
  once, scatter records to bucket files (pass 1), then build (pass 2). Accepts one or more
  dirs and reconciles their accession registries.

Calling conventions (partition path, all avoiding shell ARG\_MAX):

```bash
hoki merge-shard -t 192 out.lhg out.lhgi parts/        # parent dir — opendir internally
cat manifest.txt | hoki merge-shard -t 192 out.lhg out.lhgi -   # stdin, one dir per line
hoki merge-shard -t 4 out.lhg out.lhgi parts/A/ parts/B/        # explicit list
```

Both paths bound peak RSS by `--flush`. HOG payloads are written in completion order; a final
sort orders index entries by HOG ID before writing LHGI. `--hog-range START END` (partition
path) restricts processing to `[START, END]` for cluster jobs splitting the HOG list.

## `hoki merge-index`

```
hoki merge-index out.lhgx shard1.lhg [shard2.lhg ...|shard_dir/]
```

Combines N per-shard `.lhg` indexes into one `.lhgx` — **the recommended way to query across
shards.** It reads only each shard's index (not the data), so it is fast and cheap regardless
of shard size. `saav`/`freq` on the `.lhgx` read blocks from each shard's `.lhg` directly and
union the results; accessions are disjoint across shards, so the union is a concat. This
replaces materializing one giant merged `.lhg` with `merge`.

```bash
hoki merge-index global.lhgx shards/        # dir of per-shard .lhg
hoki freq global.lhgx N0.HOG0001527         # query across all shards
```

## `hoki merge` (deprecated — prefer `merge-index`)

Materializes one merged `.lhg` from N inputs. At scale prefer `merge-index` (no data merge).
Retained for `.lhb` consolidation and small-N cases.

```
hoki merge [-t N] [-z LVL=3] [-zo LVL=3] [--buckets N=1]
           [--acc-registry file.acc] [--hog-range START END]
           [--profile] [--hot-threshold N=100]
           out.lhg out.lhgi input1 [input2 ...]
```

Inputs may be `.lhb` or `.lhg`, mixed (a single dir arg scans for `*.lhg`). K-way heap-merges
the per-input sorted HOG lists, re-inverts each HOG group, and writes v9 blocks largest-first.
`--buckets N` writes N sharded `.lhg`/`.lhgi` pairs; `--hot-threshold N` sets the min block
count for parallel decode of one HOG. Peak RSS bounded by `--flush`.

## `hoki saav`

```
hoki saav lhg lhgi HOG_ID POS [AA] [--min-pident N]   # single index
hoki saav index.lhgx HOG_ID POS [AA] [--min-pident N] # cross-shard (union over shards)
```

TSV stdout: `acc_id \t unitig_id \t hog_pos \t obs_aa \t codon \t pident`

`POS` is 0-based `sstart`. `AA` is a single-letter filter.

`unitig_id` is `acc_id_N` where `N` is the numeric field immediately after the first
`_` in `qseqid` (e.g. `DRR000001_42_ka:f:...` → `unitig_id = DRR000001_42`).

## `hoki freq`

```
hoki freq lhg lhgi HOG_ID [--min-pident N]    # single index
hoki freq index.lhgx HOG_ID [--min-pident N]  # cross-shard (union over shards)
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

### `tN.lhp` — per-thread partition file (LHP_VERSION 4, no file header)

One file per worker thread (`partition`/`convert`). A sequence of zstd frames; each frame
batches many records column-major. Frame boundaries and per-HOG extents live in
`partition.idx` (a HOG's records may span the file).

```
Repeated until EOF:
  [4]  csz  uint32 LE             // compressed frame size
  [csz]  zstd frame

Decompressed frame (column-major):
  [4]  other_sec_len  uint32 LE
  [4]  bmp_sec_len    uint32 LE
  [other_sec_len]  other section   // contig dict + VarNT header + non-bitmap/codon columns
  [bmp_sec_len]    bitmap section
  [remainder]      codon section
```

### `partition.idx` — HOG→extent index (LHPI_VERSION 4)

```
[4]  "LHPI"
[1]  version = 4
[4]  n_threads  uint32 LE
[4]  n_hogs     uint32 LE         // total (thread, hog) entries (a HOG appears once per
                                  // thread that holds it; load accumulates extents per HOG)
per entry:
  [2]  hog_id_len  uint16 LE
  [hog_id_len]  hog_id
  [4]  n_extents   uint32 LE
  n_extents × PartitionIndexExtent (40 bytes):
    [4]  thread_idx   uint32 LE
    [8]  frame_off    uint64 LE    // compressed frame offset in tN.lhp (after its 4-byte csz)
    [4]  frame_csz    uint32 LE
    [4]  other_off    uint32 LE    // byte offsets into the decompressed frame's sections
    [4]  other_len    uint32 LE
    [4]  bmp_off      uint32 LE
    [4]  cdn_off      uint32 LE
    [4]  acc_idx      uint32 LE
    [4]  n_records    uint32 LE
```

### `acc.registry` — standalone accession registry

```
[4]  "LHGA"
[4]  n_accs  uint32 LE
per accession (position = acc_idx):
  varint(len) + acc_id
[4]  adler32  uint32 LE          // over all preceding bytes; checked on load
```

Written by `partition`/`convert` with accessions sorted; written by `ingest` in
first-seen order (`merge-shard` sorts on load for a canonical global acc id). The
per-worker `tN.hog.registry` in a spill dir uses this exact layout (position = local hog id).

### Fused spill dir — `ingest` output consumed by `merge-shard`

A spill dir holds `B` buckets × `N` workers plus registries and `spill.meta`. No
partition index — `merge-shard` builds straight from the buckets.

`spill.meta` — bucket/worker counts:

```
[4]  "LHGM"
[4]  B  uint32 LE     // bucket count
[4]  N  uint32 LE     // worker count
[4]  adler32(prev 12 bytes)
```

`tN.hog.registry` — worker `N`'s local-HOG-id → name table (LHGA byte layout; position =
local hog id). `merge-shard` unions the `N` tables into global ids and a per-worker remap.

`tN.bucket.b` — worker `N`'s spill for bucket `b`, a sequence of independently
decompressable zstd frames (one per ingest flush-drain):

```
Repeated until EOF:
  [4]  csz  uint32 LE          // compressed frame size
  [4]  usz  uint32 LE          // uncompressed size
  [csz]  zstd frame
```

Each decompressed frame concatenates spill records (a HOG appears in exactly one bucket
across all workers, since `b = hash(hog_name) % B`):

```
Repeated:
  [4]  hog_id   uint32 LE      // WORKER-LOCAL hog id (→ global via tN.hog.registry)
  [4]  acc_idx  uint32 LE      // global accession index
  [4]  cnum     uint32 LE      // unitig number (field after first '_' in qseqid)
  [4]  sstart   uint32 LE
  [4]  send     uint32 LE
  [1]  pident_u8                // min(100, pident + 0.5)
  [4]  n_obs    uint32 LE
  n_obs × (uint32 hog_offset, uint8 codon_idx)   // offset from sstart; codon 0-63
```

### `.lhg` — global inverted index (LHG_VERSION 8; v9 blocks)

```
File header (16 bytes):
  [4]  "LHGG"
  [1]  version = 8
  [1]  flags   = 0
  [2]  pad
  [8]  index_offset   uint64 LE  // byte offset of LHGI section from file start

Per-HOG entry (sorted lex by HOG ID):
  [4]  block magic   "LHH2" = v9 (current writer) | "LHHE" = v8 (legacy; reader handles both)
  [4]  stored_sz  uint32 LE
         bit 31 = 1  → raw (uncompressed); payload = stored_sz & 0x7FFFFFFF bytes
         bit 30 = 1  → zstd frame
  payload          // true byte length is data_length in the index (may exceed 30 bits)
```

A unitig — one assembly contig `(acc_idx, cnum)` — is constant across its ~150 positions, so
v9 replaces the v8 per-observation acc and cnum columns with a **unitig trailer** + a single
**Δ-uid column** (the uid is an ordinal into the trailer). ~60% smaller blocks, lossless.

```
Decompressed v9 payload:
  [local acc dict]
    varint n_local
    n_local × (varint delta_gacc, uint8 pident_u8)   // sorted global acc_idx; min pident

  [coverage]
    varint hog_length                // protein length in AA (0 = unknown)
    n_local × varint covered_aa      // AA positions covered per local acc

  [unitig trailer]                   // distinct (local_acc, cnum), sorted
    varint n_uids
    n_uids × (varint delta_local_acc, varint cnum)   // delta on local_acc; cnum = unitig number

  [position headers]
    varint n_positions
    n_positions × (varint pos_delta, varint n_obs)   // pos_delta from previous; first absolute

  [uid column — all positions concatenated]
    total_obs × varint delta_uid_ordinal   // ordinal into the trailer; delta within each
                                            // position (sorted by uid), reset at each position

  [codon column — per position]
    for each position:
      uint8 consensus_codon_idx        // 0-63, the most frequent codon
      varint n_var                     // count of observations differing from consensus
      n_var × varint delta_ordinal     // ordinals (within the position) of the differing obs
      n_var × uint8 codon_idx          // their codons
```

Index section (at index_offset; also written standalone as .lhgi):
  [4]  "LHGI"
  [4]  n_hogs  uint32 LE
  per HOG (sorted lex by hog_id):
    varint(len) + hog_id
    [8]  data_offset   uint64 LE
    [8]  data_length   uint64 LE   // true on-disk block length (8 + payload)
    [4]  n_accessions  uint32 LE
  [4]  "LHGA"
  [4]  n_accs  uint32 LE
  per accession (sorted; position = acc_idx):
    varint(len) + acc_id
  [4]  adler32  uint32 LE          // over the whole index section; checked on load
```

### `.lhgi` — standalone index

Same bytes as the LHGI section embedded at `index_offset` in `.lhg`. Load once;
range-GET `.lhg[data_offset : data_offset+data_length]` per HOG.

### `.lhgx` — cross-shard combined index (LHGX_VERSION 1)

Maps a HOG to whichever shards contain it, so queries read blocks from each shard's `.lhg`
directly and union the results (accessions are disjoint across shards → union = concat). No
merged `.lhg` is materialized. Written by `merge-index`.

```
[4]  "LHGX"
[1]  version = 1
[4]  n_shards  uint32 LE
per shard:
  [4]  path_len  uint32 LE
  [path_len]  shard .lhg path (string)
  [index bytes]  // that shard's full LHGI section (the same bytes as its .lhgi)
```

---

## Varint

Unsigned LEB128: 7 bits/byte, LSB first, high bit = continuation.
