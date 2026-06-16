# hoki

Position-centric inverted index over diamond blastx codon observations, keyed by OMA HOG.

**Deps**: C++17, zstd ≥1.4, lz4, cmake ≥3.16.

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

### Large-scale (S3 / NFS — two-phase)

```
# Phase 1: partition — one batch per accession, write indexed per-thread files
hoki convert -a ACC in.tsv ACC.lhb
hoki partition out_dir/ *.lhb           # writes t0.lhp…tN.lhp + partition.idx + acc.registry

# Phase 2: merge-shard — parallel inversion directly from the index
hoki merge-shard [-t N] [--hog-range START END] out_dir/ out.lhg out.lhgi
```

S3 example (stdin supported via `-`):
```bash
aws s3 cp s3://bucket/ACC.diamond.tsv - | hoki convert -a ACC - ACC.lhb
```

`--hog-range START END` lets cluster jobs process HOG ranges in parallel (split the
sorted HOG list from `partition.idx` across array jobs, then `hoki accregistry` to merge).

---

## `hoki convert`

```
hoki convert -a ACC [-z LVL=3] [-p MINPID=0] [-e MAXEV=1e-3] [-v] in.tsv out.lhb
```

Input: 14-col diamond blastx TSV (`qseqid qstart qend qlen qstrand sseqid sstart
send slen pident evalue cigar qseq_translated full_qseq`).

`sseqid` must be an OMA HOG ID (`bpgv2|N0.HOG...` or `panbarley.bpgv2|N0.HOG...`).
`qseqid` identifies the source assembly unitig (e.g. Logan: `ACC_N_ka:f:...`).
`full_qseq` is consumed for ACGT validation and discarded.

Skips: non-ACGT codon triplets; round-trip failures `codon_to_aa(pack(nt3)) != observed_aa`.

Writes one ShardBlock per HOG into `.lhb`, zstd-compressed at level `LVL`.

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

## `hoki merge-shard`

```
hoki merge-shard [-t N=nproc] [-z LVL=3] [--hog-range START END] \
                 partition_dir/ out.lhg out.lhgi
```

Reads `partition_dir/partition.idx` to enumerate HOGs, opens all `tN.lhp` files,
then inverts them in parallel using worker-stealing (`N` threads, default = hardware
concurrency).

HOG payloads are written in completion order to `out.lhg`; a final sort pass orders
the index entries by HOG ID before writing the LHGI section.

`--hog-range START END` restricts processing to HOGs in `[START, END)` of the sorted
index (for parallel cluster jobs splitting the HOG list).

## `hoki merge`

```
hoki merge [-t N] [-z LVL=3] [-zo LVL=3] [--buckets N=1]
           [--acc-registry file.acc] [--hog-range START END]
           [--profile] [--hot-threshold N=100]
           out.lhg out.lhgi input1 [input2 ...]
```

Inputs may be `.lhb` or `.lhg`, mixed.

Scan: all inputs scanned in parallel across `-t` threads (worker-stealing). Each input
produces a sorted list of `(hog_id, block_ref)` entries; `.lhb` files are buffered in
RAM (≤ 64 MB each).

Merge: K-way heap merge of the per-file sorted lists into a single HOG-ordered sequence.

Invert + write: for each HOG group — decompress blocks, accumulate observations,
counting-sort by position, serialize column-oriented inverted block, compress, write
`LHHE` entry. Groups are dispatched largest-first to balance threads. A dedicated
writer thread serializes output in HOG-ID order.

`-zo -1` selects LZ4 instead of zstd for output (faster build, ~15% larger file).
`--buckets N` writes N sharded `.lhg`/`.lhgi` pairs.
`--hot-threshold N` sets the minimum block count for parallel decode of a single HOG.

## `hoki accregistry`

```
hoki accregistry out.acc shard.lhgi [shard.lhgi ...]
```

Merges the accession registries from multiple `.lhgi` files into a single `.acc` file.

## `hoki saav`

```
hoki saav lhg lhgi HOG_ID POS [AA] [--min-pident N]
```

TSV stdout: `acc_id \t contig_num \t hog_pos \t obs_aa \t codon`

`POS` is 0-based `sstart`. `AA` is a single-letter filter.

`contig_num` is the numeric field immediately after the first `_` in `qseqid`
(e.g. `DRR000001_42_ka:f:...` → `42`).

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

### `tN.lhp` — per-thread partition file (no file header)

One file per worker thread. Records are appended in arrival order (no HOG sorting
within the file; the index handles lookup).

```
Repeated until EOF:
  [2]  hog_id_len  uint16 LE
  [hog_id_len bytes]  hog_id string
  PartitionEntry (24 bytes, packed):
    [4]  acc_idx        uint32 LE   // global accession index
    [4]  compressed_sz  uint32 LE
    [4]  raw_sz         uint32 LE
    [4]  n_records      uint32 LE
    [4]  min_sstart     uint32 LE
    [4]  max_send       uint32 LE
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
         else        → LZ4 frame (4-byte raw_sz LE prefix)
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
