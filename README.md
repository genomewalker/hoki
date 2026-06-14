# hoki

Position-centric inverted index over diamond blastx codon observations, keyed by
OMA HOG. Answers `(HOG, hog_pos[, AA]) → {acc_id, codon, unitig}` in O(log n_hogs).

**Deps**: C++17, zstd ≥1.4, cmake ≥3.16.

```
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)
```

---

## Pipeline

```
hoki convert -a ACC in.tsv out.lhb     # one per accession, fully parallel
hoki merge out.lhg out.lhgi *.lhb      # single inversion pass
hoki saav merged.lhg merged.lhgi HOG POS [AA]
hoki freq merged.lhg merged.lhgi HOG
hoki stat file.lhb | merged.lhg [merged.lhgi]
```

---

## `hoki convert`

```
hoki convert -a ACC [-z LVL=3] [-p MINPID=0] [-e MAXEV=1e-3] [-v] in.tsv out.lhb
```

Input: 14-col diamond blastx TSV (`qseqid qstart qend qlen qstrand sseqid sstart
send slen pident evalue cigar qseq_translated full_qseq`). `sseqid` must be an OMA
HOG ID. `full_qseq` is consumed for ACGT validation and discarded.

Skips: pident=100 (no variant), non-ACGT codon triplets, round-trip failures
`codon_to_aa(pack(nt3)) != observed_aa`.

Writes one ShardBlock per HOG into `.lhb`, compressed at level `LVL`.

## `hoki merge`

```
hoki merge out.lhg out.lhgi batch1.lhb [batch2.lhb ...]
```

Pass 1: scan all `.lhb` → collect `BatchBlockRef` (offsets, no payload).  
Between passes: sort all unique acc\_ids → `acc_idx` (alphabetical rank).  
Pass 2: per HOG group — decompress each block, invert to
`hog_pos → [(acc_idx, codon_idx, contig_num)]`, sort positions/obs, serialize
inverted block, zstd-19 compress, write `LHHE` entry.  
Pass 3: write HOG index + accession registry → `.lhgi`; backfill `.lhg` header.

## `hoki saav`

```
hoki saav lhg lhgi HOG_ID POS [AA]
```

TSV stdout: `acc_id \t unitig_id \t hog_pos \t obs_aa \t codon`

`POS` is 0-based `sstart`. `AA` is a single-letter filter. `unitig_id` = source
contig (`qseqid`) for flanking-residue verification.

## `hoki freq`

```
hoki freq lhg lhgi HOG_ID
```

TSV stdout: `hog_pos \t codon \t obs_aa \t n_accessions`

Counts distinct `acc_idx` per codon (multiple unitigs from the same accession at
the same position count once).

---

## Wire formats

### `.lhb` — batch file (SHARD\_BLOCK\_VERSION 4)

```
File header:
  [4]  "LHBB"
  [1]  version = 4
  [1]  flags   = 0
  [2]  pad
  varint(len) + acc_id

Per-HOG ShardBlock:
  ShardBlockHeader (28 bytes, packed):
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
    n_contigs × (varint(len) + bytes)    // qseqid strings
    VarNT block:
      header: 6 varints
        n_records / contig_col_bytes / sstart_col_bytes /
        span_col_bytes / bmp_bytes(=n_records) / n_codons
      contig_idx  [n_records × varint]   // index into local contig list
      sstart      [n_records × varint]   // delta-encoded; first absolute
      span        [n_records × varint]   // send − sstart + 1
      bitmap      [n_records × 1 byte]   // per-record obs bitmask
      codons      [n_codons  × 1 byte]   // codon_idx 0-63 = (nt0<<4)|(nt1<<2)|nt2
```

Codon encoding: A=0 C=1 G=2 T=3; `codon_idx = (nt0 << 4) | (nt1 << 2) | nt2`.
`obs_aa` is derived at query time via `codon_to_aa(codon_idx << 2)`.

### `.lhg` — global inverted index (LHG\_VERSION 4)

```
File header (16 bytes):
  [4]  "LHGG"
  [1]  version = 4
  [1]  flags   = 0
  [2]  pad
  [8]  index_offset   uint64 LE  // byte offset of LHGI section

Per-HOG entry (sorted lex by HOG ID):
  [4]  "LHHE"
  [4]  stored_sz  uint32 LE
         bit 31 = 0  → zstd frame;  payload size = stored_sz
         bit 31 = 1  → raw;         payload size = stored_sz & 0x7FFFFFFF
  payload (size bytes)

Decompressed HOG payload:
  varint n_unitigs
  n_unitigs × varint(contig_num)        // uint32; reconstruct as acc_id + "_" + contig_num
  varint n_positions
  per position (ascending hog_pos):
    varint hog_pos_delta                // delta from prev; first = absolute
    varint n_obs
    acc_idx    [n_obs × delta-varint]   // sorted ascending; first absolute
    codon_idx  [n_obs × 1 byte]         // 0-63
    unitig_idx [n_obs × varint]         // index into local unitig dict above

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

Raw fallback: singleton HOGs whose payload would expand under zstd are stored
uncompressed (bit 31 of `stored_sz` set). Saves the ~10-byte zstd frame overhead
on the ~90% singleton HOGs typical at <10 accessions.

### `.lhgi` — standalone index

Same bytes as the LHGI section embedded at `index_offset` in `.lhg`. Load once;
range-GET `.lhg[data_offset : data_offset+data_length]` per HOG.

---

## Varint

Unsigned LEB128: 7 bits/byte, LSB first, high bit = continuation.

---

## Compression (5 barley SRA accessions, 39 357 HOGs)

| | bytes |
|---|---|
| TSV raw | 85 MB |
| TSV −full\_qseq | 40 MB |
| `.lhb` sum | 15 MB |
| `.lhg` | **14 MB** |
| `.lhgi` | 1.3 MB |

vs TSV raw −83.5%; vs TSV−full\_qseq −64.5%.  
At ≥16k accessions: cross-accession codon redundancy in HOG-level zstd-19 drives
`codon_idx` column to near-zero entropy; expected >90% vs TSV.
