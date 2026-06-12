# hoki — HOG codon Index

Converts population-scale diamond blastx output (Logan/SRA assemblies × OMA HOG proteins)
into a compact, immediately-queryable container directory. Designed for 38 M parallel writers
with no merge step: each HOG shard is an independent append-only file.

---

## Quick start

```bash
# Build
cmake -B build -S . && cmake --build build -j

# Convert one accession's diamond TSV to a container
hoki convert DRR000001.diamond.tsv /scratch/logan.lhc

# Many accessions writing concurrently to the same container — safe with flock
parallel hoki convert {}.diamond.tsv /scratch/logan.lhc ::: $(cat accessions.txt)

# Query all codon-level observations for a HOG
hoki varnt /scratch/logan.lhc chr1H01393

# Query a specific position
hoki varnt /scratch/logan.lhc chr1H01393 42

# Find candidate intron gaps across accessions
hoki tile /scratch/logan.lhc chr1H01393 -n 10
```

---

## Diamond input format

hoki expects **diamond blastx outfmt 6** with these columns (0-based):

| # | Name | Description |
|---|------|-------------|
| 0 | qseqid | query contig ID (encodes SRA accession prefix) |
| 1 | qstart | 1-based nt start in contig |
| 2 | qend | nt end |
| 3 | qlen | total contig nt length |
| 4 | qstrand | `+` or `-` |
| 5 | sseqid | subject protein ID — HOG extracted after last `\|` |
| 6 | sstart | 1-based HOG protein position (coverage start) |
| 7 | send | HOG protein position (coverage end) |
| 8 | slen | HOG protein length |
| 9 | pident | percent identity |
| 10 | evalue | E-value |
| 11 | cigar | CIGAR string (M/I/D operations) |
| 12 | qseq_translated | translated query sequence (AA) |
| 13 | full_qseq | full nucleotide contig (used to extract codons) |

Alignments at **100% pident** are skipped — at exact identity every codon encodes
the reference amino acid and carries no variant information.

---

## Container format (`.lhc`)

A container is a **directory** containing one shard file per HOG:

```
logan.lhc/
  shards/
    chr1H00001.lhs    ← all observations for HOG chr1H00001
    chr1H00002.lhs
    ...
    HOG0044123.lhs
```

The shard path is derived deterministically from the HOG ID string
(characters `/`, `:`, ` ` replaced with `_`, `.lhs` appended), so any writer
can locate the correct shard without consulting a global index.

### Shard file (`.lhs`)

A shard is a **sequence of independently-appended shard blocks**:

```
[ShardBlock 0][ShardBlock 1][ShardBlock 2]...
```

Each block is written atomically under `flock(LOCK_EX)`, so multiple
concurrent writers targeting the same HOG shard are safe on NFS and local FS.

### Shard block wire format

```
Offset  Size  Field
0       4     magic: "LHSB"
4       4     compressed_sz   (uint32 LE)  — zstd payload size in bytes
8       4     raw_sz          (uint32 LE)  — decompressed size
12      4     n_records       (uint32 LE)
16      4     min_sstart      (uint32 LE)  — smallest HOG coverage start in block
20      4     max_send        (uint32 LE)  — largest HOG coverage end
24      …     zstd-compressed payload
```

Total header: **24 bytes**.

### Shard block payload (after decompression)

```
varint   n_contigs
for each contig:
  varint  len
  len bytes  (UTF-8 contig ID string, NOT null-terminated)

VarNT columnar block:
  varint   n_records
  varint   hdr_section_bytes
  varint   bmp_section_bytes
  varint   cdn_section_bytes
  [hdr section]  — n_records × record header (see below)
  [bmp section]  — n_records × presence bitmap
  [cdn section]  — packed codon bytes for all M-positions
```

`contig_idx` values inside the VarNT block reference the **local** contig list
prepended in this payload (0-based). Each block is self-describing.

### VarNT record header (per record in hdr section)

```
varint   contig_idx   — index into local contig list
varint   hog_idx      — always 0 for single-HOG shards (present for format symmetry)
uint24   sstart       — 1-based HOG position, coverage start
uint24   send         — 1-based HOG position, coverage end
uint32   qstart       — nt start in query contig
uint32   qend
uint32   qlen
int8     qframe       — ±1..±3 (sign=strand, magnitude=reading frame)
float32  pident
float64  evalue
varint   n_M          — number of M-position (non-gap, non-unk) codons in this record
```

### VarNT presence bitmap (per record in bmp section)

`ceil((send − sstart + 1) / 8)` bytes, LSB-first.  
Bit `i` = 1 → HOG position `(sstart + i)` has an observed codon.

### VarNT codon stream (cdn section)

One byte per set bit in the bitmaps, in ascending bit-position order across
all records. Each byte encodes one codon:

```
bits [7:6]  nt0  (A=0, C=1, G=2, T=3)
bits [5:4]  nt1
bits [3:2]  nt2
bits [1:0]  unused
```

The observed amino acid is **not stored**; it is derived at query time via the
standard genetic code lookup `CODON_TO_AA[packed_codon >> 2]` (see `format.hpp`).

---

## Querying

### `hoki varnt`

```
hoki varnt [-p MINPID] [-e MAXEV] <container.lhc> <HOG_ID> [POS [OBS_AA]]
```

Scans all shard blocks for the given HOG and emits one TSV row per codon observation:

```
contig_id   hog_pos   obs_aa   codon   pident   evalue
```

`POS` (1-based HOG position) and `OBS_AA` (single letter) filter the output.

### `hoki tile`

```
hoki tile [-p MINPID] [-e MAXEV] [-n MIN_ACC] <container.lhc> <HOG_ID>
```

Detects candidate intron/structural gaps: positions in a HOG where a single accession's
contigs cover flanking regions but leave a gap in between.  `-n` sets the minimum number
of accessions showing the same gap to report it.

Output:
1. Per-accession gap table: `accession  gap_start  gap_end  gap_aa_size  n_exon_records`
2. Gap summary: positions seen in ≥ N accessions

---

## Reading from C++

```cpp
#include "container.hpp"

// Iterate all blocks for one HOG
lhi::read_shard_file("logan.lhc/shards/chr1H01393.lhs",
    [](const std::vector<std::string>& contigs,
       const std::vector<lhi::VarNTRecord>& recs) {
        for (auto& r : recs) {
            for (auto& v : r.vars) {
                uint32_t hog_pos = r.sstart + v.hog_offset;
                char aa = lhi::AA_ALPHA[v.obs_aa];
                char codon[3]; lhi::unpack_codon(v.packed_codon, codon);
                // contig string: contigs[r.contig_idx]
            }
        }
    });
```

The shard path for any HOG ID is:
```cpp
auto shard = std::filesystem::path(container_dir) / "shards"
           / lhi::hog_to_filename(hog_id);
```

---

## Scale

| Dataset | Records | Shards | Size |
|---------|---------|--------|------|
| 1 accession (DRR000001) | ~270 | ~230 HOGs | ~300 KB |
| Logan full (38 M accessions) | ~5 × 10¹⁰ | ~50 K HOGs | ~2.5 TB est. |

At full scale each HOG shard is ~50 MB on average, and a HOG query reads only
that one file regardless of how many accessions are in the container.

---

## Build dependencies

- C++20 compiler (GCC ≥ 11 or Clang ≥ 14)
- zstd ≥ 1.4
- xxhash (for CIGAR codon extraction, header-only path)

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```
