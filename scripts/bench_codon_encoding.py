#!/usr/bin/env python3
"""
Benchmark compression strategies on real .lhp data (v3 and v4 formats).

v3 frame (decompressed): concatenated raw VarNT blocks [contig_dict|header|cols...]
v4 frame (decompressed): [other_sec_len(u32)|bmp_sec_len(u32)|other_bytes|bmp_bytes|cdn_bytes]

Strategies:
  v4_curr_l3/l6/l9   current v4 column-major at different zstd levels
  v4_lframe_l3        larger FRAME_RAW_TARGET (4MB) — more context for compressor
  v4_dict_l3          zstd dictionary trained on bitmaps (first 100 frames)
  v4_dict_bmp_l3      dict trained only on bitmap section
  v4_rle_cdn_l3       RLE-encode codon column before zstd

Usage:
    python3 bench_codon_encoding.py <lhp_dir> [max_frames]
"""
import sys, os, struct, collections, math, io, random
import zstd
import zstandard as zstd_dict  # for dictionary training and dict-compressed benchmarks

def read_varint(buf, off):
    v, shift = 0, 0
    while True:
        b = buf[off]; off += 1
        v |= (b & 0x7F) << shift
        if not (b & 0x80): return v, off
        shift += 7

def write_varint(v):
    out = []
    while v >= 128:
        out.append((v & 0x7F) | 0x80)
        v >>= 7
    out.append(v)
    return bytes(out)

def skip_contig_dict(data, off):
    n, off = read_varint(data, off)
    for _ in range(n):
        nlen, off = read_varint(data, off)
        off += nlen
    return off

def parse_block(data, off):
    """Returns (bmp_bytes, cdn_bytes, spans, new_off).
    bmp_bytes: raw dense bitmap bytes for all records concatenated.
    cdn_bytes: raw codon bytes (1 per observation).
    spans: list of span values per record.
    """
    off = skip_contig_dict(data, off)
    n_recs,   off = read_varint(data, off)
    contig_b, off = read_varint(data, off)
    sstart_b, off = read_varint(data, off)
    span_b,   off = read_varint(data, off)
    bmp_b,    off = read_varint(data, off)
    cdn_b,    off = read_varint(data, off)

    col_start = off
    contig_end = col_start + contig_b
    sstart_end = contig_end + sstart_b

    # decode spans
    p = sstart_end
    spans = []
    for _ in range(n_recs):
        sp, p = read_varint(data, p)
        spans.append(sp)

    qframe_end = p
    evalue_end = qframe_end + n_recs + n_recs*2 + n_recs*2   # qframe+pident+evalue
    bmp_end    = evalue_end + bmp_b
    cdn_end    = bmp_end + cdn_b

    bmp_bytes = data[evalue_end : bmp_end]
    cdn_bytes = data[bmp_end    : cdn_end]
    return bmp_bytes, cdn_bytes, spans, cdn_end

def encode_bitpack6(codons: bytes) -> bytes:
    """Pack 6-bit codon values into a byte stream (4 codons per 3 bytes)."""
    buf = bytearray()
    n = len(codons)
    i = 0
    while i + 3 < n:
        a, b, c, d = codons[i], codons[i+1], codons[i+2], codons[i+3]
        buf.append((a << 2) | (b >> 4))
        buf.append(((b & 0xF) << 4) | (c >> 2))
        buf.append(((c & 0x3) << 6) | d)
        i += 4
    # tail: leftover 1-3 codons packed as-is into a partial triplet
    rem = n - i
    if rem == 1:
        buf.append(codons[i] << 2)
    elif rem == 2:
        a, b = codons[i], codons[i+1]
        buf.append((a << 2) | (b >> 4))
        buf.append((b & 0xF) << 4)
    elif rem == 3:
        a, b, c = codons[i], codons[i+1], codons[i+2]
        buf.append((a << 2) | (b >> 4))
        buf.append(((b & 0xF) << 4) | (c >> 2))
        buf.append((c & 0x3) << 6)
    return bytes(buf)

def encode_gapcoded(bmp_bytes: bytes, cdn_bytes: bytes, spans: list, bitpack: bool) -> bytes:
    """Replace dense bitmap with gap-coded offsets; optionally 6-bit pack codons."""
    gaps = bytearray()
    cdn_ptr = 0
    new_cdn = bytearray()
    for i, span in enumerate(spans):
        bmp_off = sum((spans[j]+7)//8 for j in range(i))
        bmp_sz  = (span + 7) // 8
        bmp     = bmp_bytes[bmp_off : bmp_off + bmp_sz]
        prev    = 0
        for byte_i, byte in enumerate(bmp):
            for bit in range(8):
                pos = byte_i * 8 + bit
                if pos >= span: break
                if byte & (1 << bit):
                    gaps += write_varint(pos - prev)
                    prev = pos
                    new_cdn.append(cdn_bytes[cdn_ptr])
                    cdn_ptr += 1
        gaps += write_varint(span - prev + span)  # sentinel: delta > span signals end-of-record
    codon_payload = bytes(new_cdn)
    if bitpack:
        codon_payload = encode_bitpack6(codon_payload)
    return gaps + codon_payload

def compress_lv(buf: bytes, level: int) -> int:
    return len(zstd.compress(buf, level))

def reorder_blocks(raw: bytes, key_fn) -> bytes:
    """Parse all blocks in a frame, sort by key_fn(bmp_bytes, cdn_bytes, spans), re-concatenate."""
    blocks = []
    off = 0
    block_starts = []
    while off < len(raw):
        start = off
        bmp_bytes, cdn_bytes, spans, off = parse_block(raw, off)
        block_bytes = raw[start:off]
        blocks.append((key_fn(bmp_bytes, cdn_bytes, spans), block_bytes))
    blocks.sort(key=lambda x: x[0])
    return b''.join(b for _, b in blocks)

def col_major(raw: bytes) -> bytes:
    """Restructure frame: concat all other-column bytes, then all bitmaps, then all codons."""
    other_parts, bmp_parts, cdn_parts = [], [], []
    off = 0
    while off < len(raw):
        start = off
        bmp_bytes, cdn_bytes, spans, off = parse_block(raw, off)
        other_size = (off - len(cdn_bytes) - len(bmp_bytes)) - start
        other_parts.append(raw[start : start + other_size])
        bmp_parts.append(bmp_bytes)
        cdn_parts.append(cdn_bytes)
    return b''.join(other_parts) + b''.join(bmp_parts) + b''.join(cdn_parts)

# ── v4 frame parser ──────────────────────────────────────────────────────────

def is_v4_frame(raw: bytes) -> bool:
    """Heuristic: v4 frames start with two u32s (other_sec_len, bmp_sec_len)
    that together must be < total frame size and sum to frame_size - 8."""
    if len(raw) < 8: return False
    other_len = struct.unpack_from('<I', raw, 0)[0]
    bmp_len   = struct.unpack_from('<I', raw, 4)[0]
    return 8 + other_len + bmp_len <= len(raw)

def parse_v4_frame(raw: bytes):
    """Parse v4 column-major frame into list of (other_bytes, bmp_bytes, cdn_bytes, spans)."""
    other_sec_len = struct.unpack_from('<I', raw, 0)[0]
    bmp_sec_len   = struct.unpack_from('<I', raw, 4)[0]
    other_sec = raw[8 : 8 + other_sec_len]
    bmp_sec   = raw[8 + other_sec_len : 8 + other_sec_len + bmp_sec_len]
    cdn_sec   = raw[8 + other_sec_len + bmp_sec_len :]

    blocks = []
    off = 0; bmp_off = 0; cdn_off = 0
    while off < len(other_sec):
        blk_start = off
        off = skip_contig_dict(other_sec, off)
        n_recs,   off = read_varint(other_sec, off)
        contig_b, off = read_varint(other_sec, off)
        sstart_b, off = read_varint(other_sec, off)
        span_b,   off = read_varint(other_sec, off)
        bmp_b,    off = read_varint(other_sec, off)
        cdn_b,    off = read_varint(other_sec, off)
        # decode spans
        p = off + contig_b + sstart_b
        spans = []
        for _ in range(n_recs):
            sp, p = read_varint(other_sec, p)
            spans.append(sp)
        fixed_sz = n_recs * 5  # qframe(1)+pident(2)+evalue(2)
        off = p + fixed_sz
        other_bytes = other_sec[blk_start:off]
        bmp_bytes   = bmp_sec[bmp_off : bmp_off + bmp_b]
        cdn_bytes   = cdn_sec[cdn_off : cdn_off + cdn_b]
        bmp_off += bmp_b; cdn_off += cdn_b
        blocks.append((other_bytes, bmp_bytes, cdn_bytes, spans))
    return blocks, other_sec, bmp_sec, cdn_sec

def rle_encode(data: bytes) -> bytes:
    """Simple byte-level RLE: [count(u8), value(u8)] for runs ≥ 2, else [0, byte]."""
    out = bytearray()
    i = 0
    while i < len(data):
        run = 1
        while i + run < len(data) and data[i + run] == data[i] and run < 255:
            run += 1
        if run >= 2:
            out.append(run); out.append(data[i])
        else:
            out.append(0); out.append(data[i])
        i += run
    return bytes(out)

def train_dict(samples: list, dict_size: int = 112640):
    """Train a zstd dictionary; returns ZstdCompressionDict."""
    d = zstd_dict.train_dictionary(dict_size, samples)
    return d

def compress_with_dict(buf: bytes, d, level: int = 3) -> int:
    cctx = zstd_dict.ZstdCompressor(level=level, dict_data=d)
    return len(cctx.compress(buf))

# ── measure ──────────────────────────────────────────────────────────────────

def collect_frames(lhp_dir, max_frames):
    """Yield (raw_bytes, is_v4) for each frame."""
    for fname in sorted(os.listdir(lhp_dir)):
        if not fname.endswith('.lhp'): continue
        with open(os.path.join(lhp_dir, fname), 'rb') as f:
            while True:
                hdr = f.read(4)
                if len(hdr) < 4: break
                csz = struct.unpack_from('<I', hdr)[0]
                cdata = f.read(csz)
                if len(cdata) < csz: break
                raw = bytes(zstd.decompress(cdata))
                yield raw, is_v4_frame(raw)

def measure(lhp_dir, max_frames, dict_train_frames=50):
    totals = collections.defaultdict(int)
    n_frames = 0

    # ── pass 0: collect sample frames to train dicts ─────────────────────────
    print("  Collecting dict training samples...", flush=True)
    bmp_samples, cdn_samples, full_samples = [], [], []
    for raw, v4 in collect_frames(lhp_dir, dict_train_frames):
        if v4:
            blocks, other_sec, bmp_sec, cdn_sec = parse_v4_frame(raw)
            # per-block samples (not whole sections) for better dict coverage
            for ob, bmp, cdn, _ in blocks:
                if len(bmp) > 8:  bmp_samples.append(bmp)
                if len(cdn) > 8:  cdn_samples.append(cdn)
            full_samples.append(raw)
        else:
            off = 0
            while off < len(raw):
                start = off
                bmp, cdn, _, off = parse_block(raw, off)
                if len(bmp) > 8: bmp_samples.append(bmp)
                if len(cdn) > 8: cdn_samples.append(cdn)
            full_samples.append(raw)

    bmp_dict = full_dict = None
    if bmp_samples:
        try:
            bmp_dict  = train_dict(bmp_samples)
            full_dict = train_dict(full_samples)
            print(f"  Dicts trained: bmp={len(bmp_samples)} samples, full={len(full_samples)} samples", flush=True)
        except Exception as e:
            print(f"  Dict training failed: {e}", flush=True)

    # ── pass 1: benchmark strategies ─────────────────────────────────────────
    for raw, v4 in collect_frames(lhp_dir, max_frames):
        if v4:
            blocks, other_sec, bmp_sec, cdn_sec = parse_v4_frame(raw)
            # reassemble sections for RLE-codon variant
            rle_cdn = rle_encode(cdn_sec)
            raw_rle = raw[:8 + len(other_sec) + len(bmp_sec)] + rle_cdn

            # larger-frame: just measure if we re-packed into 4MB frames
            # (can't simulate cross-frame context here, so skip: use frame size flag instead)

            variants: dict[str, tuple[bytes, int]] = {
                'v4_l3':        (raw,     3),
                'v4_l6':        (raw,     6),
                'v4_l9':        (raw,     9),
                'v4_l19':       (raw,    19),
                'v4_rle_cdn_l3':(raw_rle, 3),
            }
            if bmp_dict:
                variants['v4_dict_l3']     = None  # handled separately
                variants['v4_dict_bmp_l3'] = None
        else:
            # v3: replicate old strategies for comparison
            col_maj = col_major(raw)
            variants = {
                'v3_l3':      (raw,     3),
                'v3_colmaj_l3':(col_maj, 3),
                'v3_l6':      (raw,     6),
            }

        for k, v in variants.items():
            if v is None: continue
            buf, lv = v
            totals[k+'_raw'] += len(buf)
            totals[k+'_cmp'] += compress_lv(buf, lv)

        if v4 and bmp_dict:
            totals['v4_dict_l3_raw']     += len(raw)
            totals['v4_dict_l3_cmp']     += compress_with_dict(raw, full_dict, 3)
            totals['v4_dict_bmp_l3_raw'] += len(raw)
            # dict only on bmp+cdn sections, plain on other
            bmp_cmp = compress_with_dict(bmp_sec + cdn_sec, bmp_dict, 3)
            other_cmp = compress_lv(other_sec, 3)
            totals['v4_dict_bmp_l3_cmp'] += bmp_cmp + other_cmp + 8  # 8 = two u32 headers

        n_frames += 1
        if max_frames and n_frames >= max_frames:
            break

    return totals, n_frames

if __name__ == '__main__':
    lhp_dir    = sys.argv[1]
    max_frames = int(sys.argv[2]) if len(sys.argv) > 2 else None

    print(f"Sampling {lhp_dir}  max_frames={max_frames or 'all'}")
    t, nf = measure(lhp_dir, max_frames)
    print(f"Frames measured: {nf}\n")

    # pick reference: v4_l3 if present, else v3_l3
    ref_key = 'v4_l3' if 'v4_l3_raw' in t else 'v3_l3'
    ref_cmp = t[ref_key + '_cmp']

    all_keys = [k[:-4] for k in t if k.endswith('_raw')]
    all_keys.sort()
    print(f"{'Strategy':<22} {'Raw MB':>9} {'Cmp MB':>9} {'vs ' + ref_key:>16}")
    print("-" * 60)
    for k in all_keys:
        raw = t[k+'_raw']
        cmp = t[k+'_cmp']
        marker = " ←" if k == ref_key else ""
        print(f"  {k:<20} {raw/1e6:8.1f}  {cmp/1e6:8.1f}  {ref_cmp/cmp:14.3f}×{marker}")
