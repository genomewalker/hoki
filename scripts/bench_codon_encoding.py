#!/usr/bin/env python3
"""
Benchmark codon + bitmap encoding strategies on real .lhp data.

Strategies measured per VarNT block's bitmap+codon columns:
  current   dense bitmap (1bit/pos) + 1byte/codon, flat
  bitpack6  dense bitmap + 6-bit packed codons (no padding bits)
  gapcoded  gap-coded set-position offsets (varint) + 1byte/codon
  gapcoded6 gap-coded set-position offsets (varint) + 6-bit packed codons

All strategies then zstd-compressed (level 3) at the frame level (2MB),
matching what hoki does in production.

Usage:
    python3 bench_codon_encoding.py <lhp_dir> [max_frames]
"""
import sys, os, struct, collections, math, io
import zstd

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
        # "other" = everything in the raw block before the bitmap column
        other_size = (off - len(cdn_bytes) - len(bmp_bytes)) - start
        other_parts.append(raw[start : start + other_size])
        bmp_parts.append(bmp_bytes)
        cdn_parts.append(cdn_bytes)
    return b''.join(other_parts) + b''.join(bmp_parts) + b''.join(cdn_parts)

def measure(lhp_dir, max_frames):
    totals = collections.defaultdict(int)
    n_frames = 0

    for fname in sorted(os.listdir(lhp_dir)):
        if not fname.endswith('.lhp'): continue
        path = os.path.join(lhp_dir, fname)
        nf = 0
        with open(path, 'rb') as f:
            while True:
                hdr = f.read(4)
                if len(hdr) < 4: break
                csz = struct.unpack_from('<I', hdr)[0]
                cdata = f.read(csz)
                if len(cdata) < csz: break

                raw = zstd.decompress(cdata)
                raw_bytes = bytes(raw)

                # ── sort variants ──────────────────────────────────────────
                # bmp_prefix: sort by first 8 bytes of bitmap (same-HOG proxy)
                bmp_sorted = reorder_blocks(raw_bytes,
                    lambda bmp, cdn, sp: bmp[:8].ljust(8, b'\x00'))
                # random shuffle (ablation control)
                import random
                blocks_list = []
                off2 = 0
                while off2 < len(raw_bytes):
                    s = off2
                    _, _, _, off2 = parse_block(raw_bytes, off2)
                    blocks_list.append(raw_bytes[s:off2])
                random.shuffle(blocks_list)
                shuffled = b''.join(blocks_list)

                col_maj = col_major(raw_bytes)

                variants = {
                    'curr_l3':  (raw_bytes,  3),
                    'curr_l19': (raw_bytes, 19),
                    'bmp_sort_l3':  (bmp_sorted,  3),
                    'bmp_sort_l19': (bmp_sorted, 19),
                    'shuffle_l3':   (shuffled,    3),
                    'colmaj_l3':    (col_maj,     3),
                    'colmaj_l19':   (col_maj,    19),
                }
                for k, (buf, lv) in variants.items():
                    totals[k+'_raw'] += len(buf)
                    totals[k+'_cmp'] += compress_lv(buf, lv)

                n_frames += 1
                if max_frames and n_frames >= max_frames:
                    break
        if max_frames and n_frames >= max_frames:
            break

    return totals, n_frames

if __name__ == '__main__':
    lhp_dir    = sys.argv[1]
    max_frames = int(sys.argv[2]) if len(sys.argv) > 2 else None

    print(f"Sampling {lhp_dir}  max_frames={max_frames or 'all'}")
    t, nf = measure(lhp_dir, max_frames)
    print(f"Frames: {nf}")

    ref_raw = t['curr_l3_raw']
    ref_cmp = t['curr_l3_cmp']

    keys = ['curr_l3','curr_l19','bmp_sort_l3','bmp_sort_l19','shuffle_l3','colmaj_l3','colmaj_l19']
    print(f"\n{'Strategy':<16} {'Raw MB':>9} {'Cmp MB':>9} {'vs curr_l3_cmp':>16}")
    print("-" * 54)
    for k in keys:
        raw = t[k+'_raw']
        cmp = t[k+'_cmp']
        print(f"  {k:<14} {raw/1e6:8.1f}  {cmp/1e6:8.1f}  {ref_cmp/cmp:14.3f}×")
