#!/usr/bin/env python3
"""
Measure per-column byte shares in a grouped-frames .lhp directory (v3 format).

Usage:
    python3 measure_lhp.py <out_dir> [max_frames_per_thread]

Scans each tN.lhp frame-by-frame. Each frame is [frame_csz(u32) | zstd_frame].
The decompressed frame contains concatenated VarNT blocks; blocks are self-delimiting
(header encodes column byte counts) so the index is not needed.
"""
import sys, os, struct, collections
import zstd

def read_u32le(data, off): return struct.unpack_from('<I', data, off)[0]

def read_varint(buf, off):
    v, shift = 0, 0
    while off < len(buf):
        b = buf[off]; off += 1
        v |= (b & 0x7F) << shift
        if not (b & 0x80): return v, off
        shift += 7
    raise ValueError("truncated varint")

def skip_contig_dict(data, off):
    """Skip local contig dict: n_contigs (varint) + each (len varint + name bytes)."""
    n_contigs, off = read_varint(data, off)
    for _ in range(n_contigs):
        nlen, off = read_varint(data, off)
        off += nlen
    return off

def block_size(data, off):
    """Parse one raw HOG payload (contig_dict + VarNT block). Returns (cols, new_off).
    VarNT header: n_recs, contig_b, sstart_b, span_b, bmp_b, cdn_b
    Columns on disk: contig | sstart | span | qframe | pident | evalue | bitmap | codons
    qframe/pident/evalue sizes not in header — derived from n_recs.
    """
    off = skip_contig_dict(data, off)
    n_recs,   off = read_varint(data, off)
    contig_b, off = read_varint(data, off)
    sstart_b, off = read_varint(data, off)
    span_b,   off = read_varint(data, off)
    bmp_b,    off = read_varint(data, off)
    cdn_b,    off = read_varint(data, off)
    payload_b = contig_b + sstart_b + span_b + n_recs + n_recs*2 + n_recs*2 + bmp_b + cdn_b
    cols = {
        'n_recs':  n_recs,
        'contig':  contig_b,
        'sstart':  sstart_b,
        'span':    span_b,
        'qframe':  n_recs,
        'pident':  n_recs * 2,
        'evalue':  n_recs * 2,
        'bitmap':  bmp_b,
        'codons':  cdn_b,
    }
    return cols, off + payload_b

def measure(out_dir, max_frames=None):
    totals = collections.defaultdict(int)
    total_compressed = 0
    total_frames = 0

    for fname in sorted(os.listdir(out_dir)):
        if not fname.endswith('.lhp'):
            continue
        path = os.path.join(out_dir, fname)
        nframes = 0
        with open(path, 'rb') as f:
            while True:
                hdr = f.read(4)
                if len(hdr) < 4:
                    break
                csz = struct.unpack_from('<I', hdr)[0]
                cdata = f.read(csz)
                if len(cdata) < csz:
                    break
                total_compressed += csz
                total_frames += 1
                nframes += 1
                raw = zstd.decompress(cdata)
                off = 0
                while off < len(raw):
                    cols, off = block_size(raw, off)
                    for k in ('contig','sstart','span','qframe','pident','evalue','bitmap','codons'):
                        totals[k] += cols[k]
                    totals['n_recs'] += cols['n_recs']
                if max_frames and nframes >= max_frames:
                    break

    return totals, total_compressed, total_frames

def fmt_mb(b): return f"{b/1e6:9.1f} MB"
def pct(b, t): return f"{100*b/t:.1f}%" if t else "—"

if __name__ == '__main__':
    out_dir    = sys.argv[1]
    max_frames = int(sys.argv[2]) if len(sys.argv) > 2 else None

    print(f"Reading {out_dir}  (max_frames_per_thread={max_frames or 'all'})")
    totals, total_compressed, total_frames = measure(out_dir, max_frames)

    cols = ['contig','sstart','span','qframe','pident','evalue','bitmap','codons']
    col_total = sum(totals[c] for c in cols)

    print(f"\n{'Column':<10} {'Raw bytes':>12} {'% of raw':>9}")
    print("-" * 34)
    for c in cols:
        print(f"  {c:<8} {fmt_mb(totals[c])}  {pct(totals[c], col_total)}")
    print("-" * 34)
    print(f"  {'total':<8} {fmt_mb(col_total)}")

    print(f"\nFrames:           {total_frames:,}")
    print(f"Records:          {totals['n_recs']:,}")
    print(f"Raw (col total):  {fmt_mb(col_total)}")
    print(f"Compressed:       {fmt_mb(total_compressed)}")
    if total_compressed:
        print(f"Ratio (raw→zstd): {col_total/total_compressed:.2f}×")
