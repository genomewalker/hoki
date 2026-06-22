#!/usr/bin/env python3
"""Integration tests for hoki ingest correctness and error handling."""
import subprocess, sys, os, tempfile, shutil

HOKI = os.path.join(os.path.dirname(__file__), "../build/hoki")

def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)

def hoki(*args):
    r = run([HOKI, *args])
    return r.stdout + r.stderr

# ── helpers ──────────────────────────────────────────────────────────────────

VALID_HDR = "qseqid\tqstart\tqend\tqlen\tqstrand\tsseqid\tsstart\tsend\tslen\tpident\tevalue\tcigar\tqseq_aa\tfull_qseq\n"

def make_tsv(rows, tmp):
    p = os.path.join(tmp, "test.tsv")
    with open(p, "w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    return p

# One valid row: pident=85, 3M cigar, AA=MKG, full_nt encodes those codons.
# Codons: ATG(M)=0x22, AAA(K)=0x00, GGT(G)=0x2A → packed per format.hpp
VALID_ROW = [
    "SRR0000001_1_ka:f:17_L:+:550:+",
    4, 12, 15, "+",
    "bpgv2|N0.HOG0000001",
    0, 2, 10,
    85.0, 1e-20, "3M", "MKG",
    "AAGATGAAAGGT"
]

def ingest_and_merge(rows, tmp, flush="2G"):
    tsv = make_tsv(rows, tmp)
    part = os.path.join(tmp, "part")
    lhg  = os.path.join(tmp, "out.lhg")
    lhgi = os.path.join(tmp, "out.lhgi")
    r = run([HOKI, "ingest", "-a", "auto", "--flush", flush, tsv, part])
    if r.returncode != 0:
        return None, r.stderr
    r = run([HOKI, "merge-shard", lhg, lhgi, part])
    if r.returncode != 0:
        return None, r.stderr
    return (lhg, lhgi), ""

# ── test 1: bad pident is skipped, not silently processed ────────────────────

def test_bad_float():
    tmp = tempfile.mkdtemp()
    try:
        bad = list(VALID_ROW); bad[9] = "not_a_float"
        good = list(VALID_ROW)
        tsv = make_tsv([bad, good], tmp)
        part = os.path.join(tmp, "part")
        out = hoki("ingest", "-a", "auto", tsv, part)
        assert "1 skipped" in out, f"expected 1 skipped, got: {out!r}"
        assert "1 records" in out, f"expected 1 records, got: {out!r}"
        print("test_bad_float: PASS")
    finally:
        shutil.rmtree(tmp)

# ── test 2: NT-absent rows don't inflate n_written ──────────────────────────

def test_nt_absent_row():
    tmp = tempfile.mkdtemp()
    try:
        # 13-column row (no full_qseq) — passes filter but produces no codon obs
        no_nt = list(VALID_ROW[:13])  # drop full_qseq
        good  = list(VALID_ROW)
        rows = [no_nt, good]
        tsv = make_tsv(rows, tmp)
        part = os.path.join(tmp, "part")
        out = hoki("ingest", "-a", "auto", tsv, part)
        # Only the row with NT data should count as a written record
        assert "1 records" in out, f"expected 1 records (with obs), got: {out!r}"
        print("test_nt_absent_row: PASS")
    finally:
        shutil.rmtree(tmp)

# ── test 3: flush + no-flush produce identical freq output ───────────────────

def test_flush_vs_nofflush():
    tmp = tempfile.mkdtemp()
    try:
        rows = [list(VALID_ROW) for _ in range(20)]
        for i, r in enumerate(rows):
            r[0] = f"SRR000000{(i%5)+1}_c{i}_ka:f:1_L:+:10:+"
            r[9] = 85.0 + i * 0.1
        d1 = os.path.join(tmp, "a"); os.makedirs(d1)
        d2 = os.path.join(tmp, "b"); os.makedirs(d2)
        (lhg1, lhgi1), err1 = ingest_and_merge(rows, d1, flush="1K")
        (lhg2, lhgi2), err2 = ingest_and_merge(rows, d2, flush="2G")
        assert lhg1 and lhg2, f"ingest failed: {err1} {err2}"
        f1 = hoki("freq", lhg1, lhgi1, "N0.HOG0000001")
        f2 = hoki("freq", lhg2, lhgi2, "N0.HOG0000001")
        assert f1 == f2, f"freq mismatch:\n1K: {f1}\n2G: {f2}"
        print("test_flush_vs_nofflush: PASS")
    finally:
        shutil.rmtree(tmp)

# ── test 4: pident==100 is skipped ──────────────────────────────────────────

def test_perfect_match_skipped():
    tmp = tempfile.mkdtemp()
    try:
        perfect = list(VALID_ROW); perfect[9] = 100.0
        good    = list(VALID_ROW)
        tsv = make_tsv([perfect, good], tmp)
        part = os.path.join(tmp, "part")
        out = hoki("ingest", "-a", "auto", tsv, part)
        assert "1 skipped" in out, f"expected 1 skipped, got: {out!r}"
        print("test_perfect_match_skipped: PASS")
    finally:
        shutil.rmtree(tmp)

# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    failed = 0
    for fn in [test_bad_float, test_nt_absent_row, test_flush_vs_nofflush, test_perfect_match_skipped]:
        try:
            fn()
        except Exception as e:
            print(f"{fn.__name__}: FAIL — {e}")
            failed += 1
    sys.exit(failed)
