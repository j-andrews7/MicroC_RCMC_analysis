#!/usr/bin/env python3
"""Compare outputs of the legacy Snakemake pipeline with outputs of this Nextflow pipeline.

Usage:
    compare_outputs.py LEGACY_OUTDIR NEXTFLOW_OUTDIR [--weight-rtol 1e-6] [--weight-atol 1e-12]

Compares every file the legacy pipeline produces (beds/, dedup/, repmerge/, repcools/, conditioncools/):

* ``*.pairs.gz``  - records (header lines starting with '#' are ignored because they embed command lines,
                    file paths and thread counts). Records are compared as multisets within each block of identical
                    sort keys (chrom1, chrom2, pos1, pos2, pair_type), reporting records found in only one file and,
                    separately (as a warning), records that only differ in order between tied sort keys.
* ``*.stats``     - every key/value pair.
* ``*.bed``       - exact text.
* ``*.px2``       - presence only (binary index, depends on the bgzip block layout).
* ``*.cool``/``*.mcool`` - for every resolution: chromosomes and bins exactly; pixels as a sparse join (pixels found
                    in only one file, shared pixels with different counts, total contacts moved); balancing weights
                    within tolerance; and key attributes (creation date, generator, etc. are ignored as they always
                    differ between runs).

Exit status is 0 if everything matches and 1 otherwise. Requires numpy and h5py (both are available in the
pipeline's tools container, e.g. ``singularity exec <microc_tools.sif> python compare_outputs.py ...``).
"""

import argparse
import gzip
import itertools
import shutil
import subprocess
import sys
from array import array
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

LEGACY_DIRS = ["beds", "dedup", "repmerge", "repcools", "conditioncools"]
COOL_ATTRS = ["bin-size", "bin-type", "nbins", "nchroms", "nnz", "storage-mode", "sum", "genome-assembly", "format"]
CHUNK = 10_000_000


class Report:
    def __init__(self):
        self.failures = 0
        self.warnings = 0

    def ok(self, msg):
        print(f"  [PASS] {msg}")

    def warn(self, msg):
        self.warnings += 1
        print(f"  [WARN] {msg}")

    def fail(self, msg):
        self.failures += 1
        print(f"  [FAIL] {msg}")

    def info(self, msg):
        print(f"  [INFO] {msg}")


def open_text(path):
    """Return a binary line iterator for a (b)gzipped or plain file, using a fast external decompressor if possible."""
    path = str(path)
    if path.endswith(".gz"):
        for cmd in (["bgzip", "-dc", "-@", "4", path], ["gzip", "-dc", path]):
            if shutil.which(cmd[0]):
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 20)
                return proc.stdout
        return gzip.open(path, "rb")
    return open(path, "rb")


def split_header(lines):
    header = []
    first_body = None
    for line in lines:
        if line.startswith(b"#"):
            header.append(line)
        else:
            first_body = line
            break
    body = itertools.chain([first_body], lines) if first_body is not None else iter(())
    return header, body


PAIRS_SORT_COLS = ["chrom1", "chrom2", "pos1", "pos2", "pair_type"]
DEFAULT_PAIRS_COLS = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "pair_type"]
MAX_EXAMPLE_BLOCK = 1000


def pairs_key_function(header):
    """Key matching the `pairtools sort` order: chrom1, chrom2, pos1 (numeric), pos2 (numeric), pair_type."""
    cols_line = next((h for h in header if h.startswith(b"#columns")), None)
    names = cols_line.split(b":", 1)[1].split() if cols_line else [c.encode() for c in DEFAULT_PAIRS_COLS]
    names = [n.decode() for n in names]
    fields = [(names.index(c), c.startswith("pos")) for c in PAIRS_SORT_COLS if c in names]
    max_idx = max(i for i, _ in fields)

    def key(line):
        f = line.split(b"\t", max_idx + 1)
        return tuple(int(f[i]) if numeric else f[i].rstrip(b"\n") for i, numeric in fields)

    return key


def iter_pair_blocks(body, key_fn):
    """Yield (key, line hashes, lines-or-None) for runs of records sharing the same sort key."""
    key, hashes, lines = None, array("q"), []
    for line in body:
        k = key_fn(line)
        if k != key:
            if key is not None:
                yield key, hashes, lines
            key, hashes, lines = k, array("q"), []
        hashes.append(hash(line))
        if lines is not None:
            lines.append(line)
            if len(lines) > MAX_EXAMPLE_BLOCK:
                lines = None
    if key is not None:
        yield key, hashes, lines


def multiset_difference(hashes_a, hashes_b):
    xa = np.frombuffer(hashes_a, dtype=np.int64)
    xb = np.frombuffer(hashes_b, dtype=np.int64)
    ua, ca = np.unique(xa, return_counts=True)
    ub, cb = np.unique(xb, return_counts=True)
    _, ia, ib = np.intersect1d(ua, ub, assume_unique=True, return_indices=True)
    matched = int(np.minimum(ca[ia], cb[ib]).sum())
    return len(xa) - matched, len(xb) - matched


def compare_pairs(a, b, rep):
    """Compare pairs bodies block by block (records sharing the same sort key), as multisets and as sequences."""
    ha, body_a = split_header(open_text(a))
    hb, body_b = split_header(open_text(b))
    cols_a = [h for h in ha if h.startswith(b"#columns")]
    cols_b = [h for h in hb if h.startswith(b"#columns")]
    if cols_a != cols_b:
        rep.fail(f"{a.name}: different #columns header: {cols_a} vs {cols_b}")
        return
    key_fn = pairs_key_function(ha)
    blocks_a, blocks_b = iter_pair_blocks(body_a, key_fn), iter_pair_blocks(body_b, key_fn)
    ba, bb = next(blocks_a, None), next(blocks_b, None)
    n_a = n_b = only_a = only_b = reordered = 0
    examples_a, examples_b = [], []
    while ba is not None or bb is not None:
        if bb is None or (ba is not None and ba[0] < bb[0]):
            n_a += len(ba[1])
            only_a += len(ba[1])
            if ba[2] and len(examples_a) < 3:
                examples_a.extend(ba[2][: 3 - len(examples_a)])
            ba = next(blocks_a, None)
        elif ba is None or bb[0] < ba[0]:
            n_b += len(bb[1])
            only_b += len(bb[1])
            if bb[2] and len(examples_b) < 3:
                examples_b.extend(bb[2][: 3 - len(examples_b)])
            bb = next(blocks_b, None)
        else:
            n_a += len(ba[1])
            n_b += len(bb[1])
            if ba[1] != bb[1]:
                da, db = multiset_difference(ba[1], bb[1])
                if da == 0 and db == 0:
                    reordered += len(ba[1])
                else:
                    only_a += da
                    only_b += db
                    if ba[2] is not None and bb[2] is not None:
                        ca, cb = Counter(ba[2]), Counter(bb[2])
                        examples_a.extend(list((ca - cb).elements())[: max(0, 3 - len(examples_a))])
                        examples_b.extend(list((cb - ca).elements())[: max(0, 3 - len(examples_b))])
            ba, bb = next(blocks_a, None), next(blocks_b, None)

    if only_a == 0 and only_b == 0:
        if reordered == 0:
            rep.ok(f"{a.name}: {n_a:,} records identical")
        else:
            rep.warn(
                f"{a.name}: same {n_a:,} records, but {reordered:,} records with tied sort keys are in a different order"
            )
        return
    rep.fail(
        f"{a.name}: {n_a:,} (legacy) vs {n_b:,} (nextflow) records; "
        f"{only_a:,} only in legacy ({100.0 * only_a / max(n_a, 1):.4f}%), "
        f"{only_b:,} only in nextflow ({100.0 * only_b / max(n_b, 1):.4f}%)"
        + (f"; {reordered:,} further records only reordered" if reordered else "")
    )
    for label, examples in (("legacy only", examples_a), ("nextflow only", examples_b)):
        for line in examples:
            rep.info(f"{label}: {line!r}")


def parse_stats(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            key, _, value = line.partition("\t")
            out[key] = value
    return out


def compare_stats(a, b, rep):
    sa, sb = parse_stats(a), parse_stats(b)
    keys = sorted(set(sa) | set(sb))
    diff = [k for k in keys if sa.get(k) != sb.get(k)]
    if not diff:
        rep.ok(f"{a.name}: {len(keys)} statistics identical")
    else:
        rep.fail(f"{a.name}: {len(diff)} of {len(keys)} statistics differ")
        for k in diff[:15]:
            rep.info(f"{k}: legacy={sa.get(k)} nextflow={sb.get(k)}")


def compare_text(a, b, rep):
    if Path(a).read_bytes().replace(b"\r\n", b"\n") == Path(b).read_bytes().replace(b"\r\n", b"\n"):
        rep.ok(f"{a.name}: identical")
    else:
        rep.fail(f"{a.name}: content differs")


def cool_groups(h5):
    if "resolutions" in h5:
        return {f"resolutions/{r}": h5["resolutions"][r] for r in sorted(h5["resolutions"], key=int)}
    return {"/": h5}


def compare_arrays_chunked(da, db, name, rep, label):
    if da.shape != db.shape:
        rep.fail(f"{label}: {name} length differs ({da.shape[0]:,} vs {db.shape[0]:,})")
        return False
    n_diff = 0
    for start in range(0, da.shape[0], CHUNK):
        n_diff += int(np.count_nonzero(da[start : start + CHUNK] != db[start : start + CHUNK]))
    if n_diff:
        rep.fail(f"{label}: {name} has {n_diff:,} differing values out of {da.shape[0]:,}")
        return False
    return True


def compare_weights(xa, xb, rep, label, rtol, atol):
    has_wa, has_wb = "weight" in xa["bins"], "weight" in xb["bins"]
    if has_wa != has_wb:
        rep.fail(f"{label}: weight column present in only one file")
        return None
    if not has_wa:
        return "no weights"
    wa, wb = xa["bins"]["weight"][:], xb["bins"]["weight"][:]
    nan_a, nan_b = np.isnan(wa), np.isnan(wb)
    if not np.array_equal(nan_a, nan_b):
        rep.fail(
            f"{label}: weight NaN masks differ ({int(nan_a.sum()):,} vs {int(nan_b.sum()):,} NaN bins, "
            f"{int(np.count_nonzero(nan_a != nan_b)):,} bins differ)"
        )
        nan_a = nan_a | nan_b
    fin = ~nan_a
    if np.array_equal(wa[fin], wb[fin]):
        return "weights identical"
    rel = np.abs(wa[fin] - wb[fin]) / np.maximum(np.abs(wa[fin]), 1e-300)
    close = np.isclose(wa[fin], wb[fin], rtol=rtol, atol=atol)
    summary = f"max rel diff {float(rel.max()):.2e}, median {float(np.median(rel)):.2e}"
    if close.all():
        return f"weights equal within tolerance ({summary})"
    rep.fail(
        f"{label}: {int((~close).sum()):,} of {int(fin.sum()):,} weights differ beyond tolerance ({summary}, "
        f"Pearson r = {float(np.corrcoef(wa[fin], wb[fin])[0, 1]):.6f})"
    )
    return None


def compare_pixels(xa, xb, rep, label):
    """Sparse comparison of pixel tables, processed in ranges of bin1 using the bin1_offset index."""
    off_a = xa["indexes"]["bin1_offset"][:]
    off_b = xb["indexes"]["bin1_offset"][:]
    nbins = len(off_a) - 1
    if len(off_b) - 1 != nbins:
        rep.fail(f"{label}: different number of bins")
        return False
    pa, pb = xa["pixels"], xb["pixels"]
    only_a = only_b = changed = 0
    moved = 0
    total = 0
    start = 0
    while start < nbins:
        end = int(np.searchsorted(off_a, off_a[start] + CHUNK, side="right")) - 1
        end = min(max(end, start + 1), nbins)
        sa, sb = slice(int(off_a[start]), int(off_a[end])), slice(int(off_b[start]), int(off_b[end]))
        ka = pa["bin1_id"][sa].astype(np.int64) * nbins + pa["bin2_id"][sa]
        kb = pb["bin1_id"][sb].astype(np.int64) * nbins + pb["bin2_id"][sb]
        ca, cb = pa["count"][sa].astype(np.int64), pb["count"][sb].astype(np.int64)
        _, ia, ib = np.intersect1d(ka, kb, assume_unique=True, return_indices=True)
        mask_a = np.ones(len(ka), bool)
        mask_a[ia] = False
        mask_b = np.ones(len(kb), bool)
        mask_b[ib] = False
        only_a += int(mask_a.sum())
        only_b += int(mask_b.sum())
        changed += int(np.count_nonzero(ca[ia] != cb[ib]))
        moved += int(np.abs(ca[ia] - cb[ib]).sum() + ca[mask_a].sum() + cb[mask_b].sum())
        total += int(ca.sum())
        start = end
    if only_a == only_b == changed == 0:
        return True
    rep.fail(
        f"{label}: pixels only in legacy: {only_a:,}, only in nextflow: {only_b:,}, shared with different counts: "
        f"{changed:,}; summed |count difference| = {moved:,} ({100.0 * moved / max(total, 1):.4f}% of {total:,} contacts)"
    )
    return False


def compare_cool(a, b, rep, rtol, atol):
    with h5py.File(a, "r") as ha, h5py.File(b, "r") as hb:
        ga, gb = cool_groups(ha), cool_groups(hb)
        if set(ga) != set(gb):
            rep.fail(f"{a.name}: different resolutions: {sorted(ga)} vs {sorted(gb)}")
        for key in [k for k in ga if k in gb]:
            label = f"{a.name}::{key}"
            xa, xb = ga[key], gb[key]
            problems = 0
            for attr in COOL_ATTRS:
                va, vb = xa.attrs.get(attr), xb.attrs.get(attr)
                if (va is not None or vb is not None) and not np.array_equal(np.asarray(va), np.asarray(vb)):
                    rep.fail(f"{label}: attribute {attr} differs ({va} vs {vb})")
                    problems += 1
            for col in ["name", "length"]:
                if not np.array_equal(xa["chroms"][col][:], xb["chroms"][col][:]):
                    rep.fail(f"{label}: chroms/{col} differs")
                    problems += 1
            bins_ok = True
            for col in ["chrom", "start", "end"]:
                if not compare_arrays_chunked(xa["bins"][col], xb["bins"][col], f"bins/{col}", rep, label):
                    problems += 1
                    bins_ok = False
            if bins_ok and not compare_pixels(xa, xb, rep, label):
                problems += 1
            weight_msg = compare_weights(xa, xb, rep, label, rtol, atol)
            if weight_msg is None:
                problems += 1
            if problems == 0:
                rep.ok(f"{label}: {xa['pixels']['count'].shape[0]:,} pixels identical, {weight_msg}")


def list_outputs(root):
    return sorted(
        p.relative_to(root)
        for d in LEGACY_DIRS
        if (root / d).is_dir()
        for p in (root / d).iterdir()
        if p.is_file()
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("legacy", type=Path, help="Output directory of the legacy Snakemake pipeline")
    parser.add_argument("nextflow", type=Path, help="Output directory (--outdir) of the Nextflow pipeline")
    parser.add_argument("--weight-rtol", type=float, default=1e-6, help="Relative tolerance for balancing weights")
    parser.add_argument("--weight-atol", type=float, default=1e-12, help="Absolute tolerance for balancing weights")
    parser.add_argument("--skip-pairs", action="store_true", help="Do not compare *.pairs.gz files (slow for big data)")
    args = parser.parse_args()

    rep = Report()
    legacy_files = list_outputs(args.legacy)
    nf_files = list_outputs(args.nextflow)
    if not legacy_files:
        sys.exit(f"No legacy outputs found under {args.legacy} ({', '.join(LEGACY_DIRS)})")

    print(f"Legacy:   {args.legacy}\nNextflow: {args.nextflow}\n")
    missing = [f for f in legacy_files if f not in nf_files]
    print("== File inventory")
    if missing:
        for f in missing:
            rep.fail(f"missing from Nextflow outputs: {f}")
    else:
        rep.ok(f"all {len(legacy_files)} legacy output files are present")
    for f in nf_files:
        if f not in legacy_files:
            rep.info(f"only in Nextflow outputs: {f}")

    for rel in legacy_files:
        if rel in missing:
            continue
        a, b = args.legacy / rel, args.nextflow / rel
        name = rel.name
        print(f"== {rel}")
        if name.endswith(".pairs.gz"):
            if args.skip_pairs:
                rep.info("skipped (--skip-pairs)")
            else:
                compare_pairs(a, b, rep)
        elif name.endswith(".stats"):
            compare_stats(a, b, rep)
        elif name.endswith(".bed"):
            compare_text(a, b, rep)
        elif name.endswith(".px2"):
            rep.ok("index present (binary content not compared)")
        elif name.endswith(".cool") or name.endswith(".mcool"):
            compare_cool(a, b, rep, args.weight_rtol, args.weight_atol)
        else:
            rep.info("not compared")

    print()
    if rep.failures:
        print(f"RESULT: {rep.failures} difference(s) found" + (f", {rep.warnings} warning(s)" if rep.warnings else ""))
        sys.exit(1)
    if rep.warnings:
        print(f"RESULT: all compared outputs match ({rep.warnings} warning(s): record order differs within tied sort keys)")
    else:
        print("RESULT: all compared outputs match")


if __name__ == "__main__":
    main()
