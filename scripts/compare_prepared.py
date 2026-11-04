#!/usr/bin/env python3
"""Compare a PREPARE output against the archived reference h5ad.

Usage:
    compare_prepared.py REFERENCE CANDIDATE [--skip-doublets]

REFERENCE is the h5ad that produced the published Cell Explorer store; CANDIDATE
is the PREPARE output. Prints one PASS/FAIL line per check and exits 0 only if
every check passes. X is never loaded whole: its CSR indptr/indices/data are read
from h5py in blocks. obs and var (hundreds of MB) are loaded.

The two files encode their indexes differently (the reference has a named,
nullable-string obs index; the candidate has flat unnamed ones), so obs and var
are read with anndata's element reader and compared by values, not codes.

--skip-doublets leaves out is_doublet, scdbl_score and doublet_status, which the
merge-only PREPARE does not produce.
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

BLOCK = 50_000_000  # values compared at a time
IDENTITY = ["library", "site", "patient", "tcr_status"]
DOUBLETS = {"is_doublet", "scdbl_score", "doublet_status"}


def load_prepare():
    path = Path(__file__).resolve().parent.parent / "bin" / "prepare_h5ad.py"
    spec = importlib.util.spec_from_file_location("prepare_h5ad", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prepare_h5ad"] = mod
    spec.loader.exec_module(mod)
    return mod


def read_frame(f, key):
    return ad.io.read_elem(f[key])


def same_values(a: pd.Series, b: pd.Series) -> bool:
    # Categoricals are compared by value, never by code; missing equals missing.
    a, b = a.astype(object), b.astype(object)
    na_a, na_b = a.isna().to_numpy(), b.isna().to_numpy()
    if not np.array_equal(na_a, na_b):
        return False
    return bool((a[~na_a].to_numpy() == b[~na_b].to_numpy()).all())


def compare_x(fr, fc):
    xr, xc = fr["X"], fc["X"]
    if not isinstance(xr, h5py.Group) or not isinstance(xc, h5py.Group):
        return "X is not stored as CSR in both files"
    for key in ("indptr", "indices", "data"):
        a, b = xr[key], xc[key]
        if a.shape != b.shape:
            return f"{key} length {a.shape[0]:,} vs {b.shape[0]:,}"
        if a.dtype != b.dtype:
            return f"{key} dtype {a.dtype} vs {b.dtype}"
        for start in range(0, a.shape[0], BLOCK):
            x, y = a[start : start + BLOCK], b[start : start + BLOCK]
            bad = np.flatnonzero(x != y)
            if bad.size:
                return f"first mismatch at {key}[{start + int(bad[0])}]"
    if tuple(xr.attrs["shape"]) != tuple(xc.attrs["shape"]):
        return f"shape {tuple(xr.attrs['shape'])} vs {tuple(xc.attrs['shape'])}"
    return None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("reference")
    p.add_argument("candidate")
    p.add_argument("--skip-doublets", action="store_true")
    args = p.parse_args(argv)

    failed = []

    def report(name, problem):
        print(f"{'PASS' if problem is None else 'FAIL'} {name}" + (f": {problem}" if problem else ""))
        if problem is not None:
            failed.append(name)

    with h5py.File(args.reference, "r") as fr, h5py.File(args.candidate, "r") as fc:
        ro, co = read_frame(fr, "obs"), read_frame(fc, "obs")
        rv, cv = read_frame(fr, "var"), read_frame(fc, "var")

        n = min(len(ro), len(co))
        bad = np.flatnonzero(ro.index.to_numpy(dtype=object)[:n] != co.index.to_numpy(dtype=object)[:n])
        report(
            "obs_names",
            None if len(ro) == len(co) and not bad.size
            else f"{len(ro):,} vs {len(co):,} cells, first difference at {int(bad[0]) if bad.size else n}",
        )
        report("var_names", None if list(rv.index) == list(cv.index) else "differ")
        report("X", compare_x(fr, fc))

        skip = DOUBLETS if args.skip_doublets else set()
        for col in IDENTITY:
            if col not in ro or col not in co:
                report(col, "missing in " + ("reference" if col not in ro else "candidate"))
            elif len(ro) != len(co):
                report(col, "cell counts differ")
            else:
                report(col, None if same_values(ro[col], co[col]) else "values differ")

        # Candidate modality = reference modality filled from Modality, by the
        # same rule PREPARE used to merge the case-duplicate columns.
        if "modality" in co and "modality" in ro and len(ro) == len(co):
            prep = load_prepare()
            pair = ro[["modality"] + (["Modality"] if "Modality" in ro else [])].copy()
            try:
                prep.merge_case_duplicate_columns(pair)
                problem = None if same_values(pair["modality"], co["modality"]) else "values differ"
            except SystemExit as e:
                problem = f"reference modality columns conflict: {e}"
            report("modality", problem)
        else:
            report("modality", "column missing or cell counts differ")

        expected = set(ro.columns) - {"Modality"} - skip
        actual = set(co.columns) - skip
        report(
            "obs columns",
            None if expected == actual
            else f"missing {sorted(expected - actual)}, extra {sorted(actual - expected)}",
        )

        for col in sorted(skip & set(ro.columns) & set(co.columns)):
            print(f"SKIP {col}")
        if "leiden" in co:
            print("INFO leiden clusters: %d" % co["leiden"].nunique())
            print("INFO leiden counts:", co["leiden"].value_counts().sort_index().to_dict())

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
