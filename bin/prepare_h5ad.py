#!/usr/bin/env python3
"""Prepare a merged h5ad for cell2zarr / Cell Explorer.

Ported from the on-cluster prepare_gbm.py. Adds the identity columns (library,
site, patient), the TCR status, merges obs columns whose names differ only by
case, computes the PCA / UMAP / Leiden embedding on a normalized copy (X keeps
raw counts), checks the result and writes the h5ad plus provenance.json.
"""

import argparse
import datetime
import hashlib
import importlib.metadata
import json
import logging
import re
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

log = logging.getLogger("prepare_h5ad")

MISSING = ["", "NA", "nan", "None"]


def site_of(library: str, site_regex: re.Pattern) -> str:
    m = site_regex.search(library)
    return m.group(1).replace("-", "") if m else "unknown"


def patient_of(library: str, site: str, patient_regex: re.Pattern, patient_prefix: str) -> str:
    # Libraries that encode the patient get prefix + group 1; the rest come from
    # a single patient per site, so the site stands in for the patient.
    m = patient_regex.search(library)
    return f"{patient_prefix}{m.group(1)}" if m else site


def is_missing(s: pd.Series) -> pd.Series:
    return s.isna() | s.astype(str).isin(MISSING)


def add_identity_columns(
    adata,
    library_column: str,
    site_regex: re.Pattern,
    patient_regex: re.Pattern,
    patient_prefix: str,
) -> int:
    lib = adata.obs[library_column].astype(str)
    site = lib.map(lambda l: site_of(l, site_regex))
    patient = [patient_of(l, s, patient_regex, patient_prefix) for l, s in zip(lib, site)]
    adata.obs["library"] = pd.Categorical(lib)
    adata.obs["site"] = pd.Categorical(site)
    adata.obs["patient"] = pd.Categorical(patient)
    n_unknown = int((site == "unknown").sum())
    if n_unknown:
        log.warning("%s cells have a library that matches no site pattern", f"{n_unknown:,}")
    else:
        log.info("all libraries matched a site pattern")
    return n_unknown


def add_tcr_status(adata, tcr_column: str) -> None:
    if not tcr_column:
        return
    # Compare as text: astype(bool) would turn the string "False" into True.
    positive = adata.obs[tcr_column].astype(str) == "True"
    n_pos = int(positive.sum())
    log.info("tcr_status: %s TCR+ of %s cells", f"{n_pos:,}", f"{adata.n_obs:,}")
    if n_pos == 0:
        log.warning("0 TCR+ cells: column %r holds no value equal to the text 'True' "
                    "(numeric or lower-case encodings are not recognised)", tcr_column)
    adata.obs["tcr_status"] = pd.Categorical(
        ["TCR+" if p else "TCR-" for p in positive], categories=["TCR-", "TCR+"]
    )


def merge_case_duplicate_columns(obs: pd.DataFrame) -> list[str]:
    # Zarr stores each column as a directory, so names that differ only by case
    # collide on case-insensitive filesystems (macOS, Windows) and a copy there
    # silently loses files. Keep the first column of each group, fill its missing
    # values from the others, and stop if any cell has two different values.
    groups = {}
    for col in obs.columns:
        groups.setdefault(col.lower(), []).append(col)
    dropped = []
    for keep, *drop in (g for g in groups.values() if len(g) > 1):
        merged = obs[keep].astype(object)
        for col in drop:
            other = obs[col].astype(object)
            both = ~is_missing(merged) & ~is_missing(other)
            conflicts = both & (merged.astype(str) != other.astype(str))
            if conflicts.any():
                raise SystemExit(
                    f"obs columns {keep!r} and {col!r} disagree on {int(conflicts.sum()):,} cells; "
                    "merge them by hand."
                )
            merged = merged.where(~is_missing(merged) | is_missing(other), other)
        # Keep numeric columns numeric; anndata turns string columns back into
        # categoricals on write.
        obs[keep] = pd.to_numeric(merged) if pd.api.types.is_numeric_dtype(obs[keep]) else merged
        for col in drop:
            del obs[col]
        dropped.extend(drop)
        log.info("merged obs columns %s into %r", drop, keep)
    return dropped


BLOCK = 10_000_000  # values checked at a time, to bound temporary memory


def check_input(adata, normalize: bool) -> list[str]:
    # Normalizing an object that is already normalized silently gives a wrong
    # embedding, so require raw counts. A merge can hold per-library pieces, so
    # every stored value is checked, in blocks.
    if not normalize:
        return []
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).reshape(-1)
    for start in range(0, values.shape[0], BLOCK):
        block = np.asarray(values[start:start + BLOCK])
        if not np.all(block == np.round(block)) or (block < 0).any():
            return ["X is not non-negative integer counts; --normalize needs raw counts "
                    "(this object looks already normalized)."]
    return []


@contextmanager
def step(name: str):
    t0 = time.monotonic()
    log.info("%s: start", name)
    yield
    log.info("%s: done in %.1fs", name, time.monotonic() - t0)


def check_config(adata, library_column: str, tcr_column: str, atomic_axes: list[str]) -> list[str]:
    # Everything a bad parameter would break, reported together before any
    # expensive work. Mirrors what the later steps read.
    problems = []
    if library_column not in adata.obs:
        problems.append(f"library_column {library_column!r} is not an obs column.")
    if tcr_column and tcr_column not in adata.obs:
        problems.append(f"tcr_column {tcr_column!r} is not an obs column.")
    created = {"leiden", "library", "site", "patient"} | ({"tcr_status"} if tcr_column else set())
    for axis in atomic_axes:
        if axis not in created and axis not in adata.obs:
            problems.append(f"atomic axis {axis!r} is neither created by PREPARE nor an obs column.")
    return problems


def embed(adata, n_hvg: int, n_pcs: int, n_neighbors: int, leiden_resolution: float) -> None:
    work = adata.copy()
    with step("normalize/log1p"):
        sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)
    with step("HVG"):
        sc.pp.highly_variable_genes(work, n_top_genes=n_hvg, flavor="seurat")
        work = work[:, work.var["highly_variable"]].copy()
    with step("scale"):
        # Densify the HVG subset (n_cells x n_hvg) so scale does not warn about it.
        if sparse.issparse(work.X):
            work.X = work.X.toarray()
        sc.pp.scale(work, max_value=10)
    with step("PCA"):
        sc.pp.pca(work, n_comps=n_pcs, svd_solver="randomized", random_state=0)
    with step("neighbors"):
        sc.pp.neighbors(work, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=0)
    with step("UMAP"):
        sc.tl.umap(work, random_state=0)
    with step("Leiden"):
        sc.tl.leiden(work, resolution=leiden_resolution, flavor="igraph", n_iterations=2,
                     random_state=0)
    adata.obsm["X_pca"] = work.obsm["X_pca"].astype(np.float32)
    adata.obsm["X_umap"] = work.obsm["X_umap"].astype(np.float32)
    # cell2zarr atomic axes must be categorical. Keep scanpy's own categorical
    # (natural-sorted: "2" before "10"); re-wrapping would sort as text.
    adata.obs["leiden"] = work.obs["leiden"].values
    adata.var["highly_variable"] = adata.var_names.isin(work.var_names)


def finalize_indexes(adata) -> None:
    # An index name shared with a column breaks zarr writing once cell2zarr makes
    # the column categorical. Object dtype keeps anndata writing `_index` as a flat
    # array: pandas 3 infers StringDtype, which anndata writes as a values+mask
    # group that cell2zarr cannot read as an array.
    for df in (adata.obs, adata.var):
        df.index.name = None
        df.index = df.index.astype(object)


def check_prepared(adata, atomic_axes: list[str]) -> list[str]:
    problems = []
    for label, df in (("obs", adata.obs), ("var", adata.var)):
        name = df.index.name
        if name is not None and name in df.columns:
            problems.append(f"{label} index name {name!r} is also a {label} column.")
    groups = {}
    for col in adata.obs.columns:
        groups.setdefault(col.lower(), []).append(col)
    for cols in groups.values():
        if len(cols) > 1:
            problems.append(f"obs columns {cols} differ only by case.")
    for key in ("X_pca", "X_umap"):
        if key not in adata.obsm:
            problems.append(f"obsm[{key!r}] is missing.")
    for axis in atomic_axes:
        if axis not in adata.obs:
            problems.append(f"atomic axis {axis!r} is not an obs column.")
        elif not isinstance(adata.obs[axis].dtype, pd.CategoricalDtype):
            problems.append(f"atomic axis {axis!r} is {adata.obs[axis].dtype}, not categorical.")
    return problems


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(16 * 1024 * 1024):
            h.update(block)
    return h.hexdigest()


# Names only; the pins live in containers/prepare/requirements.txt. Kept here
# because Nextflow stages bin/ scripts without the rest of the repo.
PINNED_PACKAGES = ("scanpy", "anndata", "pandas", "numpy", "scipy", "h5py", "igraph",
                   "leidenalg", "umap-learn", "numba", "scikit-learn")


def pinned_versions() -> dict[str, str]:
    return {n: importlib.metadata.version(n) for n in PINNED_PACKAGES}


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--provenance", required=True)
    ap.add_argument("--n-hvg", type=int, default=3000)
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--leiden-resolution", type=float, default=1.0)
    ap.add_argument("--library-column", default="orig.ident")
    ap.add_argument("--site-regex", default=r"(DFCI-?\d+|JHU\d+|MSK\d+|GBM-001)")
    ap.add_argument("--patient-regex", default=r"GBM-001-(\d{3})")
    ap.add_argument("--patient-prefix", default="GBM001-P")
    ap.add_argument("--tcr-column", default="tcr_has_tcr_match", help="empty to skip")
    ap.add_argument("--atomic-axes", default="", help="comma-separated; empty for none")
    ap.add_argument("--source-uri", default="", help="original location of --input")
    ap.add_argument("--normalize", action="store_true")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--pipeline-revision", default="")
    ap.add_argument("--container", default="")
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    sc.settings.n_jobs = args.threads
    steps = []
    environment = pinned_versions()  # fail now, not after the write

    adata = sc.read_h5ad(args.input)
    steps.append(f"read {adata.n_obs:,} x {adata.n_vars:,}")

    axes = [a for a in args.atomic_axes.split(",") if a]
    problems = check_input(adata, args.normalize) + check_config(
        adata, args.library_column, args.tcr_column, axes)
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1

    adata.X = adata.X.astype(np.float32)
    steps.append("cast X to float32")

    dropped = merge_case_duplicate_columns(adata.obs)
    steps.append(f"merged case-duplicate obs columns, dropped {dropped}")
    n_unknown = add_identity_columns(
        adata, args.library_column, re.compile(args.site_regex),
        re.compile(args.patient_regex), args.patient_prefix,
    )
    steps.append(f"added library, site, patient; {n_unknown:,} cells with unknown site")
    add_tcr_status(adata, args.tcr_column)
    steps.append(f"tcr_status from {args.tcr_column!r}" if args.tcr_column else "tcr_status skipped")

    embed(adata, args.n_hvg, args.n_pcs, args.n_neighbors, args.leiden_resolution)
    steps.append(f"embedding: {args.n_hvg} HVGs, {args.n_pcs} PCs, "
                 f"{args.n_neighbors} neighbors, UMAP, Leiden resolution {args.leiden_resolution}")
    finalize_indexes(adata)
    no_doublets = "doublet flags not added: no singlets companion (merge-only path)"
    log.warning(no_doublets)
    steps.append(no_doublets)

    problems = check_prepared(adata, axes)
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1

    out = Path(args.output)
    adata.write_h5ad(out)
    steps.append(f"wrote {out.name}")

    src = Path(args.input)
    provenance = {
        "file": {"name": out.name, "size_bytes": out.stat().st_size, "sha256": sha256_of(out),
                 "n_obs": adata.n_obs, "n_vars": adata.n_vars},
        "source": {"uri": args.source_uri, "path": args.input, "size_bytes": src.stat().st_size,
                   "sha256": sha256_of(src)},
        "parameters": vars(args),
        "pipeline": {"revision": args.pipeline_revision, "container": args.container},
        "environment": environment,
        "steps": steps,
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    Path(args.provenance).write_text(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
