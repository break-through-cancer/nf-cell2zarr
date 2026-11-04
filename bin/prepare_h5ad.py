#!/usr/bin/env python3
"""Prepare a merged h5ad for cell2zarr / Cell Explorer.

Ported from the on-cluster prepare_gbm.py. This part adds the identity columns
(library, site, patient), the TCR status and merges obs columns whose names
differ only by case.
"""

import logging
import re

import pandas as pd

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
