"""Unit tests for bin/prepare_h5ad.py."""

import importlib.util
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

SPEC = importlib.util.spec_from_file_location(
    "prepare_h5ad", Path(__file__).parent.parent / "bin" / "prepare_h5ad.py"
)
prep = importlib.util.module_from_spec(SPEC)
sys.modules["prepare_h5ad"] = prep
SPEC.loader.exec_module(prep)

SITE = re.compile(r"(DFCI-?\d+|JHU\d+|MSK\d+|GBM-001)")
PAT = re.compile(r"GBM-001-(\d{3})")


def make_adata(**cols):
    n = len(next(iter(cols.values())))
    obs = pd.DataFrame(cols, index=[f"c{i}" for i in range(n)])
    return ad.AnnData(X=np.zeros((n, 2), dtype=np.float32), obs=obs)


def test_site_and_patient_from_library_names():
    cases = {
        "BTC-GBM-001-002-GEX_S2_L002": ("GBM001", "GBM001-P002"),
        "GBM1.DFCI1.S1": ("DFCI1", "DFCI1"),
        "BTC-DFCI-7-S4-GEX-LIB": ("DFCI7", "DFCI7"),
        "JHU1-S1": ("JHU1", "JHU1"),
        "MSK1-TCR": ("MSK1", "MSK1"),
    }
    for lib, (site, patient) in cases.items():
        assert prep.site_of(lib, SITE) == site
        assert prep.patient_of(lib, site, PAT, "GBM001-P") == patient


def test_unknown_library_is_counted_not_fatal():
    adata = make_adata(**{"orig.ident": ["GBM1.DFCI1.S1", "NEWSITE-S1", "NEWSITE-S2"]})
    n_unknown = prep.add_identity_columns(adata, "orig.ident", SITE, PAT, "GBM001-P")
    assert n_unknown == 2
    assert list(adata.obs["site"]) == ["DFCI1", "unknown", "unknown"]
    assert list(adata.obs["patient"]) == ["DFCI1", "unknown", "unknown"]
    assert isinstance(adata.obs["site"].dtype, pd.CategoricalDtype)


def test_tcr_status_from_bools_and_strings():
    adata = make_adata(
        tcr_has_tcr_match=pd.Series([True, False, "True", "False", "", None], dtype=object).tolist()
    )
    prep.add_tcr_status(adata, "tcr_has_tcr_match")
    assert list(adata.obs["tcr_status"]) == ["TCR+", "TCR-", "TCR+", "TCR-", "TCR-", "TCR-"]
    assert list(adata.obs["tcr_status"].cat.categories) == ["TCR-", "TCR+"]


def test_tcr_status_skipped_when_column_blank():
    adata = make_adata(x=[1, 2])
    prep.add_tcr_status(adata, "")
    assert "tcr_status" not in adata.obs


def test_case_duplicates_merge_without_inventing_missing_values():
    obs = pd.DataFrame(
        {"modality": ["GEX", "NA", "NA", "NA"], "Modality": ["NA", "GEX", None, "NA"]}
    )
    dropped = prep.merge_case_duplicate_columns(obs)
    assert dropped == ["Modality"]
    assert list(obs["modality"]) == ["GEX", "GEX", "NA", "NA"]


def test_case_duplicates_that_conflict_stop_the_run():
    obs = pd.DataFrame({"a": ["x", "y"], "A": ["x", "z"]})
    with pytest.raises(SystemExit, match="'a' and 'A' disagree on 1 cells"):
        prep.merge_case_duplicate_columns(obs)


def test_case_duplicates_keep_numeric_dtype():
    nan = np.nan
    obs = pd.DataFrame(
        {"score": [1.0, nan, 3.0], "Score": [1.0, 2.0, nan], "SCORE": [nan, nan, 3.0]}
    )
    prep.merge_case_duplicate_columns(obs)
    assert list(obs.columns) == ["score"]
    assert obs["score"].tolist() == [1.0, 2.0, 3.0]
    assert obs["score"].dtype.kind == "f"
