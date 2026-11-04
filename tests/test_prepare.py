"""Unit tests for bin/prepare_h5ad.py."""

import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

SPEC = importlib.util.spec_from_file_location(
    "prepare_h5ad", Path(__file__).parent.parent / "bin" / "prepare_h5ad.py"
)
prep = importlib.util.module_from_spec(SPEC)
sys.modules["prepare_h5ad"] = prep
SPEC.loader.exec_module(prep)

FIXTURE = str(Path(__file__).parent / "data" / "tiny_merge.h5ad")
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


def test_check_input_rejects_non_integer_x_when_normalizing():
    adata = ad.AnnData(X=np.array([[0.5, 1.0], [2.0, 0.0]], dtype=np.float32))
    assert any("integer" in p for p in prep.check_input(adata, normalize=True))
    assert prep.check_input(adata, normalize=False) == []


def test_check_input_checks_every_block(monkeypatch):
    monkeypatch.setattr(prep, "BLOCK", 4)
    x = np.ones((3, 4), dtype=np.float32)
    x[2, 3] = 0.5  # beyond the first block
    assert prep.check_input(ad.AnnData(X=x), normalize=True)
    x[2, 3] = -1.0
    assert prep.check_input(ad.AnnData(X=x), normalize=True)
    assert prep.check_input(ad.AnnData(X=np.ones((3, 4), dtype=np.float32)), normalize=True) == []


def test_pinned_versions_match_requirements_names():
    reqs = Path(__file__).parent.parent / "containers" / "prepare" / "requirements.txt"
    names = {re.split(r"[=<>]", l)[0] for l in reqs.read_text().splitlines() if l.strip()}
    assert set(prep.pinned_versions()) == names


def test_check_prepared_names_every_problem():
    obs = pd.DataFrame(
        {"lib": ["a", "b"], "x": [1, 2], "X": [1, 2], "flag": [True, False]},
        index=pd.Index(["c0", "c1"], name="lib"),
    )
    adata = ad.AnnData(X=np.zeros((2, 2), dtype=np.float32), obs=obs)
    text = "\n".join(prep.check_prepared(adata, ["leiden", "flag"]))
    for needle in ("'lib'", "'x'", "X_pca", "X_umap", "'leiden'", "'flag'"):
        assert needle in text


def cli_args(tmp_path, *extra):
    return ["--input", FIXTURE, "--output", str(tmp_path / "out.h5ad"),
            "--provenance", str(tmp_path / "provenance.json"),
            "--n-hvg", "20", "--n-pcs", "10", "--n-neighbors", "10", *extra]


def test_cli_end_to_end_on_fixture(tmp_path):
    rc = prep.main(["--input", FIXTURE, "--output", str(tmp_path / "out.h5ad"),
                    "--provenance", str(tmp_path / "provenance.json"),
                    "--n-hvg", "20", "--n-pcs", "10", "--n-neighbors", "10",
                    "--atomic-axes", "leiden,library,patient,site,tcr_status",
                    "--normalize", "--threads", "2",
                    "--source-uri", "s3://bucket/merged.h5ad",
                    "--pipeline-revision", "test", "--container", "test"])
    assert rc == 0
    out = ad.read_h5ad(tmp_path / "out.h5ad")
    assert out.X.dtype == np.float32
    src = ad.read_h5ad(FIXTURE)
    assert np.array_equal(out.X.toarray(), src.X.astype(np.float32).toarray())
    assert list(out.obs["leiden"].cat.categories) == sorted(
        out.obs["leiden"].cat.categories, key=int)
    assert "Modality" not in out.obs and "modality" in out.obs
    assert {"library", "site", "patient", "tcr_status", "leiden"} <= set(out.obs)
    assert out.obsm["X_pca"].shape == (200, 10) and out.obsm["X_umap"].shape == (200, 2)
    assert out.obs.index.name is None and out.var.index.name is None
    with h5py.File(tmp_path / "out.h5ad") as f:  # Review Focus 1
        assert isinstance(f["obs"]["_index"], h5py.Dataset)
        assert isinstance(f["var"]["_index"], h5py.Dataset)
    prov = json.loads((tmp_path / "provenance.json").read_text())
    assert prov["file"]["sha256"] == hashlib.sha256((tmp_path / "out.h5ad").read_bytes()).hexdigest()
    assert prov["file"]["n_obs"] == 200
    assert prov["environment"]["scanpy"] == "1.12.5"
    assert prov["source"]["uri"] == "s3://bucket/merged.h5ad"
    assert prov["source"]["path"] == FIXTURE
    assert "doublet flags not added: no singlets companion (merge-only path)" in prov["steps"]


def test_cli_fails_fast_on_bad_input(tmp_path):
    src = ad.read_h5ad(FIXTURE)
    src.X = src.X * 0.5
    src.write_h5ad(tmp_path / "half.h5ad")
    rc = prep.main(["--input", str(tmp_path / "half.h5ad"),
                    "--output", str(tmp_path / "out.h5ad"),
                    "--provenance", str(tmp_path / "provenance.json"), "--normalize"])
    assert rc != 0
    assert not (tmp_path / "out.h5ad").exists()
    assert not (tmp_path / "provenance.json").exists()


def test_check_config_reports_everything_at_once():
    adata = make_adata(x=["a", "b"])
    problems = prep.check_config(adata, "orig.ident", "tcr_col",
                                 ["doublet_status", "leiden", "x", "tcr_status", "nope"])
    text = "\n".join(problems)
    for needle in ("'orig.ident'", "'tcr_col'", "'doublet_status'", "'nope'"):
        assert needle in text
    assert "'x'" not in text and "'leiden'" not in text
    # tcr_status is only created when a tcr column is named
    assert any("'tcr_status'" in p for p in prep.check_config(adata, "x", "", ["tcr_status"]))
    assert prep.check_config(adata, "x", "x", ["tcr_status", "library", "site", "patient"]) == []


def test_cli_fails_before_embedding_on_unknown_axis(tmp_path, monkeypatch, capsys):
    def boom(*a, **k):
        raise AssertionError("embedding started")
    monkeypatch.setattr(prep, "embed", boom)
    rc = prep.main(cli_args(tmp_path, "--atomic-axes", "leiden,doublet_status"))
    assert rc != 0
    assert "doublet_status" in capsys.readouterr().err
    assert not (tmp_path / "out.h5ad").exists()
    assert not (tmp_path / "provenance.json").exists()


def test_cli_fails_on_missing_library_and_tcr_columns(tmp_path, capsys):
    rc = prep.main(cli_args(tmp_path, "--library-column", "nolib", "--tcr-column", "notcr"))
    err = capsys.readouterr().err
    assert rc != 0 and "'nolib'" in err and "'notcr'" in err


def test_doublet_step_is_warned(tmp_path, caplog):
    with caplog.at_level("WARNING", logger="prepare_h5ad"):
        assert prep.main(cli_args(tmp_path)) == 0
    assert "doublet flags not added" in caplog.text


def test_tcr_status_warns_when_nothing_is_positive(caplog):
    adata = make_adata(tcr=["false", "0", "no"])
    with caplog.at_level("INFO", logger="prepare_h5ad"):
        prep.add_tcr_status(adata, "tcr")
    assert any(r.levelname == "WARNING" and "0 TCR+" in r.message for r in caplog.records)


def test_tcr_status_logs_positive_count(caplog):
    adata = make_adata(tcr=[True, False, True])
    with caplog.at_level("INFO", logger="prepare_h5ad"):
        prep.add_tcr_status(adata, "tcr")
    assert "2 TCR+ of 3" in caplog.text
    assert not any(r.levelname == "WARNING" for r in caplog.records)


def test_embed_logs_each_step_with_elapsed_seconds(tmp_path, caplog):
    with caplog.at_level("INFO", logger="prepare_h5ad"):
        assert prep.main(cli_args(tmp_path)) == 0
    for step in ("normalize/log1p", "HVG", "scale", "PCA", "neighbors", "UMAP", "Leiden"):
        assert re.search(rf"{re.escape(step)}.*\d+(\.\d+)?s", caplog.text), step
