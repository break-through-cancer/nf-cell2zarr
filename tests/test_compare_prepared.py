"""Tests for scripts/compare_prepared.py."""

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

SPEC = importlib.util.spec_from_file_location(
    "compare_prepared", Path(__file__).parent.parent / "scripts" / "compare_prepared.py"
)
cmp = importlib.util.module_from_spec(SPEC)
sys.modules["compare_prepared"] = cmp
SPEC.loader.exec_module(cmp)

N = 6


def make_obs(reference: bool, **overrides):
    names = [f"cell{i}" for i in range(N)]
    obs = pd.DataFrame(
        {
            "library": pd.Categorical(["L1", "L1", "L2", "L2", "L3", "L3"]),
            "site": pd.Categorical(["DFCI1", "DFCI1", "JHU2", "JHU2", "MSK3", "MSK3"]),
            "patient": pd.Categorical(["P1", "P1", "P2", "P2", "P3", "P3"]),
            "tcr_status": pd.Categorical(["TCR-", "TCR+"] * 3, categories=["TCR-", "TCR+"]),
            "leiden": pd.Categorical(["0", "1", "0", "1", "0", "2"]),
        }
    )
    if reference:
        # The older script's quirks: a named nullable-string index (written as a
        # values+mask group) and a capitalised duplicate of modality.
        obs["modality"] = ["RNA", "", "RNA", "RNA", np.nan, "RNA"]
        obs["Modality"] = ["", "ATAC", "", "", "ATAC", ""]
        obs.index = pd.Index(names, dtype="string", name="bpcells_name")
    else:
        obs["modality"] = ["RNA", "ATAC", "RNA", "RNA", "ATAC", "RNA"]
        obs.index = pd.Index(names, dtype=object)
    for k, v in overrides.items():
        obs[k] = v
    return obs


def write(path, obs, reference: bool):
    rng = np.random.default_rng(0)
    X = sparse.random(N, 5, density=0.6, format="csr", random_state=rng, dtype=np.float32)
    var = pd.DataFrame(
        index=pd.Index([f"g{i}" for i in range(5)], dtype="string" if reference else object)
    )
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(path)
    return X


@pytest.fixture
def pair(tmp_path):
    ref, cand = tmp_path / "ref.h5ad", tmp_path / "cand.h5ad"
    ref_obs = make_obs(True)
    write(ref, ref_obs, True)
    # Same X: rewrite candidate from the reference's matrix.
    a = ad.read_h5ad(ref)
    cand_obs = make_obs(False)
    ad.AnnData(
        X=a.X, obs=cand_obs, var=pd.DataFrame(index=pd.Index(list(a.var_names), dtype=object))
    ).write_h5ad(cand)
    return ref, cand


def run(capsys, *args):
    code = cmp.main([str(a) for a in args])
    return code, capsys.readouterr().out


def test_identical_files_pass(pair, capsys):
    code, out = run(capsys, *pair)
    assert code == 0, out
    assert "FAIL" not in out


def test_reference_uses_named_nullable_index(pair):
    import h5py

    with h5py.File(pair[0]) as f, h5py.File(pair[1]) as g:
        assert f["obs"].attrs["_index"] == "bpcells_name"
        assert isinstance(f["obs"]["bpcells_name"], h5py.Group)
        assert isinstance(g["obs"]["_index"], h5py.Dataset)


def test_changed_site_value_fails_naming_site(pair, capsys):
    ref, cand = pair
    a = ad.read_h5ad(cand)
    a.obs["site"] = pd.Categorical(["DFCI1", "DFCI1", "JHU2", "JHU2", "MSK3", "DFCI1"])
    a.write_h5ad(cand)
    code, out = run(capsys, ref, cand)
    assert code == 1
    assert "site" in out


def test_changed_x_value_fails_reporting_location(pair, capsys):
    ref, cand = pair
    a = ad.read_h5ad(cand)
    a.X.data[2] += 1
    a.write_h5ad(cand)
    code, out = run(capsys, ref, cand)
    assert code == 1
    assert "FAIL X" in out and "data[2]" in out


def test_modality_mismatch_fails(pair, capsys):
    ref, cand = pair
    a = ad.read_h5ad(cand)
    a.obs["modality"] = ["RNA", "RNA", "RNA", "RNA", "ATAC", "RNA"]
    a.write_h5ad(cand)
    code, out = run(capsys, ref, cand)
    assert code == 1
    assert "modality" in out


def test_doublet_columns_ignored_with_skip_flag(pair, capsys):
    ref, cand = pair
    a = ad.read_h5ad(ref)
    a.obs["is_doublet"] = [False] * N
    a.write_h5ad(ref)
    code, _ = run(capsys, ref, cand, "--skip-doublets")
    assert code == 0
    code, out = run(capsys, ref, cand)
    assert code == 1
    assert "is_doublet" in out


def test_x_compared_across_blocks(pair, capsys, monkeypatch):
    monkeypatch.setattr(cmp, "BLOCK", 2)
    code, out = run(capsys, *pair)
    assert code == 0, out


def test_inlined_merge_rule_agrees_with_prepare_h5ad():
    spec = importlib.util.spec_from_file_location(
        "prepare_h5ad", Path(__file__).parent.parent / "bin" / "prepare_h5ad.py"
    )
    prep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prep)

    frames = [
        pd.DataFrame({"modality": ["GEX", "NA", "NA", "NA"], "Modality": ["NA", "GEX", None, "NA"]}),
        pd.DataFrame({"modality": ["RNA", "", np.nan, "None", "nan", "RNA"],
                      "Modality": ["", "ATAC", "ATAC", "ATAC", "ATAC", "RNA"]}),
        make_obs(True)[["modality", "Modality"]],
        ad.read_h5ad(Path(__file__).parent / "data" / "tiny_merge.h5ad").obs[["modality", "Modality"]],
    ]
    for frame in frames:
        expected = frame.copy()
        prep.merge_case_duplicate_columns(expected)
        got = cmp.merge_modality(frame["modality"], frame["Modality"])
        assert list(got.astype(str)) == list(expected["modality"].astype(str))


def test_script_does_not_depend_on_the_repo_layout(tmp_path):
    import shutil
    import subprocess

    alone = tmp_path / "compare_prepared.py"
    shutil.copy(Path(cmp.__file__), alone)
    out = subprocess.run([sys.executable, str(alone), "--help"], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "importlib" not in alone.read_text()  # no loading of ../bin by path
