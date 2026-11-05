"""Unit tests for the Cirro preprocess script."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SPEC = importlib.util.spec_from_file_location(
    "preprocess", Path(__file__).parent.parent / ".cirro" / "preprocess.py"
)
preprocess = importlib.util.module_from_spec(SPEC)
sys.modules["preprocess"] = preprocess
SPEC.loader.exec_module(preprocess)


class FakeLogger:
    def __init__(self):
        self.messages = []

    def info(self, msg):
        self.messages.append(msg)


class FakeDataset:
    def __init__(self, files, params):
        self.files = pd.DataFrame({"file": files})
        self.params = dict(params)
        self.logger = FakeLogger()

    # Mirrors cirro.helpers.preprocess_dataset.PreprocessDataset (cirro 1.10.1),
    # including its asserts, so a call the real SDK would reject fails here too.
    def add_param(self, name, value, overwrite=False, log=True):
        assert overwrite or name not in self.params, \
            f"Cannot add parameter {name}, already exists (and overwrite=False)"
        self.params[name] = value

    def remove_param(self, name, force=False):
        assert force or name in self.params, \
            f"Cannot remove parameter {name}, does not exist (and force=False)"
        if name in self.params:
            del self.params[name]


def test_resolves_the_single_h5ad():
    ds = FakeDataset(["s3://bucket/project/atlas.h5ad"], {"two_phase": True})
    preprocess.run(ds)
    assert ds.params["input"] == "s3://bucket/project/atlas.h5ad"


def test_rejects_loom_input():
    ds = FakeDataset(["s3://bucket/project/atlas.loom"], {"two_phase": True})
    with pytest.raises(ValueError, match="no .h5ad file"):
        preprocess.run(ds)


def test_rejects_multiple_h5ad_files():
    ds = FakeDataset(
        ["s3://bucket/a.h5ad", "s3://bucket/b.h5ad"], {"two_phase": True}
    )
    with pytest.raises(ValueError, match="exactly one"):
        preprocess.run(ds)


def test_rejects_two_phase_only_options_under_single_phase():
    ds = FakeDataset(
        ["s3://bucket/a.h5ad"],
        {"two_phase": False, "shard_size": 30, "dtype": "float16", "normalize": True},
    )
    with pytest.raises(ValueError) as excinfo:
        preprocess.run(ds)

    message = str(excinfo.value)
    assert "shard_size" in message
    assert "dtype" in message
    assert "normalize" in message


def test_single_phase_without_those_options_is_fine():
    ds = FakeDataset(["s3://bucket/a.h5ad"], {"two_phase": False})
    preprocess.run(ds)
    assert ds.params["input"] == "s3://bucket/a.h5ad"


def test_blank_optional_params_are_removed():
    ds = FakeDataset(
        ["s3://bucket/a.h5ad"],
        {"two_phase": True, "atomic_axes": "", "add_layers": "   ", "coarse": None},
    )
    preprocess.run(ds)

    assert "atomic_axes" not in ds.params
    assert "add_layers" not in ds.params
    assert "coarse" not in ds.params


def test_strata_axes_are_preserved_when_set():
    ds = FakeDataset(
        ["s3://bucket/a.h5ad"],
        {"two_phase": True, "atomic_axes": " cell_type , treatment "},
    )
    preprocess.run(ds)
    assert ds.params["atomic_axes"] == "cell_type,treatment"


def test_prepare_selects_h5ad_by_name():
    ds = FakeDataset(
        ["s3://b/d/merged.h5ad", "s3://b/d/merged.tcr_only.h5ad"],
        {"two_phase": True, "prepare": True, "h5ad_name": "merged.h5ad"},
    )
    preprocess.run(ds)
    assert ds.params["input"] == "s3://b/d/merged.h5ad"
    assert "h5ad_name" not in ds.params


def test_prepare_with_missing_name_lists_what_exists():
    ds = FakeDataset(
        ["s3://b/d/a.h5ad", "s3://b/d/b.h5ad"],
        {"two_phase": True, "prepare": True, "h5ad_name": "merged.h5ad"},
    )
    with pytest.raises(ValueError, match=r"merged\.h5ad.*a\.h5ad.*b\.h5ad"):
        preprocess.run(ds)


def test_without_prepare_multiple_h5ads_still_rejected():
    ds = FakeDataset(
        ["s3://b/a.h5ad", "s3://b/b.h5ad"],
        {"two_phase": True, "prepare": False, "h5ad_name": "a.h5ad"},
    )
    with pytest.raises(ValueError, match="exactly one"):
        preprocess.run(ds)


def test_h5ad_name_is_removed_even_without_prepare():
    ds = FakeDataset(
        ["s3://b/a.h5ad"],
        {"two_phase": True, "prepare": False, "h5ad_name": "merged.h5ad"},
    )
    preprocess.run(ds)
    assert "h5ad_name" not in ds.params


def test_blank_tcr_column_passes_through_as_empty_string():
    # "" means "skip tcr_status" to Nextflow, so it must not be removed.
    ds = FakeDataset(
        ["s3://b/a.h5ad"], {"two_phase": True, "tcr_column": ""}
    )
    preprocess.run(ds)
    assert ds.params["tcr_column"] == ""


def test_regex_params_are_not_mangled():
    regex = r"GBM-001-(\d{1,3}) | x"
    ds = FakeDataset(
        ["s3://b/a.h5ad"],
        {"two_phase": True, "site_regex": regex, "patient_regex": regex},
    )
    preprocess.run(ds)
    assert ds.params["site_regex"] == regex
    assert ds.params["patient_regex"] == regex


def test_blank_prepare_params_are_removed_so_defaults_apply():
    ds = FakeDataset(
        ["s3://b/a.h5ad"],
        {
            "two_phase": True,
            "n_hvg": None,
            "n_pcs": "",
            "library_column": "  ",
            "site_regex": "",
            "n_neighbors": 15,
        },
    )
    preprocess.run(ds)
    for name in ("n_hvg", "n_pcs", "library_column", "site_regex"):
        assert name not in ds.params
    assert ds.params["n_neighbors"] == 15


def test_prepare_false_with_no_h5ad_name_key_at_all():
    ds = FakeDataset(["s3://b/a.h5ad"], {"two_phase": True, "prepare": False})
    preprocess.run(ds)
    assert ds.params["input"] == "s3://b/a.h5ad"


def test_run_with_atomic_axes_set_succeeds():
    ds = FakeDataset(
        ["s3://b/a.h5ad"],
        {"two_phase": True, "atomic_axes": "cell_type", "add_layers": "counts", "coarse": "a;b"},
    )
    preprocess.run(ds)
    assert ds.params["atomic_axes"] == "cell_type"


def test_run_overwrites_a_preexisting_input_param():
    ds = FakeDataset(["s3://b/a.h5ad"], {"two_phase": True, "input": "s3://stale"})
    preprocess.run(ds)
    assert ds.params["input"] == "s3://b/a.h5ad"


def test_prepare_rejects_ambiguous_h5ad_name():
    ds = FakeDataset(
        ["s3://b/x/merged.h5ad", "s3://b/y/merged.h5ad"],
        {"two_phase": True, "prepare": True, "h5ad_name": "merged.h5ad"},
    )
    with pytest.raises(ValueError, match=r"more than one.*x/merged\.h5ad.*y/merged\.h5ad"):
        preprocess.run(ds)


# A dataset produced by another pipeline (e.g. sazaidi-scratch-tcrmerge, which
# has no file-mapping rules) reaches preprocess with an empty files.csv. The
# input dataset folders arrive instead as input_paths ($.inputs[*].dataPath).
TCRMERGE = "s3://project-x/datasets/cee37219/data"
TCRMERGE_FILES = [
    f"{TCRMERGE}/merged.h5ad",
    f"{TCRMERGE}/merged.tcr_only.h5ad",
    f"{TCRMERGE}/tables/merged_obs.tsv",
]


@pytest.fixture
def s3_listing(monkeypatch):
    calls = []

    def fake_list(prefix):
        calls.append(prefix)
        return [f for f in TCRMERGE_FILES if f.startswith(prefix.rstrip("/") + "/")]

    monkeypatch.setattr(preprocess, "_list_s3", fake_list)
    return calls


def test_pipeline_output_input_is_listed_from_input_paths(s3_listing):
    ds = FakeDataset([], {"two_phase": True, "prepare": True, "h5ad_name": "merged.h5ad",
                          "input_paths": [TCRMERGE]})
    preprocess.run(ds)
    assert ds.params["input"] == f"{TCRMERGE}/merged.h5ad"
    assert s3_listing == [TCRMERGE]


def test_input_paths_given_as_a_single_string(s3_listing):
    ds = FakeDataset([], {"two_phase": True, "prepare": True, "h5ad_name": "merged.h5ad",
                          "input_paths": TCRMERGE})
    preprocess.run(ds)
    assert ds.params["input"] == f"{TCRMERGE}/merged.h5ad"


def test_input_paths_never_reaches_nextflow(s3_listing):
    listed = FakeDataset([], {"two_phase": True, "prepare": True, "h5ad_name": "merged.h5ad",
                              "input_paths": [TCRMERGE]})
    uploaded = FakeDataset(["s3://b/a.h5ad"], {"two_phase": True, "input_paths": ["s3://b"]})
    for ds in (listed, uploaded):
        preprocess.run(ds)
        assert "input_paths" not in ds.params


def test_files_csv_wins_over_input_paths(s3_listing):
    ds = FakeDataset(["s3://b/a.h5ad"], {"two_phase": True, "input_paths": [TCRMERGE]})
    preprocess.run(ds)
    assert ds.params["input"] == "s3://b/a.h5ad"
    assert s3_listing == []


def test_pipeline_output_without_prepare_still_needs_exactly_one_h5ad(s3_listing):
    ds = FakeDataset([], {"two_phase": True, "input_paths": [TCRMERGE]})
    with pytest.raises(ValueError, match="2 .h5ad files"):
        preprocess.run(ds)


def test_empty_files_and_no_input_paths_reports_no_h5ad(s3_listing):
    ds = FakeDataset([], {"two_phase": True})
    with pytest.raises(ValueError, match="no .h5ad file"):
        preprocess.run(ds)
