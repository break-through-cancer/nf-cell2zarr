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

    def add_param(self, key, value, overwrite=True):
        self.params[key] = value

    def remove_param(self, key, force=False):
        self.params.pop(key, None)


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
