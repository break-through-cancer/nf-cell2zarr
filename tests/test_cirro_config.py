"""Checks on the Cirro configuration that only a real Cirro run would otherwise catch.

The first Cirro run (2026-11-05) failed on both of the first two: the input
mapping resolved to an empty params.json, and AWS Batch rejected the
tag-plus-digest image names.
"""

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
CIRRO = ROOT / ".cirro"


def _form_paths(schema, prefix=()):
    for name, prop in schema.get("properties", {}).items():
        path = prefix + (name,)
        yield path
        yield from _form_paths(prop, path)


def _input_mapping():
    return json.loads((CIRRO / "process-input.json").read_text())


def _form():
    form = json.loads((CIRRO / "process-form.json").read_text())
    return form.get("form", form)


# Keys Cirro supplies from the run itself rather than from the form.
NOT_FROM_FORM = {"outdir", "input_paths"}


def test_every_form_param_maps_from_dataset_params():
    """Cirro resolves form values under $.dataset.params; any other root is empty."""
    form_paths = set(_form_paths(_form()))
    mapping = _input_mapping()
    bad = {}
    for key, expr in mapping.items():
        if key in NOT_FROM_FORM:
            continue
        prefix = "$.dataset.params."
        if not expr.startswith(prefix):
            bad[key] = f"{expr} does not start with {prefix}"
        elif tuple(expr[len(prefix):].split(".")) not in form_paths:
            bad[key] = f"{expr} names no field in process-form.json"
    assert not bad, bad


def test_every_form_field_is_mapped():
    """A form field nobody maps is silently ignored on Cirro."""
    sections = {p for p in _form_paths(_form()) if len(p) == 1}
    fields = {p for p in _form_paths(_form()) if len(p) == 2}
    mapped = {
        tuple(expr.removeprefix("$.dataset.params.").split("."))
        for expr in _input_mapping().values()
    }
    assert sections and fields
    assert fields - mapped == set()


def test_outdir_is_the_output_dataset():
    """Without this, outputs publish to ./results on the headnode and are lost."""
    assert _input_mapping().get("outdir") == "$.dataset.dataPath"


def test_input_paths_are_the_input_datasets():
    """files.csv is empty for a dataset written by another pipeline (the GBM merge);
    preprocess then lists these folders."""
    assert _input_mapping().get("input_paths") == "$.inputs[*].dataPath"


# Copied from Nextflow 25.10.8, modules/nextflow/src/main/groovy/nextflow/container/
# ContainerNameValidator.groovy, which the AWS Batch executor applies to every task's
# image: a tag OR a digest, never both.
_TAG = r"[a-zA-Z0-9-_.]+"
_SHA = r"[a-fA-F0-9]{64}"
_COMPONENT = r"[a-z0-9-_.]+"
_REGISTRY = r"[a-zA-Z0-9-.]+(:[0-9]+)?"
_IMAGE = rf"({_REGISTRY}/)?{_COMPONENT}(/{_COMPONENT})*"
AWS_BATCH_IMAGE = re.compile(rf"^{_IMAGE}((:{_TAG})|(@sha256:{_SHA}))?$")


def _config_images():
    text = (ROOT / "nextflow.config").read_text()
    return dict(re.findall(r"^\s*(container|prepare_container)\s*=\s*'([^']+)'", text, re.M))


def test_both_images_are_configured():
    assert set(_config_images()) == {"container", "prepare_container"}


@pytest.mark.parametrize("name", ["container", "prepare_container"])
def test_image_names_pass_the_aws_batch_validator(name):
    image = _config_images()[name]
    assert AWS_BATCH_IMAGE.match(image), image


@pytest.mark.parametrize("name", ["container", "prepare_container"])
def test_images_are_pinned_by_digest(name):
    assert re.search(r"@sha256:[a-f0-9]{64}$", _config_images()[name])
