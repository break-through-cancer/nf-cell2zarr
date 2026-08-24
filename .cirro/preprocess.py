#!/usr/bin/env python3
"""Cirro preprocess script for nf-cell2zarr.

Runs on the headnode after the user submits the run form and before Nextflow
launches. Resolves the input h5ad and rejects parameter combinations that
cell2zarr would otherwise silently ignore.
"""

# Options that cell2zarr's CLI accepts but silently drops under single-phase
# conversion. Rather than reproduce that silence, reject them up front and name
# every offender at once.
TWO_PHASE_ONLY = (
    "cell_chunk_size",
    "shard_size",
    "dtype",
    "obsm_cell_chunk_size",
    "encoding_config",
    "normalize",
)

# Blank strings from an untouched optional form field must not reach Nextflow as
# the literal "", which would read as "the user asked for this".
OPTIONAL = (
    "add_layers",
    "atomic_axes",
    "coarse",
    "max_atomic_cardinality",
    "n_top_genes",
    "encoding_config",
)


def _resolve_input(ds):
    """Return the single .h5ad path in the input dataset."""
    paths = [str(p) for p in ds.files["file"].tolist()]
    h5ads = [p for p in paths if p.endswith(".h5ad")]

    if not h5ads:
        raise ValueError(
            "The input dataset contains no .h5ad file. nf-cell2zarr accepts "
            "AnnData only; Loom (.loom) and Seurat (.h5seurat) files in a "
            f"Single-Cell Dataset are not supported. Found: {paths}"
        )

    if len(h5ads) > 1:
        raise ValueError(
            "nf-cell2zarr converts exactly one dataset per run, but the input "
            f"contains {len(h5ads)} .h5ad files: {h5ads}. Run it once per file."
        )

    return h5ads[0]


def _check_conversion_mode(ds):
    """Reject two-phase-only options under a single-phase run."""
    if ds.params.get("two_phase", True):
        return

    offenders = [
        name
        for name in TWO_PHASE_ONLY
        if ds.params.get(name) not in (None, "", False)
    ]
    if offenders:
        raise ValueError(
            "These options only apply to two-phase conversion and would be "
            f"silently ignored: {', '.join(offenders)}. Either clear them or "
            "enable two-phase conversion."
        )


def _clean_optional(ds):
    """Drop blank optional params; normalise the ones that survive."""
    for name in OPTIONAL:
        value = ds.params.get(name)
        if value is None or (isinstance(value, str) and not value.strip()):
            ds.remove_param(name)
        elif isinstance(value, str):
            # " cell_type , treatment " -> "cell_type,treatment"
            ds.add_param(name, ",".join(p.strip() for p in value.split(",") if p.strip()))


def run(ds):
    """Entry point, factored out so it can be unit tested."""
    _check_conversion_mode(ds)
    _clean_optional(ds)

    h5ad = _resolve_input(ds)
    ds.add_param("input", h5ad)

    ds.logger.info(f"nf-cell2zarr input: {h5ad}")
    ds.logger.info(f"nf-cell2zarr params: {ds.params}")

    return ds


if __name__ == "__main__":
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    run(PreprocessDataset.from_running())
