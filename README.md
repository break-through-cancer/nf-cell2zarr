# nf-cell2zarr

Nextflow pipeline that converts single-cell `h5ad` datasets to Zarr v3 stores for
[cell-explorer](https://github.com/cBioPortal/cell-explorer-py), designed to run on the
[Cirro](https://cirro.bio) platform (Nextflow on AWS Batch).

For the on-prem *iris* SLURM cluster, use `cell-explorer-scripts` instead. The two are
independent by design.

## Run locally

```bash
uv run --with anndata --with scipy python tests/data/make_fixture.py
nextflow run . -profile test,docker
```

The container image is `linux/amd64` only. On Apple Silicon, Docker runs it under
emulation — fine for tests, and irrelevant on AWS Batch.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `input` | — | Path to an `.h5ad` file. Required. |
| `outdir` | `results` | Where the finished store is published. |
| `prepare` | `false` | Run the optional PREPARE stage (embed and annotate a raw merge) before conversion. |
| `prepare_container` | `ghcr.io/break-through-cancer/nf-cell2zarr-prepare@sha256:be57110e…` (tag `sha-b1923f2`) | PREPARE image, pinned by digest (AWS Batch rejects tag@digest). Published by `.github/workflows/prepare-image.yml`. |
| `n_hvg` | `3000` | Highly variable genes used for PCA (prepare). |
| `n_pcs` | `50` | Principal components computed (prepare). |
| `n_neighbors` | `15` | Neighbours in the kNN graph behind UMAP and Leiden (prepare). |
| `leiden_resolution` | `1.0` | Leiden clustering resolution (prepare). |
| `library_column` | `orig.ident` | obs column holding the library name; `library`, `site` and `patient` derive from it (prepare). |
| `site_regex` | `(DFCI-?\d+\|JHU\d+\|MSK\d+\|GBM-001)` | Regex whose group 1 is the site, matched against the library (prepare). |
| `patient_regex` | `GBM-001-(\d{3})` | Regex whose group 1 is the patient number; libraries that do not match use the site (prepare). |
| `patient_prefix` | `GBM001-P` | Prefix put before the matched patient number (prepare). |
| `tcr_column` | `tcr_has_tcr_match` | obs column holding TCR presence, giving `tcr_status`; empty skips it (prepare). |
| `two_phase` | `true` | Two-phase streaming conversion. Keeps memory bounded. |
| `cell_chunk_size` | `10000` | Cells per chunk in phase 1. Does not affect output layout. |
| `var_chunk_size` | `1` | Genes per chunk. One gene per chunk means a single-gene fetch reads one chunk. |
| `shard_size` | `30` | Genes per shard. Keeps the object count down. |
| `dtype` | `float16` | Output dtype. |
| `add_layers` | — | Comma-separated layer names, one task each. `convert` writes every other key already; `layers/` is the only thing left to add. |
| `atomic_axes` | — | Comma-separated obs columns. Strata are built only when set. |
| `coarse` | — | Semicolon-separated coarse table specs, e.g. `cell_type;cell_type,treatment`. |
| `max_atomic_cardinality` | — | Reject the strata build above this many atomic strata. |

Defaults match the shipped SPECTRUM store.

## Testing

```bash
nf-test test --tag stub                                # stubs, no image pull
docker build --platform linux/amd64 -t nf-cell2zarr-prepare:dev containers/prepare  # once, and after any containers/prepare change
nf-test test --tag container --profile test,docker     # real conversion and real PREPARE
```

The PREPARE container test runs the local `nf-cell2zarr-prepare:dev` image, so build it
first. Unit tests for the prepare script and the Cirro preprocess script run with
`uv run --with-requirements containers/prepare/requirements.txt --group dev pytest tests/ -q`.

`nf-test --tag` is include-only — there is no negation — so every test carries an
explicit tag and the two suites are selected by name.
