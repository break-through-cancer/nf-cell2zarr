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
nf-test test --tag container --profile test,docker     # real conversion
```

`nf-test --tag` is include-only — there is no negation — so every test carries an
explicit tag and the two suites are selected by name.
