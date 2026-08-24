# Convert h5ad to Zarr (cell2zarr)

Converts a single-cell `h5ad` dataset into a Zarr v3 store for cell-explorer.

**Input:** a Single-Cell Dataset containing exactly one `.h5ad` file. Loom and
Seurat files are not supported.

**Output:** a `.zarr` store in this project, with any layers you selected and,
optionally, precomputed strata tables.

## What gets written

The conversion writes every key present in your `h5ad` — `X`, `obs`, `var`,
`obsm`, `obsp` and `uns` — on its own. Layers are the exception: name the ones
you want and each is added as a separate parallel task.

## Choosing strata axes

Strata are built only if you name the `obs` columns to use as atomic axes.
These are specific to your dataset — if you do not know them, leave the field
blank, run the conversion, and inspect the resulting store's `obs` columns
before running again with the axes filled in.

## Defaults

The conversion defaults match the shipped SPECTRUM store: one gene per chunk so
a single-gene fetch reads exactly one chunk, and 30-gene shards to keep the
object count manageable.

## Where the output goes

The store is written into this Cirro project. Publishing it to the
cell-explorer serving bucket is a separate step and is not performed by this
pipeline.
