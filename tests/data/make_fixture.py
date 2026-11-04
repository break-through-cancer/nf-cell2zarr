"""Generate the tiny h5ad fixture used by the pipeline tests.

Run: uv run --with anndata --with scipy python tests/data/make_fixture.py
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

N_OBS = 200
N_VARS = 50


def INDEX(names: list[str]) -> pd.Index:
    """Build a row/column index the way real single-cell data stores one.

    Must be object dtype, not pandas StringDtype. anndata writes a nullable
    string index as a *group* (values + validity mask) under the
    "nullable-string-array" encoding, and cell2zarr's `_write_metadata` reads
    `obs/_index` as a flat array — it slices it with `[:]`, which on a group is
    a path lookup and raises TypeError. Real datasets carry an object index and
    take the flat-array path, so the fixture must too.

    Passing a bare list here would infer StringDtype on current pandas.
    """
    return pd.Index(names, dtype=object)


def main() -> None:
    rng = np.random.default_rng(0)

    X = rng.standard_normal((N_OBS, N_VARS)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["A", "B", "C", "D"] * (N_OBS // 4)),
            "treatment": pd.Categorical(["ctrl", "drug"] * (N_OBS // 2)),
            "n_counts": rng.integers(100, 10000, N_OBS).astype(np.float32),
        },
        index=INDEX([f"cell_{i}" for i in range(N_OBS)]),
    )
    var = pd.DataFrame(
        {"gene_name": [f"gene_{i}" for i in range(N_VARS)]},
        index=INDEX([f"gene_{i}" for i in range(N_VARS)]),
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.standard_normal((N_OBS, 2)).astype(np.float32)
    adata.layers["counts"] = rng.integers(0, 100, (N_OBS, N_VARS)).astype(np.float32)
    # obsp exists so the fixture covers `cell2zarr convert`'s sparse obsp path.
    adata.obsp["connectivities"] = sparse.random(
        N_OBS, N_OBS, density=0.02, format="csr", dtype=np.float32, random_state=0
    )
    adata.uns["method"] = "fixture"

    out = Path(__file__).parent / "tiny.h5ad"
    adata.write_h5ad(out)
    print(f"wrote {out} ({out.stat().st_size} bytes)")

    make_merge_fixture()


def make_merge_fixture() -> None:
    """A merge-shaped object for tests/test_prepare.py: integer counts, both
    indexes named like the column they duplicate, case-duplicate obs columns."""
    rng = np.random.default_rng(0)
    n_vars = 60
    X = sparse.csr_matrix(rng.poisson(1, (N_OBS, n_vars)).astype(np.float64))
    obs_index = INDEX([f"cell_{i}" for i in range(N_OBS)])
    obs_index.name = "bpcells_name"
    var_index = INDEX([f"gene_{i}" for i in range(n_vars)])
    var_index.name = "bpcells_name"
    libs = ["BTC-GBM-001-002-GEX_S2_L002", "GBM1.DFCI1.S1", "JHU1-S1"]
    obs = pd.DataFrame(
        {
            "orig.ident": [libs[i % 3] for i in range(N_OBS)],
            "tcr_has_tcr_match": [i % 2 == 0 for i in range(N_OBS)],
            "modality": ["GEX" if i % 2 == 0 else "NA" for i in range(N_OBS)],
            "Modality": ["NA" if i % 2 == 0 else "GEX" for i in range(N_OBS)],
        },
        index=obs_index,
    )
    # Same dtype as the index, or anndata calls the column's values "different".
    var = pd.DataFrame({"bpcells_name": var_index.to_series().astype(object)}, index=var_index)
    out = Path(__file__).parent / "tiny_merge.h5ad"
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(out)
    print(f"wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
