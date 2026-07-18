# Spatialsnake Bioconda Packaging

This directory contains the focused Bioconda contribution draft for
`spatialsnake`.

The final PR should add only:

```text
recipes/spatialsnake/meta.yaml
```

The recipe installs Spatialsnake source code and conda-managed dependencies
that build quickly. Heavy packages such as the torch/scvi chain are completed
by the built-in helper after installation. This draft does not package
third-party projects such as `napari-spatialdata`, `cell2location`,
`cellphonedb`, `cellcharter`, `ktplotspy`, `torchgmm`, or `sknw` as new
Bioconda recipes.

## Why Only One Recipe

Adding third-party recipes in the same PR is not automatically a copyright
violation if the upstream licenses allow redistribution and the recipes use
stable sources with hashes. But it increases review scope and maintenance
responsibility. For this project, the cleaner plan is:

- Bioconda recipe: Spatialsnake source plus dependencies already available from
  Bioconda/conda-forge/defaults.
- Built-in command: exact-version pip/R packages that are unavailable from
  conda or conflict with the conda-managed base environment.

## Built-In Package Completion

These remain outside mandatory Bioconda `run` dependencies and are installed by
`spatialsnake install-packages` when needed:

- `spatialdata==0.5.0`, `spatialdata-io==0.3.0`,
  `spatialdata-plot==0.2.11`, and `squidpy==1.6.5`
- `napari-spatialdata==0.5.7`
- `cell2location==0.1.5`
- `scvi-tools==1.4.0`, `torch==2.8.0`,
  `pytorch-lightning==2.5.5`, `pyro-ppl==1.9.1`,
  `pyro-api==0.1.1`, and `opt-einsum==3.3.0`
- `pyarrow==21.0.0`, `geopandas==1.1.1`, and `pyogrio==0.12.1`
- `igraph==0.11.9`, `louvain==0.8.2`, and `leidenalg==0.10.2`
- `dask==2024.11.2`, `dask-expr==1.1.19`, and
  `dask-image==2024.5.3`
- `opencv-python==4.13.0.92`
- `numba==0.62.1` and `llvmlite==0.45.1` are present in the recipe
  dependency solve, then verified/reinstalled as PyPI wheels by the helper so
  the conda-native path matches the original PyPI/test environment after
  `pyarrow==21.0.0` is installed from PyPI
- optional downstream Python packages such as `cellphonedb`, `cellcharter`,
  `ktplotspy`, `torchgmm`, `sknw`, `liana`, `loompy`, and `pyscenic`
- conflict-managed `pybanksy==1.3.4 --no-deps`
- R/GitHub packages such as `Seurat 5.5.0`, `ggsankey`, `schard`,
  `CellChat`, `spacexr`, `NMF 0.26`, and `circlize 0.4.15`

Conda/Bioconda users get a minimal runnable environment with:

```bash
spatialsnake install-packages
```

Full downstream support is explicit:

```bash
spatialsnake install-packages --extended
```

## User Install Path

After the Bioconda PR is reviewed, merged, and propagated:

```bash
conda create -n spatialsnake_env -c bioconda -c conda-forge spatialsnake
conda activate spatialsnake_env
spatialsnake --version
spatialsnake install-packages
```

For full downstream support:

```bash
spatialsnake install-packages --extended
```

See `CONTRIBUTING_RECIPES.md` for the detailed official workflow.
