# 11-genome-hgt-publication

Notebooks to reproduce the figures from the 11-genome pioneer-panel paper.

## Setup

1. Create the conda environment and register the ipykernel for use with Jupyter:

   ```
   conda env create -f environment.yml
   conda activate 11-genome-hgt-publication
   ### Register jupyter kernel for running jupyter notebooks
   python -m ipykernel install --user --name 11-genome-hgt-publication --display-name "11-genome-hgt-publication"

   ```

2. Populate `data/` (see **Data** below) -- this folder is gitignored, you create it locally.
3. Run the notebooks in the order listed under **Notebooks** below.

## System Requirements

Most of this repo is plain Python/Jupyter. The exception is `gene_linear_model_code/` (see
below), which fits per-gene linear mixed models via R's `lme4` through
[pymer4](https://github.com/ejolly/pymer4)/`rpy2`, and a Bayesian LMM via PyTorch/Pyro:

- **R is not something you install separately.** `environment.yml`'s `ejolly` + `conda-forge`
  channels resolve `pymer4`, `rpy2`, `r-base`, `r-lme4`, `r-lmerTest`, and `r-emmeans` together as
  conda packages -- `conda env create -f environment.yml` is the only setup step needed.
- **OS support follows pymer4's own constraints:** Windows is not officially supported (use WSL
  instead); macOS on Intel can hit `rpy2` compiler errors during environment creation, fixable via
  `brew install gcc` and setting the `CC`/`CFLAGS` environment variables (see pymer4's
  [installation docs](https://eshinjolly.com/pymer4/pages/installation.html)).
- **Multi-core CPU recommended** for `run_lmer_lmm_lb_salt.py` and
  `run_lmer_lmm_lb_salt_direction.py` -- they fit one `lme4` model per gene (~48K genes),
  parallelized via `multiprocessing`. Tune `N_PROCESSES` at the top of each script for your
  machine's core count.
- **Runtime expectations:** the Bayesian LMM (`run_bayesian_lmm_lb_salt.ipynb`) trains 3
  independent seeds at 1000 SVI iterations each; on CPU this took roughly 30 minutes per seed in
  the original run (a GPU, if available, is used automatically). The two `lme4`-based scripts fit tens of thousands of
  per-gene models and will also take a while, scaling down with more cores.

## Data

This repo does not bundle data. Create a local `data/` folder at the repo root and place the
following files in it before running any notebook:

| File | Used by | Source |
|---|---|---|
| `data/gene_fitness_lmm_results.parquet` | figure3, figure4, sup_figure2, sup_figure3 | `<FILL_IN: public release link>` |
| `data/library_insert_data.parquet` | figure1 | `<FILL_IN: public release link>` |
| `data/selection_insert_data.parquet` | figure2, sup_figure1, sup_figure2 | `<FILL_IN: public release link>` |
| `data/gene_annotation.parquet` | figure1 | `<FILL_IN: public release link -- not yet in the internal Figshare-prep pipeline, needs to be added>` |
| `data/strand_direction_stats.parquet` | sup_figure2 | `<FILL_IN: public release link -- not yet in the internal Figshare-prep pipeline, needs to be added>` |
| `data/clonal_validation_data.csv` | figure2, sup_figure1 | `<FILL_IN: public release link>` |
| `data/fit_organism_Keio.tab` | figure4 | [E coli FitnessBrowser Data](https://fit.genomics.lbl.gov/cgi-bin/createFitData.cgi?orgId=Keio) (see note below) |
| `data/gtdb/bac120_taxonomy_r232.tsv` | figure1, figure3 | [GTDB r232 taxonomy](https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) |
| `data/gtdb/sp_clusters_r232.tsv` | figure1 | [GTDB r232 sp clusters](https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/sp_clusters_r232.tsv) |
| `data/gtdb/bac120_r232.tree` | figure1 | [GTDB r232 phylogenetic tree](https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/bac120_r232.tree) |
| `data/GCF_023521615.1_ASM2352161v1_genomic.gff.gz` | figure2 | [Bacillus subtilis PY79 GFF](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/521/615/GCF_023521615.1_ASM2352161v1/GCF_023521615.1_ASM2352161v1_genomic.gff.gz) |

Column descriptions for the three main data tables are in `gene_fitness_lmm_results_columns.md`,
`library_insert_columns.md`, and `selection_insert_columns.md`.

`data/strand_direction_stats.parquet` are not yet part of the
public dataset release -- they'll need to be added there before an external reader can populate
`data/` and actually run `sup_figure2.ipynb`.

GFF files are standard NCBI RefSeq GFF3 annotations -- place them directly in `data/` under
their unmodified downloaded name (only `GCF_023521615.1_ASM2352161v1_genomic.gff.gz` is needed
here, for `Bacillus_subtilis_PY79_GCF_023521615.1`). `utils.load_local_gff_databases` builds a
local `gffutils` database from each one the first time it's called and caches it in
`data/.cache/`.

The KEIO fitness data (`data/fit_organism_Keio.tsv`) is served behind Cloudflare bot
protection, so `wget`/`curl` get a 403 instead of the file -- it has to be downloaded through
an actual browser:

1. Open [https://fit.genomics.lbl.gov/cgi-bin/createFitData.cgi?orgId=Keio](https://fit.genomics.lbl.gov/cgi-bin/createFitData.cgi?orgId=Keio) in a browser.
2. Save the download as `fit_organism_Keio.tsv` into `data/` (this is the default filename the browser gives it, so usually no rename needed).

## Notebooks

Run `figure1.ipynb` before `figure3.ipynb` and `figure4.ipynb` -- both read an intermediate
file (`results/pioneer_genome_gtdb_phylogeny.tsv`) that `figure1.ipynb` writes. The rest have no
cross-notebook dependencies.

| Notebook | Produces | Depends on |
|---|---|---|
| `figure1.ipynb` | Figure 1a, 1c, 1d, 1e, 1f | -- (writes `results/pioneer_genome_gtdb_phylogeny.tsv`) |
| `figure2.ipynb` | Figure 2b, 2c, 2d, gene-analysis Sankey | -- |
| `figure3.ipynb` | Figure 3a-3f | `figure1.ipynb` |
| `figure4.ipynb` | Figure 4a-4e | `figure1.ipynb` |
| `sup_figure1.ipynb` | Sup Figure 1a | -- |
| `sup_figure2.ipynb` | Sup Figure 2a-2c | -- |
| `sup_figure3.ipynb` | Sup Figure 3a-3c | -- |

Each notebook writes its output figures to a local `results/` folder (gitignored --
regenerated, not committed).

## Gene-level linear models

`gene_linear_model_code/` contains the three models that produced the `bayes_*`/`lmer_*` columns
already baked into `data/gene_fitness_results_with_annotations.parquet`, plus the
`data/directional_lmer_stats.parquet` table and a per-gene insert-coverage/VIF metrics notebook
(also folded into `gene_fitness_results_with_annotations.parquet`'s `n_inserts`/`n_sense_inserts`/
`n_antisense_inserts`/`VIF` columns) -- so you can reproduce (or rerun with different parameters)
the underlying statistics rather than just the figures built from them. No additional downloads
are needed -- both inputs they read (`data/selection_experiment_insert_data.parquet` and
`data/gene_fitness_results_with_annotations.parquet`, the latter reused here purely as a
gene-coordinate lookup) are already listed in the **Data** table above. See **System
Requirements** above before running the two `lme4`-based scripts -- `run_lmer_lmm_lb_salt.py`/
`run_lmer_lmm_lb_salt_direction.py` need R/`lme4` via `pymer4`.

| Script | Computes | Produces |
|---|---|---|
| `run_bayesian_lmm_lb_salt.ipynb` | Bayesian LMM (Pyro SVI, `pioneer_hgt_core.bayesian_lmm`) | `results/gene_bayesian_lmm_stats.parquet` + diagnostic plots (ELBO convergence, seed stability, volcano plot) |
| `run_lmer_lmm_lb_salt.py` | Frequentist coordinate-window LMM (`lme4` via `pymer4`) | `results/gene_lmer_lmm_stats.parquet` |
| `run_lmer_lmm_lb_salt_direction.py` | Direction-aware extension of `run_lmer_lmm_lb_salt.py` (adds a sense/antisense term per gene) | `results/gene_lmer_lmm_direction_stats.parquet` |
| `gene_insert_metrics.ipynb` | Per-gene insert coverage (total/sense/antisense) + Variance Inflation Factor within a 5kb coordinate window | `results/gene_insert_metrics.parquet` |

All four are independent -- there's no run order between them. Run from inside
`gene_linear_model_code/` (paths inside each script/notebook are relative, matching the
`figure_code/` notebooks' convention). `run_bayesian_lmm_lb_salt.ipynb` and
`gene_insert_metrics.ipynb` are notebooks -- open them in Jupyter and run all cells, same as any
`figure_code/` notebook (this also lets you see the Bayesian LMM's diagnostic plots inline as
each cell executes). The other two are plain scripts:

```
cd gene_linear_model_code
python run_lmer_lmm_lb_salt.py
python run_lmer_lmm_lb_salt_direction.py
```

(If you'd rather run a notebook headlessly, e.g. `jupyter nbconvert --to notebook --execute
--inplace run_bayesian_lmm_lb_salt.ipynb`, that works too, but needs `nbconvert` installed
separately -- it isn't part of `environment.yml`, matching the rest of this repo's notebooks.)
`pymer4_lmm.py` is a shared helper module imported by the two `lme4` scripts (gene-level `lme4`
fitting + FDR correction) -- not meant to be run directly.

`run_lmer_lmm_lb_salt_direction.py`'s output schema already matches the public
`directional_lmer_stats.parquet` exactly (documented in `directional_lmer_stats_columns.md`); the
other scripts' outputs use the same `bayes_*`/`lmer_*`/`n_inserts`/`VIF` column names as
`gene_fitness_results_with_annotations.parquet`, but re-running them won't reproduce that file
bit-for-bit -- it's a per-gene merge of all these outputs plus several other annotation joins not
included here.

## License

MIT -- see `LICENSE`.
