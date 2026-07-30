# 11-genome-paper-figures

Notebooks to reproduce the figures from the 11-genome pioneer-panel paper.

## Setup

1. Create the conda environment:

   ```
   conda env create -f environment.yml
   conda activate 11-genome-paper-figures
   ```

2. Populate `data/` (see **Data** below) -- this folder is gitignored, you create it locally.
3. Run the notebooks in the order listed under **Notebooks** below.

## Data

This repo does not bundle data. Create a local `data/` folder at the repo root and place the
following files in it before running any notebook:

| File | Used by | Source |
|---|---|---|
| `data/gene_fitness_lmm_results.parquet` | figure3, figure4, sup_figure1, sup_figure2 | `<FILL_IN: public release link>` |
| `data/library_insert_data.parquet` | figure1 | `<FILL_IN: public release link>` |
| `data/selection_insert_data.parquet` | figure2, sup_figure1 | `<FILL_IN: public release link>` |
| `data/gene_annotation.parquet` | figure1 | `<FILL_IN: public release link -- not yet in the internal Figshare-prep pipeline, needs to be added>` |
| `data/strand_direction_stats.parquet` | sup_figure1 | `<FILL_IN: public release link -- not yet in the internal Figshare-prep pipeline, needs to be added>` |
| `data/keio_fitness.tsv` | figure4 | `<FILL_IN: Fitness Browser Keio dataset URL (E. coli BW25113)>` |
| `data/gtdb/bac120_taxonomy_r232.tsv` | figure1, figure3 | `<FILL_IN: GTDB release r232 download URL>` |
| `data/gtdb/sp_clusters_r232.tsv` | figure1 | `<FILL_IN: GTDB release r232 download URL>` |
| `data/gtdb/bac120_r232.tree` | figure1 | `<FILL_IN: GTDB release r232 download URL>` |
| `data/gff/Bacillus_subtilis_PY79_GCF_023521615.1.gff` | figure2 | `<FILL_IN: NCBI RefSeq GFF3 for GCF_023521615.1>` |

Column descriptions for the three main data tables are in `gene_fitness_lmm_results_columns.md`,
`library_insert_columns.md`, and `selection_insert_columns.md`.

`data/gene_annotation.parquet` and `data/strand_direction_stats.parquet` are not yet part of the
public dataset release -- they'll need to be added there before an external reader can populate
`data/` and actually run `figure1.ipynb` / `sup_figure1.ipynb`.

`data/gff/*.gff` files are standard NCBI RefSeq GFF3 annotations, downloaded by GCF accession
(only `Bacillus_subtilis_PY79_GCF_023521615.1` is needed here) -- `utils.load_local_gff_databases`
builds a local `gffutils` database from each one the first time it's called and caches it in
`data/gff/.cache/`.

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
| `sup_figure1.ipynb` | Sup Figure 1a-1c | -- |
| `sup_figure2.ipynb` | Sup Figure 2a-2c | -- |

Each notebook writes its output figures to a local `results/` folder (gitignored --
regenerated, not committed).

## License

MIT -- see `LICENSE`.
