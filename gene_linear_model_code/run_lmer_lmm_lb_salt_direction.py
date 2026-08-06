# %%
import os
import subprocess
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
import time
from pymer4_lmm import run_coord_lmm_direction_for_gene
from pymer4_lmm import fdr_correction

# Ensure R_HOME is visible to rpy2 — R lives in the active conda env but may not be on PATH
if 'R_HOME' not in os.environ:
    _conda_r = os.path.join(sys.prefix, 'bin', 'R')
    try:
        os.environ['R_HOME'] = subprocess.check_output([_conda_r, 'RHOME'], text=True).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass  # set R_HOME manually if this fails

# --- data directory (populate yourself -- see README's Data section) ---
DATA_DIR = Path("../data")
RESULTS_DIR = Path("../results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# --- paths ---
INPUT_PATH        = DATA_DIR / "selection_experiment_insert_data.parquet"
GENE_COORDS_PATH  = DATA_DIR / "gene_fitness_results_with_annotations.parquet"
OUTPUT_STATS_PATH = RESULTS_DIR / "gene_lmer_lmm_direction_stats.parquet"

# --- parameters ---
N_PROCESSES          = min(12, os.cpu_count() or 4)
N_GENES_TO_TEST       = None   # set to an integer to test a subset
COORD_WINDOW_BP       = 5000
ENVIRONMENT           = "LB_4_salt"
MAX_EMPTY_INSERT_IDS  = 1000

### Gene coordinate lookup -- derived from the public gene-level fitness table, which carries
### one row per gene with locus_tag/gene_chrom/gene_start/gene_end (same columns the original
### script read from a private S3 gene-annotation parquet)
gene_coords_df = pd.read_parquet(GENE_COORDS_PATH)
locus_tag_coords: dict[str, tuple[str, int, int]] = {
    row.locus_tag: (row.gene_chrom, int(row.gene_start), int(row.gene_end))
    for row in gene_coords_df.itertuples(index=False)
}

# --- fitness data ---
fitness_data     = pd.read_parquet(INPUT_PATH)
fitness_data_filt = fitness_data[fitness_data['environment'] == ENVIRONMENT]

# --- reshape to standard fields ---
prepped_data = fitness_data_filt[[
    "bc_sequence", "fitness", "gene_ids_fully_covered",
    "replicate", "insert_type",
    "reference_name", "reference_start", "reference_end",
]].rename(columns={"bc_sequence": "insert_id"})

# --- gene–insert relationships ---
gene_matrix_long = (
    prepped_data[prepped_data["insert_type"] == "aligned"][["insert_id", "gene_ids_fully_covered"]]
    .explode("gene_ids_fully_covered")
    .drop_duplicates()
)

gene_matrix_long.fillna("Empty", inplace=True)
gene_matrix_long['locus_tags_fully_covered'] = (
    gene_matrix_long['gene_ids_fully_covered'].str.replace('gene-', '', regex=False)
)

gene_matrix_long_with_genes = gene_matrix_long[gene_matrix_long['locus_tags_fully_covered'] != "Empty"]

insert_to_genes = (
    gene_matrix_long_with_genes
    .groupby('insert_id')['locus_tags_fully_covered']
    .apply(frozenset).to_dict()
)

# --- empty insert IDs (null group) ---
empty_insert_ids = list(
    prepped_data[prepped_data["insert_type"] == "empty_insert"]["insert_id"].unique()
)[:MAX_EMPTY_INSERT_IDS]

# --- per-insert coordinate lookup ---
aligned_coords = prepped_data[prepped_data["insert_type"] == "aligned"][
    ["insert_id", "reference_name", "reference_start", "reference_end"]
].drop_duplicates(subset="insert_id")

insert_coords: dict[str, tuple[str, int, int]] = {
    row.insert_id: (row.reference_name, int(row.reference_start), int(row.reference_end))
    for row in aligned_coords.itertuples(index=False)
}

# ---
# Build antisense gene lookup: insert_id -> frozenset of locus_tags (prefix stripped)
# that are both fully covered AND in the antisense orientation on that insert.
# Only inserts with at least one antisense fully-covered gene are stored; all others
# return frozenset() via dict.get() default.
# ---
print("Building antisense gene lookup...")
antisense_cols = fitness_data_filt[
    (fitness_data_filt['insert_type'] == 'aligned') &
    (fitness_data_filt['gene_ids_fully_covered'].notna())
][['bc_sequence', 'gene_ids_fully_covered', 'antisense_gene_ids']].drop_duplicates(subset='bc_sequence')

insert_to_antisense_genes: dict[str, frozenset[str]] = {}
for row in tqdm(antisense_cols.itertuples(index=False), total=len(antisense_cols),
                desc="Antisense lookup"):
    if len(row.antisense_gene_ids) == 0:
        continue
    fully_covered_set = set(row.gene_ids_fully_covered)
    antisense_covered = frozenset(
        gid.replace('gene-', '')
        for gid in row.antisense_gene_ids
        if gid in fully_covered_set
    )
    if antisense_covered:
        insert_to_antisense_genes[row.bc_sequence] = antisense_covered

print(f"  {len(insert_to_antisense_genes):,} inserts have ≥1 antisense fully-covered gene")

# ---
# Module-level globals for shared data — set before Pool() so forked workers
# inherit them via copy-on-write without any pickling overhead.
# ---
_PREPPED_DATA          = None
_LOCUS_TAG_COORDS      = None
_INSERT_COORDS         = None
_INSERT_TO_GENES       = None
_EMPTY_INSERT_IDS      = None
_COORD_WINDOW_BP       = None
_ANTISENSE_GENES       = None


def _run_gene(locus_tag):
    return run_coord_lmm_direction_for_gene(
        locus_tag, _PREPPED_DATA, _LOCUS_TAG_COORDS,
        _INSERT_COORDS, _INSERT_TO_GENES, _EMPTY_INSERT_IDS,
        _COORD_WINDOW_BP, _ANTISENSE_GENES,
    )


genes_to_run = (
    gene_matrix_long_with_genes["locus_tags_fully_covered"].unique()[:N_GENES_TO_TEST]
    if N_GENES_TO_TEST is not None
    else gene_matrix_long_with_genes["locus_tags_fully_covered"].unique()
)

# Assign shared data to globals before creating the pool
_PREPPED_DATA     = prepped_data
_LOCUS_TAG_COORDS = locus_tag_coords
_INSERT_COORDS    = insert_coords
_INSERT_TO_GENES  = insert_to_genes
_EMPTY_INSERT_IDS = empty_insert_ids
_COORD_WINDOW_BP  = COORD_WINDOW_BP
_ANTISENSE_GENES  = insert_to_antisense_genes

if __name__ == "__main__":
    lm_result_list = []
    t_start = time.perf_counter()
    with mp.Pool(processes=N_PROCESSES) as pool:
        for result in tqdm(pool.imap_unordered(_run_gene, genes_to_run, chunksize=10),
                           total=len(genes_to_run)):
            lm_result_list.append(result)
    elapsed = time.perf_counter() - t_start

    print(f"\ndirection model: {len(genes_to_run)} genes in {elapsed:.1f}s  ({elapsed / len(genes_to_run):.2f}s per gene)")

    # ---
    # FDR correction — applied separately to each hypothesis:
    #   sense effect, sense-vs-antisense delta, and total antisense effect
    # ---
    lm_stats_df = pd.DataFrame(lm_result_list)
    lm_stats_df["lmer_q_value"]            = fdr_correction(lm_stats_df["lmer_p_value"].values)
    lm_stats_df["lmer_rev_q_value"]        = fdr_correction(lm_stats_df["lmer_rev_p_value"].values)
    lm_stats_df["lmer_antisense_q_value"]  = fdr_correction(lm_stats_df["lmer_antisense_p"].values)

    lm_stats_df.to_parquet(OUTPUT_STATS_PATH)
    print(f"Saved {len(lm_stats_df):,} rows to {OUTPUT_STATS_PATH}")
