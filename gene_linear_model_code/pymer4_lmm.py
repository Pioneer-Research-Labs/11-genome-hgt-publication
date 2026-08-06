import numpy as np
import pandas as pd
import polars as pl
from pymer4.models import lmer
from polars import col as pcol
from scipy.stats import shapiro, norm
import pymer4.tidystats as ts
from statsmodels.stats.multitest import multipletests

def run_coord_lmm_for_gene(locus_tag, 
                           barcode_fitness_table, 
                           locus_tag_coords, 
                           insert_coords, 
                           insert_to_genes, 
                           empty_insert_ids, 
                           COORD_WINDOW_BP):
    '''
    Run a linear mixed model for a given gene.
    
    Args:
        locus_tag: The gene ID to run the model for. Note cannot contain '-' in locus_tag.
        barcode_fitness_table: A pandas DataFrame with columns 'insert_id', 'fitness', and 'replicate'.
        gene_coords: A dictionary with gene IDs as keys and tuples of (chromosome, start, end) as values.
        insert_coords: A dictionary with insert IDs as keys and tuples of (chromosome, start, end) as values.
        insert_to_genes: A dictionary with insert IDs as keys and sets of gene IDs as values.
        empty_insert_ids: A list of the IDs of the empty inserts.
        COORD_WINDOW_BP: The window size in base pairs around the gene to include in the model.

    Returns:
        A dictionary containing model results.
    '''
    # Imported inside the function so rpy2 is not initialised in the parent process before fork

    chrom, gene_start, gene_stop = locus_tag_coords[locus_tag]
    window_start = gene_start - COORD_WINDOW_BP
    window_end   = gene_stop  + COORD_WINDOW_BP

    # --- Select aligned inserts with at least one endpoint inside the padded window ---
    selected_inserts = [
        ins_id for ins_id, (c, s, e) in insert_coords.items()
        if c == chrom and (window_start <= s <= window_end or window_start <= e <= window_end)
    ]
    
    # --- Build binary indicator matrix for genes fully covered by selected inserts ---
    genes_in_window = sorted(frozenset().union(*(insert_to_genes.get(i, frozenset()) for i in selected_inserts)))
    gene_idx = {g: j for j, g in enumerate(genes_in_window)}

    indicator = np.zeros((len(selected_inserts), len(genes_in_window)), dtype=np.int8)
    genes_set = set(genes_in_window)
    for i, ins_id in enumerate(selected_inserts):
        for gene in (insert_to_genes.get(ins_id, frozenset()) & genes_set):
            indicator[i, gene_idx[gene]] = 1

    gene_insert_df = pd.DataFrame(indicator, index=selected_inserts, columns=genes_in_window)

    # Add empty inserts as all-zero rows (the null group)
    empty_rows = pd.DataFrame(0, index=empty_insert_ids, columns=genes_in_window, dtype=np.int8)
    neighbor_data_wide = pd.concat([gene_insert_df, empty_rows])

    ### Join back to relevant fitness info
    filt_fitness = pd.merge(
        barcode_fitness_table[["insert_id", "fitness", "replicate"]],
        neighbor_data_wide,
        left_on="insert_id",
        right_index=True,
    )
        
    # Cast to str so R treats replicate as categorical
    filt_fitness['replicate'] = filt_fitness['replicate'].astype(str)

    all_predictors = " + ".join([col for col in filt_fitness.columns
                                  if col not in ['insert_id', 'fitness', 'replicate']])
    formula = f'fitness ~ {all_predictors} + replicate + (1|insert_id)'

    try:
        model = lmer(formula, data=pl.from_pandas(filt_fitness))
        model.fit(verbose=False)
        shapiro_mask = filt_fitness[locus_tag] == 1
        shapiro_stat, shapiro_p = shapiro(ts.resid(model)[shapiro_mask])
        row = model.result_fit.filter(pcol("term") == locus_tag)
        insert_var = model.ranef_var.filter(pcol("group") == "insert_id")["estimate"].item()

        return {'locus_tag':           locus_tag,
                'lmer_estimate':       row['estimate'].item(),
                'lmer_standard_error': row['std_error'].item(),
                'lmer_shapiro_stat':   shapiro_stat,
                'lmer_shapiro_p':      shapiro_p,
                'lmer_insert_var':     insert_var,
                'lmer_z_value':        row['t_stat'].item(),
                'lmer_p_value':        row['p_value'].item()}

    except Exception:
        return {'locus_tag':           locus_tag,
                'lmer_estimate':       np.nan,
                'lmer_standard_error': np.nan,
                'lmer_shapiro_stat':   np.nan,
                'lmer_shapiro_p':      np.nan,
                'lmer_insert_var':     np.nan,
                'lmer_z_value':        np.nan,
                'lmer_p_value':        np.nan}

def run_coord_lmm_direction_for_gene(locus_tag,
                                     barcode_fitness_table,
                                     locus_tag_coords,
                                     insert_coords,
                                     insert_to_genes,
                                     empty_insert_ids,
                                     COORD_WINDOW_BP,
                                     insert_to_antisense_genes):
    '''
    Extended version of run_coord_lmm_for_gene that adds a direction term for the
    target gene only.

    The model is:
        fitness ~ locus_tag + locus_tag_rev + [neighbor genes] + replicate + (1|insert_id)

    where locus_tag_rev = 1 when the insert fully covers the target gene AND the gene is
    in the antisense orientation on that insert.

    This yields:
      - lmer_estimate:           effect when covered sense (locus_tag coeff)
      - lmer_rev_estimate:       antisense - sense delta  (locus_tag_rev coeff)
      - lmer_antisense_estimate: total antisense effect = sense + delta, with SE from the
                                 variance-covariance matrix (contrast via delta method)

    Args:
        insert_to_antisense_genes: dict mapping insert_id -> frozenset of locus_tags
            (prefix stripped) that are antisense AND fully covered on that insert.
            Inserts absent from the dict are treated as having no antisense genes.
    '''
    import rpy2.robjects as robjects

    col_rev = f'{locus_tag}_rev'

    chrom, gene_start, gene_stop = locus_tag_coords[locus_tag]
    window_start = gene_start - COORD_WINDOW_BP
    window_end   = gene_stop  + COORD_WINDOW_BP

    selected_inserts = [
        ins_id for ins_id, (c, s, e) in insert_coords.items()
        if c == chrom and (window_start <= s <= window_end or window_start <= e <= window_end)
    ]

    genes_in_window = sorted(frozenset().union(*(insert_to_genes.get(i, frozenset()) for i in selected_inserts)))
    gene_idx = {g: j for j, g in enumerate(genes_in_window)}

    indicator = np.zeros((len(selected_inserts), len(genes_in_window)), dtype=np.int8)
    genes_set = set(genes_in_window)
    for i, ins_id in enumerate(selected_inserts):
        for gene in (insert_to_genes.get(ins_id, frozenset()) & genes_set):
            indicator[i, gene_idx[gene]] = 1

    gene_insert_df = pd.DataFrame(indicator, index=selected_inserts, columns=genes_in_window)

    # --- Direction indicator: 1 if insert covers target gene AND target is antisense ---
    antisense_set = insert_to_antisense_genes  # local alias for brevity
    gene_insert_df[col_rev] = np.array([
        1 if (locus_tag in insert_to_genes.get(ins_id, frozenset()) and
              locus_tag in antisense_set.get(ins_id, frozenset()))
        else 0
        for ins_id in selected_inserts
    ], dtype=np.int8)

    n_sense     = int(gene_insert_df[locus_tag].sum()) - int(gene_insert_df[col_rev].sum()) \
                  if locus_tag in gene_insert_df.columns else 0
    n_antisense = int(gene_insert_df[col_rev].sum())

    empty_rows = pd.DataFrame(0, index=empty_insert_ids,
                              columns=gene_insert_df.columns, dtype=np.int8)
    neighbor_data_wide = pd.concat([gene_insert_df, empty_rows])

    filt_fitness = pd.merge(
        barcode_fitness_table[["insert_id", "fitness", "replicate"]],
        neighbor_data_wide,
        left_on="insert_id",
        right_index=True,
    )
    filt_fitness['replicate'] = filt_fitness['replicate'].astype(str)

    # Put locus_tag then col_rev first; exclude both from the general predictor sweep
    other_genes = [c for c in filt_fitness.columns
                   if c not in {'insert_id', 'fitness', 'replicate', locus_tag, col_rev}]
    predictors = ' + '.join([locus_tag, col_rev] + other_genes)
    formula = f'fitness ~ {predictors} + replicate + (1|insert_id)'

    _nan_row = {
        'locus_tag':               locus_tag,
        'lmer_estimate':           np.nan,
        'lmer_standard_error':     np.nan,
        'lmer_z_value':            np.nan,
        'lmer_p_value':            np.nan,
        'lmer_rev_estimate':       np.nan,
        'lmer_rev_standard_error': np.nan,
        'lmer_rev_z_value':        np.nan,
        'lmer_rev_p_value':        np.nan,
        'lmer_antisense_estimate': np.nan,
        'lmer_antisense_se':       np.nan,
        'lmer_antisense_z':        np.nan,
        'lmer_antisense_p':        np.nan,
        'lmer_n_sense_inserts':    n_sense,
        'lmer_n_antisense_inserts':n_antisense,
        'lmer_shapiro_stat':       np.nan,
        'lmer_shapiro_p':          np.nan,
        'lmer_insert_var':         np.nan,
    }

    # Pre-flight: either degenerate case → rank-deficient, model will fail.
    #   n_antisense == 0: col_rev is all-zero (no antisense coverage in window)
    #   n_sense == 0:     locus_tag == col_rev everywhere (all coverage is antisense)
    if n_antisense == 0 or n_sense == 0:
        return _nan_row

    try:
        model = lmer(formula, data=pl.from_pandas(filt_fitness))
        model.fit(verbose=False)

        shapiro_mask = filt_fitness[locus_tag] == 1
        shapiro_stat, shapiro_p = shapiro(ts.resid(model)[shapiro_mask])

        row_fwd = model.result_fit.filter(pcol('term') == locus_tag)
        row_rev = model.result_fit.filter(pcol('term') == col_rev)
        insert_var = model.ranef_var.filter(pcol('group') == 'insert_id')['estimate'].item()

        beta_fwd = row_fwd['estimate'].item()
        beta_rev = row_rev['estimate'].item()

        # --- VCV-based contrast for total antisense effect (β_fwd + β_rev) ---
        vcov_r   = robjects.r['as.matrix'](robjects.r['vcov'](model.r_model))
        vcov_np  = np.array(vcov_r)
        fe_names = list(robjects.r['rownames'](vcov_r))
        tag_idx  = fe_names.index(locus_tag)
        rev_idx  = fe_names.index(col_rev)

        var_fwd     = float(vcov_np[tag_idx, tag_idx])
        var_rev     = float(vcov_np[rev_idx, rev_idx])
        cov_fwd_rev = float(vcov_np[tag_idx, rev_idx])

        antisense_estimate = beta_fwd + beta_rev
        antisense_var      = var_fwd + var_rev + 2.0 * cov_fwd_rev
        antisense_se       = float(np.sqrt(max(antisense_var, 0.0)))
        antisense_z        = antisense_estimate / antisense_se if antisense_se > 0 else np.nan
        antisense_p        = float(2.0 * norm.cdf(-abs(antisense_z))) if not np.isnan(antisense_z) else np.nan

        return {
            'locus_tag':               locus_tag,
            'lmer_estimate':           float(row_fwd['estimate'].item()),
            'lmer_standard_error':     float(row_fwd['std_error'].item()),
            'lmer_z_value':            float(row_fwd['t_stat'].item()),
            'lmer_p_value':            float(row_fwd['p_value'].item()),
            'lmer_rev_estimate':       float(row_rev['estimate'].item()),
            'lmer_rev_standard_error': float(row_rev['std_error'].item()),
            'lmer_rev_z_value':        float(row_rev['t_stat'].item()),
            'lmer_rev_p_value':        float(row_rev['p_value'].item()),
            'lmer_antisense_estimate': antisense_estimate,
            'lmer_antisense_se':       antisense_se,
            'lmer_antisense_z':        antisense_z,
            'lmer_antisense_p':        antisense_p,
            'lmer_n_sense_inserts':    n_sense,
            'lmer_n_antisense_inserts':n_antisense,
            'lmer_shapiro_stat':       float(shapiro_stat),
            'lmer_shapiro_p':          float(shapiro_p),
            'lmer_insert_var':         float(insert_var),
        }

    except Exception as e:
        import sys, traceback
        print(f'[direction_lmm:{locus_tag}] {type(e).__name__}: {e}', file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return _nan_row


def fdr_correction(p_values, alpha=0.05):
    """Applies Benjamini-Hochberg FDR correction."""
    mask = ~np.isnan(p_values)
    p_adj = np.full(p_values.shape, np.nan)
    if mask.sum() > 0:
        p_adj[mask] = multipletests(p_values[mask], alpha=alpha, method='fdr_bh')[1]
    return p_adj