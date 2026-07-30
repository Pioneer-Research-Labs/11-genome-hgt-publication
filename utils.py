"""
Plotting utilities for reproducing the 11-genome paper figures.

Includes a small set of gene-region plotting helpers vendored from the
internal (private) `seq2anno.visualization` module -- see the "Vendored from
seq2anno" section below -- plus a local GFF-loading helper that replaces
`seq2anno.load_genome_resources_from_s3` (which streamed private S3 data).

Call plt.show() or fig.savefig() after adding any additional layers to the
axes returned by plot_gene_region_plain().
"""

import logging
from pathlib import Path

import gffutils
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

logger = logging.getLogger(__name__)


# =============================================================================
# Vendored from seq2anno.visualization (Pioneer-Research-Labs/Seq2Anno, private)
# Copied rather than imported so this repo has no dependency on private infra.
# =============================================================================

def _draw_directional_marker(ax, start, end, y, strand, color, linewidth=3, marker_height=0.15):
    """
    Draw a directional marker (|-> or <-|) for a genomic feature or insert.

    Args:
        ax: matplotlib Axes object
        start: genomic start position
        end: genomic end position
        y: y-axis position
        strand: '+' for forward, '-' for reverse, or bool (True=reverse)
        color: line color
        linewidth: line width
        marker_height: height of directional markers
    """
    length = end - start
    marker_width = length * 0.02

    # Normalize strand to boolean is_reverse
    if isinstance(strand, bool):
        is_reverse = strand
    else:
        is_reverse = (strand == '-')

    if is_reverse:
        # Reverse/antisense: <-|
        ax.plot([start, start + marker_width],
                [y, y - marker_height / 2],
                color=color, linewidth=linewidth, zorder=6)
        ax.plot([start, start + marker_width],
                [y, y + marker_height / 2],
                color=color, linewidth=linewidth, zorder=6)
        ax.plot([end, end],
                [y - marker_height / 2, y + marker_height / 2],
                color=color, linewidth=linewidth, zorder=6)
    else:
        # Forward/sense: |->
        ax.plot([start, start],
                [y - marker_height / 2, y + marker_height / 2],
                color=color, linewidth=linewidth, zorder=6)
        ax.plot([end - marker_width, end],
                [y - marker_height / 2, y],
                color=color, linewidth=linewidth, zorder=6)
        ax.plot([end - marker_width, end],
                [y + marker_height / 2, y],
                color=color, linewidth=linewidth, zorder=6)


def _resolve_fitness(df, fitness_mode, barcode_col):
    """
    Resolve which fitness column to use, computing it if necessary.
    Prefer corrected mean fitness > mean fitness > fitness.

    Returns:
        (use_fitness, fitness_column_name, df) - df may be a copy with computed column
    """
    if fitness_mode == 'none':
        return False, None, df

    if 'mean_fitness_corrected' in df.columns:
        print("Using mean_fitness_corrected")
        return True, 'mean_fitness_corrected', df

    if 'mean_fitness' in df.columns:
        print("Using mean_fitness")
        return True, 'mean_fitness', df

    if 'fitness' in df.columns:
        print("Using fitness")
        df = df.copy()
        df['calculated_fitness'] = df.groupby(barcode_col)['fitness'].transform(fitness_mode)
        return True, 'calculated_fitness', df

    logger.warning("No fitness columns found. Setting fitness_mode to 'none'.")
    return False, None, df


def _compute_y_layout(barcodes_df, use_fitness, fitness_column_name):
    """
    Compute y-axis positions for barcodes and genes.

    Returns:
        (min_fitness, max_fitness, y_gene_center, fixed_barcode_y)
        fixed_barcode_y is None when using fitness.
    """
    if use_fitness and fitness_column_name:
        max_fitness = barcodes_df[fitness_column_name].max() if not barcodes_df.empty else 1
        min_fitness = barcodes_df[fitness_column_name].min() if not barcodes_df.empty else 0
        return min_fitness, max_fitness, max_fitness + 2.0, None
    else:
        fixed_y = 0.0
        return fixed_y, fixed_y, fixed_y + 2.0, fixed_y


def _get_genes_in_region(species_dbs_map, species_name, chromosome, plot_start, plot_end):
    """Query gffutils database for genes in a region."""
    if species_name not in species_dbs_map:
        logger.error("No database found for species '%s'.", species_name)
        return []

    db = species_dbs_map[species_name]
    genes = []
    for feature in db.region(seqid=chromosome, start=plot_start, end=plot_end, featuretype='gene'):
        genes.append({
            'start': feature.start,
            'end': feature.end,
            'name': feature.attributes.get('Name', ['N/A'])[0],
            'strand': feature.strand
        })
    return genes


def _fitness_label(fitness_mode):
    """Return display label for a fitness mode."""
    return {'mean': 'Mean Fitness', 'min': 'Minimum Fitness', 'max': 'Maximum Fitness'}.get(
        fitness_mode, 'Fitness')


# =============================================================================
# Local replacement for seq2anno.s3_genomes.load_genome_resources_from_s3
# =============================================================================

def load_local_gff_databases(gff_dir, species_ids, cache_dir=None):
    """
    Build gffutils FeatureDB objects from local GFF3 files.

    Replaces the private seq2anno.load_genome_resources_from_s3() call, which
    streamed pre-built databases from an internal S3 bucket. Here, each
    species' GFF3 annotation is expected at `gff_dir / f"{species_id}.gff"`
    (species_id matching the pioneer_genome naming, e.g.
    "Bacillus_subtilis_PY79_GCF_023521615.1") -- download the matching NCBI
    RefSeq GFF3 file for each species listed in the README's Data section.

    Args:
        gff_dir: Directory containing one {species_id}.gff file per species.
        species_ids: List of species_id strings to load.
        cache_dir: Where to cache built .gffdb files (default: gff_dir/.cache).

    Returns:
        Dict mapping species_id -> gffutils.FeatureDB.
    """
    gff_dir = Path(gff_dir)
    cache_dir = Path(cache_dir) if cache_dir else gff_dir / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    species_dbs = {}
    for species_id in species_ids:
        db_path = cache_dir / f"{species_id}.gffdb"
        if db_path.exists():
            species_dbs[species_id] = gffutils.FeatureDB(str(db_path))
            continue

        gff_path = gff_dir / f"{species_id}.gff"
        if not gff_path.exists():
            raise FileNotFoundError(
                f"Missing GFF3 annotation for '{species_id}' at {gff_path}. "
                "See the README's Data section for the NCBI RefSeq download link."
            )
        species_dbs[species_id] = gffutils.create_db(
            str(gff_path),
            str(db_path),
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    return species_dbs


# =============================================================================
# Gene-region plotting (plain / non-directional variant)
# =============================================================================

# --- colors for plain (non-directional) variant ---
_COLOR_COVERING    = "forestgreen"    # insert fully covers gene of interest
_COLOR_NONCOVERING = "#A0A0A0"  # insert does not fully cover
_COLOR_TARGET_GENE = "forestgreen"    # the focal gene block arrow
_COLOR_OTHER_GENE  = "#A0A0A0"  # all other gene block arrows

# vertical offsets for staggered gene labels (cycles through these in genomic order)
_TEXT_OFFSETS = [-1, 0.5]


def _draw_block_arrow(ax, start, end, y, strand, color, body_height=0.25):
    """Draw a filled block arrow for one gene. Points right for +, left for -."""
    length = end - start
    if length <= 0:
        return
    head_len = max(50, length * 0.18)
    kwargs = dict(
        width=body_height,
        head_width=body_height * 1.8,
        head_length=head_len,
        length_includes_head=True,
        color=color,
        zorder=3,
    )
    if strand == "-":
        ax.add_patch(FancyArrow(end, y, -length, 0, **kwargs))
    else:
        ax.add_patch(FancyArrow(start, y, length, 0, **kwargs))


def _plot_genes_plain(ax, genes, y_gene_center, target_gene_name, covered_gene_names=None):
    """Block-arrow gene track: black for target gene, gray for all others. No strand coloring."""
    sorted_genes = sorted(genes, key=lambda g: (g["start"] + g["end"]) / 2)
    for i, gene in enumerate(sorted_genes):
        is_target = gene["name"] == target_gene_name
        color  = _COLOR_TARGET_GENE if is_target else _COLOR_OTHER_GENE
        weight = "bold" if is_target else "normal"

        _draw_block_arrow(ax, gene["start"], gene["end"], y_gene_center, gene["strand"], color)

        show_label = is_target or (covered_gene_names is None) or (gene["name"] in covered_gene_names)
        if show_label:
            text_y = y_gene_center + _TEXT_OFFSETS[i % len(_TEXT_OFFSETS)]
            ax.text(
                (gene["start"] + gene["end"]) / 2, text_y, gene["name"],
                ha="center", va="bottom", fontsize=10, color=color, weight=weight,
            )

    if not genes:
        return y_gene_center
    return y_gene_center + max(_TEXT_OFFSETS)


def plot_gene_region_plain(
    gene_name,
    species_name,
    annotated_df,
    species_dbs_map,
    window_kb=1,
    fitness_mode="mean",
    barcode_col="bc_sequence",
    figsize=(15, 8),
):
    """
    Simplified variant of plot_gene_region with no directional coloring.

    Inserts are black when they fully span the gene of interest, gray otherwise.
    The gene track uses filled block arrows: black for the focal gene, gray for
    all others. No strand orientation is encoded by color anywhere.

    Returns (fig, ax) without displaying.
    """
    if species_name not in species_dbs_map:
        logger.error("No database found for species '%s'.", species_name)
        return None, None

    species_df = annotated_df[annotated_df["species_full"] == species_name].copy()
    if species_df.empty:
        logger.error("No data found for species '%s'.", species_name)
        return None, None

    try:
        gene_feature = species_dbs_map[species_name][gene_name]
    except Exception:
        logger.error("Gene '%s' not found in database for species '%s'.", gene_name, species_name)
        return None, None

    gene_start = gene_feature.start
    gene_end   = gene_feature.end
    gene_chrom = gene_feature.seqid
    # _get_genes_in_region returns the GFF Name attribute, not the feature ID.
    # Resolve the display name here so the target-gene check in _plot_genes_plain matches.
    target_display_name = next(iter(gene_feature.attributes.get("Name", [])), gene_name)

    window_bp  = window_kb * 1000
    plot_start = max(0, gene_start - window_bp)
    plot_end   = gene_end + window_bp

    barcodes_in_window = species_df[
        (species_df["reference_name"] == gene_chrom) &
        (species_df["reference_start"] <= plot_end) &
        (species_df["reference_end"] >= plot_start)
    ].copy()

    if barcodes_in_window.empty:
        logger.warning("No inserts found in window around gene '%s'.", gene_name)
        return None, None

    use_fitness, fitness_col, barcodes_in_window = _resolve_fitness(
        barcodes_in_window, fitness_mode, barcode_col)

    unique_barcodes = barcodes_in_window.groupby(barcode_col).first().reset_index()

    covering_barcodes = set(
        unique_barcodes[
            (unique_barcodes["reference_start"] <= gene_start) &
            (unique_barcodes["reference_end"] >= gene_end)
        ][barcode_col]
    )

    genes = _get_genes_in_region(
        species_dbs_map, species_name, gene_chrom, plot_start, plot_end)

    fig, ax = plt.subplots(figsize=figsize)
    min_fit, max_fit, y_gene_center, fixed_y = _compute_y_layout(
        unique_barcodes, use_fitness, fitness_col)

    covered_gene_names = {
        gene["name"]
        for gene in genes
        if plot_start <= (gene["start"] + gene["end"]) / 2 <= plot_end
    }

    max_text_y = _plot_genes_plain(
        ax, genes, y_gene_center, target_display_name, covered_gene_names)

    for _, barcode in unique_barcodes.iterrows():
        barcode_y   = barcode[fitness_col] if (use_fitness and fitness_col) else fixed_y
        is_covering = barcode[barcode_col] in covering_barcodes
        is_reverse  = barcode.get("is_reverse", None)
        line_color  = _COLOR_COVERING if is_covering else _COLOR_NONCOVERING

        ax.plot(
            [barcode["reference_start"], barcode["reference_end"]],
            [barcode_y, barcode_y],
            color=line_color, linewidth=2, solid_capstyle="butt", zorder=5,
        )
        if is_reverse is not None:
            _draw_directional_marker(
                ax, barcode["reference_start"], barcode["reference_end"], barcode_y,
                is_reverse, line_color, linewidth=2,
            )

    ax.set_xlim(plot_start, plot_end)
    ax.set_ylim(min_fit - 0.3, max_text_y + 0.8)
    ax.set_xlabel(f"Genomic Position on {gene_chrom} (bp)", fontsize=13)
    ax.set_ylabel(_fitness_label(fitness_mode) if use_fitness else "Barcode Position", fontsize=13)
    ax.tick_params(labelsize=11)
    ax.set_title(
        f'Inserts Near Gene "{gene_name}" in {species_name}\n'
        f'({len(unique_barcodes)} inserts in window, {len(covering_barcodes)} fully covering gene)'
    )
    plt.tight_layout()
    return fig, ax
