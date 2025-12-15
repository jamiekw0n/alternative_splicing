import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def _distinct_colors(cmap_name: str, n: int):
    """
    Return n visually distinct colors from the given colormap.
    - Uses discrete palette entries if available (e.g., tab20).
    - Otherwise samples evenly from the continuous colormap.
    """
    cmap = plt.get_cmap(cmap_name)
    # If it's a ListedColormap, try to use its discrete list
    base = getattr(cmap, "colors", None)
    if isinstance(base, (list, tuple)) and len(base) > 0:
        if n <= len(base):
            # pick n roughly evenly-spaced indices across the palette
            idxs = np.linspace(0, len(base) - 1, n)
            idxs = np.round(idxs).astype(int)
            # ensure we have exactly n (fallback if rounding collided)
            if len(np.unique(idxs)) == n:
                return [base[i] for i in idxs]
        # if n > len(base) (or rounding collided), fall back to continuous sampling
    # Sample the continuous map evenly in [0,1)
    return [cmap(t) for t in np.linspace(0, 1, n, endpoint=False)]


def plot_count_isoform_usage(
    adata,
    gene,
    cell_types=None,   # <- optional: auto-picks CTs expressing the gene
    ct_col='gpt-3.5-turbo_simplified_ai_cell_type',
    transcript_id_col='transcript_id',
    min_frac=0.1,
    figsize=(36, 18),
    cmap='tab20',
    bar_width=0.8,
    title_fs=16,
    label_fs=10,
    tick_fs=12,
    legend_fs=12,
    *,
    layer=None,                    # e.g. "counts" if your raw counts are in a layer
    annotate_cells=True,
    annotate_counts=True,
    cells_label='expressing',      # 'expressing' (default) or 'total'
    min_frac_label=0.03,           # don't label segments smaller than this fraction
    counts_label_fs=10,
    cells_label_fs=10,
    return_table=False             # return a DataFrame of what was plotted
):
    """
    One stacked vertical bar per cell type showing isoform usage of `gene`.
    Adds:
      - # cells above each bar (default = cells with >0 gene counts in that CT)
      - raw counts on each stacked segment (hidden if < min_frac_label of bar)
    """

    # --- helpers ---
    def get_mat(A):
        mat = A.layers[layer] if layer is not None else A.X
        # Convert sparse to dense if needed for indexing
        if hasattr(mat, 'toarray'):
            return mat
        return mat

    def sum_over_axis0(A):
        mat = get_mat(A)
        if hasattr(mat, 'sum'):
            result = mat.sum(axis=0)
            if hasattr(result, 'A1'):  # sparse matrix
                return result.A1
            return np.asarray(result).ravel()
        return np.asarray(mat.sum(axis=0)).ravel()

    def sum_over_axis1(A):
        mat = get_mat(A)
        if hasattr(mat, 'sum'):
            result = mat.sum(axis=1)
            if hasattr(result, 'A1'):  # sparse matrix
                return result.A1
            return np.asarray(result).ravel()
        return np.asarray(mat.sum(axis=1)).ravel()

    # sanity checks
    if transcript_id_col not in adata.var.columns:
        raise ValueError(f"adata.var must contain '{transcript_id_col}'")
    if 'gene_name' not in adata.var.columns:
        raise ValueError("adata.var must contain 'gene_name'")
    if ct_col not in adata.obs.columns:
        raise ValueError(f"adata.obs must contain '{ct_col}'")

    # find transcripts for this gene
    gene_mask = (adata.var['gene_name'] == gene)
    if gene_mask.sum() == 0:
        raise ValueError(f"Gene {gene!r} not found in adata.var['gene_name']")
    
    # Subset adata to gene transcripts for easier indexing
    adata_gene = adata[:, gene_mask].copy()
    
    # Get transcript IDs in the same order as adata_gene.var (to match counts order)
    iso_all = adata_gene.var[transcript_id_col].values

    # auto-pick cell types if not provided: those with >0 total counts for this gene
    if cell_types is None:
        total_per_cell = sum_over_axis1(adata_gene)
        express_mask = total_per_cell > 0
        cell_types = adata.obs.loc[express_mask, ct_col].unique().tolist()

    # only those cell-types present in the data
    present = [ct for ct in cell_types if (adata.obs[ct_col] == ct).any()]
    if not present:
        raise ValueError("None of the requested cell-types are present")

    # pick the single most abundant transcript overall (for x-axis label)
    global_counts = sum_over_axis0(adata_gene)
    most_idx = int(global_counts.argmax())
    iso_most = iso_all[most_idx]

    # compute total gene counts & fraction of iso_most per CT (to sort bars)
    total_by_ct = {}
    frac_most_by_ct = {}
    for ct in present:
        ct_mask = adata.obs[ct_col] == ct
        ad_ct_gene = adata_gene[ct_mask]
        counts = sum_over_axis0(ad_ct_gene)
        total = float(counts.sum())
        total_by_ct[ct] = total
        frac_most_by_ct[ct] = 0.0 if total == 0 else counts[most_idx] / total

    # drop CTs with zero expression
    present = [ct for ct in present if total_by_ct[ct] > 0]
    if not present:
        raise ValueError(f"No cells express {gene!r} in any of your chosen cell-types")

    # sort remaining CTs by fraction of iso_most (desc)
    present.sort(key=lambda ct: frac_most_by_ct[ct], reverse=True)

    # determine which isoforms get their own color (>= min_frac in at least one CT)
    shown_isos = set()
    for ct in present:
        ct_mask = adata.obs[ct_col] == ct
        ad_ct_gene = adata_gene[ct_mask]
        counts = sum_over_axis0(ad_ct_gene)
        if counts.sum() == 0:
            continue
        fracs = counts / counts.sum()
        order = np.argsort(fracs)[::-1]
        iso_ord = iso_all[order]
        keep = fracs[order] >= min_frac
        shown_isos.update(iso_ord[keep])
    shown_isos = sorted(shown_isos)

    # ---- DISTINCT COLORS (new) ----
    n_colors = max(1, len(shown_isos))
    base_colors = _distinct_colors(cmap, n_colors)
    iso2col = {iso: base_colors[i] for i, iso in enumerate(shown_isos)}
    other_col = (0.8, 0.8, 0.8, 1.0)
    # --------------------------------

    # plot
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(present))
    legend_patches = {}
    y_max = 1.0

    # optional table rows
    table_rows = []

    for i, ct in enumerate(present):
        ct_mask = adata.obs[ct_col] == ct
        ad_ct_gene = adata_gene[ct_mask]

        # counts per transcript for this gene & CT
        counts = sum_over_axis0(ad_ct_gene)
        total_ct = counts.sum()
        if total_ct == 0:
            continue  # Skip if no expression
        fracs = counts / total_ct
        order = np.argsort(fracs)[::-1]
        iso_ord, fracs_ord, counts_ord = iso_all[order], fracs[order], counts[order]

        # keep majors, lump tail
        keep = fracs_ord >= min_frac
        iso_k = list(iso_ord[keep])
        frac_k = list(fracs_ord[keep])
        count_k = list(counts_ord[keep])

        tail_frac = float(fracs_ord[~keep].sum())
        tail_count = float(counts_ord[~keep].sum())
        if tail_frac > 0:
            iso_k.append('Other')
            frac_k.append(tail_frac)
            count_k.append(tail_count)
        
        # Normalize fractions to ensure they sum to exactly 1.0 (fix floating point issues)
        total_frac = sum(frac_k)
        if total_frac > 0:
            frac_k = [f / total_frac for f in frac_k]

        # cells counts
        gene_counts_by_cell = sum_over_axis1(ad_ct_gene)
        n_cells_expressing = int(np.count_nonzero(gene_counts_by_cell > 0))
        n_cells_total = int(ad_ct_gene.n_obs)

        # stacked bar
        bottom = 0.0
        for iso, frac, cnt in zip(iso_k, frac_k, count_k):
            if iso == 'Other':
                col = other_col
            elif iso in iso2col:
                col = iso2col[iso]
            else:
                # Fallback: isoform not in shown_isos but >= min_frac for this CT
                # Use a default color (gray)
                col = other_col
            ax.bar(i, frac, bottom=bottom, width=bar_width, color=col, edgecolor='white')

            # segment label: raw counts
            if annotate_counts and frac >= min_frac_label and cnt > 0:
                ax.text(
                    i, bottom + frac/2.0, f"{int(round(cnt)):,}",
                    ha='center', va='center', fontsize=counts_label_fs
                )

            # record for table
            if return_table:
                table_rows.append({
                    'gene': gene,
                    'cell_type': ct,
                    'isoform': iso,
                    'count': int(round(cnt)),
                    'fraction': float(frac),
                    'n_cells_expressing': n_cells_expressing,
                    'n_cells_total': n_cells_total
                })

            bottom += frac
            if iso not in legend_patches:
                legend_patches[iso] = Patch(color=col, label=iso)

        # cells label above the bar
        if annotate_cells:
            n_label = n_cells_total if cells_label == 'total' else n_cells_expressing
            ax.text(
                i, y_max + 0.02, f"n={n_label:,}",
                ha='center', va='bottom', fontsize=cells_label_fs, clip_on=False
            )

    # formatting
    ax.set_xticks(x)
    ax.set_xticklabels(present, rotation=45, ha='right', fontsize=tick_fs)
    ax.set_ylabel("Fraction of transcripts", fontsize=label_fs)
    ax.set_xlabel(f"Cell types sorted by fraction of {iso_most}", fontsize=label_fs)
    ax.set_title(f"Isoform usage for {gene}", fontsize=title_fs)
    ax.set_ylim(0, 1.10)  # headroom for n= labels

    ax.legend(
        handles=list(legend_patches.values()),
        bbox_to_anchor=(1.02, 1), loc='upper left',
        frameon=False, title='Transcript ID',
        fontsize=legend_fs, title_fontsize=legend_fs
    )
    # Adjust layout to accommodate legend
    plt.subplots_adjust(right=0.85)  # Leave space for legend
    try:
        plt.tight_layout()
    except:
        pass  # If tight_layout fails, subplots_adjust should handle it

    if return_table:
        df = pd.DataFrame(table_rows)
        df = df[['gene','cell_type','isoform','count','fraction','n_cells_expressing','n_cells_total']]
        return fig, df
    else:
        return fig

