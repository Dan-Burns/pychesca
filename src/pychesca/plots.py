"""
Plotting functions for pychesca CHESCA and CHESPA analyses.
"""

from __future__ import annotations

import warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram


# ---------------------------------------------------------------------------
# CHESCA — correlation heatmap
# ---------------------------------------------------------------------------

def plot_corr(hac, cutoff: float | None = None, save_file: str | None = None) -> None:
    """
    Plot a binary correlation heatmap showing residue pairs above *cutoff*.

    Parameters
    ----------
    hac : pychesca.cluster.HAC
    cutoff : float, optional
        Override the cutoff stored on *hac*. Same units as ``hac.cutoff``.
    save_file : str, optional
        If provided, save the figure to this path.
    """
    sns.set(font_scale=2)
    effective_cutoff = hac.cutoff if cutoff is None else cutoff
    corr = hac.get_corr_above(effective_cutoff)
    fig, ax = plt.subplots(figsize=(25, 20))
    sns.heatmap(corr, cbar=False, ax=ax)
    ax.invert_yaxis()
    if save_file is not None:
        fig.savefig(save_file)


# ---------------------------------------------------------------------------
# CHESCA — dendrogram
# ---------------------------------------------------------------------------

def _plot_dendrogram(model, **kwargs) -> None:
    """
    Wrap a scikit-learn ``AgglomerativeClustering`` model for scipy's dendrogram.

    Adapted from:
    https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
    """
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)
    dendrogram(linkage_matrix, **kwargs)


def _get_cluster_annotation_positions(clusters, threshold, ax, orientation, max_cor):
    """
    Compute x/y positions for cluster ID annotations on a dendrogram.

    Parameters
    ----------
    clusters : pd.DataFrame
        ``HAC.clusters`` — 'cluster' column indexed by residue/state.
    threshold : float
        The distance threshold used for coloring branches.
    ax : matplotlib.axes.Axes
    orientation : {'top', 'right'}
    max_cor : float
        Maximum correlation distance in the data (used when all residues
        cluster below the threshold).

    Returns
    -------
    xs, ys : list, list
    """
    ordered_indices = []
    x_positions = []

    tick_labels = ax.get_xticklabels() if orientation == "top" else ax.get_yticklabels()
    for t in tick_labels:
        text = t.get_text()
        try:
            ordered_indices.append(float(text))
        except ValueError:
            ordered_indices.append(text)
        x_positions.append(t.get_position())

    x_positions = np.array(x_positions)

    xs = []
    for cluster_id in range(clusters["cluster"].max() + 1):
        mask = clusters.loc[ordered_indices]["cluster"] == cluster_id
        if orientation == "top":
            vals = x_positions[np.where(mask)][:, 0]
        else:
            vals = x_positions[np.where(mask)][:, 1]
        xs.append(float(np.median(vals)))

    y_val = max_cor if max_cor < threshold else threshold
    ys = [y_val] * len(xs)
    return xs, ys


def show_dendrogram(
    hac,
    save_file: str | None = None,
    orientation: str = "right",
    leaf_rotation: float | None = None,
    annotate_clusters: bool = True,
    sub_cluster: int | None = None,
) -> None:
    """
    Plot a CHESCA dendrogram.

    Parameters
    ----------
    hac : pychesca.cluster.HAC
    save_file : str, optional
        Path to save the figure.
    orientation : {'right', 'top'}
        Direction the dendrogram grows. 'right' → tall figure; 'top' → wide.
    leaf_rotation : float, optional
        Rotation angle for leaf axis tick labels.
    annotate_clusters : bool
        Whether to add numeric cluster ID annotations.
    sub_cluster : int, optional
        If set, the title will read "Cluster *N* States".
    """
    max_cor = float(hac.corr_distance.max().max())
    threshold = 100 - hac.cutoff

    figsize = (8, 5) if orientation == "top" else (15, 40)
    fig, ax = plt.subplots(figsize=figsize)

    _plot_dendrogram(
        hac.model,
        color_threshold=threshold,
        labels=hac.corr_distance.index,
        ax=ax,
        orientation=orientation,
        leaf_rotation=leaf_rotation,
        leaf_font_size=14,
    )

    if annotate_clusters:
        xs, ys = _get_cluster_annotation_positions(
            hac.clusters, threshold, ax, orientation, max_cor
        )
        if orientation == "right":
            xs, ys = ys, xs
        for i, (x, y) in enumerate(zip(xs, ys)):
            ax.text(x, y, str(i + 1))

    ax.set_title(f"Cluster {sub_cluster} States" if sub_cluster is not None else "CHESCA Clusters")
    ax.grid(visible=False)
    ax.set_facecolor("white")
    ax.xaxis.set_tick_params(labelsize=15)

    if save_file is not None:
        fig.savefig(save_file)


# ---------------------------------------------------------------------------
# SVD biplot
# ---------------------------------------------------------------------------

def plot_svd(
    svd,
    centering: str = "column",
    save_file: str | None = None,
) -> None:
    """
    Biplot of the first two SVD components (residues + states).

    Parameters
    ----------
    svd : pychesca.svd.SVD
    centering : {'column', 'row'}
        Which centered SVD to plot.
    save_file : str, optional
        Path to save the figure.
    """
    if centering == "column":
        data = svd.column_svd
    elif centering == "row":
        data = svd.row_svd
    else:
        raise ValueError("centering must be 'column' or 'row'")

    uds = data["uds"]
    v = data["V"]

    x_min, x_max = uds[:, 0].min(), uds[:, 0].max()
    y_min, y_max = uds[:, 1].min(), uds[:, 1].max()

    fig, ax = plt.subplots()
    ax.scatter(uds[:, 0], uds[:, 1], marker="o", facecolors="none", color="black")
    ax.vlines(0, y_min, y_max, color="blue")
    ax.hlines(0, x_min, x_max, color="red")
    ax.scatter(v[:, 0], v[:, 1], marker="D", color="magenta")
    for i, label in enumerate(svd.cols):
        ax.annotate(label, (v[i, 0], v[i, 1]), color="magenta")

    if save_file is not None:
        fig.savefig(save_file)


# ---------------------------------------------------------------------------
# Correlation cutoff heatmap
# ---------------------------------------------------------------------------

def heatmap_correlation_cutoffs(
    df,
    min_corr: float = 94.0,
    save_file: str | None = None,
) -> None:
    """
    Heatmap of pairwise absolute correlation coefficients above *min_corr*.

    Useful for choosing the CHESCA clustering cutoff (see Protein NMR book,
    chapter 18, p. 403 Note 5).

    Parameters
    ----------
    df : pd.DataFrame
        Combined chemical shift DataFrame.
    min_corr : float
        Minimum correlation to display (percentage or fraction).
    save_file : str, optional
        Path to save the figure.
    """
    corr = df.T.corr().abs().fillna(0)
    threshold = min_corr / 100 if min_corr > 1 else min_corr

    mask = np.triu(np.ones_like(corr, dtype=bool)) | (corr.abs() < threshold)

    fig, ax = plt.subplots(figsize=(25, 20))
    sns.heatmap(corr, mask=mask, cmap="plasma", square=True, linewidths=0.5, ax=ax)
    ax.set_title(f"Correlation Coefficients Above {min_corr}")

    if save_file is not None:
        fig.savefig(save_file)


# ---------------------------------------------------------------------------
# CHESPA bar plots
# ---------------------------------------------------------------------------

def plot_chespa(chespa, save_file: str | None = None) -> None:
    """
    Three-panel bar plot of CHESPA results: Δδ, cos θ, and fractional activation X.

    Parameters
    ----------
    chespa : pychesca.chespa.Chespa
    save_file : str, optional
        Path to save the figure.
    """
    D = "\u0394"   # Δ
    d = "\u03B4"   # δ
    theta = "\u03B8"  # θ

    resis = [str(r) for r in chespa.resis]
    fig, axs = plt.subplots(3, 1, figsize=(40, 20))

    data = [
        (chespa.ref_to_A, f"{D}{d}{chespa.het_nuc.upper()}Hcomb (ppm)"),
        (chespa.cos_theta, f"cos{theta}"),
        (chespa.X, "X (fractional activation)"),
    ]

    for ax, (values, title) in zip(axs, data):
        ax.bar(resis, values)
        ax.grid(visible=False)
        ax.set_title(title)
        ax.set_xticklabels(resis, fontdict={"size": 12}, rotation=90)

    fig.tight_layout()
    if save_file is not None:
        fig.savefig(save_file)
