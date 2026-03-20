"""
Command-line interface for pychesca.

Entry point: ``chesca``
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

from pychesca.cluster import HAC
from pychesca.plots import heatmap_correlation_cutoffs, plot_corr, plot_svd, show_dendrogram
from pychesca.svd import SVD
from pychesca.visualize import clusters_to_pymol


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chesca",
        description=(
            "Chemical Shift Covariance Analysis (CHESCA). "
            "Clusters NMR residues by chemical shift correlation and optionally "
            "runs SVD and sub-cluster state analysis."
        ),
    )
    parser.add_argument(
        "-file",
        required=True,
        metavar="FILE",
        help="Path to the combined chemical shift CSV (index column must be named 'RESI').",
    )
    parser.add_argument(
        "-output",
        required=True,
        metavar="DIR",
        help="Output directory (will be created if it does not exist).",
    )
    parser.add_argument(
        "-cutoff",
        type=float,
        default=98.0,
        metavar="CUTOFF",
        help="Absolute correlation cutoff as a percentage (default: 98.0).",
    )
    parser.add_argument(
        "-minimum_cluster",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Minimum number of residues in a cluster to use for sub-cluster "
            "state analysis. Omit to skip sub-clustering."
        ),
    )
    parser.add_argument(
        "-linkage",
        default="single",
        choices=["complete", "single", "average", "ward"],
        metavar="METHOD",
        help=(
            "Linkage method for hierarchical clustering: "
            "'single' (default), 'complete', 'average', or 'ward'. "
            "Note: 'ward' requires metric='euclidean'."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """
    Run a full CHESCA analysis pipeline.

    Parameters
    ----------
    argv : list of str, optional
        Override ``sys.argv`` — useful for testing.

    Returns
    -------
    int
        Exit code (0 = success).
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    # --- Load data ---
    try:
        df = pd.read_csv(args.file, index_col="RESI")
    except FileNotFoundError:
        print(f"ERROR: File not found: {args.file!r}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: Failed to read '{args.file}': {e}", file=sys.stderr)
        return 1

    # --- Output directory ---
    os.makedirs(args.output, exist_ok=True)
    out = args.output
    cutoff = args.cutoff
    print(f"Correlation cutoff: {cutoff}%")

    # --- Clustering ---
    print(f"Linkage method    : {args.linkage}")
    if args.minimum_cluster is not None:
        hac = HAC(df, cutoff, method=args.linkage, sub_cluster_cutoff=args.minimum_cluster)
    else:
        hac = HAC(df, cutoff, method=args.linkage)

    # --- Plots ---
    plot_corr(hac, cutoff, save_file=f"{out}/correlation_matrix.pdf")
    show_dendrogram(hac, save_file=f"{out}/dendrogram.pdf")
    hac.clusters.sort_values("cluster").to_csv(f"{out}/cluster_assignments.csv")

    # --- SVD ---
    dims = SVD(df)
    plot_svd(dims, centering="column", save_file=f"{out}/svd_plot.pdf")

    # --- Sub-cluster state dendrograms ---
    if args.minimum_cluster is not None:
        for cluster_id in hac.sub_cluster_ids:
            sub_resis = hac.clusters[hac.clusters["cluster"] == cluster_id].index
            state_corr = df.loc[sub_resis].corr().abs()
            hac_states = HAC(state_corr, cluster_states=True)
            show_dendrogram(
                hac_states,
                orientation="top",
                annotate_clusters=False,
                sub_cluster=cluster_id,
                save_file=f"{out}/sub_cluster_{cluster_id}.pdf",
            )

    # --- PyMOL output ---
    if args.minimum_cluster is not None:
        pml_df = hac.clusters[hac.clusters["cluster"].isin(hac.sub_cluster_ids)]
    else:
        pml_df = hac.clusters

    clusters_to_pymol(pml_df, output=f"{out}/pymol_selections.pml")

    print(f"Done. Results written to {out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
