"""
PyMOL visualization output for CHESCA cluster results.

Generates ``.pml`` script files that color residues in PyMOL according
to their cluster assignments.
"""

from __future__ import annotations

import pandas as pd

# PyMOL named colors, cycled for cluster assignments
_COLORS = [
    "red", "blue", "yellow", "black", "white", "purpleblue", "magenta",
    "orange", "pink", "purple", "lightblue", "teal", "marine",
    "brightorange", "chocolate", "sand", "slate", "neptunium", "zinc", "radium",
]


def clusters_to_pymol(df: pd.DataFrame, output: str) -> None:
    """
    Write a PyMOL ``.pml`` selection script from a cluster assignment DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with a ``'cluster'`` column, indexed by residue identifier
        (e.g. from :attr:`pychesca.cluster.HAC.clusters`).
    output : str
        File path for the output ``.pml`` file.
    """
    clusters: dict[int, list] = {}
    for cluster_label in sorted(df["cluster"].unique()):
        clusters[int(cluster_label)] = list(df.loc[df["cluster"] == cluster_label].index)

    with open(output, "w") as f:
        f.write("set sphere_scale, 0.8\n")
        for i, (cluster_id, residues) in enumerate(clusters.items()):
            color = _COLORS[i % len(_COLORS)]
            # int(float(r)) handles decimal naming conventions (e.g. ILE methyl pseudo-residues)
            sel = "+".join(str(int(float(r))) for r in residues)
            f.write(f"select c_{cluster_id}, resi {sel}\n")
            f.write(f"color {color}, c_{cluster_id}\n")
            f.write(f"show spheres, c_{cluster_id}\n")
