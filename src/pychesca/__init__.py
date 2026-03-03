"""
pychesca — Chemical Shift Covariance Analysis in Python

Implements the CHESCA framework from Boulton et al. (2014),
Scientific Reports, 4, 7306.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pychesca")
except PackageNotFoundError:
    __version__ = "unknown"

from pychesca.cluster import HAC
from pychesca.svd import SVD
from pychesca.chespa import Chespa
from pychesca.plots import (
    plot_corr,
    show_dendrogram,
    plot_svd,
    plot_chespa,
    heatmap_correlation_cutoffs,
)
from pychesca.visualize import clusters_to_pymol
from pychesca.utils import open_file

__all__ = [
    "HAC",
    "SVD",
    "Chespa",
    "plot_corr",
    "show_dendrogram",
    "plot_svd",
    "plot_chespa",
    "heatmap_correlation_cutoffs",
    "clusters_to_pymol",
    "open_file",
]
