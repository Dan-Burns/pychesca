"""
Hierarchical Agglomerative Clustering (HAC) for CHESCA.

Reference
---------
Boulton, S., Akimoto, M., Selvaratnam, R., Bashiri, A., & Melacini, G. (2014).
A tool set to map allosteric networks through the NMR chemical shift covariance
analysis. Scientific Reports, 4, 7306.
"""

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering


class HAC:
    """
    Hierarchical Agglomerative Clustering on NMR combined chemical shift data.

    Parameters
    ----------
    df : pd.DataFrame
        Combined chemical shift DataFrame. Rows = residues (indexed by RESI),
        columns = perturbation states.
    cutoff : float
        Correlation cutoff as a percentage (e.g. 98.0) or fraction (e.g. 0.98).
        Residues with absolute correlation above this threshold are grouped.
    method : str
        Linkage method passed to ``AgglomerativeClustering``. Default 'complete'.
    metric : str
        Distance metric. Default 'euclidean'.
    cluster_states : bool
        If ``True``, cluster the perturbation *states* instead of residues.
        In this mode, ``df`` should already be a correlation matrix.
    sub_cluster_cutoff : int or None
        If provided, compute which clusters have more than this many residues,
        stored in :attr:`sub_cluster_ids`.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        cutoff: float = 98.0,
        method: str = "complete",
        metric: str = "euclidean",
        cluster_states: bool = False,
        sub_cluster_cutoff: int | None = None,
    ):
        self.df = df
        # Support cutoff expressed as fraction (e.g. 0.98) or percentage (e.g. 98)
        self.cutoff = cutoff * 100 if cutoff < 1 else cutoff

        if cluster_states:
            self.absolute_corr = df.corr().abs()
        else:
            self.absolute_corr = df.T.corr().abs().fillna(0)

        self.corr_distance = 1 - self.absolute_corr
        self.method = method
        self.metric = metric

        if self.df.index.duplicated().any():
            raise ValueError(
                "Duplicate indices found in the DataFrame. Please remove them before proceeding."
            )

        self.model: AgglomerativeClustering = self._fit_model(method=method, metric=metric)
        self.clusters: pd.DataFrame = self._get_clusters()
        self.n_clusters: int = int(self.clusters["cluster"].max()) + 1

        if sub_cluster_cutoff is not None:
            self.sub_cluster_ids: set[int] = self.get_clusters_above_cutoff(cutoff=sub_cluster_cutoff)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _fit_model(self, method: str = "complete", metric: str = "euclidean") -> AgglomerativeClustering:
        """Fit AgglomerativeClustering on the correlation-distance matrix."""
        model = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=100 - self.cutoff,
            metric=metric,
            linkage=method,
        )
        model.fit(self.corr_distance)
        return model

    def _get_clusters(self) -> pd.DataFrame:
        """Return a DataFrame mapping each residue/state to its cluster label."""
        dfc = pd.DataFrame(index=self.corr_distance.index)
        dfc["cluster"] = self.model.labels_
        return dfc

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_corr_above(self, cutoff: float) -> pd.DataFrame:
        """
        Return a boolean DataFrame of entries with absolute correlation above *cutoff*.

        Parameters
        ----------
        cutoff : float
            Threshold as a percentage (>1) or fraction (≤1).
        """
        threshold = cutoff / 100 if cutoff > 1 else cutoff
        return self.absolute_corr > threshold

    def get_clusters_above_cutoff(self, cutoff: int = 3) -> set[int]:
        """
        Return the set of cluster IDs that contain more than *cutoff* residues.

        Parameters
        ----------
        cutoff : int
            Minimum number of members for a cluster to be included.
        """
        return {
            c_id
            for c_id in self.clusters["cluster"].unique()
            if len(self.clusters[self.clusters["cluster"] == c_id]) > cutoff
        }
