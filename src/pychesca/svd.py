"""
Singular Value Decomposition of NMR combined chemical shift data.

Reference
---------
Boulton, S., Akimoto, M., Selvaratnam, R., Bashiri, A., & Melacini, G. (2014).
A tool set to map allosteric networks through the NMR chemical shift covariance
analysis. Scientific Reports, 4, 7306. (Supplementary, SVD section)

"If an antagonist or reverse agonist state is not available, we suggest trying
multiple reference states or the row mean centering approach."
"""

from __future__ import annotations

import pandas as pd
import scipy.linalg


class SVD:
    """
    Perform SVD on a combined chemical shift DataFrame.

    Two centerings are computed automatically:

    * **Row-centered** — subtract the per-residue mean across states.
    * **Column-centered** — subtract the per-state mean across residues.

    Both results are stored as dictionaries under :attr:`row_svd` and
    :attr:`column_svd` with keys: ``U``, ``S``, ``V``, ``sdiag``, ``uds``
    (the score matrix U·S), ``s_variance`` (fraction of variance explained),
    and ``vdf`` (V as a labeled DataFrame).

    The first two PCs should explain ≥90% of variance; if not, consider
    adding additional reference states.

    Parameters
    ----------
    df : pd.DataFrame
        Combined chemical shift data. Rows = residues, columns = states.
    """

    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df
        self.cols = df.columns

        self.row_svd = self._compute_svd(df.subtract(df.mean(axis=1), axis=0))
        self.column_svd = self._compute_svd(df.subtract(df.mean(axis=0), axis=1))

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _compute_svd(self, centered: pd.DataFrame) -> dict:
        u, s, v = scipy.linalg.svd(centered)
        sdiag = scipy.linalg.diagsvd(s, *centered.shape)
        n_cols = len(self.cols)
        return {
            "U": u,
            "S": s,
            "V": v,
            "sdiag": sdiag,
            "uds": u @ sdiag,
            "s_variance": s**2 / (s**2).sum(),
            "vdf": pd.DataFrame(
                v,
                columns=[f"PC{i + 1}" for i in range(n_cols)],
                index=self.cols,
            ),
        }
