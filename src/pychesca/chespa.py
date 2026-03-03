"""
Chemical Shift Projection Analysis (CHESPA).

Identifies false-positive CHESCA clusters by projecting chemical shift
perturbations onto a reference axis and computing fractional activation.

Reference
---------
Boulton, S., Akimoto, M., Selvaratnam, R., Bashiri, A., & Melacini, G. (2014).
A tool set to map allosteric networks through the NMR chemical shift covariance
analysis. Scientific Reports, 4, 7306.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Heteronucleus-to-scaling-coefficient mapping (SI eq. 5)
_HET_COEF: dict[str, float] = {
    "N": 0.2,
    "C": 0.25,
}


class Chespa:
    """
    Chemical Shift Projection Analysis (CHESPA).

    Computes the combined chemical shift distance (Δδ), the cosine of the
    angle between perturbation vectors (cos θ), and the fractional activation
    (X) for each residue.

    Parameters
    ----------
    file_path : str or Path
        Path to a CSV with columns ``RESI``, ``refw1``, ``refw2``,
        ``Aw1``, ``Aw2``, ``Bw1``, ``Bw2``.
        ``w1`` is the heteronucleus dimension; ``w2`` is the proton dimension.
        ``ref`` = reference state, ``A`` = antagonist/reverse-agonist,
        ``B`` = agonist/activated.
    het_nuc : str
        Heteronucleus type — ``'N'`` (nitrogen) or ``'C'`` (carbon).

    Attributes
    ----------
    resis : list
        Residue identifiers (index of the input file).
    ref_to_A : np.ndarray
        Combined chemical shift distance from ref to state A (SI eq. 5).
    ref_to_B : np.ndarray
        Combined chemical shift distance from ref to state B.
    cos_theta : np.ndarray
        Per-residue cosine of the angle between Δref→A and Δref→B (SI eq. 6).
    X : np.ndarray
        Fractional activation (SI eq. 7).  Positive → towards activation.
    """

    def __init__(self, file_path: str, het_nuc: str = "N") -> None:
        het_nuc = het_nuc.upper()
        if het_nuc not in _HET_COEF:
            raise ValueError(f"het_nuc must be 'N' or 'C', got {het_nuc!r}")
        self.het_nuc = het_nuc
        self.het_coef = _HET_COEF[het_nuc]

        df = self._load(file_path)
        self.df = df
        self.resis = df.index.tolist()

        self.states = {
            "ref": df[["refw1", "refw2"]].values,
            "A": df[["Aw1", "Aw2"]].values,
            "B": df[["Bw1", "Bw2"]].values,
        }

        self.ref_to_B = self._state_distance("B")
        self.ref_to_A = self._state_distance("A")
        self.cos_theta = self._state_angle()
        self.X = self._fractional_activation()

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _load(file_path: str) -> pd.DataFrame:
        """Load and validate the CSV, returning a DataFrame indexed by RESI."""
        try:
            df = pd.read_csv(file_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"CHESPA data file not found: {file_path!r}")

        # Find the RESI column case-insensitively — it may not be the very first
        # column if pandas reads an unnamed auto-index column first.
        resi_col = next(
            (col for col in df.columns if col.upper() == "RESI"),
            None,
        )
        if resi_col is None:
            raise ValueError(
                "No 'RESI' column found in the CSV. "
                "Ensure one column is named 'RESI' (case-insensitive)."
            )

        df = df.set_index(resi_col)
        df.index.name = "RESI"
        return df

    def _state_distance(self, state: str) -> np.ndarray:
        """
        Combined chemical shift distance from ref → *state* (SI eq. 5).

        Δδ_comb = sqrt((ΔδH)² + (α·ΔδX)²)
        where α is the heteronucleus scaling coefficient.
        """
        ref = self.states["ref"]
        s = self.states[state]
        delta_h = (ref[:, 1] - s[:, 1]) ** 2
        delta_het = (self.het_coef * (ref[:, 0] - s[:, 0])) ** 2
        return np.sqrt(delta_h + delta_het)

    def _state_angle(self) -> np.ndarray:
        """Cosine of angle between per-residue perturbation vectors (SI eq. 6)."""
        ref, A, B = self.states["ref"], self.states["A"], self.states["B"]
        vb = (B - ref) / np.linalg.norm(B - ref)
        va = (A - ref) / np.linalg.norm(A - ref)
        dot = (va * vb).sum(axis=1)
        denominator = self.ref_to_A * self.ref_to_B
        cos_theta = dot / denominator
        return np.nan_to_num(cos_theta, nan=0.0)

    def _fractional_activation(self) -> np.ndarray:
        """
        Fractional activation X (SI eq. 7).

        Normalized projection of Δref→A onto Δref→B.
        X > 0 → towards activation; X < 0 → towards inactivation.
        """
        ref, A, B = self.states["ref"], self.states["A"], self.states["B"]
        vb = (B - ref) / np.linalg.norm(B - ref)
        va = (A - ref) / np.linalg.norm(A - ref)
        dot = (va * vb).sum(axis=1)
        return dot / (vb * vb).sum(axis=1)
