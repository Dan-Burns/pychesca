"""Tests for pychesca.chespa.Chespa."""

import pytest
import numpy as np
from pathlib import Path

from pychesca.chespa import Chespa

DATA = Path(__file__).parent / "data" / "chespa_cs.csv"


@pytest.fixture
def chespa_n():
    return Chespa(str(DATA), het_nuc="N")


def test_chespa_constructs(chespa_n):
    assert chespa_n is not None


def test_chespa_resis_match_index(chespa_n):
    # resis should be the DataFrame index, not a column
    assert len(chespa_n.resis) > 0
    # All resis should also exist in the underlying df index
    assert set(chespa_n.resis) == set(chespa_n.df.index.tolist())


def test_cos_theta_range(chespa_n):
    cos = chespa_n.cos_theta
    assert np.all(cos >= -1.0 - 1e-6), "cos_theta below -1"
    assert np.all(cos <= 1.0 + 1e-6), "cos_theta above 1"


def test_ref_to_a_nonnegative(chespa_n):
    assert np.all(chespa_n.ref_to_A >= 0)


def test_ref_to_b_nonnegative(chespa_n):
    assert np.all(chespa_n.ref_to_B >= 0)


def test_file_not_found_raises():
    with pytest.raises(FileNotFoundError):
        Chespa("definitely_not_a_real_file.csv")


def test_invalid_het_nuc_raises():
    with pytest.raises(ValueError, match="het_nuc"):
        Chespa(str(DATA), het_nuc="X")


def test_het_nuc_normalized_lowercase():
    # lowercase 'n' should work the same as 'N'
    c = Chespa(str(DATA), het_nuc="n")
    assert c.het_nuc == "N"
    assert c.het_coef == 0.2
