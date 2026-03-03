"""Tests for pychesca.svd.SVD."""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path

from pychesca.svd import SVD

DATA = Path(__file__).parent / "data" / "test_data.csv"


@pytest.fixture
def df():
    return pd.read_csv(DATA, index_col="RESI")


@pytest.fixture
def svd_obj(df):
    return SVD(df)


def test_svd_constructs(svd_obj):
    assert svd_obj is not None


@pytest.mark.parametrize("key", ["U", "S", "V", "sdiag", "uds", "s_variance", "vdf"])
def test_row_svd_keys(svd_obj, key):
    assert key in svd_obj.row_svd


@pytest.mark.parametrize("key", ["U", "S", "V", "sdiag", "uds", "s_variance", "vdf"])
def test_column_svd_keys(svd_obj, key):
    assert key in svd_obj.column_svd


def test_s_variance_sums_to_one(svd_obj):
    for centering in ("row_svd", "column_svd"):
        var = getattr(svd_obj, centering)["s_variance"]
        assert abs(var.sum() - 1.0) < 1e-6, f"{centering} s_variance does not sum to 1"


def test_uds_shape(svd_obj, df):
    n_residues, n_states = df.shape
    uds = svd_obj.column_svd["uds"]
    assert uds.shape[0] == n_residues


def test_vdf_index_matches_columns(svd_obj, df):
    vdf = svd_obj.column_svd["vdf"]
    assert list(vdf.index) == list(df.columns)
