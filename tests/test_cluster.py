"""Tests for pychesca.cluster.HAC."""

import pytest
import pandas as pd
from pathlib import Path

from pychesca.cluster import HAC

DATA = Path(__file__).parent / "data" / "test_data.csv"


@pytest.fixture
def df():
    return pd.read_csv(DATA, index_col="RESI")


def test_hac_constructs(df):
    hac = HAC(df, cutoff=98.0)
    assert hac is not None


def test_hac_cluster_df_shape(df):
    hac = HAC(df, cutoff=98.0)
    assert "cluster" in hac.clusters.columns
    assert set(hac.clusters.index) == set(df.index)


def test_hac_n_clusters_positive(df):
    hac = HAC(df, cutoff=98.0)
    assert hac.n_clusters >= 1


def test_get_corr_above_boolean_df(df):
    hac = HAC(df, cutoff=98.0)
    result = hac.get_corr_above(98.0)
    assert result.dtypes.unique().tolist() == [bool]


def test_get_corr_above_fraction_equiv(df):
    hac = HAC(df)
    pct = hac.get_corr_above(98)
    frac = hac.get_corr_above(0.98)
    pd.testing.assert_frame_equal(pct, frac)


def test_get_clusters_above_cutoff(df):
    hac = HAC(df, cutoff=90.0, sub_cluster_cutoff=2)
    # sub_cluster_ids should be a subset of all cluster ids
    assert hac.sub_cluster_ids.issubset(set(hac.clusters["cluster"].unique()))


def test_duplicate_index_raises(df):
    duped = pd.concat([df, df.iloc[:1]])
    with pytest.raises(ValueError, match="Duplicate"):
        HAC(duped)


def test_cutoff_fraction_normalised(df):
    hac_pct = HAC(df, cutoff=98.0)
    hac_frac = HAC(df, cutoff=0.98)
    assert hac_pct.cutoff == hac_frac.cutoff == 98.0
