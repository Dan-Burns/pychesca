"""
File I/O utilities for pychesca.
"""

from __future__ import annotations

import pandas as pd


def open_file(file: str) -> pd.DataFrame:
    """
    Open a spreadsheet or delimited text file into a :class:`pandas.DataFrame`.

    Supported formats: ``.xlsx``, ``.xls``, ``.tsv``, ``.csv`` (default).

    Parameters
    ----------
    file : str
        Path to the input file.

    Returns
    -------
    pd.DataFrame
    """
    ext = file.rsplit(".", 1)[-1].lower()
    if ext in ("xlsx", "xls"):
        return pd.read_excel(file)
    elif ext == "tsv":
        return pd.read_csv(file, sep="\t")
    else:
        return pd.read_csv(file)
