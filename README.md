# pychesca

**Chemical Shift Covariance Analysis (CHESCA) in Python**

Implements the CHESCA framework from:

> Boulton, S., Akimoto, M., Selvaratnam, R., Bashiri, A., & Melacini, G. (2014).
> A tool set to map allosteric networks through the NMR chemical shift covariance analysis.
> *Scientific Reports*, 4, 7306.

## Run pychesca in the cloud

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Dan-Burns/pychesca/blob/main/pychesca_colab.ipynb)

## Installation

```bash
pip install git+https://github.com/Dan-Burns/pychesca.git
```

Or install from source in editable mode:

```bash
git clone https://github.com/Dan-Burns/pychesca.git
cd pychesca
pip install -e ".[dev]"
```

---

## Data Format

### CHESCA / SVD

A CSV file with a `RESI` index column and one column per perturbation state
(ligand, mutation, etc.):

```
RESI,apo,state1,state2,state3
10,0.012,0.045,0.078,0.031
11,0.003,0.009,0.011,0.007
...
```

### CHESPA

A CSV with `RESI` as the index and columns `refw1`, `refw2`, `Aw1`, `Aw2`,
`Bw1`, `Bw2` (`w1` = heteronucleus dimension, `w2` = proton):

```
RESI,refw1,refw2,Aw1,Aw2,Bw1,Bw2
10,120.1,8.21,120.4,8.35,119.9,8.10
...
```

---

## Python API

```python
import pandas as pd
from pychesca import HAC, SVD, Chespa
from pychesca import plot_corr, show_dendrogram, plot_svd, plot_chespa

# Load combined chemical shift data
df = pd.read_csv("my_shifts.csv", index_col="RESI")

# --- CHESCA clustering ---
hac = HAC(df, cutoff=98.0)
print(hac.clusters)         # DataFrame: residue → cluster ID
print(hac.n_clusters)       # number of clusters

plot_corr(hac, save_file="correlation.pdf")
show_dendrogram(hac, save_file="dendrogram.pdf")
hac.clusters.sort_values("cluster").to_csv("cluster_assignments.csv")

# --- SVD ---
dims = SVD(df)
plot_svd(dims, centering="column", save_file="svd.pdf")

# --- CHESPA (requires ref + A + B states) ---
chespa = Chespa("chespa_data.csv", het_nuc="N")
print(chespa.cos_theta)      # per-residue cos θ
print(chespa.X)              # fractional activation
plot_chespa(chespa, save_file="chespa.pdf")

# --- PyMOL export ---
from pychesca import clusters_to_pymol
clusters_to_pymol(hac.clusters, output="selections.pml")
```

---

## Command-Line Interface

```
chesca -file my_shifts.csv -output results/ -cutoff 98.0
```

Options:

| Flag | Default | Description |
|---|---|---|
| `-file` | *(required)* | Path to the combined chemical shift CSV |
| `-output` | *(required)* | Output directory |
| `-cutoff` | `98.0` | Correlation cutoff (%) for clustering |
| `-minimum_cluster` | `None` | Min cluster size to trigger sub-cluster state analysis |

Outputs saved to the output directory:

- `correlation_matrix.pdf`
- `dendrogram.pdf`
- `cluster_assignments.csv`
- `svd_plot.pdf`
- `pymol_selections.pml`
- `sub_cluster_<id>.pdf` *(if `-minimum_cluster` is set)*

---

## Development

```bash
# Run tests
pytest tests/ -v

# Run tests with coverage
pytest tests/ --cov=pychesca
```
