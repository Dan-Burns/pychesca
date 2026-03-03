import numpy as np
import pandas as pd 
from sklearn.cluster import AgglomerativeClustering as AG
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

# probably not going to need this.
def determine_filetype(file):
    return file.split('.')[-1]

def open_file(file):
    
    ext = determine_filetype(file)
    if ext == 'xlsx' or ext == 'xls':
        return pd.read_excel(file)
    elif ext == 'tsv':
        return pd.read_csv(file, sep='\t')
    else:
        return pd.read_csv(file)

