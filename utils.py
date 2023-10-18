import numpy as np
import pandas as pd 
from sklearn.cluster import AgglomerativeClustering as AG
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage


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

def get_corr(df, correlation_cutoff):
    return df.T.corr().abs() > correlation_cutoff

def get_corr_distance(df):
    return 1-df.T.corr().abs().fillna(0)


def get_distance_matrix(corr_dist, method='complete', metric='euclidean'):
    '''
    df : pd.DataFrame of combined chemical shifts

    

    https://predictivehacks.com/hierarchical-clustering-in-python/
    '''
 
    return linkage(corr_dist, method=method, metric=metric)

def get_dendrogram(distance_matrix, cutoff):
    '''
    cutoff : int of float
        distance (1-100) to consider in the same cluster
    '''
    # Create a dendrogram (change color threshhold to see clusters above cutoff 2 means 98% and up)
    return dendrogram(distance_matrix, color_threshold=100-cutoff)

def get_clusters(distance_matrix, corr_dist, cutoff):
    dfc = pd.DataFrame(index=corr_dist.index)
    # Assign cluster labels (2 is 98%)
    dfc['cluster_labels'] = fcluster(distance_matrix, 100-cutoff, criterion='distance')
    return dfc