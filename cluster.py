import numpy as np
import pandas as pd 
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
import scipy
import matplotlib.pyplot as plt


class HAG:
    '''
    Handle the Hiearchical Agglomerative Clustering of covarying chemical shifts

    TODO: handle reference state column
    '''

    def __init__(self,
                 df):
        self.df = df
        self.absolute_corr = df.T.corr().abs().fillna(0)
        self.corr_distance = 1 - self.absolute_corr
        self.distance_matrix = self.get_distance_matrix()
        self.linkage = None


    def get_corr_above(self, cutoff):
        'Return a df with correlation above cutoff'
        return self.absolute_corr > cutoff/100


    def get_distance_matrix(self, method='complete', metric='euclidean'):
        '''
        https://predictivehacks.com/hierarchical-clustering-in-python/
        Returns 
        -------
        scipy.cluster.hierarchy.linkage
        '''

        self.linkage = linkage(self.corr_distance, method=method, metric=metric)

    def get_dendrogram(self, cutoff=98, method='complete', metric='euclidean',
                       ax=None):
        '''
        cutoff : int of float
            distance (1-100) to consider in the same cluster
        '''
        # TODO
        # If you use this function and specify different method or metric
        # the dendrogram will not reflect the changes
        # Have to run get_distance_matrix with the method and metric specified
        # there first
        if (self.linkage is None) or ((method !='complete') and\
                                       (metric != 'euclidean')):
            self.get_distance_matrix(method, metric)
        # Create a dendrogram (change color threshhold to see clusters above cutoff 2 means 98% and up)
        return dendrogram(self.linkage, 
                          color_threshold=100-cutoff, ax=ax)

    def get_clusters(self, cutoff=98, criterion='distance',
                    method='complete', 
                    metric='euclidean'):
        

        self.dfc = pd.DataFrame(index=self.corr_distance.index)
        # Assign cluster labels
        if self.linkage == None:
            self.linkage = self.get_distance_matrix(method, metric)

        self.dfc['cluster_labels'] = fcluster(self.linkage, 100-cutoff, 
                                              criterion)

# clustering of clusters is recommended.


# complete-linkage AC and fragment clustering
# Boulton, S., Akimoto, M., Selvaratnam, R., Bashiri, A., & Melacini, G. (2014). 
# A tool set to map allosteric networks through the NMR chemical shift covariance analysis. 
# Scientific Reports, 4, 7306.