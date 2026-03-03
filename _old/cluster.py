import numpy as np
import pandas as pd 
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import AgglomerativeClustering
import scipy
import matplotlib.pyplot as plt


class HAC:
    '''
    Handle the Hiearchical Agglomerative Clustering 

    TODO: handle reference state column
    
    '''

    def __init__(self,
                 df,
                 cutoff=98,
                 method='complete', 
                 metric='euclidean',
                 cluster_states=False,
                 sub_cluster_cutoff=None):
        self.df = df
        # dealing with cutoff as a percentage since the clustering algorithms
        # don't treat it as a decimal 
        if cutoff < 1:
            self.cutoff = cutoff * 100
        else:
            self.cutoff = cutoff
        if cluster_states == False:
            self.absolute_corr = df.T.corr().abs().fillna(0)
        else:
            self.absolute_corr = df.corr().abs()
        self.corr_distance = 1 - self.absolute_corr
        self.method = method
        self.metric = metric
        self.linkage = self.get_distance_matrix(method=self.method,metric=self.metric)
        
        self.clusters = self.get_clusters()
        self.n_clusters = self.clusters['cluster'].max()
        if sub_cluster_cutoff is not None:
            self.sub_cluster_ids = self.get_clusters_above_cutoff(cutoff=sub_cluster_cutoff)
        check_duplicate = self.check_duplicate_indices()
        if check_duplicate:
            raise ValueError('Duplicate indices found in the dataframe. Please fix this before proceeding.')
        
            
    def check_duplicate_indices(self):
        'Check if there are duplicate indices in the dataframe'
        if self.df.index.duplicated().any():
            print('Duplicate indices found in the dataframe. Please fix this before proceeding.')
            return True
        else:
            return False

    def get_corr_above(self, cutoff):
        'Return a df with correlation above cutoff'
        if cutoff > 1:
            return self.absolute_corr > cutoff/100
        else:
            return self.absolute_corr > cutoff


    def get_distance_matrix(self, method='complete', metric='euclidean'):
        '''
        https://predictivehacks.com/hierarchical-clustering-in-python/
        Returns 
        -------
        AgglomerativeClustering
        '''
        self.linkage = AgglomerativeClustering(n_clusters=None, 
                                               distance_threshold=100-self.cutoff, 
                                               metric=metric, 
                                               linkage=method)
        return self.linkage.fit(self.corr_distance)

    def get_clusters(self, criterion='distance',
                    method='complete', 
                    metric='euclidean'):
    
        dfc = pd.DataFrame(index=self.corr_distance.index)
        # Assign cluster labels
        if self.linkage is None:
            self.get_distance_matrix(method, metric)

        dfc['cluster'] = self.linkage.labels_
        
        return dfc

    def get_clusters_above_cutoff(self, cutoff=3):
        '''
        Return cluster ids for clusters that have more than {cutoff} residues.

        '''
        real_clusters = set()
        for c_id in self.clusters['cluster'].unique():
            if len(self.clusters[self.clusters['cluster']==c_id]) > cutoff:
                real_clusters.add(c_id)

        return real_clusters


# complete-linkage AC and fragment clustering
# Boulton, S., Akimoto, M., Selvaratnam, R., Bashiri, A., & Melacini, G. (2014). 
# A tool set to map allosteric networks through the NMR chemical shift covariance analysis. 
# Scientific Reports, 4, 7306.
