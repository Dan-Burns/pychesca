import numpy as np
import pandas as pd 
from scipy.cluster.hierarchy import fcluster, linkage
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
            # can loop through all the indices and make all the sub cluster dendrograms available here.
            self.sub_clusters = 
            for clust_id in self.sub_clusters

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
        scipy.cluster.hierarchy.linkage
        '''

        self.linkage = linkage(self.corr_distance, method=method, metric=metric)

    # def get_dendrogram(self, method='complete', metric='euclidean',
    #                    ax=None):
    #     '''
    #     cutoff : int of float
    #         distance (1-100) to consider in the same cluster
    #     '''
    #     # TODO
    #     # If you use this function and specify different method or metric
    #     # the dendrogram will not reflect the changes
    #     # Have to run get_distance_matrix with the method and metric specified
    #     # there first
        
    #     if (self.linkage is None) or ((method !='complete') or\
    #                                    (metric != 'euclidean')):
    #         self.get_distance_matrix(method, metric)
    #     # Create a dendrogram (change color threshhold to see clusters above cutoff 2 means 98% and up)
    #     return dendrogram(self.linkage, 
    #                       color_threshold=100-self.cutoff, 
    #                       labels=self.df.index,
    #                       ax=ax)

    def get_clusters(self,criterion='distance',
                    method='complete', 
                    metric='euclidean'):
    

        dfc = pd.DataFrame(index=self.corr_distance.index)
        # Assign cluster labels
        if self.linkage is None:
            self.get_distance_matrix(method, metric)

        dfc['cluster'] = fcluster(self.linkage, 100-self.cutoff, 
                                              criterion)
        
        return dfc

# clustering of clusters is recommended.

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