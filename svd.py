#https://towardsdatascience.com/singular-value-decomposition-and-its-applications-in-principal-component-analysis-5b7a5f08d0bd
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html
#https://stats.stackexchange.com/questions/435338/what-do-the-matrix-s-u-v-returned-by-singular-value-decomposition-represent

# Protein NMR Book SVD section page 405

import scipy
import pandas as pd

class SVD:
    '''
    Handles SVD of chemical shift data

    TODO : deal with reference state - add another attribute
    '''

    def __init__(self,
                 df):
        self.df = df
        self.cols = df.columns
        self.row_mean = df.mean(axis=1)
        self.row_centered = df[self.cols].subtract(self.row_mean.values, axis=0)
        self.column_mean = df.mean(axis=0)
        self.column_centered = df[self.cols].subtract(self.column_mean.values, axis=1)
        # row SVD
        u, s, v = scipy.linalg.svd(self.row_centered)
        sdiag = scipy.linalg.diagsvd(s,*self.row_centered.shape)
        ### s_variance sum of first two PCs should be greater than 90%
        ### or another reference state is needed
        ### TODO: make a check of row and mean centered s_variance[:2].sum()
        # SVD on row-centered data
        self.row_svd = {'U':u, 'S':s, 'V':v, 
                        'sdiag': sdiag,
                        # u dot s is score matrix of each res on each pc
                        'uds': u@sdiag, 
                        's_variance': (s**2/(sum(s**2))), # explained_variance
                        # rows are perturbation states
                        'vdf':pd.DataFrame(v, 
                        columns=[f'PC{i+1}' for i in range(len(self.cols))], 
                        index=self.cols)
                        }
        # SVD on column-centered data
        u, s, v = scipy.linalg.svd(self.column_centered)
        sdiag = scipy.linalg.diagsvd(s,*self.column_centered.shape)
        self.column_svd = {'U':u, 'S':s, 'V':v, 
                        'sdiag': sdiag,
                        # u dot s is score matrix of each res on each pc
                        'uds': u@sdiag,
                        's_variance': (s**2/(sum(s**2))), # explained_variance
                        # rows are perturbation states
                        'vdf':pd.DataFrame(v, 
                        columns=[f'PC{i+1}' for i in range(len(self.cols))], 
                        index=self.cols)
                        }
        # SVD on reference state-centered data
        ## Not implemented

# "If an antagonist or reverse agonist stat is not available, we suggest trying
# multiple reference states or the row mean centering approach."