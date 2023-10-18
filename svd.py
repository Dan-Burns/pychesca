#https://towardsdatascience.com/singular-value-decomposition-and-its-applications-in-principal-component-analysis-5b7a5f08d0bd
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html
#https://stats.stackexchange.com/questions/435338/what-do-the-matrix-s-u-v-returned-by-singular-value-decomposition-represent

# Protein NMR Book SVD section page 405

import scipy


# calculate the row mean ref to use as reference state
# wild type reference will probably be better
chnh_shifts['mean'] = chnh_shifts[['H189G', 'H189N', 'H189Q']].mean(axis=1)

# CCS difference from mean 
mean_dif = chnh_shifts[['H189G', 'H189N', 'H189Q']].subtract(chnh_shifts['mean'], axis=0)

# column center the data
centered_difs = mean_dif.subtract(mean_dif.mean(axis=0),axis=1)

# svd
U, s, Vh = scipy.linalg.svd(centered_difs)

# protein nmr svd section
s_square = s**2
s_sum = s_square.sum()
explained_variance = s_square/s_sum