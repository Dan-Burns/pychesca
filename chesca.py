import numpy as np
import pandas as pd 
from sklearn.cluster import AgglomerativeClustering as AG
from scipy.cluster.hierarchy import dendrogram 
import scipy
import os
import sys

cs_file = sys.argv[1]
correlation_cutoff = sys.argv[2]
output_folder = sys.argv[3]


#save dendrograms
#save cluster label spreadsheet "dfc"
