import numpy as np
import pandas as pd 
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
import seaborn as sns
import scipy
import os
import sys
import argparse
import matplotlib.pyplot as plt
from cluster import HAC
from plots import *
from svd import SVD
from visualize import clusters_to_pymol
# Was the sandbox script
# this will end up importing the finished functions and run the whole analysis in one go


# Create the parser
parser = argparse.ArgumentParser(description='Load a csv into pandas')


parser.add_argument('-file', type=str, help='Path to the spreadsheet')
parser.add_argument('-output', type=str, help='Output folder path')
parser.add_argument('-cutoff', type=float, default=98.0, help='Correlation cutoff (%) for clustering')
parser.add_argument('-minimum_cluster', type=int, default=None, help='The minimum number of residues in an original cluster to us for subclustering.')
# Parse the command line arguments
args = parser.parse_args()

# Read the Excel file into a pandas DataFrame
try:
    df = pd.read_csv(args.file,index_col='RESI')
except Exception as e:
    print(f"Failed to read the spreadsheet file. Case Sensitive. Error: {e}")

# Set up folders
output_folder = args.output
os.makedirs(output_folder)
cutoff = args.cutoff
print(cutoff)

# Create a clustering object
if args.minimum_cluster is not None:
    hac = HAC(df, args.cutoff, sub_cluster_cutoff=args.minimum_cluster)
else:
    hac = HAC(df, args.cutoff)

# make the correlation plot
plot_corr(hac,args.cutoff, save_file=f'{output_folder}/correlation_matrix.pdf')


# dendrogram
show_dendrogram(hac, save_file=f'{output_folder}/dendrogram.pdf')


# SVD
dims = SVD(df)
plot_svd(dims, centering='column', save_file=f'{output_folder}/svd_plot.pdf')

if args.minimum_cluster is not None:
    for cluster_id in hac.sub_cluster_ids:
        sub_cluster_resis = hac.clusters[hac.clusters['cluster']==cluster_id].index
        state_corr = df.loc[sub_cluster_resis].corr().abs()
        hac_states = HAC(state_corr, cluster_states=True)
        show_dendrogram(hac_states, orientation='top', annotate_clusters=False,
                        sub_cluster=cluster_id,
                        save_file=f'{output_folder}/sub_cluster_{cluster_id}.pdf')



#Pymol output
if args.minimum_cluster is not None:
    clusters_to_pymol(hac.clusters[hac.clusters['cluster'].isin(hac.sub_cluster_ids)],
                  output=f'{output_folder}/pymol_selections.pml')
else:
    clusters_to_pymol(hac.clusters['cluster'],
                  output=f'{output_folder}/pymol_selections.pml')
            

# if __name__ == "__main__":
#     main()

# cs_file = sys.argv[1]
# correlation_cutoff = sys.argv[2]
# output_folder = sys.argv[3]


#save dendrograms
#save cluster label spreadsheet "dfc"
