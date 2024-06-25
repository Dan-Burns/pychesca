import numpy as np
import pandas as pd 
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
import seaborn as sns
import scipy
import os
import sys
import argparse
import matplotlib.pyplot as plt
# Was the sandbox script
# this will end up importing the finished functions and run the whole analysis in one go


# Create the parser
parser = argparse.ArgumentParser(description='Load a csv into pandas')


parser.add_argument('-file', type=str, help='Path to the spreadsheet')
parser.add_argument('-output', type=str, help='Output folder path')
parser.add_argument('-cutoff', type=float, default=98.0, help='Correlation cutoff (%) for clustering')
# Parse the command line arguments
args = parser.parse_args()

# Read the Excel file into a pandas DataFrame
try:
    df = pd.read_csv(args.file,index_col='Resi')
except Exception as e:
    print(f"Failed to read the spreadsheet file. Case Sensitive. Error: {e}")

output_folder = args.output
os.makedirs(output_folder)
cutoff = args.cutoff
print(cutoff)


# make the correlation plot
sns.set(font_scale=2)
fig, ax = plt.subplots(figsize=(25,20))
corr = df.T.corr().abs() > cutoff/100
sns.heatmap(corr,cbar=False,ax=ax)
ax.invert_yaxis()
ax.figure.savefig(f'{output_folder}/correlation_matrix.pdf')

# clustering
# define the distance as 1 - correlation (0 was 1 and means it's so close it's on top of each other)
cor_dist = 1-df.T.corr().abs().fillna(0) # fills the Na's with 0's before 1 - correlation so they turn to 1 after (far distance)
cutoff = args.cutoff
threshold = 100-cutoff
# https://predictivehacks.com/hierarchical-clustering-in-python/
# Import the fcluster and linkage functions
data = cor_dist
# Use the linkage() function
distance_matrix = linkage(data, method = 'complete', metric = 'euclidean')

dfc = pd.DataFrame(index=data.index)
# Assign cluster labels (2 is 98%)
dfc['cluster_labels'] = fcluster(distance_matrix, threshold, criterion='distance')
dfc.to_csv(f'{output_folder}/cluster_labels_{cutoff}.csv')
# make dendrogram figure
# to add cluster labels refer (possibly) to
# https://stackoverflow.com/questions/35873273/display-cluster-labels-for-a-scipy-dendrogram
fig, ax = plt.subplots(figsize=(30,10))
dn = dendrogram(distance_matrix, color_threshold=threshold,ax=ax)
ax.set_title('test')
ax.grid(visible=False)
ax.set_facecolor('white')
ax.xaxis.set_tick_params(labelsize=15)
# Display the dendogram
fig.savefig(f'{output_folder}/dendrogram_{cutoff}.pdf')

# SVD
########### WORKING ON THIS ################
# calculate the row mean ref to use as reference state
# wild type reference will probably be better
cols = df.columns
row_mean = df.mean(axis=1)
# CCS difference from mean 
row_centered = df[cols].subtract(row_mean.values, axis=0)
# column center the data
column_mean = row_centered.mean(axis=0)
column_centered = row_centered[cols].subtract(column_mean.values, axis=1)


# SVD
u, s, v = scipy.linalg.svd(column_centered)
# make diagonal s matrix of correct shape
sdiag = scipy.linalg.diagsvd(s,*column_centered.shape)
# u dot s is score matrix of each res on each pc
uds = u@sdiag
# explained_variance
s_variance = (s**2/(sum(s**2)))
# V is loading score of each state
vdf = pd.DataFrame(Vh, columns=[f'PC{i+1}' for i in range(len(df.columns))], index=df.columns)

fig, ax = plt.subplots()
ax.scatter(uds[:,0],uds[:,1])
ax.scatter(v[:,0],v[:,1], marker="D")
ax.vlines(0,1.5,-1.5, color='r')
ax.hlines(0,1.5,-1.5, color='black')
for i, txt in enumerate(list(df.columns)):
    ax.annotate(txt, (v[i,0], v[i,1]))
#ax.figure.savefig('/Users/dburns/Desktop/unfiltered_SVD_96.png')

#### Filtering out residues that don't fall in large clusters before doing SVD again #######
# Filter out residues that do not fall into large groups 
above_cutoff = dfc['cluster_labels'].value_counts()>9
clusters_to_keep = above_cutoff.index[above_cutoff].values
resis_to_keep = dfc.index[dfc['cluster_labels'].isin(clusters_to_keep)]
filtered = column_centered.loc[resis_to_keep]


########################################################################










#Pymol output
def clusters_to_pymol(df,output):
    '''
    take a df of residue indices and cluster labels
    and make pymol selections
    '''
    colors = ['red','blue','yellow','black','white','purpleblue','magenta','orange','neptunium','zinc','radium','pink',
             'purple','lighblue','teal','marine','brightorange','chocolate','sand','slate']
    
    clusters = {}
    # Get a set of the cluster ids {1,2,3,...}
    for cluster_label in set(df['cluster_labels']):
        # Add the residue ids to the dictionary that correspond to the cluster id
        clusters[int(cluster_label)] = list(df.loc[df['cluster_labels'] == cluster_label].index)
    
    with open(output,'w') as f:
        f.write('set sphere_scale, 0.8\n')
        for cluster in clusters:
            sel = ''
            for resid in clusters[cluster]:
                # int(float) to catch Aayushi's residue decimal naming scheme to differentiate CH/NH/ILE methyls etc.
                sel += f'{int(float(resid))}+'
                #sel += f'{(resid)}+'

            f.write(f'select c_{cluster}, resi {sel}\n') 
            f.write(f'color {colors[cluster-1]}, c_{cluster}\n show spheres, c_{cluster}\n')

pymol_output = f'{output_folder}/pymol_{cutoff}_clusters.pml'
clusters_to_pymol(dfc,pymol_output)
            

# if __name__ == "__main__":
#     main()

# cs_file = sys.argv[1]
# correlation_cutoff = sys.argv[2]
# output_folder = sys.argv[3]


#save dendrograms
#save cluster label spreadsheet "dfc"
