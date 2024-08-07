import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.cluster.hierarchy import dendrogram

from matplotlib.colors import Normalize, LinearSegmentedColormap, ListedColormap
import warnings
warnings.filterwarnings("ignore")

def plot_corr(HAC, cutoff=None, save_file=None):
    sns.set(font_scale=2)
    if cutoff == None:
        corr = HAC.get_corr_above(HAC.cutoff)
    else:
        corr = HAC.get_corr_above(cutoff)
    fig, ax = plt.subplots(figsize=(25,20))
    sns.heatmap(corr,cbar=False, ax=ax)
    ax.invert_yaxis()
    if save_file is not None:
        fig.savefig(save_file)

def get_cluster_annotation_positions(clusters, threshold):
    '''
    find the x and y values to add cluster labels to the dendrogram.
    
    Parameters
    ----------
    clusters : pd.DataFrame
        From HAC.clusters. 'cluster' column with id that the resid in index
        is assigned to.

    threshold : int or float
        The threshold used for coloring the dendrogram

    Returns
    -------
    x : list
    list of x positions at the center of each cluster
    y : list
    list of y positions where the cluster annotation will appear.
    This is just the threshold value repeated for each x.
    
    '''
    
    # clusters go from left to right in dendrogram

    # track the farthest right x value because it will serve as the first point
    # for the next cluster position 
    r_pos = 0
    xs = []
    ys = []
    for cluster in range(1,clusters['cluster'].max()+1):
        # indices of residues that fall in same cluster
        c = np.where((clusters['cluster']==cluster))[0] # [0] to unpack tuple
        # size of cluster * 10 tells you how wide it is in the dendrogram
        # (assuming it's always a multiple of 10...)
        # TODO: can check against dn['icoords'].max()/len(clusters.index)-1
        # fyi: icoords are x1a,x2a,x1b,x2b for the two vertical line endpoints in dcoords.
        width = c.shape[0]*10
        if cluster == 1:
            midpoint = width/2
            xs.append(midpoint)
            ys.append(threshold)
            r_pos = width
        else:
            midpoint = (width/2)+r_pos
            xs.append(midpoint)
            ys.append(threshold)
            r_pos += width
    return xs, ys









def show_dendrogram(HAC,
                     save_file=None,
                     orientation='right',
                     leaf_rotation=None,
                     cutoff_line=False,
                     annotate_clusters=True,
                     sub_cluster=None):
    if orientation == 'top':
        fig, ax = plt.subplots(figsize=(8,5)) # 40, 15
    elif orientation == 'right':
        fig, ax = plt.subplots(figsize=(15,40))
    # dn is available for plt
    #dn = HAC.get_dendrogram(method, metric ,ax=ax)
    threshold = 100-HAC.cutoff
    dn = dendrogram(HAC.linkage, 
                          color_threshold=threshold, 
                          labels=HAC.df.index,
                          ax=ax,
                          orientation=orientation,
                          leaf_rotation=leaf_rotation,
                          leaf_font_size=14
                          )
    clusters = HAC.clusters
    
    # Annotate the first occurrence of each cluster label
    xs, ys = get_cluster_annotation_positions(clusters, threshold)
    if orientation == 'top':
        pass
    elif orientation == 'right':
        xs, ys = ys, xs
    if annotate_clusters == True:
        for i, (x, y) in enumerate(zip(xs,ys)):
            # TODO: option to rotate these ids by 90 so you can view the dendrogram
            # enlarged somewhere else and rotated to more easily see which residues
            # fall in a cluster -- or just orient the dn vertically with sp option
            # BUUUUT - then annotation xs and ys need to be flipped.....
            ax.text(x,y,f'{i+1}')
    if sub_cluster is not None:
        ax.set_title(f"Cluster {sub_cluster} States")
    else:
        ax.set_title('CHESCA Clusters')
    ax.grid(visible=False)
    ax.set_facecolor('white')
    # if cutoff_line == True:
    #     ax.axhline(y=threshold, linestyle='dashed')
    ax.xaxis.set_tick_params(labelsize=15)
    #fig.show()
    if save_file is not None:
        fig.savefig(save_file)



def plot_svd(SVD, centering='column', save_file=None):
    '''
    Plots the first two components from the SVD, the residues
    and the states all on the same axis system.

    TODO: add options for all the figure elements/colors etc

    Parameters
    ----------
    SVD : pychesca.svd.SVD
        The SVD object 
    
    centering : str
        valid options are 'column' for column centered 
        or 'row' for row centered

    TODO : add additional PC projections
    '''
    if centering == 'column':
        data = SVD.column_svd
    elif centering == 'row':
        data = SVD.row_svd
    else:
        raise ValueError("centering must be 'column' or 'row'")
    
    try:
        uds = data['uds']
        v = data['V']
    except KeyError as e:
        raise KeyError(f"Missing expected data in SVD results: {e}")
    
    uds = data['uds']
    v = data['V']
    x_min, x_max = uds[:,0].min(), uds[:,0].max()
    y_min, y_max = uds[:,1].min(), uds[:,1].max()
    fig, ax = plt.subplots()
    ax.scatter(uds[:,0],uds[:,1],marker='o',facecolors='none', color='black')
    
    ax.vlines(0,y_min,y_max, color='blue')
    ax.hlines(0,x_min,x_max, color='red')
    ax.scatter(v[:,0],v[:,1], marker="D", color='magenta')
    # might need to get some arrows involved and move the letters
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.annotate.html
    # https://matplotlib.org/stable/gallery/text_labels_and_annotations/text_commands.html#sphx-glr-gallery-text-labels-and-annotations-text-commands-py
    for i, txt in enumerate(list(SVD.cols)):
        ax.annotate(txt, (v[i,0], v[i,1]), color='magenta')
    if save_file is not None:
        fig.savefig(save_file)



def heatmap_correlation_cutoffs(df, min_corr=94, save_file='None'):
    # Protein NMR chesca chapter 18 p. 403 Note 5
    # suggests making heatmaps of correlation coefficent cutoffs

    corr = df.T.corr().abs().fillna(0)
    if min_corr > 1:
        min_corr = min_corr/100

    # Mask out correlations below the minimum and diagonal
    mask = np.triu(np.ones_like(corr, dtype=bool)) | (np.abs(corr) < min_corr)
    
    fig, ax = plt.subplots(figsize=(25,20))

    
    
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, mask=mask, cmap='plasma',
                square=True, linewidths=.5, ax=ax)

    ax.set_title(f'Heatmap of Correlation Coefficients Above {min_corr}')

    if save_file is not None:
        fig.savefig(save_file)



def plot_chespa(Chespa):

    chsp = Chespa
    # symbols for plot titles
    D = '\u0394'
    d = '\u03B4'
    theta = '\u03B8'

    fig, axs = plt.subplots(3,1,figsize=(40,20))
    resis = [str(i) for i in chsp.resis]
    for i, ax in enumerate(axs):
        if i == 0:
            ax.bar(resis,chsp.ref_to_A)
            ax.grid(visible=False)
            ax.set_title(f'{D}{d}{chsp.het_nuc.upper()}Hcomb(ppm)')
            ax.set_xticklabels(resis, fontdict={'size':12}, rotation=90);
        elif i == 1:
            ax.bar(resis, chsp.cos_theta)
            ax.grid(visible=False)
            ax.set_title(f'cos{theta}')
            ax.set_xticklabels(resis, fontdict={'size':12}, rotation=90);
        else:
            ax.bar(resis, chsp.X)
            ax.grid(visible=False)
            ax.set_title(f'X (fractional activation)')
            ax.set_xticklabels(resis, fontdict={'size':12}, rotation=90);
    fig.tight_layout()


def plot_everything():
    '''
    Have all of the plots shown in one set of plt.subplots with dendrogram
    oriented vertically
    '''
    pass
