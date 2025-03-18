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

def get_cluster_annotation_positions(clusters, threshold, ax, orientation):
    '''
    find the x and y values to add cluster labels to the dendrogram.

    
    Parameters
    ----------
    clusters : pd.DataFrame
        From HAC.clusters. 'cluster' column with id that the resid in index
        is assigned to.

    threshold : int or float
        The threshold used for coloring the dendrogram

    ax : matplotlib.axes.Axes
        The axes object of the dendrogram plot

    Returns
    -------
    x : list
    list of x positions at the center of each cluster
    y : list
    list of y positions where the cluster annotation will appear.
    This is just the threshold value repeated for each x.
    
    '''
    ordered_indices = []
    x_positions = []
    if orientation == 'top':
        for t in ax.get_xticklabels():
            ordered_indices.append(float(t.get_text())) # TODO: the df indices need to be strings so this will work for clustering states
            x_positions.append(t.get_position())
        x_positions = np.array(x_positions)
    else:
        for t in ax.get_yticklabels():
            ordered_indices.append(float(t.get_text()))
            x_positions.append(t.get_position())
        x_positions = np.array(x_positions)

    xs = []
    for cluster in range(clusters['cluster'].max() + 1):
        # index the x_positions by first putting the clusters df in the same order
        # as the dendrogram and then using the indices to get the x positions for
        # each residue in the cluster (taken from ax.get_xticklabels)
        # and then us the midpoint of these values 
        # the list of values will be in the order that the clusters appear in the dendrogram from
        # left to right
        if orientation == 'top':
            vals = x_positions[np.where(clusters.loc[ordered_indices]['cluster']==cluster)][:,0]
        else:
            vals = x_positions[np.where(clusters.loc[ordered_indices]['cluster']==cluster)][:,1]
        xs.append(np.median(vals))

    ys = [threshold] * len(xs)

    return xs, ys

def plot_dendrogram(model, **kwargs):
    # https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
    # need this to plot the sklearn version
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)



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
    kwargs = {'color_threshold':threshold,
              'labels':HAC.corr_distance.index,
              'ax':ax,
              'orientation':orientation,
              'leaf_rotation':leaf_rotation,
              'leaf_font_size':14
            }
    dn = plot_dendrogram(HAC.linkage, 
                          **kwargs
                          )
    clusters = HAC.clusters
    
    # Annotate the first occurrence of each cluster label
    # treats the label as being on the x axis and then flips it after if it's oriented "right"
    xs, ys = get_cluster_annotation_positions(clusters, threshold, ax, orientation)
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
