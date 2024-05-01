import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from matplotlib.colors import Normalize, LinearSegmentedColormap, ListedColormap


def plot_corr(HAC, cutoff=98, save_file=None):
    sns.set(font_scale=2)
    if cutoff == None:
        corr = HAC.absolute_corr
    else:
        corr = HAC.get_corr_above(cutoff)
    fig, ax = plt.subplots(figsize=(25,20))
    sns.heatmap(corr,cbar=False, ax=ax)
    ax.invert_yaxis()
    if save_file is not None:
        fig.savefig(save_file)

def show_dendrogram(HAC, cutoff=98,
                     method='complete', metric='euclidean',
                     save_file=None):

    fig, ax = plt.subplots(figsize=(30,10))
    # dn is available for plt
    dn = HAC.get_dendrogram(cutoff, method, metric ,ax=ax)
    ax.set_title('CHESCA Clusters')
    ax.grid(visible=False)
    ax.set_facecolor('white')
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