import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def plot_corr(corr, output):
    sns.set(font_scale=2)
    fig, ax = plt.subplots(figsize=(25,20))
    sns.heatmap(corr,cbar=False)
    fig.savefig(output)

def show_dendrogram(distance_matrix,cutoff,save_file=None):

    fig, ax = plt.subplots(figsize=(30,10))
    # dn is available for plt.show
    dn = dendrogram(distance_matrix, color_threshold=100-cutoff,ax=ax)
    ax.set_title('CHESCA Clusters')
    ax.grid(visible=False)
    ax.set_facecolor('white')
    ax.xaxis.set_tick_params(labelsize=15)
    if save_file is not None:
        fig.savefig(save_file)