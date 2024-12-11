import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


def plot_umap(adata, embedding, cmap, use_rep, legend_list,
              figsize=(5, 5), markersize=8, fontsize=8, legendtitlesize=10,
              legend_title=None, loc='upper left', titlesize=12, title=None, save_path=None, save=False, show=True):

    n_spots = adata.shape[0]
    size = 10000 / n_spots
    order = np.arange(n_spots)

    colors = adata.obs[use_rep].astype('str').map(cmap)

    plt.figure(figsize=figsize)
    plt.scatter(embedding[order, 0], embedding[order, 1], s=size, c=colors)
    plt.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                    labelleft=False, labelbottom=False, grid_alpha=0)

    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markersize=markersize,
               markerfacecolor=cmap[str(i)], label=str(i)) for i in legend_list
    ]
    plt.legend(handles=legend_handles, fontsize=fontsize, title=legend_title,
               title_fontsize=legendtitlesize, loc=loc)
    if title:
        plt.title(title, fontsize=titlesize)

    plt.tight_layout()

    if save:
        plt.savefig(save_path)
    if show:
        plt.show()


def plot_gene_expression(adata, coord_rep, gene_name, log_norm=True, y_inverse=False,
                         cmap='viridis', figsize=(6, 4), s=100,
                         titlesize=12, norm0to1=True, title=None, save_path=None, save=False, show=True):

    adata_copy = adata.copy()

    if log_norm:
        sc.pp.normalize_total(adata_copy, target_sum=1e4)
        sc.pp.log1p(adata_copy)

    plt.subplots(figsize=figsize)

    gene_expression = adata_copy[:, gene_name].X.toarray()

    if norm0to1:
        gene_expression = (gene_expression - np.min(gene_expression)) / (np.max(gene_expression) - np.min(gene_expression))

    plt.scatter(adata_copy.obsm[coord_rep][:, 0], adata_copy.obsm[coord_rep][:, 1], linewidth=0, s=s, marker=".",
                c=gene_expression, cmap=cmap)
    if y_inverse:
        plt.gca().invert_yaxis()
    plt.axis('off')

    plt.colorbar(label='Color')

    if title:
        plt.title(title, fontsize=titlesize)

    plt.tight_layout()

    if save:
        plt.savefig(save_path)
    if show:
        plt.show()

