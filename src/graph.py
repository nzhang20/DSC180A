'''
graph.py contains functions used to run causal discovery algorithms on the cleaned dataset
'''
import pandas as pd
import numpy as np
import pydot
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import networkx as nx
import graphviz
from causallearn.search.ConstraintBased.PC import pc
from causallearn.utils.cit import kci, fastkci
from causallearn.search.ConstraintBased.FCI import fci
from causallearn.search.ScoreBased.GES import ges
from causallearn.utils.GraphUtils import GraphUtils

from causallearn.utils.cit import CIT
import itertools


def run_pc(data, fp, indep_test):
    '''
    Run PC on the data with the specified independence test and save the plot of the resulting causal graph.

    :param: data: dataframe to run PC on
    :param: fp: filepath name of causal graph result file
    :param: ind_test: string of the choice of independence test available for running PC from the causal-learn package ("fisherz", "chisq", "gsq", "kci", "mv_fisherz") 

    :return: the causallearn object, CausalGraph
    '''
    cg = pc(data.values, alpha=0.05, indep_test=indep_test)

    pyd = GraphUtils.to_pydot(cg.G, labels=data.columns)
    pyd.write_png(f"graphs/{fp}.png")
    return cg


def run_fci(data, fp):
    '''
    Run FCI on the data and save the plot of the resulting causal graph.
    
    :param: data: dataframe to run FCI on
    :param: fp: filepath name of causal graph result file 

    :return: the causallearn object, CausalGraph
    '''
    g, edges = fci(data.values, independence_test_method="fisherz")

    pyd = GraphUtils.to_pydot(g, labels=data.columns)
    pyd.write_png(f"graphs/{fp}.png")
    return g


def run_ges(data, fp):
    '''
    Run GES on the data and save the plot of the resulting causal graph.
    
    :param: data: dataframe to run GES on
    :param: fp: filepath name of causal graph result file 

    :return: the causallearn object, CausalGraph
    '''
    Record = ges(data.values)
    
    pyd = GraphUtils.to_pydot(Record["G"], labels=data.columns)
    pyd.write_png(f"graphs/{fp}.png")
    
    return Record


def fisherz(genus_A, genus_B, genus_C, data):
    '''
    Run the Fisher Z (conditional) independence test between genus_A and genus_B conditioned on genus_c.

    :param: genus_A: index of the first genus to test independence
    :param: genus_B: index of the second genus to test independence
    :param: genus_C: index/list of index of the conditioning set
    :param: data: data containing the genus values
    '''
    fisherz_obj = CIT(data.values, "fisherz")
    pValue = fisherz_obj(genus_A, genus_B, genus_C)
    return pValue


def run_ouralg(adj_matrix, data, genus_lst, indep_test):
    '''
    Run our algorithm on the data and save the plot of the resulting causal graph. The methodology can be found in the associated paper. Returns the resulting adjacency matrix.

    :param: adj_matrix: adjacency matrix of the correlation graph from the paper
    :param: data: data for the corresponding cohort
    :param: genus_lst: list of genus corresponding to the adjacency matrix as node labels
    :param: indep_test: choice of independence test

    :return: the updated adjacency matrix
    '''
    adj_df = pd.DataFrame(adj_matrix, columns=genus_lst, index=genus_lst)

    
    def create_tuples(genus_lst, genus_A, genus_B):
        '''
        Creates tuples of 2 combinations for all items in genus_lst except genus_a and genus_b.
        
        :param: genus_lst: list of genus names
        :param: genus_A: index of the first genus to exclude
        :param: genus_B: index of the second genus to exclude
    
        :return: list of tuples, where each tuple contains two distinct genus names from genus_lst, excluding genus_a and genus_b.
        '''
        filtered_genus_lst = [x for x in range(len(genus_lst)) if (x != genus_A) and (x != genus_B)]
        combinations = list(itertools.combinations(filtered_genus_lst, 2))
        
        return combinations
        
    
    removed_edges = {}

    for i in range(len(genus_lst)):
        for j in range(len(genus_lst)):
            if adj_df.iloc[i, j] == 1.0:
                genus_A = i
                genus_B = j
                for k in range(len(genus_lst)):
                    if k != i and k != j:
                        genus_C = [k]
                        p_value = indep_test(genus_A, genus_B, genus_C, data)
    
                        if p_value < 0.05: # no multiple testing correction
                            adj_df.iloc[i, j] = 0
                            removed_edges[(i, j)] = [genus_C]
                            break
    
    for i in range(len(genus_lst)):
        for j in range(len(genus_lst)):
            if adj_df.iloc[i, j] == 1.0:
                genus_A = i
                genus_B = j
                combinations = create_tuples(genus_lst, genus_A, genus_B)
    
                for k in combinations:
                    genus_C = k
                    p_value = indep_test(genus_A, genus_B, list(genus_C), data)
    
                    if p_value < 0.05:
                        adj_df.iloc[i, j] = 0
                        removed_edges[(i, j)] = [genus_C]
                        break
                    else:
                        adj_df.iloc[i, j] = 2
                        
    return adj_df.values


def graph_networkx(adj_matrix, genus_lst, fp):
    '''
    Create a NetworkX graph of the adjacency matrix and corresponding genus as node labels. Currently uses the circular layout for optimal viewing.

    :param: adj_matrix: adjacency matrix of size 45 x 45 of the genus to genus correlations
    :param: genus_lst: list of genus corresponding to the adjacency matrix as node labels
    '''
    if 1 in np.unique(adj_matrix):
        edge_val = 1
        edge_color = "gray"

    if 2 in np.unique(adj_matrix):
        edge_val = 2
        edge_color = "black"
        
    rows, cols = np.where(adj_matrix == edge_val)
    edges = zip(rows.tolist(), cols.tolist())
    edges = zip(rows.tolist(), cols.tolist())
    edge_list = list(edges)
    graph_labels = {i: genus_lst[i]for i in range(len(genus_lst))}
    
    plt.figure(figsize=(20,15))
    G = nx.Graph()
    G.add_nodes_from(range(len(genus_lst)))
    G.add_edges_from(edge_list)
    pos_spaced = nx.circular_layout(G)
    nx.draw_networkx(G, pos_spaced, node_size=500, font_size=6, labels=graph_labels, with_labels=True, edge_color=edge_color, node_color="#b9d1f0")
    plt.savefig(f"graphs/{fp}.png")


# the following two functions are taken from matplotlib.pyplot's guide/documentation
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current Axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom", fontsize=12)

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels, fontsize=12)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels, fontsize=12)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-60, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="white", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def build_adjacency_matrix_heatmap_of_edges(IR_graph, IS_graph, genera, alg, fp):
    """
    Given the two IR and IS causal graphs, overlay the two and plot their comparative edges. 

    :param: IR_graph: causallearn CausalGraph object for the IR cohort
    :param: IS_graph: causallearn CausalGraph object for the IS cohort
    :param: genera: list of genera in the same order of the columns of the graphs
    :param: alg: string of algorithm type from the following list, ["pc", "fci", "ges", "ouralg"]
    :param: fp: filepath name of resulting heatmap
    """
    if alg == "pc":
        IR_graph.to_nx_graph()
        IS_graph.to_nx_graph()

        IR_adjmat = nx.to_numpy_array(IR_graph.nx_graph)
        IS_adjmat = nx.to_numpy_array(IS_graph.nx_graph)

        tail = 0

    if alg == "fci" or alg == "ges":
        IR_adjmat = IR_graph.graph
        IS_adjmat = IS_graph.graph

        tail = -1

        if alg == "fci": 
            circle = 2

    if alg == "ouralg":
        IR_adjmat = IR_graph
        IS_adjmat = IS_graph

    overlay = np.zeros(IR_adjmat.shape)
    overlay_directed = np.zeros(IR_adjmat.shape)

    for i in range(45):
        for j in range(45):
            IR_ij = IR_adjmat[i, j]
            IR_ji = IR_adjmat[j, i]
            IS_ij = IS_adjmat[i, j]
            IS_ji = IS_adjmat[j, i]

            # our algorithm edges
            if alg == "ouralg":
                if IR_ij == 2 and IS_ij == 2:
                    overlay[i, j] = 3
                elif IR_ij == 0 and IS_ij == 2:
                    overlay[i, j] = 2
                elif IR_ij == 2 and IS_ij == 0:
                    overlay[i, j] = 1
                else:
                    overlay[i, j] = 0 
                    
            # PC, FCI, GES share 2 direction contradiction scenarios
            elif ((IR_ij == 1 and IR_ji == tail and IS_ij == tail and IS_ji == 1) or # directed (L) and directed (R)
                (IR_ij == tail and IR_ji == 1 and IS_ij == 1 and IS_ji == tail) # directed(R) and directed (L)
               ):
                overlay[i, j] = 4
                overlay_directed[i, j] = 1

            # additional FCI direction contradiction scenarios (4)
            elif (alg == "fci") and (
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == tail and IS_ji == 1) or # bidirected and directed (R)
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == 1 and IS_ji == tail) or # bidirected and directed (L)
                  (IR_ij == tail and IR_ji == 1 and IS_ij == 1 and IS_ji == 1) or # directed (R) and bidirected
                  (IR_ij == 1 and IR_ji == tail and IS_ij == 1 and IS_ji == 1) # directed(L) and bidirected
                 ):
                overlay[i, j] = 4
                overlay_directed[i, j] = 3

            # PC, FCI, GES share 2 both scenarios
            elif ((IR_ij == 1 and IR_ji == tail and IS_ij == 1 and IS_ji == tail) or # directed (L) and directed (L)
                  (IR_ij == tail and IR_ji == 1 and IS_ij == tail and IS_ji == 1) # directed (R) and directed(R)
                 ): 
                overlay[i, j] = 3
                overlay_directed[i, j] = 1

            # additional PC both scenarios
            elif (alg == "pc" and 
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == tail and IS_ji == 1) or # undirected and directed (R)
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == 1 and IS_ji == tail) or # undirected and directed (L)
                  (IR_ij == tail and IR_ji == 1 and IS_ij == 1 and IS_ji == 1) or # directed (R) and undirected
                  (IR_ij == 1 and IR_ji == tail and IS_ij == 1 and IS_ji == 1) or # directed (L) and undirected
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == 1 and IS_ji == 1) # undirected and undirected
                 ):
                overlay[i, j] = 3
                overlay_directed[i, j] = 1
                if (IR_ij == 1 and IR_ji == 1 and IS_ij == 1 and IS_ji == 1):
                    overlay_directed[i, j] = 0
                    
            # additional FCI both scenarios with circle or C annotation
            elif (alg == "fci") and (
                  (IR_ij == circle and IR_ji == 1 and IS_ij == tail and IS_ji == 1) or # circle directed (R) and directed (R)
                  (IR_ij == tail and IR_ji == 1 and IS_ij == circle and IS_ji == 1) or # directed (R) and circle directed (R)
                  (IR_ij == circle and IR_ji == 1 and IS_ij == circle and IS_ji == 1) or # circle directed (R) and circle directed (R)
                  (IR_ij == 1 and IR_ji == circle and IS_ij == 1 and IS_ji == tail) or # circle directed (L) and directed (L)
                  (IR_ij == 1 and IR_ji == tail and IS_ij == 1 and IS_ji == circle) or # directed (L) and circle directed (L)
                  (IR_ij == 1 and IR_ji == circle and IS_ij == 1 and IS_ji == circle) or # circle directed (L) and circle directed (L)
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == 1 and IS_ji == 1)
                 ):
                overlay[i, j] = 3
                overlay_directed[i, j] = 2
                if (IR_ij == 1 and IR_ji == 1 and IS_ij == 1 and IS_ji == 1):
                    overlay_directed[i, j] = 3

            # additional FCI both scenarios with circle and C annotation
            elif (alg == "fci") and (
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == circle and IS_ji == 1) or # bidirected and circle directed (R)
                  (IR_ij == 1 and IR_ji == 1 and IS_ij == 1 and IS_ji == circle) or # bidirected and circle directed (L)
                  (IR_ij == circle and IR_ji == 1 and IS_ij == 1 and IS_ji == 1) or # circle directed (R) and bidirected
                  (IR_ij == 1 and IR_ji == circle and IS_ij == 1 and IS_ji == 1)
                 ):
                overlay[i, j] = 3
                overlay_directed[i, j] = 4

            # IS only
            elif (IR_ij == 0 and IR_ji == 0 and IS_ij == 1):
                overlay[i, j] = 2
                if IS_ji == tail:
                    overlay_directed[i, j] = 1
                if alg == "fci":
                    if IS_ji == circle:
                        overlay_directed[i, j] = 1
                    if IS_ji == 1:
                        overlay_directed[i, j] = 3

            # IR only
            elif (IR_ij == 1 and IS_ij == 0 and IS_ji == 0):
                overlay[i, j] = 1
                if IR_ji == tail:
                    overlay_directed[i, j] = 1
                if alg == "fci": 
                    if IR_ji == circle:
                        overlay_directed[i, j] = 1
                    if IR_ji == 1:
                        overlay_directed[i, j] = 3
                    
            # no edges (universal)
            else:
                overlay[i, j] = 0
                overlay_directed[i, j] = 0

    fig, ax = plt.subplots(figsize=(15, 15))

    edge_comp = ["None", "IR only", "IS only", "both", "direction contradiction"]
    norm = matplotlib.colors.BoundaryNorm(np.linspace(-0.5, 4.5, 6), 5)
    fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: edge_comp[norm(x)])
    
    upper_tri_mask = np.tri(45, 45, k=-1).T
    
    cmap = ListedColormap(["#d3d5d3", "#217786", "#e39a33", "#567119", "#ec061d"])

    # np.ma.array(overlay, mask=upper_tri_mask).T
    if alg != "ouralg":
        overlay = overlay.T
        overlay_directed = overlay_directed.T
    im, _ = heatmap(overlay, genera, genera, ax=ax, 
                    cmap=cmap, norm=norm, 
                    cbar_kw=dict(ticks=np.arange(0, 5), format=fmt, fraction=0.046, pad=0.04), cbarlabel="Edge Comparison")
    
    def func(x, pos):
        if x == 1:
            return "▶"
        elif x == 2:
            return "●"
        elif x == 3:
            return "C"
        elif x == 4:
            return "C●"
        else:
            return ""

    # np.ma.array(overlay_directed, mask=upper_tri_mask).T
    annotate_heatmap(im, data=overlay_directed, valfmt=matplotlib.ticker.FuncFormatter(func), threshold=0)
    
    plt.tight_layout()
    plt.savefig(f"graphs/{fp}.png")
    plt.show()

    return


def cohort_all_alg_heatmap(cohort_list, genera, fp):
    """
    Given a list of the adjacency matrices for the cohort (IR or IS), add all of the entries to produce a heatmap showing the frequency of each edge represented in an algorithm.

    :param: cohort_list: list of adjacency matrices for a single cohort
    :param: genera: list of genera in the same order of the columns of the graphs
    :param: fp: filepath name of resulting heatmap
    """
    fig, ax = plt.subplots(figsize=(15, 15))

    edge_comp = np.arange(0, 4)
    norm = matplotlib.colors.BoundaryNorm(np.linspace(-0.5, 3.5, 5), 4)
    fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: edge_comp[norm(x)])
    
    cmap = ListedColormap(["#ede5e0", "#edca85", "#b9d1f0", "#3a4a68"])
    
    im, _ = heatmap(np.sum(cohort_list, 0), genera, genera, ax=ax, cmap=cmap, norm=norm, 
                        cbar_kw=dict(ticks=np.arange(0, 4), format=fmt, fraction=0.046, pad=0.04), cbarlabel="Frequency of Edge")
    
    plt.tight_layout()
    plt.savefig(f"graphs/{fp}.png")