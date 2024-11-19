'''
graph.py contains functions used to run causal discovery algorithms on the cleaned dataset
'''
import pandas as pd
import numpy as np
import pydot
import matplotlib.pyplot as plt
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
    '''
    cg = pc(data.values, alpha=0.05, indep_test=indep_test)

    pyd = GraphUtils.to_pydot(cg.G, labels=data.columns)
    pyd.write_png(f"graphs/{fp}.png")


def run_fci(data, fp):
    '''
    Run FCI on the data and save the plot of the resulting causal graph.
    
    :param: data: dataframe to run FCI on
    :param: fp: filepath name of causal graph result file 
    '''
    g, edges = fci(data.values, independence_test_method="fisherz")

    pyd = GraphUtils.to_pydot(g, labels=data.columns)
    pyd.write_png(f"graphs/{fp}.png")


def run_ges(data, fp):
    '''
    Run GES on the data and save the plot of the resulting causal graph.
    
    :param: data: dataframe to run GES on
    :param: fp: filepath name of causal graph result file 
    '''
    Record = ges(data.values)
    
    pyd = GraphUtils.to_pydot(Record["G"], labels=data.columns)
    pyd.write_png(f"graphs/{fp}.png")
    
    return


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
    