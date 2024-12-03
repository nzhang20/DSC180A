#!/usr/bin/env python

import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pylab as py
import itertools
import glob, os
from causallearn.search.ConstraintBased.FCI import fci
from causallearn.search.ConstraintBased.PC import pc
from causallearn.utils.cit import kci, fastkci
from causallearn.utils.GraphUtils import GraphUtils
from causallearn.utils.PCUtils.BackgroundKnowledge import BackgroundKnowledge
from causallearn.graph.GraphNode import GraphNode
import pydot
from IPython.display import Image, display


from src.etl import *
from src.eda import *
from src.graph import *


def main(targets):
    if "data" in targets:
        with open("config/data-params.json") as fh:
            data_params = json.load(fh)
        get_bacteria_and_covariates(merge_gut_sample_subject(), **data_params)

    if "eda" in targets:
        # subject data eda
        subject_df = get_subject_data()

        check_discrete_distribution(subject_df)
        check_gaussian_distribution(subject_df)
        check_linearity(subject_df)

        # clean dataset
        data = pd.read_csv("data/clean.csv")

        # numerically encode the clean dataset for correlation plots
        clean_num = numerical_encoding(data)
        make_corr_plot(corr(clean_num), "full_data_corr_plot")
        
        # optionally check correlation differences between IR and IS groups
        IR_clean_num, IS_clean_num = IR_IS_split(clean_num, True)
        make_corr_plot(corr(IR_clean_num), "IR_data_corr_plot")
        make_corr_plot(corr(IS_clean_num), "IS_data_corr_plot")
        make_corr_plot(corr(IR_clean_num) - corr(IS_clean_num), "IR_IS_diff_corr_plot", "Difference-in-Correlation Matrix")
        
        # check linearity of all variables (includes bacteria, numerical covariates)
        check_all = input(f"Would you like to generate {len(list(itertools.combinations(clean_num.columns, 2)))} scatter plots to check for linearity in {len(clean_num.columns)} variables? (This may take a while.) \n(Y/[N]): ")
        if check_all.lower() == 'y':
            check_linearity_large(clean_num)

    
    if "graph" in targets:
        data = pd.read_csv("data/clean.csv")
        data = numerical_encoding(data)
        IR, IS = IR_IS_split(data, True)
        non_genera = ["Ethnicity", "Gender", "Adj.age", "BMI", "SSPG"]

        # PC with genera only 
        # run_pc(data.drop(columns=non_genera + ["IRIS"]), "pc_genera_causal_graph", indep_test=fastkci)
        try: 
            IR_pc = run_pc(IR.drop(columns=non_genera), "pc_IR_genera_causal_graph", indep_test=fastkci)
            IS_pc = run_pc(IS.drop(columns=non_genera), "pc_IS_genera_causal_graph", indep_test=fastkci)
        except ValueError:
            IR_pc = run_pc(IR.drop(columns=non_genera), "pc_IR_genera_causal_graph", indep_test=kci)
            IS_pc = run_pc(IS.drop(columns=non_genera), "pc_IS_genera_causal_graph", indep_test=kci)

        # FCI with genera only
        # run_fci(data.drop(columns=non_genera + ["IRIS"]), "fci_genera_causal_graph")
        IR_fci = run_fci(IR.drop(columns=non_genera), "fci_IR_genera_causal_graph")
        IS_fci = run_fci(IS.drop(columns=non_genera), "fci_IS_genera_causal_graph")

        # GES with genera only
        # run_ges(data.drop(columns=non_genera + ["IRIS"]), "ges_genera_causal_graph")
        IR_ges = run_ges(IR.drop(columns=non_genera), "ges_IR_genera_causal_graph")
        IS_ges = run_ges(IS.drop(columns=non_genera), "ges_IS_genera_causal_graph")

        # GES with genera + covariates
        # run_ges(data, "graphs/ges_causal_graph")
        # run_ges(IR, "graphs/ges_IR_causal_graph")
        # run_ges(IS, "graphs/ges_IS_causal_graph")
        
        # FCI with genera + covariates
        # run_fci(data, "fci_causal_graph")
        # run_fci(IR, "fci_IR_causal_graph")
        # run_fci(IS, "fci_IS_causal_graph")

        # our algorithm
        with open("config/data-params.json") as fh:
            data_params = json.load(fh)

        genus_lst = data_params["genus"]

        sparcc_data = pd.read_csv("data/raw/sparcc.csv")
        rmcorr_data = pd.read_csv("data/raw/rmcorr.csv")
        matrix_IR_sparcc, matrix_IS_sparcc = create_adj_matrix(sparcc_data, "sparcc", genus_lst)
        matrix_IR_rmcorr, matrix_IS_rmcorr = create_adj_matrix(rmcorr_data, "rmcorr", genus_lst)

        graph_networkx(matrix_IR_sparcc, genus_lst, "IR_sparcc_graph")
        graph_networkx(matrix_IS_sparcc, genus_lst, "IS_sparcc_graph")
        graph_networkx(matrix_IR_rmcorr, genus_lst, "IR_rmcorr_graph")
        graph_networkx(matrix_IS_rmcorr, genus_lst, "IS_rmcorr_graph")

        IR_ouralg_sparcc = run_ouralg(matrix_IR_sparcc, IR.drop(columns=non_genera), genus_lst, fisherz)
        IS_ouralg_sparcc = run_ouralg(matrix_IS_sparcc, IS.drop(columns=non_genera), genus_lst, fisherz)
        IR_ouralg_rmcorr = run_ouralg(matrix_IR_rmcorr, IR.drop(columns=non_genera), genus_lst, fisherz)
        IS_ouralg_rmcorr = run_ouralg(matrix_IS_rmcorr, IR.drop(columns=non_genera), genus_lst, fisherz)
        
        graph_networkx(IR_ouralg_sparcc, genus_lst, "ouralg_IR_sparcc_causal_graph")
        graph_networkx(IS_ouralg_sparcc, genus_lst, "ouralg_IS_sparcc_causal_graph")
        graph_networkx(IR_ouralg_rmcorr, genus_lst, "ouralg_IR_rmcorr_causal_graph")
        graph_networkx(IS_ouralg_rmcorr, genus_lst, "ouralg_IS_rmcorr_causal_graph")

        # adjacency matrix heatmaps
        build_adjacency_matrix_heatmap_of_edges(IR_pc, IS_pc, genus_lst, "pc", "pc_heatmap")
        build_adjacency_matrix_heatmap_of_edges(IR_fci, IS_fci, genus_lst, "fci", "fci_heatmap")
        build_adjacency_matrix_heatmap_of_edges(IR_ges['G'], IS_ges['G'], genus_lst, "ges", "ges_heatmap")
        build_adjacency_matrix_heatmap_of_edges(IR_ouralg_sparcc, IS_ouralg_sparcc, genus_lst, "ouralg", "ouralg_sparcc_heatmap")
        build_adjacency_matrix_heatmap_of_edges(IR_ouralg_rmcorr, IS_ouralg_rmcorr, genus_lst, "ouralg", "ouralg_rmcorr_heatmap")

        # clean up adjacency matrices for the combined frequency matrix
        IR_sparcc_graph = np.copy(IR_ouralg_sparcc)
        IR_sparcc_graph = np.vectorize({0:0, 2:1}.get)(IR_sparcc_graph)
        IR_sparcc_graph += np.tril(IR_sparcc_graph.T, -1)
        
        IS_sparcc_graph = np.copy(IS_ouralg_sparcc)
        IS_sparcc_graph = np.vectorize({0:0, 2:1}.get)(IS_sparcc_graph)
        IS_sparcc_graph += np.tril(IS_sparcc_graph.T, -1)
        
        IR_rmcorr_graph = np.copy(IR_ouralg_rmcorr)
        IR_rmcorr_graph = np.vectorize({0:0, 2:1}.get)(IR_rmcorr_graph)
        IR_rmcorr_graph += np.tril(IR_rmcorr_graph.T, -1)
        
        IS_rmcorr_graph = np.copy(IS_ouralg_rmcorr)
        IS_rmcorr_graph = np.vectorize({0:0, 2:1}.get)(IS_rmcorr_graph)
        IS_rmcorr_graph += np.tril(IS_rmcorr_graph.T, -1)
        
        IR_adj_matrices = [nx.to_numpy_array(IR_pc.nx_graph), 
                           np.vectorize({0:0, -1:1, 1:1, 2:1}.get)(IR_fci.graph), 
                           np.vectorize({0:0, -1:1, 1:1}.get)(IR_ges['G'].graph), 
                           IR_sparcc_graph, 
                           IR_rmcorr_graph]
        
        IS_adj_matrices = [nx.to_numpy_array(IS_pc.nx_graph), 
                           np.vectorize({0:0, -1:1, 1:1, 2:1}.get)(IS_fci.graph), 
                           np.vectorize({0:0, -1:1, 1:1}.get)(IS_ges['G'].graph), 
                           IS_sparcc_graph, 
                           IS_rmcorr_graph]

        cohort_all_alg_heatmap(IR_adj_matrices, genus_lst, "IR_combined_heatmap")
        cohort_all_alg_heatmap(IS_adj_matrices, genus_lst, "IS_combined_heatmap")

        # GraphicalLASSO + PC
        
    if "clean" in targets:
        for f in glob.glob("data/*.csv"):
            os.remove(f)
        for f in glob.glob("plots/*.png"):
            os.remove(f)
        for f in glob.glob("graphs/*.png"):
            os.remove(f)
        
    
if __name__ == '__main__':
    targets = sys.argv[1:]
    if "all" in targets:
        main(["data", "eda", "graph"])
    else:
        main(targets)