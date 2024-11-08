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
       
        # PC is slow with KCI, waiting on Fast KCI
        # run_pc(data, "graphs/pc_causal_graph", indep_test="fastkci")
        # run_pc(IR, "graph/pc_IR_causal_graph", indep_test="fastkci")
        # run_pc(IS, "graph/pc_IS_causal_graph", indep_test="fastkci")
        
        run_fci(data, "fci_causal_graph")
        run_fci(IR, "fci_IR_causal_graph")
        run_fci(IS, "fci_IS_causal_graph")

        # ges slow
        # run_ges(data, "graphs/ges_causal_graph")
        # run_ges(IR, "graphs/ges_IR_causal_graph")
        # run_ges(IS, "graphs/ges_IS_causal_graph")
        
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