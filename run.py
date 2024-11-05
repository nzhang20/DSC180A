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
        get_bacteria_and_covariates(merge_gut_subject(), **data_params)

    if "eda" in targets:
        # subject info dataset
        subject_df = get_subject_data()
        # merged dataset
        merge_df = merge_gut_subject()

        check_discrete_distribution(subject_df)
        check_gaussian_distribution(subject_df)
        check_linearity(subject_df)

        # clean merged data and generate correlation plots
        X_clean = clean_merged_df(merge_df)
        make_corr_plot(X_clean, "plots/merged_corr_plot")
        
        # optionally check correlation differences between IR and IS groups
        corr_IR_IS(merge_df)

        # run check_corr_plot with IR/IS group analysis if specified
        check_corr_plot(merge_df, IR_IS=True)
        
        # check distributions of discrete variables in merged data
        check_merged_discrete(merge_df)

    
    if "graph" in targets:
        data = merge_gut_subject()
        IR, IS = IR_IS_classify(data)
       
        # PC is slow with KCI, waiting on Fast KCI
        # run_pc(data, "graphs/pc_causal_graph", indep_test="kci")
        # run_pc(IR, "graph/pc_IR_causal_graph", indep_test="kci")
        # run_pc(IS, "graph/pc_IS_causal_graph", indep_test="kci")
        
        run_fci(data, "graphs/fci_causal_graph")
        run_fci(IR, "graphs/fci_IR_causal_graph")
        run_fci(IS, "graphs/fci_IS_causal_graph")
        
        run_ges(data, "graphs/ges_causal_graph")
        run_ges(IR, "graphs/ges_IR_causal_graph")
        run_ges(IS, "graphs/ges_IS_causal_graph")
        
    if "clean" in targets:
        for f in glob.glob("data/*.csv"):
            os.remove(f)
        for f in glob.glob("plots/*.png"):
            os.remove(f)
        for f in glob.glob("graphs/"):
            os.remove(f)
        
    

if __name__ == '__main__':
    targets = sys.argv[1:]
    if "all" in targets:
        main(["data", "eda", "graph"])
    else:
        main(targets)