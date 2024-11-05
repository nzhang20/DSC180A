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

from src.etl import merge_gut_subject, get_bacteria_and_covariates
from src.eda import make_corr_plot
from src.graph import run_pc, run_fci, run_ges


def main(targets):
    if "data" in targets:
        with open("config/data-params.json") as fh:
            data_params = json.load(fh)
        get_bacteria_and_covariates(merge_gut_subject(), **data_params)

    if "eda" in targets:
        # make corr plots, scatter plots, qq plots (NOT FINISHED)
        # make_corr_plot(data, "plots/corr_plot")

    
    if "graph" in targets:
        data = pd.read_csv("data/clean.csv")
        IR = data[data["IR_IS_classification" == "IR"]]
        IS = data[data["IR_IS_classification" == "IS"]]
        
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