#!/usr/bin/env python

import sys
import json

# necessary imports for functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pylab as py
import itertools

from src.etl import merge_gut_subject, get_bacteria_and_covariates
from src.graph import make_corr_plot, run_ges, run_fci


def main(targets):
    if "data" in targets:
        with open("config/data-params.json") as fh:
            data_params = json.load(fh)
        get_bacteria_and_covariates(merge_gut_subject(), **data_params)

    if "eda" in targets:
        # make corr plots, scatter plots, qq plots

    
    if "graph" in targets:
        data = pd.read_csv("data/clean.csv")
        make_corr_plot(data, "plots/corr_plot")
        run_ges(data, "graphs/ges_causal_graph")
        run_fci(data, "graphs/fci_causal_graph")
        
    

if __name__ == '__main__':
    targets = sys.argv[1:]
    main(targets)