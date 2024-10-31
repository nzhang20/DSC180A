#!/usr/bin/env python

import sys
import json

from src.etl import merge_gut_subject, get_genus_and_covariates
from src.graph import make_corr_plot, run_ges


def main(targets):
    if "data" in targets:
        with open("config/data-params.json") as fh:
            data_params = json.load(fh)
        get_genus_and_covariates(merge_gut_subject(), **data_params)
        
    if "graph" in targets:
        data = pd.read_csv("data/clean.csv")
        make_corr_plot(data, "plots/corr_plot")
        run_ges(data, "graphs/ges_causal_graph")
        
    

if __name__ == '__main__':
    targets = sys.argv[1:]
    main(targets)