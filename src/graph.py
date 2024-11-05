'''
graph.py contains functions used to run causal discovery algorithms on the cleaned dataset
'''

from causallearn.search.ConstraintBased.PC import pc
from causallearn.utils.cit import kci
from causallearn.search.ConstraintBased.FCI import fci
from causallearn.search.ScoreBased.GES import ges
from causallearn.utils.GraphUtils import GraphUtils


def run_pc(data, fp, indep_test):
    '''
    Run PC on the data with the specified independence test and save the plot of the resulting causal graph.

    :param: data: dataframe to run PC on
    :param: fp: filepath name of causal graph result file
    :param: ind_test: string of the choice of independence test available for running PC from the causal-learn package ("fisherz", "chisq", "gsq", "kci", "mv_fisherz") 
    '''
    cg = pc(data.values, alpha=0.05, indep_test=indep_test)

    pyd = GraphUtils.to_pydot(cg.G)
    pyd.write_png(f"{fp}.png")


def run_fci(data, fp):
    '''
    Run FCI on the data and save the plot of the resulting causal graph.
    
    :param: data: dataframe to run FCI on
    :param: fp: filepath name of causal graph result file 
    '''
    data_array = np.array(data)
    g, edges = fci(data_array)

    pyd = GraphUtils.to_pydot(g, labels=data.columns)
    pyd.write_png(f"{fp}.png")


def run_ges(data, fp):
    '''
    Run GES on the data and save the plot of the resulting causal graph.
    
    :param: data: dataframe to run GES on
    :param: fp: filepath name of causal graph result file 
    '''
    Record = ges(data.values)
    
    pyd = GraphUtils.to_pydot(Record["G"], labels=data.columns)
    pyd.write_png(f"{fp}.png")
    
    return