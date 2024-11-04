'''
graph.py contains functions used to run causal discovery algorithms on the cleaned dataset
'''


def make_corr_plot(data, fp):
    '''
    Create a labeled correlation matrix of the columns in the dataframe. 
    
    :param: data: dataframe of columns to compute the correlation matrix of
    :param: fp: filepath name of correlation matrix plot file (pdf for vectorized)
    '''
    corr = data.corr()
    full_columns = list(data.columns)
    
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot()
    cax = ax.matshow(corr)
    fig.colorbar(cax)

    xaxis = np.arange(len(full_columns))
    ax.set_xticks(xaxis)
    ax.set_yticks(xaxis)
    ax.set_xticklabels(full_columns, rotation=90)
    ax.set_yticklabels(full_columns)
    ax.set_title("Correlation Matrix")

    fig.savefig(f"{fp}.pdf")
    plt.close(fig)

    return


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

def run_fci(data, fp):
    '''
    Run FCI on the data and save the plot of the resulting causal graph.
    
    :param: data: dataframe to run FCI on
    :param: fp: filepath name of causal graph result file 
    '''
    data_array = np.array(data)
    g, edges = fci(data_array)

    pdy = GraphUtils.to_pydot(g, labels=data.columns)
    pdy.write_png(f"{fp}.png")