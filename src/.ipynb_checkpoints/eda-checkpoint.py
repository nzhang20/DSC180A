'''
eda.py contains functions used to explore the linearity, Gaussianity, and correlation of variables in the cleaned dataframe
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import itertools


##### subject info graphs


def check_discrete_distribution(data):
    '''
    Plots grouped bar charts of the discrete variables in the subject dataset. Saves them in the plots folder as "subject_discrete_eda.png".

    :param: data: raw subject dataset
    '''
    fig, axs = plt.subplots(1, 2, figsize=(15, 15))
    fig.suptitle('Distribution of Discrete Variables in Subject Info Data')
    discrete_var = ['Ethnicity', 'Gender']
    groups = ['IR', 'IS', 'Unknown']

    if not set(discrete_var).issubset(data.columns): 
        raise Exception("The discrete variables in the raw subject data (Study, Race, Sex) are not in the dataset provided.")

    if 'IRIS' not in data.columns:
        raise Exception("The IRIS column is not in the dataset provided.")

    axs = axs.ravel()

    for i in range(len(discrete_var)):
        col = discrete_var[i]
        
        dfp = data.pivot_table(index=col, columns='IRIS', aggfunc='size')
        dfp.plot(kind='bar', rot=0, ax=axs[i])
        axs[i].set_title(col)
        axs[i].set_ylabel('Number of Individuals')
        for container in axs[i].containers:
            axs[i].bar_label(container)
        
    fig.savefig("plots/subject_discrete_eda.png")
    plt.close(fig)


def check_gaussian_distribution(data):
    '''
    Plots qq plots of the continuous variables in the subject dataset to check for normality. Saves them in the plots folder as "subject_qqplots.png".

    :param: data: raw subject dataset
    '''
    fig, axs = plt.subplots(1, 3, figsize=(20, 10))
    fig.suptitle('QQ Plots of Continuous Variables in Subject Info Data')
    continuous_var = ['Adj.age', 'BMI', 'SSPG']

    if not set(continuous_var).issubset(data.columns): 
        raise Exception("The continuous variables in the raw subject data (Adj.age, BMI, SSPG) are not in the dataset provided.")

    for i in range(len(continuous_var)):
        col = continuous_var[i]
        sm.qqplot(data[col], line='q', ax=axs[i])
        axs[i].set_title(col)
        
    fig.savefig("plots/subject_qqplots.png")
    plt.close(fig)


def check_linearity(data):
    '''
    Plots scatter plots of the continuous variables in the subject dataset to check for linearity. Saves them in the plots folder as "subject_scatter.png".

    :param: data: raw subject dataset
    '''
    fig, axs = plt.subplots(1, 3, figsize=(15, 10))
    fig.suptitle('Linearity of Continuous Variables in Subject Info Data')
    continuous_var = ['Adj.age', 'BMI', 'SSPG']
    pairs = list(itertools.combinations(continuous_var, 2))
    colors = {'IR': 'red', 'IS': 'green', 'Unknown': 'blue'}

    if not set(continuous_var).issubset(data.columns): 
        raise Exception("The continuous variables in the raw subject data (Age, BMI, SSPG) are not in the dataset provided.")

    for i in range(len(continuous_var)):
        x = data[pairs[i][0]]
        y = data[pairs[i][1]]
        c = [colors.get(i) for i in data['IRIS']]
        axs[i].scatter(x, y, c=c)
        axs[i].set_xlabel(pairs[i][0])
        axs[i].set_ylabel(pairs[i][1])

    fig.savefig("plots/subject_scatter.png")
    plt.close(fig)

        
##### for merged data


def numerical_encoding(data):
    '''
    Encodes categorical variables in the subject datset into numeric discrete variables. Returns the resulting dataframe. 

    :param: data: raw subject dataset
    ''' 
    X = data.copy()
    X['Ethnicity'] = X['Ethnicity'].map({'C': 0, 'A': 1, 'B': 2, 'H': 3, 'unknown': 4})
    X['Gender'] = X['Gender'].map({'M': 0, 'F': 1})
    X['IRIS'] = X['IRIS'].map({'IR': 0, 'IS': 1, 'Unknown': 2})
    return X
    

def corr(data, method='pearson'):
    '''
    Calculate the correlation of the columns in the dataframe using the specified method.
    
    :param: data: dataframe of columns to compute the correlation of
    :param: method: correlation analysis method, default is pearson
    '''
    if method == 'pearson':
        corr = data.corr()
    return corr


def make_corr_plot(corr_data, fp, title="Correlation Matrix"):
    '''
    Create a labeled correlation matrix of the columns in the dataframe. 
    
    :param: data: dataframe of columns to compute the correlation matrix of
    :param: fp: filepath of correlation matrix plot file
    '''
    full_columns = list(corr_data.columns)
    
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot()
    cax = ax.matshow(corr_data, cmap="coolwarm")
    fig.colorbar(cax)

    xaxis = np.arange(len(full_columns))
    ax.set_xticks(xaxis)
    ax.set_yticks(xaxis)
    ax.set_xticklabels(full_columns, rotation=90)
    ax.set_yticklabels(full_columns)
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(f"plots/{fp}.png")
    plt.close(fig)

    return


def check_linearity_large(data):
    '''
    Plots scatter plots of all pairs of the variables in the dataset to check for linearity. Saves them in the plots folder as "many_scatter.png".

    :param: data: clean dataset
    '''
    X = data.copy()
    pairs_X = list(itertools.combinations(X.columns, 2))
    
    num_rows = len(pairs_X) // 5
    # fig, axs = plt.subplots(num_rows, 5, figsize=(20, num_rows * 3))
    # fig.suptitle('Pairwise Scatter Plots for Linearity in Data')

    # axs = axs.ravel()

    # for i in range(len(pairs_X)):
    for i in range(len(pairs_X)):
        fig, ax = plt.subplots()
        x = pairs_X[i][0]
        y = pairs_X[i][1]
        # axs[i].scatter(x = X[x], y = X[y])
        # axs[i].set_xlabel(x)
        # axs[i].set_ylabel(y)
        ax.scatter(x = X[x], y = X[y])
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(f"Scatter plot of {x} and {y}")
        fig.savefig(f"plots/{x}_{y}_scatter.png")
        plt.close(fig)
        
    # fig.subplots_adjust(top=0.975, wspace=0.2, hspace=0.15)
    # fig.savefig("plots/many_scatter.png")
    # plt.close(fig)

    return
    
