'''
etl.py contains functions used to merge and clean the two raw dataframes
'''

import pandas as pd
import numpy as np


def get_subject_data():
    '''
    Returns the subject info dataset. 
    '''
    subjects = pd.read_csv("data/raw/S1_Subjects.csv")
    subjects_IRIS_known = subjects[subjects["IRIS"] != "Unknown"]
    return subjects_IRIS_known


def get_sample_data():
    '''
    Returns the sample info dataset.
    '''
    samplelist = pd.read_csv("data/raw/S3_SampleList.csv")
    return samplelist
    

def get_gut_data():
    '''
    Returns the gut microbe dataset. 
    '''
    gut_microbes = pd.read_csv("data/raw/gut_16s_abundance.txt", sep="\t")
    return gut_microbes


def clean_sample_data(samplelist):
    '''
    Returns the sample info dataset for samples where gut microbe data were collected and the sample was a healthy visit. 

    :param: samplelist: samplelist raw dataset
    '''
    if ("Gut_16S" not in samplelist.columns) and ("CL4" not in samplelist.columns):
        raise Exception("Gut_16S (gut microbe test) and CL4 (IR IS classification) are not in the columns of samplelist.")
        
    samplelist_gut_healthy = samplelist[(samplelist["Gut_16S"] == 1) & (samplelist["CL4"] == "Healthy")]
    return samplelist_gut_healthy
    

def merge_gut_sample_subject():
    '''
    Merges the gut microbes dataset and samplelist dataset on SampleID. 
    Merges the result of the previous merge with the subject info dataset on SubjectID. 
    '''
    merged_df = pd.merge(get_gut_data(), clean_sample_data(get_sample_data()), on="SampleID", how="inner")
    merged_df = pd.merge(merged_df, get_subject_data(), on="SubjectID", how="inner")
    
    return merged_df


def get_bacteria_and_covariates(df, fp = "data/clean.csv", **columns):
    '''
    Returns a dataframe of the selected bacteria by taxa and covariates. 
    
    :param: df: merged dataset containing microbes and subject info
    :param: columns: one or more lists of specific phylum/class/order/family/genus bacteria and covariates (must include IR_IS_classification)
    '''
    all_columns = []
    
    for val in columns.values():
        all_columns += val
        
    X = df.loc[:, all_columns]
    X.to_csv(fp, index=False)

    return


def IR_IS_split(df, numerical=False):
    '''
    Returns two dataframes where each is separated by the IRIS column.

    :param: df: dataframe containing all individuals and the 'IRIS' column
    '''
    if "IRIS" not in df.columns:
        raise Exception("IRIS is not in the columns of df.")

    if numerical:
        IR_df = df[df['IRIS'] == 0]  # Insulin-resistant group
        IS_df = df[df['IRIS'] == 1]  # Insulin-sensitive group

    else:
        IR_df = df[df['IRIS'] == 'IR']
        IS_df = df[df['IRIS'] == 'IS']

    IR_df = IR_df.drop(columns='IRIS')
    IS_df = IS_df.drop(columns='IRIS')
    
    return IR_df, IS_df


def create_adj_matrix(corr_data, corr_method, genus_lst):
    '''
    Transforms the correlation data given in the paper into an adjacency matrix representing the correlations as edges. Returns a dataframe for IR and IS separately and in this order. 

    :param: corr_data: dataframe containing the correlations of IR and IS specific cohorts, as well as their correlation values and p-values
    :param: corr_method: string of the corr method; either "sparcc" or "rmcorr"
    '''
    if corr_method == "sparcc":
        corr_IR = corr_data[corr_data["IR_rho"].notnull()]
        corr_IR = corr_IR[corr_IR["IR_p"] <= 0.05]
        corr_IR["genus_A"] = corr_IR["cor_char"].apply(lambda x: x.split(" vs ")[0])
        corr_IR["genus_B"] = corr_IR["cor_char"].apply(lambda x: x.split(" vs ")[1])
        corr_IR_subset = corr_IR[["genus_A", "genus_B"]]

        corr_IS = corr_data[corr_data["IS_rho"].notnull()]
        corr_IS = corr_IS[corr_IS["IS_p"] <= 0.05]
        corr_IS["genus_A"] = corr_IS["cor_char"].apply(lambda x: x.split(" vs ")[0])
        corr_IS["genus_B"] = corr_IS["cor_char"].apply(lambda x: x.split(" vs ")[1])
        corr_IS_subset = corr_IS[["genus_A", "genus_B"]]

    if corr_method == "rmcorr": 
        corr_IR = corr_data[corr_data["IR_adj.p"] <= 0.05]
        corr_IR["genus_A"] = corr_IR["cor_char"].apply(lambda x: x.split(" vs ")[0])
        corr_IR["genus_B"] = corr_IR["cor_char"].apply(lambda x: x.split(" vs ")[1])
        corr_IR_subset = corr_IR[["genus_A", "genus_B"]]

        corr_IS = corr_data[corr_data["IS_adj.p"] <= 0.05]
        corr_IS["genus_A"] = corr_IS["cor_char"].apply(lambda x: x.split(" vs ")[0])
        corr_IS["genus_B"] = corr_IS["cor_char"].apply(lambda x: x.split(" vs ")[1])
        corr_IS_subset = corr_IS[["genus_A", "genus_B"]]

    matrix_IR = np.zeros((len(genus_lst), len(genus_lst)))
    matrix_IS = np.zeros((len(genus_lst), len(genus_lst)))
    
    for i, genus_A in enumerate(genus_lst):
      for j, genus_B in enumerate(genus_lst):
        if any((corr_IR_subset['genus_A'] == genus_A) & (corr_IR_subset['genus_B'] == genus_B)):
          matrix_IR[i, j] = 1
        if any((corr_IS_subset['genus_A'] == genus_A) & (corr_IS_subset['genus_B'] == genus_B)):
          matrix_IS[i, j] = 1

    return matrix_IR, matrix_IS