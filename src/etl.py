'''
etl.py contains functions used to merge and clean the two raw dataframes
'''
import pandas as pd


def get_subject_data():
    '''
    Returns the subject info dataset. 
    '''
    subject_info = pd.read_csv("data/raw/subject_file.csv")
    return subject_info

def get_gut_data():
    '''
    Returns the gut microbe dataset. 
    '''
    gut_microbes = pd.read_csv("data/raw/gut_16s_abundance.txt", sep="\t")
    gut_microbes["SubjectID"] = gut_microbes["SampleID"].str.split("-").str[0]
    return gut_microbes

def merge_gut_subject():
    '''
    Merges the gut microbes dataset and subject info dataset on SubjectID. 
    '''
    gut_microbes = pd.read_csv("data/raw/gut_16s_abundance.txt", sep="\t")
    subject_info = pd.read_csv("data/raw/subject_file.csv")
    
    gut_microbes["SubjectID"] = gut_microbes["SampleID"].str.split("-").str[0]
    
    merged_df = pd.merge(gut_microbes, subject_info, on="SubjectID", how="inner")
    
    return(merged_df)


def get_bacteria_and_covariates(df, fp = "data/clean.csv", **columns):
    '''
    Returns a dataframe of the selected bacteria by taxa and covariates. 
    
    :param: df: merged dataset containing microbes and subject info
    :param: columns: one or more lists of specific phylum/class/order/family/genus bacteria and covariates (must include IR_IS_classification)
    '''
    # try catch for IR_IS_classification in columns
    
    all_columns = []
    
    for val in columns.values():
        all_columns += val
        
    X = df.loc[:, all_columns]
    X.to_csv(fp)

    return


def IR_IS_classify(df):
    '''
    Returns two dataframes where each is separated by the IR_IS_classification column.

    :param: df: dataframe containing all individuals and the 'IR_IS_classification' column
    '''
    IR_df = df[df['IR_IS_classification'] == 'IR']  # Insulin-resistant group
    IS_df = df[df['IR_IS_classification'] == 'IS']  # Insulin-sensitive group
    
    return IR_df, IS_df
