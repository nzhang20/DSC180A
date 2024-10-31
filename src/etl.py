'''
etl.py contains functions used to merge and clean the two raw dataframes
'''


def merge_gut_subject():
    '''
    Merges the gut microbes dataset and subject info dataset on SubjectID. 
    '''
    gut_microbes = pd.read_csv("data/raw/gut_16s_abundance.txt", sep="\t")
    subject_info = pd.read_csv("data/raw/subject_file.csv")
    
    gut_microbes["SubjectID"] = gut_microbes["SampleID"].str.split("-").str[0]
    
    merged_df = pd.merge(gut_microbes, subject_info, on="SubjectID", how="inner")
    
    return(merged_df)


def get_genus_and_covariates(df, genus, covariates):
    '''
    Returns a dataframe of the selected genera and covariates. 
    
    :param: df: merged dataset containing microbes and subject info
    :param: genus: a list of genera to choose
    :param: covariates: a list of covariates to choose 
    '''
    X = df.loc[:, genus + covariates]
    return(X)