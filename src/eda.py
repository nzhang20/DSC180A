'''
eda.py contains functions used to explore the linearity, Gaussianity, and correlation of variables in the cleaned dataframe
'''

# subject info graphs

def check_discrete_distribution(data):
    fig, axs = plt.subplots(2, 2, figsize=(15, 15))
    fig.suptitle('Distribution of Discrete Variables in Subject Info Data')
    discrete_var = ['Study', 'Race', 'Sex']
    groups = ['IR', 'IS', 'unknown']

    axs = axs.ravel()

    for i in range(len(discrete_var)):
        col = discrete_var[i]
        
        dfp = data.pivot_table(index=col, columns='IR_IS_classification', aggfunc='size')
        dfp.plot(kind='bar', rot=0, ax=axs[i])
        axs[i].set_title(col)
        axs[i].set_ylabel('Number of Individuals')
        for container in axs[i].containers:
            axs[i].bar_label(container)
        
    plt.show()


def check_gaussian_distribution(data):
    fig, axs = plt.subplots(1, 3, figsize=(20, 10))
    fig.suptitle('Distribution of Continuous Variables in Subject Info Data')
    continuous_var = ['Age', 'BMI', 'SSPG']

    for i in range(len(continuous_var)):
        col = continuous_var[i]
        print('Number of NaNs:', data[col].isna().sum())
        sm.qqplot(data[col], line='q', ax=axs[i])
        axs[i].set_title(col)
        
    plt.show()




def check_linearity(data):
    fig, axs = plt.subplots(1, 3, figsize=(15, 10))
    fig.suptitle('Linearity of Continuous Variables in Subject Info Data')
    pairs = list(itertools.combinations(continuous_var, 2))
    colors = {'IR': 'red', 'IS': 'green', 'Unknown': 'blue'}

    for i in range(len(continuous_var)):
        x = data[pairs[i][0]]
        y = data[pairs[i][1]]
        c = [colors.get(i) for i in data['IR_IS_classification']]
        axs[i].scatter(x, y, c=c)
        axs[i].set_xlabel(pairs[i][0])
        axs[i].set_ylabel(pairs[i][1])

# for merged data

def clean_merged_df(data):
    X = merged_df.drop(columns=['SampleID', 'SubjectID', 'Study'])
    X['Race'] = X['Race'].map({'C': 0, 'A': 1, 'B': 2, 'H': 3, 'unknown': 4})
    X['Sex'] = X['Sex'].map({'M': 0, 'F': 1})
    X['IR_IS_classification'] = X['IR_IS_classification'].map({'IR': 0, 'IS': 1, 'Unknown': 2})
    return X

def check_corr_plot(data, IR_IS=false):
    
    X = clean_merged_df(data)
    if IR_IS = true:
        X_IR = X[X['IR_IS_classification'] == 0]
        X_IS = X[X['IR_IS_classification'] == 1]
        make_corr_plot(X_IR)
        make_corr_plot(X_IS)
    make_corr_plot(X)

def corr_IR_IS(data):
    X = clean_merged_df(data)
    X_IR = X[X['IR_IS_classification'] == 0]
    X_IS = X[X['IR_IS_classification'] == 1]
    X_IR_corr = X_IR.corr()
    X_IS_corr = X_IS.corr()
    corr = X_IR_corr - X_IS_corr
    full_columns = list(X_IR.columns)

    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot()
    cax = ax.matshow(corr, cmap='coolwarm')
    fig.colorbar(cax)

    xaxis = np.arange(len(full_columns))
    ax.set_xticks(xaxis)
    ax.set_yticks(xaxis)
    ax.set_xticklabels(full_columns, rotation=90)
    ax.set_yticklabels(full_columns)

    fig.show()

def check_merged_discrete(data):
    X = clean_merged_df(data)
    fig, axs = plt.subplots(255, 5, figsize=(20, 700))
    fig.suptitle('Distribution of Discrete Variables in Subject Info Data')

    pairs_X = list(itertools.combinations(X.columns, 2))

    axs = axs.ravel()

    for i in range(len(pairs_X)):
        x = pairs_X[i][0]
        y = pairs_X[i][1]
        axs[i].scatter(x = X[x], y = X[y])
        axs[i].set_xlabel(x)
        axs[i].set_ylabel(y)
        
    fig.subplots_adjust(top=0.975, wspace=0.2, hspace=0.15)