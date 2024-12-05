# Causal Discovery in Gut Microbes for Type 2 Diabetes

This project runs traditional causal discovery methods on gut microbes and type 2 diabetes (T2D) data to find causal structures among gut microbes and causal structures among gut microbes and T2D. Our goal is to determine how variable these causal discovery methods perform in contrast with the significant correlations found by the original study. We use PC, FCI, and GES from the `causal-learn` package, and we implement our own algorithm to address statistical power issues. We visualize our results using graphs and adjacency matrix heatmaps and comparing graphs using frequency heatmaps. Methods are done separately on the insulin-resistant (IR) and insulin-sensitive (IS) cohorts.

## File Descriptions
### Data
There are a few files of relevant raw data from the Stanford Integrated Human Microbiome Project.
- `data/raw/S1_Subjects.csv`: Data on the participants of the study such as their "IRIS" classification, Gender, Ethnicity, etc (Supplementary Table 1). 
- `data/raw/S3_SampleList.csv`: Data on every sample recorded in the study such as their Health status and which omics tests it corresponded to (Supplementary Table 3). 
- `data/raw/gut_16s_abundance.txt`: Normalized abundance of each gut microbe phylum/class/order/family/genus strain for each sample
- `data/raw/rmcorr.csv`: Significant correlations found using the CLR + rmcorr method (Supplementary Table 32).
- `data/raw/sparcc.csv`: Significant correlations found using the SPARCC method (Supplementary Table 32).
- `data/raw/subject_file.csv`: Demographic information and covariates collected on each individual including their disease status (`IR_IS_classification`)

### Data Cleaning
The clean dataset in `data/clean.csv` is created by filtering the merged `S1_Subjects`, `S3_SampleList`, and `gut_16s_abundance.txt` on subjects that had a known IRIS classification and samples that were collected at a "Healthy" visit.

### EDA
In order to use causal discovery algorithms, we need to verify that our data meet certain assumptions. Depending on whether there is linearity and Gaussianity in the data, we are limited to different algorithms. We also check the Pearson correlations between genera for the IR and IS cohorts separately. 

### Causal Discovery Algorithms & Graphs
There are (currently) four methods we are using to obtain a causal graph for each of the IR and IS cohorts.
1. **Peter-Clark (PC) with (Fast) Kernel-based Conditional Independence (KCI)**. This is the most well-known/widely-used causal discovery algorithm in current literature. It is a constraint-based method.
2. **Fast Causal Inference (FCI)**. This is a variant of PC that can handle unknown confounding variables by denoting these with double-arrowed edges. 
3. **Greedy Equivalence Search (GES)**. This is a score-based method that iteratively adds and removes edges based on score improvement.
4. **Our algorithm**. We test each of the correlations found in the paper with conditional independence tests on conditioning sets of size 1 and size 2. This limits the number of tests we conduct and addresses the low power issue.

### Scripts
Create and activate the conda environment with `conda env create -f environment.yml` and `conda activate capstone`. If for some reason, you are running into error importing "fastkci" from causallearn.utils.cit, please try to install causal-learn from their github repo as such: `pip install causal-learn@git+https://github.com/py-why/causal-learn.git`. 

Alternatively, activate the virtual environment with `source venv/bin/activate` or install the dependencies with `pip install -r requirements.txt`. 

To replicate our results, run...
- `python run.py all`: runs the `data`, `eda` and `graph` targets to reproduce our entire analysis from the raw data
- `python run.py data`: recreate our "clean" dataset from the raw data files described above
- `python run.py eda`: runs our EDA to produce plots on assessing linearity, Gaussianity, and correlations between features
- `python run.py graph`: runs the four causal discovery algorithms and generates results in the form of causal graphs
- `python run.py clean`: remove all generated files (e.g. the clean dataset, graphs)