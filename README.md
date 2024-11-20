# Causal Discovery in Gut Microbes for Type 2 Diabetes

This project runs traditional causal discovery methods on gut microbes and type 2 diabetes (T2D) data to find causal structures among gut microbes and causal structures among gut microbes and T2D. Our goal is to determine how variable these causal discovery methods perform in contrast with the significant correlations found by the original study. We use PC, FCI, and GES from the `causal-learn` package, and we implement our own algorithm to address statistical power issues. We visualize our results using graphs with (TBD) and comparing graphs using classification algorithms like SVMs. 

## File Descriptions
### Data
There are a few files of relevant raw data from the Stanford Integrated Human Microbiome Project.
- `data/raw/S1_Subjects.csv`: Data on the participants of the study such as their "IRIS" classification, Gender, Ethnicity, etc (Supplementary Table 1). 
- `data/raw/S3_SampleList.csv`: Data on every sample recorded in the study such as their Health status and which omics tests it corresponded to (Supplementary Table 3). 
- `data/raw/gut_16s_abundance.txt`: Normalized abundance of each gut microbe phylum/class/order/family/genus strain for each sample
- `data/raw/rmcorr.csv`: Significant correlations found using the CLR + rmcorr method (Supplementary Table 32).
- `data/raw/sparcc.csv`: Significant correlations found using the SPARCC method (Supplementary Table 32).
- `data/raw/subject_file.csv`: Demographic information and covariates collected on each individual including their disease status (`IR_IS_classification`)

### Scripts
Activate the virtual environment with `source venv/bin/activate` or install the dependencies with `pip install -r requirements.txt`. 

To replicate our results, run...
- `python run.py all`: runs the `data`, `eda` and `graph` targets to reproduce our entire analysis from the raw data
- `python run.py data`: recreate our "clean" dataset from the raw data files described above
- `python run.py eda`: runs our EDA to produce plots on assessing linearity, Gaussianity, and correlations between features
- `python run.py graph`: runs the four (soon to be five) causal discovery algorithms and generates results in the form of causal graphs
- `python run.py clean`: remove all generated files (e.g. the clean dataset, graphs)