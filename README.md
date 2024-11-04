# Causal Discovery in Gut Microbes for Type 2 Diabetes

This project runs traditional causal discovery methods on gut microbes and type 2 diabetes (T2D) data to find causal structures among gut microbes and causal structures among gut microbes and T2D. Our goal is to determine how well these causal discovery methods perform in contrast with "ground truth" studies. We use PC, FCI, and GES from the `causal-learn` package and visualize with (TBD). 

## File Descriptions
### Data
There are two files of relevant raw data from the Stanford Integrated Human Microbiome Project.
- `data/raw/gut_16s_abundance.txt`: Normalized abundance of each gut microbe phylum/class/order/family/genus strain for each sample
- `data/raw/subject_file.csv`: Demographic information and covariates collected on each individual including their disease status (`IR_IS_classification`)

### Scripts
Activate the virtual environment with `source venv/bin/activate` or install the dependencies with `pip install -r requirements.txt`. 

To replicate our results, run...
- `python run.py all`:
- `python run.py data`: recreate our "clean" dataset from the two raw data files described above
- `python run.py graph`: run the three causal discovery algorithms and generate results in the form of causal graphs
- `python run.py clean`: remove all generated files (e.g. the clean dataset, graphs)