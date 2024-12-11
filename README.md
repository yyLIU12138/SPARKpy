# SPARKpy

Here is a Python version of the SPARK method, which was developed based on the "Statistical analysis of spatial expression patterns for spatially resolved transcriptomic studies" paper. 

It should be noted that the original SPARK method includes many manual modifications, whereas my implementation largely adheres to the approach described in the original paper with minimal manual modifications. As a result, SPARKpy does not produce the same results as SPARK and should not be used as a substitute for SPARK in formal comparisons in works.

## Installation
Clone the repository. 

```
git clone https://github.com/yyLIU12138/SPARKpy.git
cd SPARKpy
```

Create an environment.

```
conda create -n SPARKpy python=3.11
conda activate SPARKpy
```

Install the required packages.

```
pip install -r requirements.txt
```

Install SPARKpy.

```
python setup.py build
python setup.py install
```
## Toy example

```
import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings('ignore')

from SPARKpy import *

data_dir = './data/processed_data/MOB/'
save_dir = './results/Rep11_MOB/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# # read in the raw data and filter it
# adata = sc.read_h5ad("XXXXXXX")
# spark_filter(adata, min_spot_percentage=0.1, min_gene_counts=10, inplace=True)
# adata.write_h5ad(data_dir + "filtered_data.h5ad")

# read in the filtered data
adata = sc.read_h5ad(data_dir + "filtered_data.h5ad")

# Initialize the model
model = SPARKpy(adata, use_rep='spatial', cov_mtx=None, max_iter=500, tol=1e-5, fit='Poisson', verbose=True)

# Fit the model and get the p-values for each kernel
# as well as the combined p-values and adjusted p-values
p_value_df = model.fit_spark()

p_value_df.to_csv(save_dir + 'p_values_sparkpy.csv', index=True)
```




