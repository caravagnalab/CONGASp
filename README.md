# Copy number genotyping jointly from scRNA and scATAC sequencing



A set of Pyro models and functions to infer CNA from scRNA-seq and scATAC-seq data. 
It comes with a companion [R package](https://github.com/caravagnalab/rcongas) that works as an interface and provides preprocessing, simulation and visualization routines.


Currently providing:

- A mixture model on segments where CNV are modelled as Categorical random variable (LatentCategorical) 

<!--Coming soon:
- A linear model in the emission that can account for known covariates
- The equivalent of MixtureGaussian but with CNVs as Categorical random variable
- A model on genes (all the other models assume a division in segments)
-->
To install:

`$ pip install congas`

<!--
To run a simple analysis on the example data

```python
import congas as cn
from congas.models import MixtureGaussian
data_dict = cn.simulation_data
params, loss = cn.run_analysis(data_dict,MixtureGaussian, steps=200, lr=0.05)
```


[Full Documentation](https://annealpyro.readthedocs.io/en/latest/)
-->