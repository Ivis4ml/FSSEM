# FSSEM
`FSSEM`(Fused Sparse Structural Equational Model) is a package for Inference of Differential Gene Regulatory Networks Based on Gene Expression and Genetic Perturbation data.

# Installation
Install R

```r
source("./src/solver.R")
```
or

```r
source("./src/solver.min.R")
```

# Input format
For now, the fssem function `multiSML_iPALM` or `genSML_iPALM` requires data of gene expression and genetic perturbation data(eQTL data). `multiSML_iPALM`: run FSSEM if `X1=X2`; `genSML_iPALM`: run FSSEM if `X1 != X2`.

I use the `list` data for representing the info of all data. One `data` list has two elements:

* `obs`: Observation data for genes and eQTLs, it has:
  + `Y1, Y2`: a gene expression matrix under 2 conditions.
  + `X` or `X1,X2`: eQTL data for each gene under 2 conditions.
  + `sk`: nonzero eQTL index for each gene.

* `var`: Inherent variable of FSSEM, it has:
  + `N`: sample size per each condition.
  + `Ng`: total gene number. 
  + `Ne`: total eQTL number. 

# Usage
```r
source("./src/demo.R")
```

Run `demo.R` code in `src` folder. And in this demo, you can simulate gene expression and genetic perturbation
data with random generated network structures. You can randomly generate networks with parameters: 

+ `N`, sample size per each condition; 
+ `Ng`, gene number; 
+ `Ne`, eQTL number; 
+ `Ns`, sparse ratio of gene network(number of nonzero entries per gene).
+ `sigma`, noise standard deviation.


