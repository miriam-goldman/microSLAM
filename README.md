# microSLAM
### An R package to perform population structure leveraged association modeling for the microbiome.
![alt text](https://github.com/miriam-goldman/mircoSLAM/blob/main/other/logo.png)

### Below is a guide showing all operative functionality in R.
## Quick start guide
### install package in R
```
library(devtools)
install_github('miriam-goldman/microSLAM')
```
### read in example simulation data which in this case is example output from MIDAS v3
```
library(microSLAM)
library(tidyverse)
library(magrittr)
n_tau=100 ### number of permutations for tau test
exp_metadata = read.csv("../example_data/exp_metadata.csv") ### read in example metadata
exp_genedata = read.csv("../example_data/genepresabs.csv") ### read in example gene data
```
### Step 1 calculate Genetic Relatedness Matrix GRM or make your own
#### calculate the GRM (1-hamming distance between each sample) from gene data a samples by genes matrix
```
GRM = calculate_GRM(exp_genedata)
```
![alt text](https://github.com/miriam-goldman/mircoSLAM/blob/main/other/exampleGRM.png)

### Step 2: tau test for population structure (strain-trait associations)
#### fit baseline glm for starting parameters in tau test

```
glm_fit0=glm("y~age+1", data = exp_metadata, family = "binomial")
```

#### fit tau test using baseline glm and GRM calculated above
```
glmm_fit=fit_tau_test(glm_fit0, GRM,species_id="test",verbose = FALSE,log_file=NA)
```

#### test the significance of the tau that was fit
```
tautestfit=run_tau_test(glm_fit0, GRM,n_tau,species_id = "test", tau0=1, phi0=1)
```
#### calculate the pvalue from the permutation test ran on tau
```
pvalue=sum(tautestfit$t>=glmm_fit$t)/n_tau
```
![alt text](https://github.com/miriam-goldman/mircoSLAM/blob/main/other/permutation.png)

### Step 3 beta test for gene-trait associations
#### transform the gene data to a long matrix to fit each gene separately

```
gene_long = exp_genedata %>% pivot_longer(cols=starts_with("gene"),
names_to ="gene_id",
values_to = "gene_value") %>%
as.data.frame()
```
#### fit a beta for each gene adjusting for the population structure
```  
gene_test_df = fit_beta(glmm_fit,glm_fit0,GRM,gene_long,SPA=TRUE)
```
#### plot a simple volcano plot of the results
```
ggplot(gene_test_df,aes(beta,-log10(SPA_pvalue)))+geom_point()
```
![alt text](https://github.com/miriam-goldman/mircoSLAM/blob/main/other/volcano.png?raw=true)
