# microSLAM
### An R package to perform population structure leveraged association modeling for the microbiome.
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/logo.png)

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
GRM = calculate_grm(exp_genedata)
```
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/exampleGRM.png)

### Step 2: tau test for population structure (strain-trait associations)
#### fit baseline glm for starting parameters in tau test

```
glm_fit0=glm("y~age+1", data = exp_metadata, family = "binomial")
```

#### fit tau test using baseline glm and GRM calculated above
```
glmm_fit=fit_tau_test(glm_fit0, GRM,species_id="test",verbose = FALSE,log_file=NA)
summary(glmm_fit)
```
```
Species ID:  test
Formula:  y~age+1+b
family:  binomial logit
Fixed-effect covariates estimates:
 (Intercept) age
 -0.287 0.005
Converged:  TRUE
Number of iterations: 5
Tau:  2.354
Phi:  1 if logit or binomail should be 1
T value of tau: 0.624
Number of Samples: 100
```


#### test the significance of the tau that was fit
```
tautestfit=run_tau_test(glm_fit0, GRM,n_tau,species_id = "test", tau0=1, phi0=1)
```
#### calculate the pvalue from the permutation test ran on tau
```
pvalue=sum(tautestfit$t>=glmm_fit$t)/n_tau
```
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/permutation.png)

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
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/volcano.png?raw=true)


#### example output dataframe

![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/betadf.png?raw=true)
Columns are:
```
species_id: the id of the bacterial species run
tau: estimate for tau variance variable
gene_id: gene id for gene run
cor: correlation between the y variable and the gene data for that gene
cor_to_b: correlation between the random b variable and the gene data for that gene
z: z value estimated for gene from GLMM
var1: variance estimated for gene from GLMM
beta: beta estimate for gene from GLMM
se_beta: standard error for beta estimated from GLMM
t_adj: t value adjusted estimated from GLMM
SPA_pvalue: saddle point adjusted pvalue from GLMM
spa_score: saddle point adjusted t score from GLMM
SPA_zvalue:  saddle point adjusted z value from GLMM
pvalue_noadj: pvalue not adjusted from saddle point approximation
converged: did the spa algrothim converge or not
```
