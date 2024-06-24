# microSLAM
### A R package to perform population structure leveraged association models for the microbiome.
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/logo.png)

### below is a guide showing all operative function in R
## Quick start guide
### install package in R
```
library(devtools)
install_github('miriam-goldman/microSLAM')
```
### read in example data which in this case is example output from MIDAS v3
```
library(microSLAM)
library(tidyverse)
library(magrittr)
exp_metadata <- read.csv("../example_data/exp_metadata.csv")
exp_genedata <- read.csv("../example_data/genepresabs.csv")
```
### Step 1 calculate Genetic Relatedness Matrix GRM or make your own
```
GRM <-calculate_GRM(exp_genedata)
```

### Step 2 fit and test tau

```
glm_fit0=glm("y~age+1", data = exp_metadata, family = "binomial")
glmm_fit=fit_tau_test(glm_fit0, GRM,species_id="test",verbose = FALSE,log_file=NA)
tautestfit=run_tau_test(glm_fit0, GRM,n_tau=100,species_id = "test", tau0=1, phi0=1)
pvalue=sum(tautestfit$t>=glmm_fit$t)/100
print(paste("t value:",glmm_fit$t,"pvalue:",pvalue))
```

### Step 3 fit and test beta for each gene

```
gene_long<-exp_genedata %>% pivot_longer(cols=starts_with("gene"),names_to ="gene_id",values_to = "gene_value") %>% as.data.frame()
gene_test_df<-fit_beta(glmm_fit,glm_fit0,GRM,gene_long,SPA=TRUE)
ggplot(gene_test_df,aes(beta,-log10(SPA_pvalue)))+geom_point()
```
