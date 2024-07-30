# microSLAM
### An R package to perform population structure leveraged association modeling for the microbiome from metagenomics data. Data are taken from outputs of MIDAS 3.
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/logo.png)

### Below is a guide showing all operative functionality in R.
## Quick start guide
### Step 1: install package in R
```
library(devtools)
install_github('miriam-goldman/microSLAM')
```
### Step 2: read in example data which in this case is simulated example output from MIDAS v3.
Any dataframe of gene data, a phenotype of interest, and any covariates will do. The samples must match between the metadata and the gene data.
```
library(microSLAM)
library(tidyverse)
library(magrittr)
n_tau = 100 ### number of permutations for tau test
exp_metadata = read.csv("../example_data/exp_metadata.csv") ### read in example metadata
## metadata is a sample by different covariates and phenotype matrix
exp_genedata = read.csv("../example_data/genepresabs.csv") ### read in example gene data
## gene data is a sample by gene matrix for one species
head(exp_metadata)
  y strain uncorrelated_strain age sample_name
1 0      0                   0  52     sample1
2 0      0                   0  45     sample2
3 1      1                   0  41     sample3
4 1      0                   0  32     sample4
5 0      0                   0  23     sample5
6 1      1                   0  29     sample6
```
```
head(exp_genedata)
 sample_name gene1 gene2 gene3 gene4 gene5 gene6 gene7 gene8 gene9
1     sample1     0     1     1     1     1     0     0     0     0
2     sample2     0     1     1     0     0     0     0     0     0
3     sample3     0     0     0     1     0     0     0     1     0
4     sample4     0     0     0     0     0     1     0     0     1
5     sample5     0     1     0     0     0     0     0     0     1
6     sample6     1     1     0     0     0     0     0     1     0
```

### Step 3: calculate Genetic Relatedness Matrix (GRM) or make your own
#### calculate the GRM (1-hamming distance between each sample) from gene data a samples by genes matrix
Dataframe must contain a sample_name column to generate GRM correctly. This GRM represents how similar our samples are within the space of the genes of this bacterial species but other representations should be valid. We will use this to estimate our random effects for our mixed model and this is used to represent the population structure.
```
GRM = calculate_grm(exp_genedata)
```
#### Visualization of the GRM with the strain information used to generate it are labeled.
This gene data was generated with a strain that is correlated to y in half of the samples and another strain or subset of semi correlated genes that are not correlated to the hypothetical phenotype of interest. We are most interested in genes that are able to explain our phenotype more than the simulated strain, in this case 3 genes were simulated to be more related to the phenotype than the strain. 
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/exampleGRM.png)

### $\tau$ test for population structure (strain-trait associations)
#### Step 4: fit baseline glm for starting parameters in tau test, this can be done as you would do a normal glm for your data.

```
glm_fit0=glm("y~age+1", data = exp_metadata, family = "binomial")
```

#### Step 5: fit tau test using baseline glm and GRM calculated above
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
##### The outputs of this can be interpreted as the microSLAM was able to converge, the $\tau$ variable is 2.35, and we have estimated our random variable b as well as our covariates coefficients.

#### Step 6: test the significance of the tau that was fit with a permuation test
```
tautestfit=run_tau_test(glm_fit0, GRM,n_tau,species_id = "test", tau0=1, phi0=1)
```
##### Calculate the pvalue from the permutation test ran on tau
```
pvalue=sum(tautestfit$t>=glmm_fit$t)/n_tau
```
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/permutationnew.png)

##### In this case compared to 100 permutation of the tau test we have the most extreme example and a significant pvalue. Therefor it makes sense to use the microSLAM mixed model to find the betas for the gene-by-gene-trait associations

### $\beta$ test for gene-trait associations
#### Step 7: transform the gene data to a long matrix to fit each gene separately

```
gene_long = exp_genedata %>% pivot_longer(cols=starts_with("gene"),
names_to ="gene_id",
values_to = "gene_value") %>%
as.data.frame()
```
#### Step 8: fit a beta for each gene adjusting for the population structure
```  
gene_test_df = fit_beta(glmm_fit,glm_fit0,GRM,gene_long,SPA=TRUE)
```
##### Plot a simple volcano plot of the results, in this case genes 1, 2, and 3, have been generated as being more related to the phenotype than the simulated strain. These are colored in red.
```
ggplot(gene_test_df,aes(beta,-log10(SPA_pvalue)))+geom_point()
```
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/volcanoplotcolor.png?raw=true)


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

##### If you want to visualize the genes that are correlated with the strain you can filter the genes to those that are > than some cut off correlated to b. Example of this is below.


```
straingenes = gene_test_df %>% filter(abs(cor_to_b)>.4)
rownames(exp_genedata) = exp_genedata$sample_name
straingenesvalues = exp_genedata[,straingenes$gene_id]
pheatmap(straingenesvalues,show_rownames=FALSE,show_colnames=FALSE,treeheight_row=0,treeheight_col = 0,labels_row="samples",labels_col="samples",
                    main=paste("Strain Associated Genes"),border_color=NA,annotation_row = exp_metadata[,-5],color=myColor,clustering_method="average")
```
![alt text](https://github.com/miriam-goldman/microSLAM/blob/main/other/strainheatmap.png?raw=true)

##### 300 genes were modeled from the strain, and we are able to recover 288 of these through this cut off. These can be used to find which genes make up substrains related to the phenotype of interest.
