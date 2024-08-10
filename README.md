<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/logo.png" width=200>

## Introduction

MicroSLAM is an R package to perform population structure leveraged association modeling for the microbiome from metagenomics data. Microbiome association studies typically link host disease or other traits to summary statistics measured in metagenomics data, such as diversity or taxonomic composition. But identifying disease-associated species based on their relative abundance does not provide insight into why these microbes act as disease markers, and it overlooks cases where disease risk is related to specific strains with unique biological functions. microSLAM is an implementation of a  mixed-effects model that performs association tests that connect host traits to the presence/absence of genes within each microbiome species, while accounting for strain genetic relatedness across hosts. Traits can be quantitative or binary (such as case/control).

MicroSLAM is fit in three steps for each species. The first step estimates a genetic relatedness matrix (GRM) that measures population structure of the microbial species across hosts. Step two calculates the association between population structure and the trait (y), enabling detection of species for which a subset of related strains confer risk. To identify specific genes (G) whose presence/absence across diverse strains is associated with the trait after adjusting for population strucuture, step three models the trait as a function of gene occurrence plus random effects (b) estimated from step two. Steps two and three can include adjustment for covariates measured on each host. 

<p align="center">
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/Newflowchart.png" width=600>
</p>

### Below is a guide showing all operative functionality in R.
Install package in R

```
library(devtools)
install_github('miriam-goldman/microSLAM')
library(microSLAM)
```

### Inputs:

Inputs are usually two dataframes:

1) Gene data: a dataframe containing information about the genes from a species' pangenome that are present versus absent in each sample (host), typically filtered to remove core genes (e.g., omiting genes that are more than 90% present across samples). Gene presence/absence can be estimated with bioiformatics tools, such as MIDAS. This will be a samples by genes matrix and must contain the column sample_name to indicate the sample name. Users will read in gene data from a file, in this example simulated gene presence/absence. 
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/genes.png" width=400>

2) Metadata: a matrix of samples by covariates and the phenotype of interest in this case y, sample names should match gene data or be filtered to match on a species by species basis.
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/metadata.png" width=400>

Example of how to import data shown below. These data are simulated output from MIDAS 3, and associated simulated metadata. These data were modeled to have a strain that is correlated with the phenotype y, and a strain that is uncorrelated as well as 3 genes that are more associated to y than the phenotype. Age was randomly simulated to be a covariate.  

```
library(tidyverse)
library(magrittr)
n_tau = 100 ### number of permutations for tau test
exp_metadata = read.csv("../example_data/exp_metadata.csv") ### read in example metadata
## metadata is a sample by different covariates and phenotype matrix
exp_genedata = read.csv("../example_data/genepresabs.csv") ### read in example gene data
## gene data is a sample by gene matrix for one species
```

### Calculate Genetic Relatedness Matrix (GRM)

From the gene matrix imported above make the GRM. This GRM represents how similar our samples are within the space of the genes of this bacterial species but other representations should be valid. We will use this to estimate our random effects for our mixed model and this is used to represent the population structure. The gene matrix must contain the column sample_name for this to run properly.  
```
GRM = calculate_grm(exp_genedata)
```
Example of GRM:

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/GRM.png" width=500>

#### Visualization of the GRM with the strain information used to generate it are labeled.
This gene data was generated with a strain that is correlated to y in half of the samples and another strain or subset of semi correlated genes that are not correlated to the hypothetical phenotype of interest. We are most interested in genes that are able to explain our phenotype more than the simulated strain, in this case 3 genes were simulated to be more related to the phenotype than the strain.

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/exampleGRM.png" width=600>

### $\tau$ test for population structure (strain-trait associations)
Fit baseline glm for starting parameters in tau test, this can be done as you would do a normal glm for your data.

```
glm_fit0 = glm("y~age+1", data = exp_metadata, family = "binomial")
```

Fit tau test using baseline glm and GRM calculated above

```
glmm_fit = fit_tau_test(glm_fit0, GRM,species_id = "test",verbose = FALSE,log_file = NA)
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
Phi:  1 if logit or binomial should be 1
T value of tau: 0.624
Number of Samples: 100
```

The outputs of this can be interpreted as the microSLAM was able to converge, the $\tau$ variable is 2.35, and we have now estimated our random variable b as well as the coefficients of the covariates.

Test the significance of the tau that was fit with a permutation test

```
tautestfit = run_tau_test(glm_fit0, GRM,n_tau,species_id = "test", tau0=1, phi0=1)
```

Calculate the pvalue from the permutation test ran on tau

```
pvalue = sum(tautestfit$t >= glmm_fit$t)/n_tau
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/permutationnew.png" with=400>

In this case compared to 100 permutations of the tau test we have the most extreme example and a significant pvalue. Therefore it makes sense to use the microSLAM mixed model to find the betas for the gene-by-gene-trait associations

### $\beta$ test for gene-trait associations

Transform the gene data to a long matrix to fit each gene separately

```
gene_long = exp_genedata %>%
    pivot_longer(cols = starts_with("gene"),
    names_to = "gene_id",
    values_to = "gene_value") %>%
    as.data.frame()
```

Fit a beta for each gene adjusting for the population structure

```  
gene_test_df = fit_beta(glmm_fit,glm_fit0,GRM,gene_long,SPA=TRUE)
```

Plot a simple volcano plot of the results, in this case genes 1, 2, and 3, have been generated as being more related to the phenotype than the simulated strain. These are colored in red.

```
ggplot(gene_test_df,aes(beta,-log10(SPA_pvalue)))+
 geom_point(data = gene_test_df[ which(gene_test_df$SPA_pvalue >= .005),], color = 'gray80', size=2)+
 geom_point(data = gene_test_df[which(gene_test_df$SPA_pvalue <= .005),], color = 'red', size = 2)+
 theme_minimal()
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/volcanoplotcolor.png?raw=true" with=400>


Example output dataframe for $\beta$ test

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/betadf1.png?raw=true">

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/betadf2.png?raw=true">


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

If you want to visualize the genes that are correlated with the strain you can filter the genes to those that are > than some cut off correlated to b. Example of this is below.


```
strain_genes = gene_test_df %>%
 filter(abs(cor_to_b) > .4)
rownames(exp_genedata) = exp_genedata$sample_name
strain_genes_values = exp_genedata[,strain_genes$gene_id]

pheatmap(strain_genes_values,
show_rownames=FALSE,
show_colnames=FALSE,
treeheight_row=0,
treeheight_col = 0,
labels_row = "samples",
labels_col = "samples",
main = paste("Strain Associated Genes"),
border_color = NA,
annotation_row = exp_metadata[,-5],
color = myColor,
clustering_method = "average")
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/strainheatmap.png?raw=true">

300 genes were modeled from the strain, and we are able to recover 288 of these through this cut off. These can be used to find which genes make up substrains related to the phenotype of interest.
