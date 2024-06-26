% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmm_functions.R
\name{fit_beta}
\alias{fit_beta}
\title{fit_beta}
\usage{
fit_beta(obj.pop.strut, glm_fit0, GRM, gene_df, SPA = FALSE)
}
\arguments{
\item{obj.pop.strut}{output of fit_tau_test; GLMM of species with GRM accounted for}

\item{glm_fit0}{glm model. Model output with no sample relatedness accounted for}

\item{GRM}{Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness}

\item{gene_df}{long data frame with, gene_id, sample_name, and gene_value}

\item{SPA}{whether to run Saddle point approximation for pvalues (will slow down output)}
}
\value{
model output for the baseline structure glmm
}
\description{
Fit a beta for each genes including the population structure model with the random effects
}
