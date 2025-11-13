# Propensity-score-matching-with-clustered-data

This repository contains supplemental code and information for the article by B. Arpino and M. Cannas, *"Propensity score matching with clustered data. An application to the estimation of the impact of caesarean section on the Apgar score"* published in [*Statistics in Medicine*](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6880) (2016).


Framework: causal inference for a binary treatment using observational data with both individual and group-level covariates. 
In the article, we compare the performance of the following algorithms aimed at reducing bias from unobserved cluster-level covariates:

a) inclusion of fixed or random effects in the propensity score model
b) pure within-cluster matching
c) our proposal: ‘preferential’ within-cluster matching. This approach first searches for control units to be matched to treated units within the same cluster. If matching is not possible within-cluster, then the algorithm searches in other clusters. 

All algorithms can be implemented in the R environment using the [CMatching](https://cran.r-project.org/web/packages/CMatching/index.html) package.

**Code for Empirical Analysis**

The article presents a case study aimed at estimating the effect of caesarean section on the Apgar score using birth register data from Sardinian hospitals.
The code for the empirical analysis is provided in this [file](./code_for_empirical_analysis.R).
This script relies only on the widely used [Matching](https://cran.r-project.org/web/packages/Matching/index.html) package and may therefore be preferred by researchers who wish to reproduce the analysis with greater flexibility or adapt the procedures beyond the predefined functions available in CMatching.


**Data Availability Statement**

The empirical analysis presented in the paper is based on data provided by the Osservatorio Epidemiologico of Regione Sardegna (IT). We thanks the Osservatorio Eidemiologico for kindly providing data. Access to these data was granted under a confidentiality agreement, and therefore they cannot be publicly shared.



