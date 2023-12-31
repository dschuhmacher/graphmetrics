# graphmetrics

  <!-- badges: start -->  [![R-CMD-check](https://github.com/dschuhmacher/graphmetrics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dschuhmacher/graphmetrics/actions/workflows/R-CMD-check.yaml)  <!-- badges: end -->

Compute distances between graphs exactly or approximately and perform distance-based statistical analyses. 

Based on [Schuhmacher and Wirth (2023)](https://arxiv.org/abs/2308.12165).

Comments, suggestions and contributions are welcome!



## Installation

Install directly from GitHub by saying in R
```
remotes::install_github("dschuhmacher/graphmetrics")
```
For doing exact computations with `method = cplex_match`, you will have to install [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer) (there are free Academic and Community Editions) and the R packages `Rcplex` and `ROI.plugin.cplex`).



## Getting started

The workhorse function for computing distances is called `gdist`. Functions for performing statistical analyses based on distances are `anderson_anova` and `msm_levene`. More information can be obtained from the respective help pages once the package has been installed.

