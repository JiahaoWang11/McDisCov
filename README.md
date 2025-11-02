
# McDisCov

Author: Jiahao Wang (<10225000478@stu.ecnu.edu.cn>)

**McDisCov** identifies taxa whose **absolute** and **relative**
abundances are associated with covariates of interest.  
It supports both continuous and categorical variables, and provides
flexible tests for arbitrary weighted taxa groups.

## Installation

You can install the development version from GitHub:

``` r
# Install devtools if needed
install.packages("devtools")

# Install McDisCov from GitHub
devtools::install_github("JiahaoWang11/McDisCov")
```

## Quick usage

``` r
# Install devtools if needed
library(McDisCov)
McDisCov(
  taxa_data, 
  meta_data,
  variables,
  imputation = FALSE,
  zero_impute = 0.5,
  zero_ratio_threshold = 1,
  alpha = 0.05,
  Wc = 0.7,
  nonzero_cut = 3,
  lib_cut = 0,
  cov_shrinkage = "diag",
  lambda_seq = seq(0, 1, by = 0.1)
)

gclr.test(model, variables, w1, w2)
```
