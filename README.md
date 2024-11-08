# changepointsPSI
Changepoint algorithms and post-selection inference for changepoints in univariate data. Deals with piecewise constant mean, piecewise linear mean, and piecewise constant variance cases.

********************************************************************************************************************************************

## Installation

To install this package:
```
devtools::install_github("rachelcarrington/changepointsPSI")
```

Inference for the change in variance model requires the package `DirichletReg`
```
install.packages("DirichletReg")
```


To use L0 inference, you will also need to install the package `ChangepointInference` (see https://github.com/jewellsean/ChangepointInference):
```
devtools::install_github("jewellsean/ChangepointInference")
```

To use PELT, you will also need to install the package `changepoint`
```
install.packages("changepoint")
```

To use seeded binary segmentation, you will need to download `seedBS.R` from https://github.com/kovacssolt/SeedBinSeg.

********************************************************************************************************************************************
## Related papers

__Improving power by conditioning on less in post-selection inference for changepoints.__ Rachel Carrington and Paul Fearnhead (2023)
https://arxiv.org/pdf/2301.05636

__Post-selection inference for quantifying uncertainty in changes in variance.__ Rachel Carrington and Paul Fearnhead (2024)
https://arxiv.org/pdf/2405.15670

********************************************************************************************************************************************

## Changepoint algorithms
The following changepoint algorithms are included in this package:
* binary segmentation: `binary_segmentation`
* wild binary segmentation: `wild_binary_segmentation`
* seeded binary segmentation: `wild_binary_segmentation` with `seeded = TRUE`
* narrowest over threshold: `narrowest_over_threshold`

### Examples:
```
x <- c(rep(1, 100), rep(-1, 100)) + rnorm(200)
results_bs <- binary_segmentation(x, threshold=4)
results_wbs <- wild_binary_segmentation(x, threshold=4, num_rand_samples=200)
results_sbs <- wild_binary_segmentation(x, threshold=4, seeded=TRUE)
results_not <- narrowest_over_threshold(x, threshold=4, num_rand_samples=200)
```

The function `find_changepoints` is a wrapper function for these algorithms, and also PELT (as implemented in the `changepoint` package). This should be used if you want to calculate p-values for the changepoints --
the output of this function is given to the function used to calculate p-values.
* binary segmentation: `method = "bs"`
* wild (or seeded) binary segmentation: `method = "wbs"`
* narrowest over threshold: `method = "not"`
* PELT: `method = "pelt"` -- for change in mean or variance only; the `changepoint` package is required for this.

```
results_bs <- find_changepoints(x, method="bs", params=list(threshold=4))
results_wbs <- find_changepoints(x, method="wbs", params=list(threshold=4, num_rand_samples=200))
results_sbs <- find_changepoints(x, method="sbs", params=list(threshold=4, seeded=TRUE))
results_not <- find_changepoints(x, method="not", params=list(threshold=4, num_rand_samples=200))
results_pelt <- find_changepoints(x, method="pelt", params=list(penalty="Manual", pen.value=10))
```
********************************************************************************************************************************************

## Post-selection inference

To calculate p-values for detected changepoints, use:
* `calculate_pvals_all` for i.i.d. data in cases where the p-values can be calculated exactly (change in mean or slope; change in variance with CUSUM loss)
* `calculate_pvals_var` for approximate p-values in the change in variance model

e.g.
```
x <- c(rep(1, 100), rep(-1, 100)) + rnorm(200)
results_bs <- find_changepoints(x, method="bs", params=list(threshold=4))
calculate_pvals_all(results_bs, h=20, sigma2=1, return_probs=TRUE)

x <- c(rnorm(100), rnorm(100, sd=2))
results_cusum <- find_changepoints(x, method="bs", params=list(threshold=10, loss="cusum"), model="var")
calculate_pvals_all(results_cusum, h=20, return_probs=TRUE)

results_lrs <- find_changepoints(x, method="pelt", params=list(penalty="Manual", pen.value=10), model="var")
calculate_pvals_var(results_lrs, h=20, N_sample=50)
```
