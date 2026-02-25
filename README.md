# mcnemarControl

[![R-CMD-check](https://github.com/PHUoL/pairwisemcnemar/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/PHUoL/pairwisemcnemar/actions/workflows/R-CMD-check.yaml)

Control-vs-treatment McNemar tests from long-format paired binary data.

## Quick start

```r
library(mcnemarControl)

dat <- mcnemar_example_long(n = 40, seed = 1)
fit <- mcnemar_control(dat, id = "id", condition = "condition", outcome = "outcome",
                       control = "Control")
fit
```
