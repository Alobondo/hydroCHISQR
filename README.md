# hydroCHISQR
Chi - Squared test for hydrological frequency analysis

# Requirements
Dependencies:
  grDevices, ggplot2, reshape2, stats, ggpubr

# Installation
You can install the development version of **hydroCHISQR** from GitHub with this R command:
```
# install.packages("remotes")
remotes::install_github("Alobondo/hydroCHISQR")
library(hydroCHISQR)
```

# Usage
Functions | Description |
--- | --- |
```hydroCHISQR(df, nc, dist_param, alpha, c_test)``` | Apply the Chi - Squared test for hydrological frequency analysis. |
```df``` | Data frame with observations and distributions adjusted. |
```nc``` | Number of classes, default 1 for Sturges, 2 for Scott, 3 for Freedman-Diaconis. |
```dist_param``` | Number of parameters of the distributions in df. |
```alpha``` | Significance level, typical 0.05. |
```c_test``` | Statistical result or plot, default 1 for statistical and 2 for plot. |

```hydroCHISQR::example_2_params``` Example of frequency analysis result for 2 parameters distributions, with first column as observations.
```hydroCHISQR::example_3_params``` Example of frequency analysis result for 3 parameters distributions, with first column as observations.

# Reporting bugs
If you find an error in some function, want to report a typo in the documentation or submit a recommendation, you can do it [here](https://github.com/Alobondo/hydroCHISQR/issues)

# Keywords
Hydrology, R package, Frequency analysis
