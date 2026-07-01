# npktRack

`npktRack` provides nonparametric kernel-based estimators of the **rank‑tracking probability (RTP)** and the **rank‑tracking probability ratio (RTPR)** for longitudinal binary outcomes, following Wu et al. (2020). It implements:

- local estimands: RTP(t, t + δ) and RTPR(t, t + δ),
- partially global estimands: mRTP(δ) and mRTPR(δ),
- fully global estimands: gRTP and gRTPR,

using Epanechnikov kernels and subject‑level bootstrap for uncertainty.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("yourusername/npktRack")
```
