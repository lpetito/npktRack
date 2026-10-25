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
remotes::install_github("lpetito/npktRack")
```

## Example

```r
library(npktRack)

set.seed(123)
id  <- rep(1:100, each = 5)
age <- rep(seq(5, 13, by = 2), times = 100)
y   <- rbinom(500, 1, plogis(-3 + 0.3 * age))

res <- npktRack(id = id, t = age, y = y,
                tmin = 5, tmax = 9, grid = 2,
                delta = c(2, 4), b = 2, nboot = 50)
res$mrtp
res$grtp
```

Note that `tmax + max(delta)` must not exceed the oldest observed time (here 9 + 4 = 13).

## Parallel computing

The bootstrap uses [future](https://future.futureverse.org/), and by default runs sequentially. To use more cores, either set a plan yourself before calling `npktRack()`:

```r
future::plan(future::multisession, workers = 4)
res <- npktRack(...)
future::plan(future::sequential)  # optional: switch back when done
```

or pass `workers` directly, which sets a temporary plan for that call only:

```r
res <- npktRack(..., workers = 4)
```
