#' Nonparametric kernel-based tracking estimands (RTP/RTPR)
#'
#' Estimate local and global rank-tracking probabilities (RTP) and
#' rank-tracking probability ratios (RTPR) using kernel smoothing,
#' as in Wu et al. (2020).
#'
#' @param id Vector of subject identifiers (character or coercible to character).
#' @param t Numeric vector of observation times (e.g. age).
#' @param y Numeric vector of binary outcomes coded 0/1, indicating membership
#'   in the high-risk set \eqn{A(t)} at time \eqn{t}.
#' @param tmin Minimum time in the domain over which to estimate tracking (start of grid).
#' @param tmax Maximum time in the domain over which to estimate tracking (end of grid).
#' @param grid Spacing of the time grid between \code{tmin} and \code{tmax}
#'   (same units as \code{t}); RTP/RTPR will be evaluated at \code{seq(tmin, tmax, by = grid)}.
#' @param delta Numeric vector of positive lag times (same units as \code{t})
#'   at which to evaluate tracking, e.g. \code{c(1, 2, 3, 4)}.
#' @param b Kernel bandwidth (same units as \code{t}); controls smoothing in both
#'   1D (marginal) and 2D (joint) kernel estimates.
#' @param nboot Number of subject-level bootstrap resamples to use for
#'   uncertainty quantification. The first iteration uses the original data.
#'
#' @details
#' The function implements nonparametric kernel estimators of:
#' \itemize{
#'   \item Local tracking probabilities RTP(t, t + \eqn{\delta})
#'   \item Local tracking ratios RTPR(t, t + \eqn{\delta})
#'   \item Partially global tracking indices mRTP(\eqn{\delta}) and mRTPR(\eqn{\delta}),
#'         obtained by integrating RTP/RTPR over \code{t} in \code{[tmin, tmax]}
#'   \item Global tracking indices gRTP and gRTPR, obtained by integrating
#'         mRTP/mRTPR over the lag domain \code{delta}
#' }
#' Estimation uses Epanechnikov kernels with subject-level weighting, and
#' subject-level bootstrap to obtain percentile confidence intervals.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{rtp}{List of RTP estimates on the time grid:
#'     \code{$rtp} (point), \code{$rtp.lci}, \code{$rtp.uci}, and
#'     \code{$rtpci} (formatted strings). Rows = \code{delta}, columns = time grid.}
#'   \item{rtpr}{Analogous list of RTPR estimates.}
#'   \item{mrtp_raw}{Data frame with mRTP/mRTPR and bootstrap CIs for each \code{delta}.}
#'   \item{grtp_raw}{Data frame with gRTP/gRTPR and bootstrap CIs.}
#'   \item{mrtp}{Compact table with \code{delta} and formatted mRTP/mRTPR CIs.}
#'   \item{grtp}{Compact table with formatted gRTP/gRTPR CIs.}
#'   \item{n.obs}{Matrix of subject counts contributing at each (t, t + \eqn{\delta}).}
#'   \item{n.y1}{Matrix of subjects with \code{y = 1} near both t and t + \eqn{\delta}.}
#'   \item{mrtp_plot}{\code{ggplot2} object for mRTP vs. \code{delta}.}
#'   \item{mrtpr_plot}{\code{ggplot2} object for mRTPR vs. \code{delta}.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' id  <- rep(1:100, each = 5)
#' age <- rep(seq(5, 13, by = 2), times = 100)
#' y   <- rbinom(500, 1, plogis(-3 + 0.3 * age))
#'
#' res <- npktRack(
#'   id    = id,
#'   t     = age,
#'   y     = y,
#'   tmin  = 5,
#'   tmax  = 11,
#'   grid  = 2,
#'   delta = c(2, 4),
#'   b     = 2,
#'   nboot = 50
#' )
#' res$grtp
#' }
#'
#' @export
npktRack <- function(id,
                     t,
                     y,
                     tmin,
                     tmax,
                     grid,
                     delta,
                     b,
                     nboot) {

  ##------- parallel setup -------##
  workers <- max(1L, parallel::detectCores() - 1L)
  future::plan(future.callr::callr, workers = workers)

  ##------- kernels -------##
  Kh <- function(x, b) {
    u <- x / b
    w <- 0.75 * (1 - u^2) / b
    w[abs(u) >= 1] <- 0
    w
  }

  Kh2D <- function(x1, x2, b1, b2) {
    u1 <- x1 / b1
    u2 <- x2 / b2
    w1 <- 0.75 * (1 - u1^2) / b1
    w2 <- 0.75 * (1 - u2^2) / b2
    w1[abs(u1) >= 1] <- 0
    w2[abs(u2) >= 1] <- 0
    w1 * w2
  }

  NW.Kernel1D <- function(id, xvec, yvec, x0, b) {
    d <- data.frame(id = id, t = xvec, y = yvec)
    by_id <- split(d, d$id)

    num_sum <- 0
    den_sum <- 0

    for (di in by_id) {
      wi <- Kh(di$t - x0, b)
      sw <- sum(wi)
      if (sw > 0) {
        ni <- nrow(di)
        num_sum <- num_sum + sum(wi * di$y) / ni
        den_sum <- den_sum + sw / ni
      }
    }
    if (den_sum > 0) num_sum / den_sum else NA_real_
  }

  NW.Kernel2D <- function(id, xvec, yvec, t1, t2, b1, b2) {
    d <- data.frame(id = id, t = xvec, y = yvec)
    by_id <- split(d, d$id)

    num_sum <- 0
    den_sum <- 0

    for (di in by_id) {
      ni <- nrow(di)
      if (ni < 2L) next

      mi <- ni * (ni - 1) / 2

      s1 <- di[abs(di$t - t1) <= b1, c("t", "y")]
      s2 <- di[abs(di$t - t2) <= b2, c("t", "y")]
      if (nrow(s1) == 0L || nrow(s2) == 0L) next

      idx <- which(outer(s2$t, s1$t, ">"), arr.ind = TRUE)
      if (nrow(idx) == 0L) next

      u1 <- s1$t[idx[, 2]] - t1
      u2 <- s2$t[idx[, 1]] - t2
      yy <- s1$y[idx[, 2]] * s2$y[idx[, 1]]
      w  <- Kh2D(u1, u2, b1, b2)

      den_i <- sum(w)
      if (den_i <= 0) next
      num_i <- sum(w * yy)

      num_sum <- num_sum + (num_i / mi)
      den_sum <- den_sum + (den_i / mi)
    }

    if (den_sum > 0) num_sum / den_sum else NA_real_
  }

  mgAUC <- function(X, TT) {
    if (length(TT) == 1L) return(X[1])
    idx <- which(!is.na(X))
    if (length(idx) < 2L) return(NA_real_)
    X2  <- X[idx]
    TT2 <- TT[idx]
    dt  <- diff(TT2)
    avg_h <- (head(X2, -1) + tail(X2, -1)) / 2
    sum(avg_h * dt) / sum(dt)
  }

  ## ------ checks -------- ##
  if (!all(y %in% c(0, 1))) stop("Outcome y must be strictly binary (0/1).")
  if (tmin >= tmax) stop("tmin must be less than tmax.")
  if (tmin < min(t) || (tmax + max(delta)) > max(t))
    stop("t-range plus max(delta) must lie inside observed time domain.")
  if (grid > min(delta)) stop("grid must be smaller than min(delta).")
  if (any(delta <= 0)) stop("All delta values must be positive.")
  if (nboot < 2L) stop("nboot must be at least 2 to compute bootstrap CIs.")

  ## ------ organize inputs -------- ##
  id    <- as.character(id)
  t     <- as.numeric(t)
  y     <- as.numeric(y)
  tmin  <- as.numeric(tmin)
  tmax  <- as.numeric(tmax)
  delta <- sort(unique(as.numeric(delta)))
  grid  <- as.numeric(grid)
  b     <- as.numeric(b)

  alltime <- as.numeric(seq(tmin, tmax, by = grid))

  df <- data.frame(
    id = id,
    t  = t,
    y  = y,
    stringsAsFactors = FALSE
  )

  ## --- bootstrap loop --- ##
  results <- future.apply::future_lapply(
    X = seq_len(nboot),
    FUN = function(nb) {

      if (nb == 1L) {
        bd <- df
        bd$id2 <- bd$id
      } else {
        unique_ids <- unique(df$id)
        bID <- data.frame(
          id   = sample(unique_ids,
                        size = length(unique_ids),
                        replace = TRUE),
          boot = seq_along(unique_ids)
        ) |>
          dplyr::mutate(id2 = paste0(id, "_", boot)) |>
          dplyr::select(-boot)

        bd <- merge(bID, df, all.x = TRUE)
      }

      probS12 <- matrix(NA_real_, nrow = length(delta), ncol = length(alltime))
      RTP     <- matrix(NA_real_, nrow = length(delta), ncol = length(alltime))
      RTPR    <- matrix(NA_real_, nrow = length(delta), ncol = length(alltime))
      mRTP    <- rep(NA_real_, length(delta))
      mRTPR   <- rep(NA_real_, length(delta))

      prob_t1 <- sapply(alltime,
                        function(x) NW.Kernel1D(bd$id2, bd$t, bd$y, x, b))

      for (j in seq_along(delta)) {
        d_j <- delta[j]

        t2_vec <- alltime + d_j
        prob_t2 <- sapply(
          t2_vec,
          function(x) NW.Kernel1D(bd$id2, bd$t, bd$y, x, b)
        )

        for (i in seq_along(alltime)) {
          probS12[j, i] <- NW.Kernel2D(
            id   = bd$id2,
            xvec = bd$t,
            yvec = bd$y,
            t1   = alltime[i],
            t2   = alltime[i] + d_j,
            b1   = b, b2 = b
          )
        }

        RTP[j, ]  <- ifelse(prob_t1 > 0, probS12[j, ] / prob_t1, NA_real_)
        RTP[j, ]  <- pmin(pmax(RTP[j, ], 0), 1)

        RTPR[j, ] <- ifelse(prob_t2 > 0, RTP[j, ] / prob_t2, NA_real_)

        mRTP[j]  <- mgAUC(X = RTP[j, ],  TT = alltime)
        mRTPR[j] <- mgAUC(X = RTPR[j, ], TT = alltime)
      }

      list(
        RTP   = RTP,
        RTPR  = RTPR,
        mRTP  = mRTP,
        mRTPR = mRTPR
      )
    },
    future.seed = TRUE
  )

  future::plan(future::sequential)

  bootRTP   <- lapply(results, `[[`, "RTP")
  bootRTPR  <- lapply(results, `[[`, "RTPR")
  bootmRTP  <- t(sapply(results, `[[`, "mRTP"))
  bootmRTPR <- t(sapply(results, `[[`, "mRTPR"))

  ## RTP output
  rtp.array <- array(
    unlist(lapply(bootRTP, as.matrix)),
    dim = c(length(delta), length(alltime), nboot)
  )

  RTP2 <- list(
    rtp     = as.data.frame(bootRTP[[1]]),
    rtp.lci = as.data.frame(apply(rtp.array[, , -1, drop = FALSE],
                                  c(1, 2), quantile, probs = 0.025, na.rm = TRUE)),
    rtp.uci = as.data.frame(apply(rtp.array[, , -1, drop = FALSE],
                                  c(1, 2), quantile, probs = 0.975, na.rm = TRUE))
  )

  RTP2[["rtpci"]] <- as.data.frame(
    matrix(
      sprintf("%.2f (%.2f, %.2f)",
              as.matrix(RTP2$rtp),
              as.matrix(RTP2$rtp.lci),
              as.matrix(RTP2$rtp.uci)),
      nrow = nrow(RTP2$rtp),
      ncol = ncol(RTP2$rtp),
      dimnames = dimnames(RTP2$rtp)
    )
  )

  ## RTPR output
  rtpr.array <- array(
    unlist(lapply(bootRTPR, as.matrix)),
    dim = c(length(delta), length(alltime), nboot)
  )

  RTPR2 <- list(
    rtpr     = as.data.frame(bootRTPR[[1]]),
    rtpr.lci = as.data.frame(apply(rtpr.array[, , -1, drop = FALSE],
                                   c(1, 2), quantile, probs = 0.025, na.rm = TRUE)),
    rtpr.uci = as.data.frame(apply(rtpr.array[, , -1, drop = FALSE],
                                   c(1, 2), quantile, probs = 0.975, na.rm = TRUE))
  )

  RTPR2[["rtprci"]] <- as.data.frame(
    matrix(
      sprintf("%.2f (%.2f, %.2f)",
              as.matrix(RTPR2$rtpr),
              as.matrix(RTPR2$rtpr.lci),
              as.matrix(RTPR2$rtpr.uci)),
      nrow = nrow(RTPR2$rtpr),
      ncol = ncol(RTPR2$rtpr),
      dimnames = dimnames(RTPR2$rtpr)
    )
  )

  colnames(RTP2$rtp)   <- colnames(RTP2$rtp.lci)   <-
    colnames(RTP2$rtp.uci)  <- colnames(RTP2$rtpci)   <- alltime
  colnames(RTPR2$rtpr) <- colnames(RTPR2$rtpr.lci) <-
    colnames(RTPR2$rtpr.uci) <- colnames(RTPR2$rtprci) <- alltime
  rownames(RTP2$rtp)   <- rownames(RTP2$rtp.lci)   <-
    rownames(RTP2$rtp.uci)  <- rownames(RTP2$rtpci)   <- delta
  rownames(RTPR2$rtpr) <- rownames(RTPR2$rtpr.lci) <-
    rownames(RTPR2$rtpr.uci) <- rownames(RTPR2$rtprci) <- delta

  ## mRTP and mRTPR
  mRTP2 <- data.frame(
    delta      = delta,
    mRTP       = bootmRTP[1, ],
    mRTP.lci   = apply(bootmRTP[-1, , drop = FALSE],  2, quantile,
                       prob = 0.025, na.rm = TRUE),
    mRTP.uci   = apply(bootmRTP[-1, , drop = FALSE],  2, quantile,
                       prob = 0.975, na.rm = TRUE),
    mRTPR      = bootmRTPR[1, ],
    mRTPR.lci  = apply(bootmRTPR[-1, , drop = FALSE], 2, quantile,
                       prob = 0.025, na.rm = TRUE),
    mRTPR.uci  = apply(bootmRTPR[-1, , drop = FALSE], 2, quantile,
                       prob = 0.975, na.rm = TRUE)
  )

  mRTP2$mRTPci <- sprintf(
    "%.2f (%.2f, %.2f)",
    mRTP2$mRTP,
    mRTP2$mRTP.lci,
    mRTP2$mRTP.uci
  )

  mRTP2$mRTPRci <- sprintf(
    "%.2f (%.2f, %.2f)",
    mRTP2$mRTPR,
    mRTP2$mRTPR.lci,
    mRTP2$mRTPR.uci
  )

  ## gRTP and gRTPR
  gRTP  <- apply(bootmRTP,  1, function(x) mgAUC(x, delta))
  gRTPR <- apply(bootmRTPR, 1, function(x) mgAUC(x, delta))

  gRTP2 <- data.frame(
    gRTP      = gRTP[1],
    gRTP.lci  = stats::quantile(gRTP[-1],  prob = 0.025, na.rm = TRUE),
    gRTP.uci  = stats::quantile(gRTP[-1],  prob = 0.975, na.rm = TRUE),
    gRTPR     = gRTPR[1],
    gRTPR.lci = stats::quantile(gRTPR[-1], prob = 0.025, na.rm = TRUE),
    gRTPR.uci = stats::quantile(gRTPR[-1], prob = 0.975, na.rm = TRUE)
  )

  gRTP2$gRTPci <- sprintf(
    "%.2f (%.2f, %.2f)",
    gRTP2$gRTP,
    gRTP2$gRTP.lci,
    gRTP2$gRTP.uci
  )

  gRTP2$gRTPRci <- sprintf(
    "%.2f (%.2f, %.2f)",
    gRTP2$gRTPR,
    gRTP2$gRTPR.lci,
    gRTP2$gRTPR.uci
  )

  ## sample sizes
  n.obs0 <- matrix(NA_integer_, nrow = length(delta), ncol = length(alltime))
  n.y10  <- matrix(NA_integer_, nrow = length(delta), ncol = length(alltime))

  for (j in seq_along(delta)) {
    d_j <- delta[j]
    for (i in seq_along(alltime)) {
      t1 <- alltime[i]
      t2 <- alltime[i] + d_j

      ids_with_both <- df |>
        dplyr::group_by(id) |>
        dplyr::summarise(
          has_t1 = any(abs(t - t1) <= b),
          has_t2 = any(abs(t - t2) <= b),
          .groups = "drop"
        ) |>
        dplyr::filter(has_t1 & has_t2)

      if (nrow(ids_with_both) == 0L) {
        n.obs0[j, i] <- 0L
        n.y10[j, i]  <- 0L
        next
      }

      n.obs0[j, i] <- nrow(ids_with_both)

      tmp <- df |>
        dplyr::filter(
          id %in% ids_with_both$id,
          abs(t - t1) <= b | abs(t - t2) <= b
        ) |>
        dplyr::group_by(id) |>
        dplyr::summarise(
          y1_t1 = any(abs(t - t1) <= b & y == 1),
          y1_t2 = any(abs(t - t2) <= b & y == 1),
          .groups = "drop"
        )

      n.y10[j, i] <- if (nrow(tmp) == 0L) 0L else sum(tmp$y1_t1 & tmp$y1_t2)
    }
  }
  colnames(n.obs0) <- alltime
  colnames(n.y10)  <- alltime
  rownames(n.obs0) <- delta
  rownames(n.y10)  <- delta

  ## plots
  mrtp_plot <- ggplot2::ggplot(mRTP2,
                               ggplot2::aes(x = delta, y = mRTP)) +
    ggplot2::geom_line() +
    ggplot2::geom_line(
      ggplot2::aes(y = mRTP.lci),
      linetype = 2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mRTP.uci),
      linetype = 2
    ) +
    ggplot2::annotate(
      "text",
      x = min(delta) + (max(delta) - min(delta)) * 0.1,
      y = min(mRTP2$mRTP.lci, na.rm = TRUE) * 1.05,
      label = paste0(
        "gRTP=", round(gRTP2$gRTP, 2), " (",
        round(gRTP2$gRTP.lci, 2), ", ",
        round(gRTP2$gRTP.uci, 2), ")"
      ),
      fontface = "bold"
    ) +
    ggplot2::labs(
      x = expression(delta ~ "(lag time)"),
      y = expression(mRTP(delta))
    ) +
    ggplot2::theme_bw()

  mrtpr_plot <- ggplot2::ggplot(mRTP2,
                                ggplot2::aes(x = delta, y = mRTPR)) +
    ggplot2::geom_line() +
    ggplot2::geom_line(
      ggplot2::aes(y = mRTPR.lci),
      linetype = 2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mRTPR.uci),
      linetype = 2
    ) +
    ggplot2::annotate(
      "text",
      x = min(delta) + (max(delta) - min(delta)) * 0.1,
      y = min(mRTP2$mRTPR.lci, na.rm = TRUE) * 1.05,
      label = paste0(
        "gRTPR=", round(gRTP2$gRTPR, 2), " (",
        round(gRTP2$gRTPR.lci, 2), ", ",
        round(gRTP2$gRTPR.uci, 2), ")"
      ),
      fontface = "bold"
    ) +
    ggplot2::labs(
      x = expression(delta ~ "(lag time)"),
      y = expression(mRTPR(delta))
    ) +
    ggplot2::theme_bw()

  list(
    rtp        = RTP2,
    rtpr       = RTPR2,
    mrtp_raw   = mRTP2,
    grtp_raw   = gRTP2,
    mrtp       = mRTP2 |>
      dplyr::select(delta, mRTPci, mRTPRci),
    grtp       = gRTP2 |>
      dplyr::select(gRTPci, gRTPRci),
    n.obs      = n.obs0,
    n.y1       = n.y10,
    mrtp_plot  = mrtp_plot,
    mrtpr_plot = mrtpr_plot
  )
}
