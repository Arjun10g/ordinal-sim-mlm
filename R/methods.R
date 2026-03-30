source("R/helpers.R")

method_oracle <- function(dat) {
  dat$xv <- dat$x_star
  fit_mlm(dat) |> mutate(method = "Oracle")
}

method_naive <- function(dat) {
  dat$xv <- std(as.numeric(dat$x_ord))
  fit_mlm(dat) |> mutate(method = "Naive")
}

method_dummy <- function(dat) {
  cats <- sort(unique(dat$x_ord))
  K    <- length(cats)
  if (K < 2) return(tibble(bx = NA_real_, se = NA_real_, ci_lo = NA_real_,
                            ci_hi = NA_real_, p = NA_real_, converged = FALSE, method = "Dummy"))
  dat$xf <- factor(dat$x_ord)
  tryCatch(suppressWarnings(suppressMessages({
    m1 <- lmer(y ~ xf + w + (1 | cluster), data = dat, REML = TRUE)
    fe <- fixef(m1)
    ce <- c(0, fe[grep("^xf", names(fe))])
    dat$xv <- std(ce[dat$x_ord - min(cats) + 1])
    fit_mlm(dat) |> mutate(method = "Dummy")
  })), error = function(e) {
    tibble(bx = NA_real_, se = NA_real_, ci_lo = NA_real_,
           ci_hi = NA_real_, p = NA_real_, converged = FALSE, method = "Dummy")
  })
}

method_isotonic <- function(dat) {
  iso_pred <- compute_isotonic(dat$y, dat$x_ord, dat$cluster)
  dat$xv <- if (sd(iso_pred) < 1e-10) std(as.numeric(dat$x_ord)) else std(iso_pred)
  fit_mlm(dat) |> mutate(method = "Isotonic")
}

method_penord <- function(dat, lambda = 1.0) {
  cats <- sort(unique(dat$x_ord))
  K    <- length(cats)
  D    <- K - 1
  if (D < 1) return(tibble(bx = NA_real_, se = NA_real_, ci_lo = NA_real_,
                            ci_hi = NA_real_, p = NA_real_, converged = FALSE, method = "PenOrd"))

  dm <- outer(dat$x_ord, cats[-1], "==") * 1.0
  cl <- as.integer(dat$cluster)
  yc <- dat$y; dc <- dm
  for (j in seq_len(max(cl))) {
    m <- cl == j
    yc[m] <- yc[m] - mean(yc[m])
    for (d in seq_len(D)) dc[m, d] <- dc[m, d] - mean(dc[m, d])
  }

  P <- matrix(0, D, D)
  P[1, 1] <- 1
  if (D > 1) {
    for (k in 2:D) {
      P[k-1, k-1] <- P[k-1, k-1] + 1; P[k, k] <- P[k, k] + 1
      P[k-1, k] <- P[k-1, k] - 1; P[k, k-1] <- P[k, k-1] - 1
    }
  }

  b <- tryCatch(solve(crossprod(dc) + lambda * P, crossprod(dc, yc)), error = function(e) NULL)
  if (is.null(b)) return(tibble(bx = NA_real_, se = NA_real_, ci_lo = NA_real_,
                                ci_hi = NA_real_, p = NA_real_, converged = FALSE, method = "PenOrd"))

  ce <- c(0, as.numeric(b))
  dat$xv <- std(ce[dat$x_ord - min(cats) + 1])
  fit_mlm(dat) |> mutate(method = "PenOrd")
}

method_polycor <- function(dat) {
  dat$xv <- std(as.numeric(dat$x_ord))
  res <- fit_mlm(dat)
  if (!res$converged || is.na(res$bx)) return(res |> mutate(method = "PolyCor"))

  rps <- compute_polyserial(dat$x_ord)
  if (rps < 0.01 || rps > 0.999) return(res |> mutate(method = "PolyCor"))
  corr <- 1 / rps
  tibble(bx = res$bx * corr, se = res$se * corr,
         ci_lo = res$ci_lo * corr, ci_hi = res$ci_hi * corr,
         p = res$p, converged = TRUE, method = "PolyCor")
}

method_lcp <- function(dat) {
  sc <- compute_latent_scores(dat$x_ord)
  dat$xv <- as.numeric(sc[as.character(dat$x_ord)])
  fit_mlm(dat) |> mutate(method = "LCP")
}

method_ic <- function(dat) {
  iso_pred <- compute_isotonic(dat$y, dat$x_ord, dat$cluster)
  if (sd(iso_pred) < 1e-10) {
    dat$xv <- std(as.numeric(dat$x_ord))
    return(fit_mlm(dat) |> mutate(method = "IC"))
  }
  dat$xv <- std(iso_pred)
  res <- fit_mlm(dat)
  if (!res$converged || is.na(res$bx)) return(res |> mutate(method = "IC"))

  K    <- length(unique(dat$x_ord))
  nper <- nrow(dat) / K
  rps  <- compute_polyserial(dat$x_ord)
  if (rps < 0.01 || rps > 0.999) return(res |> mutate(method = "IC"))

  alpha <- (1 / K) * (1 + nper / (nper + K))
  corr  <- (1 / rps)^alpha
  tibble(bx = res$bx * corr, se = res$se * corr,
         ci_lo = res$ci_lo * corr, ci_hi = res$ci_hi * corr,
         p = res$p, converged = TRUE, method = "IC")
}

ALL_METHODS <- list(method_oracle, method_naive, method_dummy, method_isotonic,
                    method_penord, method_polycor, method_lcp, method_ic)
