library(tidyverse)
library(lme4)
library(lmerTest)
library(SimDesign)

std <- function(x) {
  s <- sd(x)
  if (s < 1e-10) return(x - mean(x))
  (x - mean(x)) / s
}

fit_mlm <- function(dat) {
  tryCatch(suppressWarnings(suppressMessages({
    m  <- lmerTest::lmer(y ~ xv + w + (1 | cluster), data = dat, REML = TRUE)
    co <- coef(summary(m))
    bx <- co["xv", "Estimate"]
    se <- co["xv", "Std. Error"]
    p  <- co["xv", "Pr(>|t|)"]
    df <- co["xv", "df"]
    tc <- qt(0.975, df)
    tibble(bx = bx, se = se,
           ci_lo = bx - tc * se, ci_hi = bx + tc * se,
           p = p, converged = TRUE)
  })), error = function(e) {
    tibble(bx = NA_real_, se = NA_real_,
           ci_lo = NA_real_, ci_hi = NA_real_,
           p = NA_real_, converged = FALSE)
  })
}

compute_latent_scores <- function(x_ord) {
  cats <- sort(unique(x_ord))
  K    <- length(cats)
  cts  <- tabulate(x_ord, nbins = max(cats))[cats]
  cp   <- cumsum(cts) / length(x_ord)
  cp   <- pmin(pmax(cp[-K], 1e-4), 1 - 1e-4)
  th   <- qnorm(cp)
  te   <- c(-Inf, th, Inf)
  sc   <- map_dbl(seq_along(cats), ~ {
    d <- pnorm(te[.x + 1]) - pnorm(te[.x])
    if (d < 1e-15) 0 else (dnorm(te[.x]) - dnorm(te[.x + 1])) / d
  })
  set_names(sc, as.character(cats))
}

compute_polyserial <- function(x_ord) {
  cats <- sort(unique(x_ord))
  K    <- length(cats)
  cts  <- tabulate(x_ord, nbins = max(cats))[cats]
  cp   <- cumsum(cts) / length(x_ord)
  cp   <- pmin(pmax(cp[-K], 1e-4), 1 - 1e-4)
  sum(dnorm(qnorm(cp))) / sd(x_ord)
}

compute_isotonic <- function(y, x_ord, cluster) {
  cl <- as.integer(cluster)
  ya <- y
  for (j in seq_len(max(cl))) {
    m <- cl == j
    ya[m] <- ya[m] - mean(ya[m])
  }
  cm  <- tapply(ya, x_ord, mean)
  iso <- isoreg(as.numeric(names(cm)), cm)
  iso$yf[match(x_ord, as.integer(names(cm)))]
}
