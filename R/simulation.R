source("R/methods.R")

TRUE_BW <- 0.30

# --- SimDesign functions ---

Generate <- function(condition, fixed_objects) {
  K      <- condition$K
  dist   <- condition$dist
  beta_x <- condition$beta_x
  J      <- condition$J
  n_j    <- condition$n_j
  ICC    <- condition$ICC
  N      <- J * n_j
  tau    <- sqrt(ICC / (1 - ICC))

  cluster <- factor(rep(seq_len(J), each = n_j))
  x_star  <- rnorm(N)
  w       <- rep(rnorm(J), each = n_j)
  u       <- rep(rnorm(J, 0, tau), each = n_j)
  e       <- rnorm(N)
  y       <- beta_x * x_star + TRUE_BW * w + u + e

  probs <- switch(dist,
    symmetric = seq(0, 1, length.out = K + 1)[2:K],
    skewed    = ((1:(K - 1)) / K)^2,
    polarized = 0.5 - 0.4 * cos(pi * (1:(K - 1)) / K))
  x_ord <- as.integer(cut(x_star, c(-Inf, quantile(x_star, probs), Inf)))

  data.frame(y = y, x_star = x_star, x_ord = x_ord, w = w, cluster = cluster)
}

Analyse <- function(condition, dat, fixed_objects) {
  bind_rows(map(ALL_METHODS, ~ .x(dat)))
}

Summarise <- function(condition, results, fixed_objects) {
  truth <- condition$beta_x

  results |>
    filter(converged) |>
    group_by(method) |>
    summarise(
      n_ok          = n(),
      mean_bx       = mean(bx),
      bias          = mean(bx) - truth,
      rel_bias      = if (abs(truth) > 1e-10) (mean(bx) - truth) / truth else NA_real_,
      rmse          = sqrt(mean((bx - truth)^2)),
      emp_se        = sd(bx),
      mean_model_se = mean(se),
      se_ratio      = mean(se) / sd(bx),
      coverage      = mean(ci_lo <= truth & ci_hi >= truth),
      rej_rate      = mean(p < 0.05),
      mean_ci_width = mean(ci_hi - ci_lo),
      .groups       = "drop"
    )
}
