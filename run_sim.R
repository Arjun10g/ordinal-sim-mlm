source("R/simulation.R")

args  <- commandArgs(trailingOnly = TRUE)
PILOT <- "--pilot" %in% args

N_REPS <- if (PILOT) 50L else 1000L

design <- crossing(
  K      = c(3L, 5L, 7L),
  dist   = c("symmetric", "skewed", "polarized"),
  beta_x = c(0, 0.20, 0.50),
  J      = c(20L, 50L, 100L),
  n_j    = c(5L, 10L, 20L),
  ICC    = c(0.05, 0.15, 0.30)
) |> as.data.frame()

cat(sprintf("Design: %d conditions x %d reps = %s fits\n",
            nrow(design), N_REPS,
            format(nrow(design) * N_REPS * 8L, big.mark = ",")))

res <- runSimulation(
  design       = design,
  replications = N_REPS,
  generate     = Generate,
  analyse      = Analyse,
  summarise    = Summarise,
  seed         = rep(20260328L, nrow(design)),
  parallel     = TRUE,
  save         = TRUE,
  save_results = TRUE,
  filename     = "results/ordinal-sim",
  packages     = c("tidyverse", "lme4", "lmerTest"),
  progress     = TRUE
)

write_csv(res, "results/sim_summary.csv")
cat("Done. Results saved to results/sim_summary.csv\n")
