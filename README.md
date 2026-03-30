# Ordinal Predictors in Multilevel Models: Simulation Study

Companion code for *"How Should Ordinal Predictors Be Handled in Multilevel Models? A Comprehensive Simulation Comparison with Practical Recommendations."*

## Overview

This simulation compares eight methods for handling ordinal predictors in two-level random-intercept models across a 3^6 = 729 factorial design (1,000 replications each).

### Methods

| Method | Description |
|--------|-------------|
| Oracle | True latent X* (upper bound) |
| Naive | Standardized integer coding |
| Dummy | K-1 indicators, implied spacing, refit |
| Isotonic | Monotone category means (Barlow et al., 1972) |
| PenOrd | Difference penalty (Gertheiss & Tutz, 2009) |
| PolyCor | Polyserial correction (Olsson et al., 1982) |
| LCP | Latent-Calibrated Predictor (novel) |
| IC | Isotonic-Calibrated Predictor (novel) |

### Design Factors

| Factor | Levels |
|--------|--------|
| Categories (K) | 3, 5, 7 |
| Threshold distribution | Symmetric, Skewed, Polarized |
| True effect (beta_x) | 0, 0.20, 0.50 |
| Clusters (J) | 20, 50, 100 |
| Cluster size (n_j) | 5, 10, 20 |
| ICC | .05, .15, .30 |

### Performance Measures

Bias, relative bias, RMSE, SE ratio, 95% CI coverage, rejection rate (Type I error / power), and relative efficiency vs. Oracle.

## Requirements

```r
install.packages(c("tidyverse", "lme4", "lmerTest", "SimDesign"))
```

## Usage

```bash
# Full run (729 conditions x 1000 reps)
Rscript run_sim.R

# Pilot run (729 conditions x 50 reps)
Rscript run_sim.R --pilot
```

## Project Structure

```
R/
  helpers.R      # Shared utilities (fit_mlm, polyserial, isotonic, LCP scores)
  methods.R      # Eight method implementations
  simulation.R   # SimDesign Generate/Analyse/Summarise functions
run_sim.R        # Entry point
results/         # Output (created at runtime)
```

## Output

- `results/sim_summary.csv` -- Aggregated performance measures by condition and method
- `results/ordinal-sim.rds` -- Full SimDesign results object

## License

MIT
