# AdaptiveMCMCWorkflow

`AdaptiveMCMCWorkflow` is a teaching-oriented R package for adaptive
random-walk MCMC. The package is designed for presentation, study, and small
research demonstrations rather than for maximum software abstraction.

The package implements:

- standard random-walk Metropolis
- Adaptive Metropolis (AM)
- Robust Adaptive Metropolis (RAM)
- Delayed Rejection Adaptive Metropolis (DRAM)

The implementation is split deliberately:

- C++ via `Rcpp` and `RcppArmadillo` for the computational kernels
- R for wrappers, examples, diagnostics, summaries, and plots

The comments are written to follow the mathematical storyline of an adaptive
MCMC presentation:

1. fixed-kernel Metropolis is sensitive to proposal scale and covariance choice
2. adaptive methods change the proposal during the run
3. that flexibility helps efficiency, but it requires theoretical care
4. AM, RAM, and DRAM solve different practical problems

## Installation

From the package source directory:

```r
devtools::load_all()
```

Or from a built source tarball:

```r
install.packages(
  "AdaptiveMCMCWorkflow_0.1.0.9000.tar.gz",
  repos = NULL,
  type = "source"
)
```

Do not source files from `R/` directly. They are package source files, not
standalone scripts.

## Simplest run from this repo

```r
source("run_adaptive_mcmc_demo.R")
```

The root script `run_adaptive_mcmc_demo.R`
is the main step-by-step demonstration file. It shows:

- how to initialize the package
- how to choose a target distribution
- how to choose the starting Gaussian proposal
- how to run the four samplers
- how to compare them with tables and plots
- how to check empirical adaptive-MCMC validity diagnostics

## Simplest package calls

```r
library(AdaptiveMCMCWorkflow)

target <- make_correlated_gaussian_target()

run <- adaptive_metropolis(
  target_density = target$density,
  initial_state = target$recommended_initial,
  n_iter = 300,
  initial_proposal_covariance = 0.25 * diag(2),
  adapt_start = 30
)

summarize_sampler_run(run, burn_in = 50)
plot_sampler_path(target, run, burn_in = 50)
```

## Teaching targets

The package includes small targets chosen to highlight specific issues:

- `make_standard_normal_target()`: a simple one-dimensional baseline
- `make_small_scale_target()`: illustrates high acceptance but tiny movement
- `make_correlated_gaussian_target()`: illustrates unknown correlation geometry
- `make_correlated_t_target()`: illustrates the same geometry on a
  non-Gaussian target
- `make_dram_rescue_target()`: illustrates delayed rejection after large rejected moves
- `make_multimodal_target()`: illustrates local trapping and mode switching diagnostics
- `make_high_dimensional_gaussian_target()`: illustrates higher-dimensional covariance tuning

## Main entry points

Sampler wrappers:

- `run_metropolis()`
- `run_adaptive_metropolis()`
- `run_ram()`
- `run_dram()`

Comparison helpers:

- `compare_adaptive_methods()`
- `summarize_chain()`
- `effective_sample_size_summary()`
- `benchmark_sampler_suite()`

Plots and diagnostics:

- `plot_trace()`
- `plot_acf_chain()`
- `plot_acceptance()`
- `plot_covariance_evolution()`
- `plot_mode_switching_behavior()`
- `plot_adaptation_diagnostic()`
- `plot_containment_diagnostic()`

## Theory posture

The package is careful about theory:

- it illustrates diminishing adaptation empirically
- it gives containment-style stability diagnostics
- it does **not** claim that simulation plots alone prove validity

That distinction is intentional and is reflected in the code comments, tests,
and vignettes.
