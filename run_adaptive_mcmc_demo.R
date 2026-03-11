# AdaptiveMCMCWorkflow demo
#
# Open this file in RStudio and run it section by section.
# This file uses one example only: a correlated Gaussian target.


# STEP 0: INITIALIZE THE PACKAGE ---------------------------------------
#
# From scratch:
# 1. setwd("/Users/gabrielalwan/Documents/AdaptiveMCMCworkflow")
# 2. install.packages("devtools")   # only once, if needed
# 3. Run this section

if (!file.exists("DESCRIPTION")) {
  stop(
    "Set the working directory to /Users/gabrielalwan/Documents/AdaptiveMCMCworkflow first.",
    call. = FALSE
  )
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop(
    "Install devtools first with install.packages('devtools'), then rerun this file.",
    call. = FALSE
  )
}

devtools::load_all(".")

if (!"target_density" %in% names(formals(get("basic_metropolis", envir = asNamespace("AdaptiveMCMCWorkflow"))))) {
  stop(
    "The current package source did not load correctly. Restart R and rerun this file.",
    call. = FALSE
  )
}

# Bind the current package functions explicitly so stale objects in the global
# environment cannot shadow them.
basic_metropolis <- get("basic_metropolis", envir = asNamespace("AdaptiveMCMCWorkflow"))
adaptive_metropolis <- get("adaptive_metropolis", envir = asNamespace("AdaptiveMCMCWorkflow"))
robust_adaptive_metropolis <- get("robust_adaptive_metropolis", envir = asNamespace("AdaptiveMCMCWorkflow"))
dram <- get("dram", envir = asNamespace("AdaptiveMCMCWorkflow"))
make_gaussian_target <- get("make_gaussian_target", envir = asNamespace("AdaptiveMCMCWorkflow"))
summarize_sampler_run <- get("summarize_sampler_run", envir = asNamespace("AdaptiveMCMCWorkflow"))
summarize_adaptive_validity <- get("summarize_adaptive_validity", envir = asNamespace("AdaptiveMCMCWorkflow"))
plot_target_contours <- get("plot_target_contours", envir = asNamespace("AdaptiveMCMCWorkflow"))
plot_sampler_path <- get("plot_sampler_path", envir = asNamespace("AdaptiveMCMCWorkflow"))
plot_trace <- get("plot_trace", envir = asNamespace("AdaptiveMCMCWorkflow"))
plot_acf_chain <- get("plot_acf_chain", envir = asNamespace("AdaptiveMCMCWorkflow"))
plot_covariance_evolution <- get("plot_covariance_evolution", envir = asNamespace("AdaptiveMCMCWorkflow"))
plot_acceptance <- get("plot_acceptance", envir = asNamespace("AdaptiveMCMCWorkflow"))
diminishing_adaptation_diagnostic <- get("diminishing_adaptation_diagnostic", envir = asNamespace("AdaptiveMCMCWorkflow"))
containment_diagnostic <- get("containment_diagnostic", envir = asNamespace("AdaptiveMCMCWorkflow"))

set.seed(1)

# When this file is run non-interactively, send plots to a temporary PDF so the
# script still runs linearly without repeated plot-condition checks.
plot_device_opened <- FALSE

if (!interactive()) {
  grDevices::pdf(file = tempfile(fileext = ".pdf"), width = 7, height = 5)
  plot_device_opened <- TRUE
}


# STEP 1: PARAMETERS YOU CAN CHANGE ------------------------------------
#
# target_mean:
#   the center of the Gaussian target.
#
# target_sd_1 and target_sd_2:
#   the marginal standard deviations of the two coordinates.
#
# target_rho:
#   the correlation between the two coordinates.
#   A value near 1 creates a narrow positively sloped ridge.
#
# target_covariance:
#   built from the standard deviations and the correlation.
#   This is the matrix used by the Gaussian target distribution.
#
# initial_state:
#   the starting point theta^(0) of the chain.
#
# S_0:
#   the starting Gaussian proposal covariance.
#   Basic Metropolis keeps this fixed.
#   AM, RAM, and DRAM all start from this proposal and then adapt.
#
# n_iter:
#   total number of MCMC iterations.
#
# burn_in:
#   number of early iterations dropped from summaries.
#
# adapt_start:
#   first iteration where AM and DRAM are allowed to adapt.
#
# target_acceptance:
#   RAM target acceptance rate alpha^*.
#
# gamma:
#   DRAM second-stage scale. It must be less than 1 so the rescue proposal is
#   smaller than the first-stage proposal.
#
target_mean <- c(0, 0)
target_sd_1 <- 1.0
target_sd_2 <- 2.0
target_rho <- 0.98

target_covariance <- matrix(
  c(
    target_sd_1^2,
    target_rho * target_sd_1 * target_sd_2,
    target_rho * target_sd_1 * target_sd_2,
    target_sd_2^2
  ),
  nrow = 2,
  byrow = TRUE
)

initial_state <- c(-3, 3)
S_0 <- 0.25 * diag(2)
n_iter <- 1000
burn_in <- 200
adapt_start <- 30
target_acceptance <- 0.234
gamma <- 0.5


# STEP 2: BUILD THE CORRELATED GAUSSIAN TARGET -------------------------

target <- make_gaussian_target(
  mean = target_mean,
  covariance = target_covariance,
  name = "Correlated Gaussian ridge",
  description = paste(
    "A two-dimensional Gaussian target with very strong positive correlation.",
    "This is the standard example for showing why covariance adaptation helps."
  ),
  recommended_initial = initial_state,
  recommended_proposal_covariance = S_0,
  recommended_adapt_start = adapt_start,
  recommended_delayed_rejection_scale = gamma,
  plot_limits = list(
    x = target_mean[1] + c(-2.5, 2.5) * target_sd_1,
    y = target_mean[2] + c(-2.5, 2.5) * target_sd_2
  )
)

target$name
target$description
target$dimension


# STEP 3: SHOW WHY THE TARGET IS CORRELATED ----------------------------

target_sd_1
target_sd_2
target_rho
target_covariance
target_correlation <- cov2cor(target_covariance)
target_correlation

# The off-diagonal entry of the covariance matrix is not the correlation.
# It equals rho * sd_1 * sd_2 = 0.98 * 1 * 2 = 1.96.
# The actual correlation matrix is shown by cov2cor(...).
# The off-diagonal correlation is 0.98.
# That means when one coordinate moves up, the other tends to move up too.
# The target density is concentrated along a narrow diagonal ridge.


x_limits <- target$plot_limits$x
y_limits <- target$plot_limits$y

x_grid <- seq(x_limits[1], x_limits[2], length.out = 160)
y_grid <- seq(y_limits[1], y_limits[2], length.out = 160)
relative_density <- outer(
  x_grid,
  y_grid,
  Vectorize(function(x1, x2) target$density(c(x1, x2)))
)
relative_density <- relative_density / max(relative_density)

graphics::contour(
  x = x_grid,
  y = y_grid,
  z = relative_density,
  xlab = expression(theta[1]),
  ylab = expression(theta[2]),
  main = "Correlated Gaussian target: contour plot",
  asp = 1,
  levels = c(0.08, 0.18, 0.32, 0.50, 0.68, 0.84),
  drawlabels = TRUE,
  lwd = 1.4,
  col = "#8c2d04"
)

graphics::points(
  target_mean[1],
  target_mean[2],
  pch = 19,
  col = "#08519c"
)

graphics::image(
  x = x_grid,
  y = y_grid,
  z = relative_density,
  col = grDevices::hcl.colors(20, "YlOrRd", rev = TRUE),
  xlab = expression(theta[1]),
  ylab = expression(theta[2]),
  main = "Correlated Gaussian target: filled top view",
  asp = 1,
  useRaster = TRUE
)

graphics::contour(
  x = x_grid,
  y = y_grid,
  z = relative_density,
  levels = c(0.08, 0.18, 0.32, 0.50, 0.68, 0.84),
  add = TRUE,
  drawlabels = FALSE,
  lwd = 1.1,
  col = "grey15"
)

graphics::persp(
  x = x_grid,
  y = y_grid,
  z = relative_density,
  theta = 40,
  phi = 25,
  expand = 0.7,
  col = "#d95f02",
  border = NA,
  shade = 0.4,
  ticktype = "detailed",
  xlab = expression(theta[1]),
  ylab = expression(theta[2]),
  zlab = "Relative density",
  main = "Correlated Gaussian target: 3D surface"
)


# INPUTS 

initial_state
S_0
n_iter
burn_in
adapt_start
target_acceptance
gamma


# RUN BASIC METROPOLIS 

basic_run <- basic_metropolis(
  target_density = target$density,
  initial_state = initial_state,
  n_iter = n_iter,
  proposal_covariance = S_0
)

basic_summary <- summarize_sampler_run(basic_run, burn_in = burn_in)
basic_summary


# DICTIONARY
#
# This target is two-dimensional, so the summary has one row for coordinate 1
# and one row for coordinate 2.
#
# acceptance_rate:
#   the fraction of post-burn-in moves that were accepted.
#   This matters, but it is not enough by itself.
#   A chain can have a high acceptance rate simply because the proposal moves
#   are tiny.
#   This summary value can differ from run$acceptance_rate because run$acceptance_rate
#   is computed over the full run, including the early adaptation period.
#
# lag1_autocorrelation:
#   the correlation between one draw and the very next draw.
#   Values close to 1 mean consecutive draws are very similar.
#   That means the chain is sticky and is moving slowly through the target.
#
# effective_sample_size:
#   the rough number of independent-equivalent draws after accounting for
#   autocorrelation.
#   If you keep 800 draws but the ESS is only about 10, the chain is
#   carrying much less information than 800 independent draws would.
#
# Why these metrics matter:
# - acceptance rate tells us how often proposals are accepted
# - autocorrelation tells us how much the chain is repeating itself
# - ESS tells us how much usable information the chain has produced
#
# For comparing samplers on this example:
# - lower lag1_autocorrelation is better
# - higher effective_sample_size is better
# - acceptance_rate is a secondary diagnostic, not the main goal

basic_ess_fraction <- basic_summary$effective_sample_size / basic_summary$iterations_used
basic_interpretation <- data.frame(
  coordinate = basic_summary$coordinate,
  acceptance_rate = basic_summary$acceptance_rate,
  lag1_autocorrelation = basic_summary$lag1_autocorrelation,
  effective_sample_size = basic_summary$effective_sample_size,
  ess_fraction_of_kept_draws = basic_ess_fraction
)

basic_interpretation

plot_sampler_path(target, basic_run, burn_in = burn_in, main = "Basic Metropolis path")
plot_trace(basic_run, burn_in = burn_in)
plot_acf_chain(basic_run, dim = 1, lag_max = 20, burn_in = burn_in)
plot_acceptance(basic_run, window = 50)


# RUN ADAPTIVE METROPOLIS 

am_run <- adaptive_metropolis(
  target_density = target$density,
  initial_state = initial_state,
  n_iter = 10000,
  initial_proposal_covariance = S_0,
  adapt_start = adapt_start
)

am_summary <- summarize_sampler_run(am_run, burn_in = burn_in)
am_summary

# AM can have a different post-burn-in acceptance rate than its full-run rate.
# That is expected because the proposal changes during the early adaptation
# period and the summary table intentionally drops those burn-in iterations.
am_run$acceptance_rate


plot_sampler_path(target, am_run, burn_in = burn_in, main = "Adaptive Metropolis path")
plot_trace(am_run, burn_in = burn_in)
plot_acf_chain(am_run, dim = 1, lag_max = 20, burn_in = burn_in)
plot_acceptance(am_run, window = 50)


# STEP 7: RUN RAM ------------------------------------------------------

ram_run <- robust_adaptive_metropolis(
  target_density = target$density,
  initial_state = initial_state,
  n_iter = 10000,
  initial_proposal_covariance = S_0,
  target_acceptance = target_acceptance
)

ram_summary <- summarize_sampler_run(ram_run, burn_in = burn_in)
ram_summary

plot_sampler_path(target, ram_run, burn_in = burn_in, main = "RAM path")
plot_trace(ram_run, burn_in = burn_in)
plot_acf_chain(ram_run, dim = 1, lag_max = 20, burn_in = burn_in)
plot_acceptance(ram_run, window = 50)


# STEP 8: RUN DRAM -----------------------------------------------------

dram_run <- dram(
  target_density = target$density,
  initial_state = initial_state,
  n_iter = 10000,
  initial_proposal_covariance = S_0,
  adapt_start = adapt_start,
  delayed_rejection_scale = gamma
)

dram_summary <- summarize_sampler_run(dram_run, burn_in = burn_in)
dram_summary

# These rates answer two different questions:
# - how often did DRAM need a second-stage rescue proposal?
# - when DRAM tried that rescue step, how often did it succeed?
mean(dram_run$second_stage_attempted[(burn_in + 1):n_iter])
mean(
  dram_run$stage_two_accepted[(burn_in + 1):n_iter][
    dram_run$second_stage_attempted[(burn_in + 1):n_iter]
  ]
)

plot_sampler_path(target, dram_run, burn_in = burn_in, main = "DRAM path")
plot_trace(dram_run, burn_in = burn_in)
plot_acf_chain(dram_run, dim = 1, lag_max = 20, burn_in = burn_in)
plot_acceptance(dram_run, window = 50)


# STEP 9: PUT THE FOUR SUMMARY TABLES TOGETHER -------------------------

summary_table <- rbind(
  basic_summary,
  am_summary,
  ram_summary,
  dram_summary
)

summary_table


# STEP 9A: COMPARE THE METHODS MORE DIRECTLY ---------------------------
#
# The table above still has one row per coordinate.
# This smaller table gives one row per algorithm by averaging the
# coordinate-specific metrics.

summary_by_algorithm <- do.call(
  rbind,
  lapply(
    split(summary_table, summary_table$algorithm),
    function(one_algorithm) {
      mean_or_na <- function(values) {
        if (all(is.na(values))) {
          return(NA_real_)
        }
        mean(values, na.rm = TRUE)
      }

      data.frame(
        algorithm = one_algorithm$algorithm[1],
        acceptance_rate = mean(one_algorithm$acceptance_rate),
        mean_jump_distance = mean(one_algorithm$mean_jump_distance),
        mean_lag1_autocorrelation = mean(one_algorithm$lag1_autocorrelation),
        mean_effective_sample_size = mean(one_algorithm$effective_sample_size),
        ess_fraction_of_kept_draws = mean(
          one_algorithm$effective_sample_size / one_algorithm$iterations_used
        ),
        stage_one_accept_rate = mean_or_na(one_algorithm$stage_one_accept_rate),
        second_stage_attempt_rate = mean_or_na(one_algorithm$second_stage_attempt_rate),
        second_stage_accept_rate = mean_or_na(one_algorithm$second_stage_accept_rate),
        second_stage_move_rate = mean_or_na(one_algorithm$second_stage_move_rate)
      )
    }
  )
)

summary_by_algorithm

# - acceptance_rate is the overall move rate after burn-in.
#   For DRAM, this includes both stage-one accepts and stage-two rescue accepts.
# - second_stage_accept_rate is conditional on stage two being attempted.
# - second_stage_move_rate is the share of kept iterations rescued at stage two.
# - If acceptance_rate is high but ESS is low, the chain is accepting many
#   moves without exploring efficiently.
# - If lag1 autocorrelation is lower and ESS is higher, the sampler is mixing
#   better.
# - The best method on this target is the one with the best overall tradeoff:
#   reasonable acceptance, lower autocorrelation, and higher ESS.


# SCHECK ADAPTIVE-MCMC ASSUMPTIONS EMPIRICALLY 
#
# These are empirical diagnostics, not proofs.
# We check:
# 1. whether adaptation becomes smaller later in the run,
# 2. whether the adaptive proposal covariance stays finite and well behaved.

validity_table <- rbind(
  summarize_adaptive_validity(am_run, early_window = 100, late_window = 100),
  summarize_adaptive_validity(ram_run, early_window = 100, late_window = 100),
  summarize_adaptive_validity(dram_run, early_window = 100, late_window = 100)
)

validity_table

am_diminishing <- diminishing_adaptation_diagnostic(am_run)
graphics::plot(
  am_diminishing$iteration,
  am_diminishing$adaptation_magnitude,
  type = "l",
  xlab = "Iteration",
  ylab = "Kernel change",
  main = "Adaptive Metropolis: adaptation magnitude"
)

am_containment <- containment_diagnostic(am_run)
graphics::plot(
  am_containment$iteration,
  am_containment$average_marginal_sd,
  type = "l",
  xlab = "Iteration",
  ylab = "Average marginal proposal SD",
  main = "Adaptive Metropolis: proposal scale"
)

graphics::plot(
  am_containment$iteration,
  am_containment$condition_number,
  type = "l",
  xlab = "Iteration",
  ylab = "Condition number",
  main = "Adaptive Metropolis: proposal condition number"
)

ram_diminishing <- diminishing_adaptation_diagnostic(ram_run)
graphics::plot(
  ram_diminishing$iteration,
  ram_diminishing$adaptation_magnitude,
  type = "l",
  xlab = "Iteration",
  ylab = "Kernel change",
  main = "RAM: adaptation magnitude"
)

ram_containment <- containment_diagnostic(ram_run)
graphics::plot(
  ram_containment$iteration,
  ram_containment$average_marginal_sd,
  type = "l",
  xlab = "Iteration",
  ylab = "Average marginal proposal SD",
  main = "RAM: proposal scale"
)

graphics::plot(
  ram_containment$iteration,
  ram_containment$condition_number,
  type = "l",
  xlab = "Iteration",
  ylab = "Condition number",
  main = "RAM: proposal condition number"
)

dram_diminishing <- diminishing_adaptation_diagnostic(dram_run)
graphics::plot(
  dram_diminishing$iteration,
  dram_diminishing$adaptation_magnitude,
  type = "l",
  xlab = "Iteration",
  ylab = "Kernel change",
  main = "DRAM: adaptation magnitude"
)

dram_containment <- containment_diagnostic(dram_run)
graphics::plot(
  dram_containment$iteration,
  dram_containment$average_marginal_sd,
  type = "l",
  xlab = "Iteration",
  ylab = "Average marginal proposal SD",
  main = "DRAM: proposal scale"
)

graphics::plot(
  dram_containment$iteration,
  dram_containment$condition_number,
  type = "l",
  xlab = "Iteration",
  ylab = "Condition number",
  main = "DRAM: proposal condition number"
)

if (plot_device_opened) {
  grDevices::dev.off()
}
