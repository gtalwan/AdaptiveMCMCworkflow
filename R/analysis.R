#' Summarise One Sampler Run
#'
#' The package summary stays intentionally simple. It reports acceptance rate,
#' mean jump distance, lag-1 autocorrelation, and a rough effective sample size.
#'
#' @param run A sampler result created by this package.
#' @param burn_in Number of initial iterations to discard before computing the
#'   summary.
#' @param ess_max_lag Maximum lag used in the rough ESS calculation.
#'
#' @return A data frame with one row per coordinate of the chain.
#'
#' @examples
#' target <- make_small_scale_target()
#' run <- adaptive_metropolis(
#'   target_density = target$density,
#'   initial_state = target$recommended_initial,
#'   n_iter = 300,
#'   initial_proposal_covariance = target$recommended_proposal_covariance,
#'   adapt_start = 40
#' )
#' summarize_sampler_run(run)
#'
#' @export
summarize_sampler_run <- function(run, burn_in = 0L, ess_max_lag = NULL) {
  run <- validate_sampler_run(run)

  # Drop the early part of the chain before computing summary diagnostics.
  # This makes the table line up with the way burn-in is explained in the demo.
  kept_draws <- draws_after_burn_in(run, burn_in = burn_in)
  kept_index <- seq.int(from = burn_in + 1L, to = run$n_iter)
  kept_accepted <- run$accepted[kept_index]

  # Autocorrelation is a coordinate-wise diagnostic, so we compute it
  # separately for each column of the chain matrix.
  lag_one_acf <- vapply(
    seq_len(ncol(kept_draws)),
    function(index) {
      values <- acf_values_internal(kept_draws[, index], lag_max = 1L)
      if (length(values)) values[1] else NA_real_
    },
    numeric(1)
  )

  # ESS is also coordinate-specific. Each component can mix at a different
  # rate, especially on anisotropic targets such as correlated Gaussian ridges.
  ess <- vapply(
    seq_len(ncol(kept_draws)),
    function(index) {
      effective_sample_size_internal(kept_draws[, index], lag_max = ess_max_lag)
    },
    numeric(1)
  )

  data.frame(
    algorithm = run$algorithm,
    coordinate = seq_len(ncol(kept_draws)),
    iterations_used = nrow(kept_draws),
    # This is the realised move rate after burn-in. For DRAM this includes
    # both ordinary stage-one accepts and rescued stage-two accepts.
    acceptance_rate = mean(kept_accepted),
    # Jump distance is Euclidean distance between successive realised states.
    # Tiny average jumps often indicate slow exploration even if acceptance is high.
    mean_jump_distance = mean(run$jump_distance[kept_index]),
    # Lag-1 autocorrelation asks how similar one draw is to the very next draw.
    # Values close to 1 mean the chain is sticky.
    lag1_autocorrelation = lag_one_acf,
    # ESS converts the autocorrelated chain into an "independent-equivalent"
    # sample size. This is usually more informative than acceptance alone.
    effective_sample_size = ess,
    final_target_value = utils::tail(run$target_values[kept_index], 1L),
    # For non-DRAM runs these stage-specific columns are NA by construction.
    stage_one_accept_rate = if ("stage_one_accepted" %in% names(run)) {
      mean(run$stage_one_accepted[kept_index])
    } else {
      NA_real_
    },
    second_stage_attempt_rate = if ("second_stage_attempted" %in% names(run)) {
      mean(run$second_stage_attempted[kept_index])
    } else {
      NA_real_
    },
    second_stage_accept_rate = if ("stage_two_accepted" %in% names(run)) {
      attempted <- run$second_stage_attempted[kept_index]
      if (any(attempted)) {
        # Conditional DRAM rescue success rate:
        # among the iterations that actually reached stage two, how many moved?
        mean(run$stage_two_accepted[kept_index][attempted])
      } else {
        NA_real_
      }
    } else {
      NA_real_
    },
    second_stage_move_rate = if ("stage_two_accepted" %in% names(run)) {
      # Unconditional DRAM stage-two contribution:
      # the share of kept iterations rescued at the second stage.
      mean(run$stage_two_accepted[kept_index])
    } else {
      NA_real_
    }
  )
}

#' Summarise a Full Sampler Comparison
#'
#' This is a convenience wrapper around [summarize_sampler_run()] for the output
#' of [run_sampler_suite()].
#'
#' @param comparison Output from [run_sampler_suite()].
#' @param burn_in Number of initial iterations discarded from each sampler.
#' @param ess_max_lag Maximum lag used in the rough ESS calculation.
#'
#' @return A combined summary data frame.
#'
#' @examples
#' comparison <- run_sampler_suite(make_correlated_t_target(), n_iter = 250, seed = 1)
#' summarize_sampler_suite(comparison)
#'
#' @export
summarize_sampler_suite <- function(comparison, burn_in = 0L, ess_max_lag = NULL) {
  comparison <- validate_sampler_comparison(comparison)
  do.call(
    rbind,
    lapply(
      comparison$runs,
      summarize_sampler_run,
      burn_in = burn_in,
      ess_max_lag = ess_max_lag
    )
  )
}

#' Summarize a Chain
#'
#' This wrapper dispatches to [summarize_sampler_run()] or
#' [summarize_sampler_suite()] depending on what you pass in.
#'
#' @param object Either an `"adaptive_mcmc_run"` object or an
#'   `"adaptive_mcmc_comparison"` object.
#' @param burn_in Number of initial iterations discarded before summarizing.
#' @param ess_max_lag Maximum lag used in the rough ESS calculation.
#'
#' @return A summary data frame.
#'
#' @examples
#' target <- make_standard_normal_target()
#' run <- run_metropolis(
#'   target$density,
#'   target$recommended_initial,
#'   200,
#'   target$recommended_proposal_covariance
#' )
#' summarize_chain(run)
#'
#' @export
summarize_chain <- function(object, burn_in = 0L, ess_max_lag = NULL) {
  if (inherits(object, "adaptive_mcmc_run")) {
    return(summarize_sampler_run(object, burn_in = burn_in, ess_max_lag = ess_max_lag))
  }

  if (inherits(object, "adaptive_mcmc_comparison")) {
    return(summarize_sampler_suite(object, burn_in = burn_in, ess_max_lag = ess_max_lag))
  }

  stop("`object` must be a sampler run or a sampler comparison.", call. = FALSE)
}

#' Effective Sample Size Summary
#'
#' This helper extracts just the rough effective sample size diagnostics from
#' the chain summary helpers.
#'
#' @param object Either a sampler run or a sampler comparison.
#' @param burn_in Number of initial iterations discarded.
#' @param ess_max_lag Maximum lag used in the rough ESS calculation.
#'
#' @return A data frame containing the algorithm name, coordinate, and rough ESS.
#'
#' @examples
#' comparison <- compare_adaptive_methods(make_correlated_t_target(), n_iter = 200, seed = 1)
#' effective_sample_size_summary(comparison$comparison)
#'
#' @export
effective_sample_size_summary <- function(object, burn_in = 0L, ess_max_lag = NULL) {
  summary_table <- summarize_chain(object, burn_in = burn_in, ess_max_lag = ess_max_lag)
  summary_table[c("algorithm", "coordinate", "effective_sample_size")]
}

#' Benchmark the Teaching Samplers
#'
#' This helper records simple wall-clock timings for the four teaching
#' algorithms on the same target.
#'
#' @param target A target object created by this package.
#' @param n_iter Number of iterations used in each timed run.
#' @param repetitions Number of repeated timings per algorithm.
#' @param initial_state Optional common starting point.
#' @param controls Optional sampler-specific controls merged into the package
#'   defaults.
#' @param seed Optional seed reused within each repetition.
#'
#' @return A list with `raw` and `summary` timing tables.
#'
#' @examples
#' benchmark <- benchmark_sampler_suite(
#'   make_correlated_t_target(),
#'   n_iter = 200,
#'   repetitions = 2,
#'   seed = 1
#' )
#' benchmark$summary
#'
#' @export
benchmark_sampler_suite <- function(target,
                                    n_iter = 1000L,
                                    repetitions = 3L,
                                    initial_state = NULL,
                                    controls = NULL,
                                    seed = NULL) {
  target <- validate_target(target)
  n_iter <- ensure_positive_integer(n_iter, "n_iter")
  repetitions <- ensure_positive_integer(repetitions, "repetitions")
  initial_state <- initial_state %||% target$recommended_initial %||% rep(0, target$dimension)
  initial_state <- ensure_numeric_state(initial_state)

  if (length(initial_state) != target$dimension) {
    stop("`initial_state` must match the target dimension.", call. = FALSE)
  }

  controls <- merge_sampler_controls(default_sampler_controls(target), controls)

  algorithms <- list(
    basic = function() {
      basic_metropolis(
        target_density = target$density,
        initial_state = initial_state,
        n_iter = n_iter,
        proposal_covariance = controls$basic$proposal_covariance
      )
    },
    adaptive = function() {
      adaptive_metropolis(
        target_density = target$density,
        initial_state = initial_state,
        n_iter = n_iter,
        initial_proposal_covariance = controls$adaptive$initial_proposal_covariance,
        adapt_start = controls$adaptive$adapt_start,
        epsilon = controls$adaptive$epsilon
      )
    },
    ram = function() {
      robust_adaptive_metropolis(
        target_density = target$density,
        initial_state = initial_state,
        n_iter = n_iter,
        initial_proposal_covariance = controls$ram$initial_proposal_covariance,
        target_acceptance = controls$ram$target_acceptance,
        adapt_exponent = controls$ram$adapt_exponent
      )
    },
    dram = function() {
      dram(
        target_density = target$density,
        initial_state = initial_state,
        n_iter = n_iter,
        initial_proposal_covariance = controls$dram$initial_proposal_covariance,
        adapt_start = controls$dram$adapt_start,
        delayed_rejection_scale = controls$dram$delayed_rejection_scale,
        epsilon = controls$dram$epsilon
      )
    }
  )

  raw_rows <- vector("list", length(algorithms) * repetitions)
  row_index <- 1L

  for (replication in seq_len(repetitions)) {
    for (algorithm_name in names(algorithms)) {
      if (!is.null(seed)) {
        set.seed(as.integer(seed) + replication - 1L)
      }

      result <- NULL
      timing <- system.time({
        result <- algorithms[[algorithm_name]]()
      })

      raw_rows[[row_index]] <- data.frame(
        algorithm = result$algorithm,
        repetition = replication,
        elapsed_seconds = unname(timing["elapsed"]),
        user_seconds = unname(timing["user.self"]),
        system_seconds = unname(timing["sys.self"]),
        acceptance_rate = result$acceptance_rate,
        mean_jump_distance = mean(result$jump_distance)
      )

      row_index <- row_index + 1L
    }
  }

  raw_table <- do.call(rbind, raw_rows)
  summary_table <- stats::aggregate(
    raw_table[c(
      "elapsed_seconds",
      "user_seconds",
      "system_seconds",
      "acceptance_rate",
      "mean_jump_distance"
    )],
    by = list(algorithm = raw_table$algorithm),
    FUN = stats::median
  )

  list(raw = raw_table, summary = summary_table)
}

#' Diminishing Adaptation Diagnostic
#'
#' This function reports a practical proxy for diminishing adaptation: the
#' Frobenius norm difference between successive proposal covariance matrices.
#'
#' @param run A sampler result created by this package.
#'
#' @return A data frame with one row per iteration.
#'
#' @examples
#' comparison <- run_sampler_suite(make_correlated_t_target(), n_iter = 250, seed = 1)
#' head(diminishing_adaptation_diagnostic(comparison$runs$adaptive))
#'
#' @export
diminishing_adaptation_diagnostic <- function(run) {
  run <- validate_sampler_run(run)
  magnitude <- as.numeric(run$adaptation_magnitude)

  # The raw adaptation magnitude records how much the proposal changed at each
  # step. The cumulative mean smooths that sequence so longer-run trends are
  # easier to inspect in a plot.
  data.frame(
    iteration = seq_along(magnitude),
    adaptation_magnitude = magnitude,
    cumulative_mean = cumsum(magnitude) / seq_along(magnitude)
  )
}

#' Containment-Style Stability Diagnostic
#'
#' A simulation run cannot prove containment in general, so this function
#' reports empirical stability summaries of the proposal covariance history.
#'
#' @param run A sampler result created by this package.
#'
#' @return A data frame with one row per iteration.
#'
#' @examples
#' comparison <- run_sampler_suite(make_dram_rescue_target(), n_iter = 250, seed = 1)
#' head(containment_diagnostic(comparison$runs$dram))
#'
#' @export
containment_diagnostic <- function(run) {
  run <- validate_sampler_run(run)
  history <- run$proposal_covariance_history
  n_iter <- dim(history)[3]
  min_eigenvalue <- numeric(n_iter)
  max_eigenvalue <- numeric(n_iter)
  trace_value <- numeric(n_iter)
  log_determinant <- numeric(n_iter)
  condition_number <- numeric(n_iter)
  average_marginal_sd <- numeric(n_iter)
  dimension <- dim(history)[1]

  for (iteration in seq_len(n_iter)) {
    # Pull out the proposal covariance used at this iteration.
    covariance <- matrix(
      history[, , iteration, drop = TRUE],
      nrow = dimension,
      ncol = dimension
    )

    # The eigenvalues describe the proposal scale in the principal directions.
    # We reuse them to form several containment-style diagnostics.
    eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    min_eigenvalue[iteration] <- min(eigenvalues)
    max_eigenvalue[iteration] <- max(eigenvalues)
    trace_value[iteration] <- sum(eigenvalues)
    log_determinant[iteration] <- as.numeric(
      determinant(covariance, logarithm = TRUE)$modulus
    )

    # The condition number measures anisotropy:
    # large / small eigenvalue. It grows when the proposal becomes much more
    # stretched in one direction than another.
    condition_number[iteration] <- max_eigenvalue[iteration] / min_eigenvalue[iteration]

    # The square roots of the diagonal entries are marginal proposal standard
    # deviations. Their average is a compact summary of overall proposal scale.
    average_marginal_sd[iteration] <- sqrt(mean(diag(covariance)))
  }

  data.frame(
    iteration = seq_len(n_iter),
    min_eigenvalue = min_eigenvalue,
    max_eigenvalue = max_eigenvalue,
    trace = trace_value,
    log_determinant = log_determinant,
    condition_number = condition_number,
    average_marginal_sd = average_marginal_sd
  )
}

#' Summarize Adaptive MCMC Validity Diagnostics
#'
#' This helper turns the adaptive diagnostics into a compact teaching table.
#' It does not prove ergodicity. It simply reports whether adaptation appears
#' to shrink and whether the proposal covariance remains numerically stable in
#' the observed run.
#'
#' @param object Either one sampler run or one sampler comparison.
#' @param early_window Number of early iterations used when averaging the
#'   adaptation magnitude.
#' @param late_window Number of late iterations used when averaging the
#'   adaptation magnitude.
#'
#' @return A data frame summarizing diminishing-adaptation and
#'   containment-style diagnostics.
#'
#' @examples
#' comparison <- run_sampler_suite(make_correlated_gaussian_target(), n_iter = 250, seed = 1)
#' summarize_adaptive_validity(comparison)
#'
#' @export
summarize_adaptive_validity <- function(object,
                                        early_window = 50L,
                                        late_window = 50L) {
  early_window <- ensure_positive_integer(early_window, "early_window")
  late_window <- ensure_positive_integer(late_window, "late_window")

  if (inherits(object, "adaptive_mcmc_comparison")) {
    adaptive_runs <- object$runs[c("adaptive", "ram", "dram")]

    return(do.call(
      rbind,
      lapply(
        adaptive_runs,
        summarize_adaptive_validity,
        early_window = early_window,
        late_window = late_window
      )
    ))
  }

  run <- validate_sampler_run(object)
  diminishing <- diminishing_adaptation_diagnostic(run)
  containment <- containment_diagnostic(run)

  # Find the first iteration where adaptation is actually nonzero. This avoids
  # comparing "early" and "late" windows that still include fixed pre-adaptation
  # iterations.
  first_adaptive_iteration <- match(TRUE, diminishing$adaptation_magnitude > 0)

  if (is.na(first_adaptive_iteration)) {
    first_adaptive_iteration <- 1L
  }

  # Compare an early adaptive window against a late window. The ratio
  # late / early is a simple empirical proxy for diminishing adaptation.
  early_window <- min(early_window, nrow(diminishing) - first_adaptive_iteration + 1L)
  late_window <- min(late_window, nrow(diminishing))
  early_slice <- seq.int(
    from = first_adaptive_iteration,
    length.out = early_window
  )
  early_mean <- mean(diminishing$adaptation_magnitude[early_slice])
  late_mean <- mean(utils::tail(diminishing$adaptation_magnitude, late_window))
  kernel_change_ratio <- if (early_mean > 0) late_mean / early_mean else NA_real_

  # This table is intentionally cautious. It summarizes whether adaptation got
  # smaller and whether the proposal covariance stayed numerically reasonable.
  # It does not prove containment or ergodicity.
  data.frame(
    algorithm = run$algorithm,
    early_mean_kernel_change = early_mean,
    late_mean_kernel_change = late_mean,
    kernel_change_ratio = kernel_change_ratio,
    adaptation_is_smaller_late = late_mean < early_mean,
    min_proposal_eigenvalue = min(containment$min_eigenvalue),
    max_proposal_condition_number = max(containment$condition_number),
    final_average_marginal_sd = utils::tail(containment$average_marginal_sd, 1L)
  )
}

#' Plot Contours of a Two-Dimensional Target
#'
#' This helper evaluates the target density on a grid and draws contour lines that
#' can be overlaid with MCMC paths.
#'
#' @param target A two-dimensional target object created by this package.
#' @param n_grid Number of grid points in each coordinate direction.
#' @param levels Number of contour levels.
#' @param xlim Optional horizontal plot limits.
#' @param ylim Optional vertical plot limits.
#' @param add Should the contour plot be added to an existing figure?
#' @param ... Additional arguments passed to [graphics::contour()].
#'
#' @return Invisibly returns the grid used to draw the contours.
#'
#' @examples
#' target <- make_correlated_t_target()
#' plot_target_contours(target, n_grid = 60)
#'
#' @export
plot_target_contours <- function(target,
                                 n_grid = 100L,
                                 levels = 10L,
                                 xlim = NULL,
                                 ylim = NULL,
                                 add = FALSE,
                                 ...) {
  target <- validate_target(target)

  if (target$dimension != 2L) {
    stop("Contour plots are only defined here for two-dimensional targets.", call. = FALSE)
  }

  n_grid <- ensure_positive_integer(n_grid, "n_grid")
  levels <- ensure_positive_integer(levels, "levels")
  xlim <- xlim %||% target$plot_limits$x
  ylim <- ylim %||% target$plot_limits$y
  x_grid <- seq(xlim[1], xlim[2], length.out = n_grid)
  y_grid <- seq(ylim[1], ylim[2], length.out = n_grid)
  density_matrix <- matrix(NA_real_, nrow = length(x_grid), ncol = length(y_grid))

  for (x_index in seq_along(x_grid)) {
    for (y_index in seq_along(y_grid)) {
      density_matrix[x_index, y_index] <- target$density(c(
        x_grid[x_index],
        y_grid[y_index]
      ))
    }
  }

  relative_density <- density_matrix / max(density_matrix)

  graphics::contour(
    x = x_grid,
    y = y_grid,
    z = relative_density,
    nlevels = levels,
    add = add,
    drawlabels = FALSE,
    xlab = if (add) "" else expression(x[1]),
    ylab = if (add) "" else expression(x[2]),
    ...
  )

  invisible(list(x = x_grid, y = y_grid, z = relative_density))
}

#' Plot a Two-Dimensional Sampler Path
#'
#' This plot overlays a realised MCMC trajectory on top of the target contours.
#'
#' @param target A two-dimensional target object.
#' @param run A sampler result created by this package.
#' @param burn_in Number of initial iterations to omit from the path.
#' @param add_contours Should target contours be drawn before the chain?
#' @param line_col Color used for the path line.
#' @param point_col Color used for the sample points.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return Invisibly returns the post-burn-in draws.
#'
#' @examples
#' comparison <- run_sampler_suite(make_correlated_t_target(), n_iter = 250, seed = 1)
#' plot_sampler_path(comparison$target, comparison$runs$adaptive, burn_in = 50)
#'
#' @export
plot_sampler_path <- function(target,
                              run,
                              burn_in = 0L,
                              add_contours = TRUE,
                              line_col = "#1f78b4",
                              point_col = "#1f78b455",
                              ...) {
  target <- validate_target(target)
  run <- validate_sampler_run(run)
  draws <- draws_after_burn_in(run, burn_in = burn_in)

  if (target$dimension != 2L || ncol(draws) != 2L) {
    stop("`plot_sampler_path()` requires a two-dimensional target and run.", call. = FALSE)
  }

  if (add_contours) {
    plot_target_contours(target, col = "grey70")
  } else {
    graphics::plot(
      draws[, 1],
      draws[, 2],
      type = "n",
      xlab = expression(x[1]),
      ylab = expression(x[2]),
      ...
    )
  }

  graphics::lines(draws[, 1], draws[, 2], col = line_col, ...)
  graphics::points(draws[, 1], draws[, 2], pch = 16, cex = 0.5, col = point_col)

  invisible(draws)
}

#' Plot a Four-Way Sampler Comparison
#'
#' This helper draws one panel for each of the four teaching samplers.
#'
#' @param comparison Output from [run_sampler_suite()].
#' @param burn_in Number of initial iterations to discard in each panel.
#'
#' @return Invisibly returns the comparison object.
#'
#' @examples
#' comparison <- run_sampler_suite(make_correlated_t_target(), n_iter = 250, seed = 1)
#' plot_sampler_comparison(comparison, burn_in = 50)
#'
#' @export
plot_sampler_comparison <- function(comparison, burn_in = 0L) {
  comparison <- validate_sampler_comparison(comparison)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  set_compact_plot_layout(4L)

  plot_sampler_path(comparison$target, comparison$runs$basic, burn_in = burn_in)
  graphics::title(main = comparison$runs$basic$algorithm)

  plot_sampler_path(comparison$target, comparison$runs$adaptive, burn_in = burn_in)
  graphics::title(main = comparison$runs$adaptive$algorithm)

  plot_sampler_path(comparison$target, comparison$runs$ram, burn_in = burn_in)
  graphics::title(main = comparison$runs$ram$algorithm)

  plot_sampler_path(comparison$target, comparison$runs$dram, burn_in = burn_in)
  graphics::title(main = comparison$runs$dram$algorithm)

  invisible(comparison)
}

#' Plot Coordinate Traces
#'
#' Trace plots are a simple way to show sticking and slow exploration.
#'
#' @param run A sampler result created by this package.
#' @param dims Coordinates to plot.
#' @param burn_in Number of initial iterations to omit.
#'
#' @return Invisibly returns the post-burn-in draws used in the plot.
#'
#' @examples
#' comparison <- run_sampler_suite(make_small_scale_target(), n_iter = 250, seed = 1)
#' plot_trace(comparison$runs$basic)
#'
#' @export
plot_trace <- function(run, dims = NULL, burn_in = 0L) {
  run <- validate_sampler_run(run)
  draws <- draws_after_burn_in(run, burn_in = burn_in)
  dims <- dims %||% seq_len(min(2L, ncol(draws)))
  dims <- as.integer(dims)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  set_compact_plot_layout(length(dims))

  for (dimension in dims) {
    graphics::plot(
      draws[, dimension],
      type = "l",
      xlab = "Iteration",
      ylab = sprintf("x[%d]", dimension),
      main = sprintf("%s trace for coordinate %d", run$algorithm, dimension)
    )
  }

  invisible(draws)
}

#' Plot Jump Distance Histogram
#'
#' Jump distances are especially useful in the poor-proposal-scale example.
#'
#' @param run A sampler result created by this package.
#' @param burn_in Number of initial iterations to discard.
#' @param breaks Histogram break specification passed to [graphics::hist()].
#'
#' @return Invisibly returns the jump distances used in the histogram.
#'
#' @examples
#' comparison <- run_sampler_suite(make_small_scale_target(), n_iter = 250, seed = 1)
#' plot_jump_distance_histogram(comparison$runs$basic)
#'
#' @export
plot_jump_distance_histogram <- function(run, burn_in = 0L, breaks = "FD") {
  run <- validate_sampler_run(run)
  kept_index <- seq.int(from = burn_in + 1L, to = run$n_iter)
  jump_distance <- run$jump_distance[kept_index]

  graphics::hist(
    jump_distance,
    breaks = breaks,
    col = "grey80",
    border = "white",
    main = sprintf("Jump distances: %s", run$algorithm),
    xlab = "Euclidean jump distance"
  )

  invisible(jump_distance)
}

#' Plot an Autocorrelation Comparison
#'
#' This figure compares the sample autocorrelation function for one coordinate
#' across the four teaching samplers.
#'
#' @param comparison Output from [run_sampler_suite()].
#' @param dim Coordinate to analyse.
#' @param lag_max Maximum lag shown in the plot.
#' @param burn_in Number of initial iterations to discard.
#'
#' @return Invisibly returns the matrix of plotted autocorrelations.
#'
#' @examples
#' comparison <- run_sampler_suite(make_small_scale_target(), n_iter = 250, seed = 1)
#' plot_acf_comparison(comparison, dim = 1, lag_max = 20)
#'
#' @export
plot_acf_comparison <- function(comparison, dim = 1L, lag_max = 30L, burn_in = 0L) {
  comparison <- validate_sampler_comparison(comparison)
  dim <- as.integer(dim)
  lag_max <- ensure_positive_integer(lag_max, "lag_max")

  acf_matrix <- vapply(
    comparison$runs,
    function(run) {
      draws <- draws_after_burn_in(run, burn_in = burn_in)

      if (dim > ncol(draws)) {
        stop("`dim` is larger than the target dimension.", call. = FALSE)
      }

      acf_values_internal(draws[, dim], lag_max = lag_max)
    },
    numeric(lag_max)
  )

  graphics::matplot(
    x = seq_len(lag_max),
    y = acf_matrix,
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "Lag",
    ylab = "Sample autocorrelation",
    main = sprintf("ACF comparison for coordinate %d", dim)
  )

  graphics::abline(h = 0, lty = 2, col = "grey60")
  graphics::legend(
    "topright",
    legend = names(comparison$runs),
    col = seq_along(comparison$runs),
    lty = 1,
    lwd = 2
  )

  invisible(acf_matrix)
}

#' Plot Diminishing Adaptation Diagnostics
#'
#' This figure shows the adaptation magnitude and the average marginal proposal
#' standard deviation over time.
#'
#' @param run A sampler result created by this package.
#'
#' @return Invisibly returns the plotted diagnostic data frames.
#'
#' @examples
#' comparison <- run_sampler_suite(make_correlated_t_target(), n_iter = 250, seed = 1)
#' plot_adaptation_diagnostic(comparison$runs$adaptive)
#'
#' @export
plot_adaptation_diagnostic <- function(run) {
  run <- validate_sampler_run(run)
  diminishing <- diminishing_adaptation_diagnostic(run)
  containment <- containment_diagnostic(run)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  set_compact_plot_layout(2L)

  graphics::plot(
    diminishing$iteration,
    diminishing$adaptation_magnitude,
    type = "l",
    xlab = "Iteration",
    ylab = "Kernel change",
    main = sprintf("Diminishing adaptation proxy: %s", run$algorithm)
  )

  graphics::plot(
    containment$iteration,
    containment$average_marginal_sd,
    type = "l",
    xlab = "Iteration",
    ylab = "Average marginal proposal SD",
    main = sprintf("Proposal scale history: %s", run$algorithm)
  )

  invisible(list(diminishing = diminishing, containment = containment))
}

#' Plot Containment-Style Stability Diagnostics
#'
#' These plots display empirical stability summaries of the proposal covariance
#' history. They are heuristic checks, not proofs of containment.
#'
#' @param run A sampler result created by this package.
#'
#' @return Invisibly returns the underlying containment diagnostic data frame.
#'
#' @examples
#' comparison <- run_sampler_suite(make_dram_rescue_target(), n_iter = 250, seed = 1)
#' plot_containment_diagnostic(comparison$runs$dram)
#'
#' @export
plot_containment_diagnostic <- function(run) {
  run <- validate_sampler_run(run)
  diagnostics <- containment_diagnostic(run)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  set_compact_plot_layout(4L)

  graphics::plot(
    diagnostics$iteration,
    diagnostics$max_eigenvalue,
    type = "l",
    xlab = "Iteration",
    ylab = "Largest eigenvalue",
    main = "Largest proposal eigenvalue"
  )

  graphics::plot(
    diagnostics$iteration,
    diagnostics$min_eigenvalue,
    type = "l",
    xlab = "Iteration",
    ylab = "Smallest eigenvalue",
    main = "Smallest proposal eigenvalue"
  )

  graphics::plot(
    diagnostics$iteration,
    diagnostics$condition_number,
    type = "l",
    xlab = "Iteration",
    ylab = "Condition number",
    main = "Proposal condition number"
  )

  graphics::plot(
    diagnostics$iteration,
    diagnostics$log_determinant,
    type = "l",
    xlab = "Iteration",
    ylab = "log det(covariance)",
    main = "Proposal log determinant"
  )

  invisible(diagnostics)
}

#' Plot the Autocorrelation Function of One Chain
#'
#' This single-chain ACF plot is useful when you want to explain one algorithm
#' at a time.
#'
#' @param run A sampler result created by this package.
#' @param dim Coordinate to analyse.
#' @param lag_max Maximum lag displayed.
#' @param burn_in Number of initial iterations discarded.
#'
#' @return Invisibly returns the plotted autocorrelation values.
#'
#' @examples
#' comparison <- run_sampler_suite(make_standard_normal_target(), n_iter = 250, seed = 1)
#' plot_acf_chain(comparison$runs$basic, dim = 1, lag_max = 20)
#'
#' @export
plot_acf_chain <- function(run, dim = 1L, lag_max = 30L, burn_in = 0L) {
  run <- validate_sampler_run(run)
  draws <- draws_after_burn_in(run, burn_in = burn_in)
  dim <- as.integer(dim)
  lag_max <- ensure_positive_integer(lag_max, "lag_max")

  if (dim < 1L || dim > ncol(draws)) {
    stop("`dim` must select an existing coordinate.", call. = FALSE)
  }

  acf_values <- acf_values_internal(draws[, dim], lag_max = lag_max)

  graphics::plot(
    seq_along(acf_values),
    acf_values,
    type = "h",
    lwd = 3,
    xlab = "Lag",
    ylab = "Sample autocorrelation",
    main = sprintf("ACF: %s, coordinate %d", run$algorithm, dim)
  )
  graphics::abline(h = 0, lty = 2, col = "grey60")
  graphics::points(seq_along(acf_values), acf_values, pch = 16)

  invisible(acf_values)
}

#' Plot Acceptance Behavior Over Time
#'
#' This plot shows both cumulative acceptance and a rolling acceptance window.
#'
#' @param run A sampler result created by this package.
#' @param window Window size for the rolling acceptance rate.
#'
#' @return Invisibly returns a data frame of plotted acceptance summaries.
#'
#' @examples
#' comparison <- run_sampler_suite(make_small_scale_target(), n_iter = 250, seed = 1)
#' plot_acceptance(comparison$runs$ram)
#'
#' @export
plot_acceptance <- function(run, window = 50L) {
  run <- validate_sampler_run(run)
  window <- ensure_positive_integer(window, "window")
  accepted <- as.numeric(run$accepted)
  cumulative_acceptance <- cumsum(accepted) / seq_along(accepted)

  rolling_acceptance <- vapply(
    seq_along(accepted),
    function(index) {
      start_index <- max(1L, index - window + 1L)
      mean(accepted[start_index:index])
    },
    numeric(1)
  )

  graphics::plot(
    cumulative_acceptance,
    type = "l",
    lwd = 2,
    ylim = c(0, 1),
    xlab = "Iteration",
    ylab = "Acceptance rate",
    main = sprintf("Acceptance behavior: %s", run$algorithm)
  )
  graphics::lines(rolling_acceptance, col = 2, lwd = 2)
  graphics::legend(
    "topright",
    legend = c("Cumulative", sprintf("Rolling (%d)", window)),
    col = c(1, 2),
    lty = 1,
    lwd = 2
  )

  invisible(data.frame(
    iteration = seq_along(accepted),
    cumulative_acceptance = cumulative_acceptance,
    rolling_acceptance = rolling_acceptance
  ))
}

#' Plot Proposal Covariance Evolution
#'
#' This plot makes proposal adaptation visible. For one-dimensional runs it
#' shows the proposal variance directly. For higher-dimensional runs it can show
#' either the diagonal entries of `S_n` or the eigenvalues of `S_n`.
#'
#' @param run A sampler result created by this package.
#' @param view Either `"diagonal"` or `"eigenvalues"`.
#'
#' @return Invisibly returns the matrix of plotted values.
#'
#' @examples
#' comparison <- run_sampler_suite(make_correlated_t_target(), n_iter = 250, seed = 1)
#' plot_covariance_evolution(comparison$runs$adaptive)
#'
#' @export
plot_covariance_evolution <- function(run, view = c("diagonal", "eigenvalues")) {
  run <- validate_sampler_run(run)
  view <- match.arg(view)
  history <- run$proposal_covariance_history
  n_iter <- dim(history)[3]
  dimension <- dim(history)[1]

  raw_values <- if (view == "diagonal") {
    vapply(
      seq_len(n_iter),
      function(index) {
        covariance <- history[, , index, drop = TRUE]

        if (dimension == 1L) {
          return(as.numeric(covariance))
        }

        diag(covariance)
      },
      numeric(dimension)
    )
  } else {
    vapply(
      seq_len(n_iter),
      function(index) {
        covariance <- history[, , index, drop = TRUE]

        if (dimension == 1L) {
          return(as.numeric(covariance))
        }

        sort(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values)
      },
      numeric(dimension)
    )
  }

  values <- if (dimension == 1L) {
    matrix(raw_values, ncol = 1L)
  } else {
    t(raw_values)
  }

  graphics::matplot(
    x = seq_len(n_iter),
    y = values,
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "Iteration",
    ylab = if (view == "diagonal") "Diagonal entry of S_n" else "Eigenvalue of S_n",
    main = sprintf("Proposal covariance evolution: %s", run$algorithm)
  )

  invisible(values)
}

#' Plot Mode-Switching Behavior
#'
#' This diagnostic is mainly for the multimodal teaching example.
#'
#' @param run A sampler result created by this package.
#' @param target Optional target object.
#' @param mode_centers Optional matrix of mode centers, one row per mode.
#' @param burn_in Number of initial iterations discarded.
#'
#' @return Invisibly returns a data frame containing the assigned mode label and
#'   cumulative switch count over time.
#'
#' @examples
#' target <- make_multimodal_target()
#' run <- run_metropolis(
#'   target_density = target$density,
#'   initial_state = target$recommended_initial,
#'   n_iter = 250,
#'   proposal_covariance = target$recommended_proposal_covariance
#' )
#' plot_mode_switching_behavior(run, target = target)
#'
#' @export
plot_mode_switching_behavior <- function(run,
                                         target = NULL,
                                         mode_centers = NULL,
                                         burn_in = 0L) {
  run <- validate_sampler_run(run)
  draws <- draws_after_burn_in(run, burn_in = burn_in)

  if (!is.null(target)) {
    target <- validate_target(target)
    mode_centers <- mode_centers %||% target$mode_centers
  }

  if (is.null(mode_centers)) {
    if (ncol(draws) != 1L) {
      stop(
        "Provide `mode_centers` for multi-dimensional runs, or pass a target that contains them.",
        call. = FALSE
      )
    }

    mode_centers <- matrix(c(-1, 1), ncol = 1)
  }

  mode_centers <- as.matrix(mode_centers)

  if (ncol(mode_centers) != ncol(draws)) {
    stop("`mode_centers` must have one column per chain coordinate.", call. = FALSE)
  }

  assigned_mode <- vapply(
    seq_len(nrow(draws)),
    function(index) {
      centered_modes <- sweep(mode_centers, 2, draws[index, ], FUN = "-")
      distances <- rowSums(centered_modes^2)
      which.min(distances)
    },
    integer(1)
  )

  switch_indicator <- c(FALSE, diff(assigned_mode) != 0)
  cumulative_switches <- cumsum(switch_indicator)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  set_compact_plot_layout(2L)

  graphics::plot(
    assigned_mode,
    type = "s",
    xlab = "Iteration",
    ylab = "Assigned mode",
    main = sprintf("Mode labels: %s", run$algorithm)
  )

  graphics::plot(
    cumulative_switches,
    type = "l",
    xlab = "Iteration",
    ylab = "Cumulative switches",
    main = sprintf("Mode switching count: %s", run$algorithm)
  )

  invisible(data.frame(
    iteration = seq_len(nrow(draws)),
    assigned_mode = assigned_mode,
    switched = switch_indicator,
    cumulative_switches = cumulative_switches
  ))
}
