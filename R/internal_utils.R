# Internal helper functions live here. They are small and explicit because the
# package is meant to be read by learners. The goal is to keep argument checks
# and common bookkeeping tidy without hiding the statistical ideas.

#' @noRd
`%||%` <- function(left, right) {
  if (is.null(left)) {
    right
  } else {
    left
  }
}

#' @noRd
ensure_numeric_state <- function(initial_state) {
  state <- as.numeric(initial_state)

  if (length(state) == 0L) {
    stop("`initial_state` must contain at least one coordinate.", call. = FALSE)
  }

  if (any(!is.finite(state))) {
    stop("`initial_state` must be finite.", call. = FALSE)
  }

  state
}

#' @noRd
ensure_positive_integer <- function(value, name) {
  integer_value <- as.integer(value)

  if (length(integer_value) != 1L || is.na(integer_value) || integer_value <= 0L) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }

  integer_value
}

#' @noRd
ensure_positive_number <- function(value, name) {
  numeric_value <- as.numeric(value)

  if (length(numeric_value) != 1L || is.na(numeric_value) || numeric_value <= 0) {
    stop(sprintf("`%s` must be a positive number.", name), call. = FALSE)
  }

  numeric_value
}

#' @noRd
ensure_probability <- function(value, name) {
  numeric_value <- as.numeric(value)

  if (length(numeric_value) != 1L || is.na(numeric_value) ||
      numeric_value <= 0 || numeric_value >= 1) {
    stop(sprintf("`%s` must be strictly between 0 and 1.", name), call. = FALSE)
  }

  numeric_value
}

#' @noRd
ensure_target_density <- function(target_density) {
  if (!is.function(target_density)) {
    stop(
      "`target_density` must be a function that returns one density value.",
      call. = FALSE
    )
  }

  target_density
}

#' @noRd
resolve_target_density <- function(target_density = NULL, log_target = NULL) {
  if (!is.null(target_density) && !is.null(log_target)) {
    stop(
      "Supply either `target_density` or `log_target`, not both.",
      call. = FALSE
    )
  }

  if (!is.null(target_density)) {
    return(ensure_target_density(target_density))
  }

  if (!is.null(log_target)) {
    if (!is.function(log_target)) {
      stop("`log_target` must be a function.", call. = FALSE)
    }

    return(function(x) exp(log_target(x)))
  }

  stop(
    "Supply `target_density`. The old `log_target` argument is only kept as a compatibility alias.",
    call. = FALSE
  )
}

#' @noRd
ensure_square_matrix <- function(value, dimension, name) {
  matrix_value <- as.matrix(value)

  if (!is.numeric(matrix_value)) {
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  }

  if (!all(dim(matrix_value) == c(dimension, dimension))) {
    stop(
      sprintf("`%s` must be a %d by %d matrix.", name, dimension, dimension),
      call. = FALSE
    )
  }

  if (any(!is.finite(matrix_value))) {
    stop(sprintf("`%s` must be finite.", name), call. = FALSE)
  }

  if (max(abs(matrix_value - t(matrix_value))) > 1e-10) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }

  matrix_value
}

#' @noRd
check_initial_density <- function(target_density, initial_state) {
  value <- target_density(initial_state)

  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
    stop(
      "The target density must return one positive finite numeric value at `initial_state`.",
      call. = FALSE
    )
  }

  as.numeric(value)
}

#' @noRd
standardise_sampler_output <- function(result, initial_state, target_name = NULL) {
  result$algorithm <- as.character(result$algorithm)
  result$draws <- as.matrix(result$draws)
  result$accepted <- as.logical(result$accepted)

  for (field in c("stage_one_accepted", "stage_two_accepted", "second_stage_attempted")) {
    if (field %in% names(result)) {
      result[[field]] <- as.logical(result[[field]])
    }
  }

  result$acceptance_rate <- as.numeric(result$acceptance_rate)
  if ("target_values" %in% names(result)) {
    result$target_values <- as.numeric(result$target_values)
  }

  result$proposal_covariance_history <- as.array(result$proposal_covariance_history)
  result$adaptation_magnitude <- as.numeric(result$adaptation_magnitude)
  result$jump_distance <- as.numeric(result$jump_distance)
  result$initial_state <- as.numeric(initial_state)
  result$n_iter <- nrow(result$draws)
  result$dimension <- ncol(result$draws)
  result$target_name <- target_name

  class(result) <- c("adaptive_mcmc_run", "list")
  result
}

#' @noRd
validate_target <- function(target) {
  if (!inherits(target, "adaptive_mcmc_target")) {
    stop("`target` must be a target object created by this package.", call. = FALSE)
  }

  target
}

#' @noRd
default_sampler_controls <- function(target) {
  target <- validate_target(target)

  proposal_covariance <-
    target$recommended_proposal_covariance %||%
    diag(target$dimension)

  adapt_start <- target$recommended_adapt_start %||% 50L
  delayed_scale <- target$recommended_delayed_rejection_scale %||% 0.5

  list(
    basic = list(
      proposal_covariance = proposal_covariance
    ),
    adaptive = list(
      initial_proposal_covariance = proposal_covariance,
      adapt_start = adapt_start,
      epsilon = 1e-8
    ),
    ram = list(
      initial_proposal_covariance = proposal_covariance,
      target_acceptance = 0.234,
      adapt_exponent = 0.67
    ),
    dram = list(
      initial_proposal_covariance = proposal_covariance,
      adapt_start = adapt_start,
      delayed_rejection_scale = delayed_scale,
      epsilon = 1e-8
    )
  )
}

#' @noRd
merge_sampler_controls <- function(default_controls, user_controls = NULL) {
  if (is.null(user_controls)) {
    return(default_controls)
  }

  utils::modifyList(default_controls, user_controls)
}

#' @noRd
validate_sampler_run <- function(run) {
  if (!inherits(run, "adaptive_mcmc_run")) {
    stop("`run` must be a sampler result created by this package.", call. = FALSE)
  }

  run
}

#' @noRd
validate_sampler_comparison <- function(comparison) {
  if (!inherits(comparison, "adaptive_mcmc_comparison")) {
    stop(
      "`comparison` must be the output of `run_sampler_suite()`.",
      call. = FALSE
    )
  }

  comparison
}

#' @noRd
draws_after_burn_in <- function(run, burn_in = 0L) {
  run <- validate_sampler_run(run)
  burn_in <- as.integer(burn_in)

  if (length(burn_in) != 1L || is.na(burn_in) || burn_in < 0L || burn_in >= run$n_iter) {
    stop("`burn_in` must be between 0 and `n_iter - 1`.", call. = FALSE)
  }

  run$draws[seq.int(from = burn_in + 1L, to = run$n_iter), , drop = FALSE]
}

#' @noRd
acf_values_internal <- function(series, lag_max) {
  stats::acf(
    as.numeric(series),
    lag.max = lag_max,
    plot = FALSE,
    demean = TRUE
  )$acf[-1L]
}

#' @noRd
effective_sample_size_internal <- function(series, lag_max = NULL) {
  series <- as.numeric(series)
  sample_size <- length(series)

  if (sample_size < 3L) {
    return(NA_real_)
  }

  lag_max <- lag_max %||% min(100L, sample_size - 1L)
  lag_max <- min(as.integer(lag_max), sample_size - 1L)

  if (lag_max < 1L) {
    return(as.numeric(sample_size))
  }

  autocorrelations <- acf_values_internal(series, lag_max = lag_max)

  if (!length(autocorrelations)) {
    return(as.numeric(sample_size))
  }

  first_non_positive <- match(TRUE, autocorrelations <= 0, nomatch = length(autocorrelations) + 1L)

  positive_part <- if (first_non_positive == 1L) {
    numeric(0)
  } else {
    autocorrelations[seq_len(first_non_positive - 1L)]
  }

  integrated_autocorrelation_time <- 1 + 2 * sum(positive_part)
  sample_size / integrated_autocorrelation_time
}

#' @noRd
choose_plot_layout <- function(n_panels) {
  n_panels <- as.integer(n_panels)

  if (length(n_panels) != 1L || is.na(n_panels) || n_panels <= 0L) {
    stop("`n_panels` must be a positive integer.", call. = FALSE)
  }

  # Base graphics can fail with "figure margins too large" when a small device
  # is split into many panels using the default margins. We inspect the current
  # device size and choose a layout that is a little more forgiving in narrow
  # RStudio panes.
  device_size <- tryCatch(grDevices::dev.size("in"), error = function(...) c(7, 7))
  width <- device_size[1]
  height <- device_size[2]

  if (n_panels == 1L) {
    return(c(1L, 1L))
  }

  if (n_panels == 2L) {
    if (width >= height) {
      return(c(1L, 2L))
    }

    return(c(2L, 1L))
  }

  if (n_panels == 4L) {
    aspect_ratio <- width / height

    if (aspect_ratio >= 1.7) {
      return(c(1L, 4L))
    }

    if (aspect_ratio <= 0.7) {
      return(c(4L, 1L))
    }

    return(c(2L, 2L))
  }

  c(n_panels, 1L)
}

#' @noRd
set_compact_plot_layout <- function(n_panels) {
  layout <- choose_plot_layout(n_panels)

  # These margins are smaller than the base defaults but still leave enough
  # room for readable axis labels and titles in typical teaching plots.
  graphics::par(
    mfrow = layout,
    mar = c(3.2, 3.4, 2.0, 0.8),
    mgp = c(2.0, 0.7, 0),
    oma = c(0, 0, 0.4, 0),
    cex.main = 0.95
  )

  invisible(layout)
}
