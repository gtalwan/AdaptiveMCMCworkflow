# Internal target constructor used by the exported teaching targets.
# It keeps the object layout consistent without forcing learners to read a
# deeply abstract class system.
#
#' @noRd
new_teaching_target <- function(density,
                                dimension,
                                name,
                                description,
                                scenario,
                                recommended_initial = NULL,
                                recommended_proposal_covariance = NULL,
                                recommended_adapt_start = 50L,
                                recommended_delayed_rejection_scale = 0.5,
                                plot_limits = NULL,
                                mean = NULL,
                                covariance = NULL,
                                target_type = "custom",
                                mode_centers = NULL) {
  if (!is.function(density)) {
    stop("`density` must be a function.", call. = FALSE)
  }

  dimension <- ensure_positive_integer(dimension, "dimension")
  recommended_adapt_start <- ensure_positive_integer(
    recommended_adapt_start,
    "recommended_adapt_start"
  )
  recommended_delayed_rejection_scale <- ensure_positive_number(
    recommended_delayed_rejection_scale,
    "recommended_delayed_rejection_scale"
  )

  if (!is.null(recommended_initial)) {
    recommended_initial <- ensure_numeric_state(recommended_initial)

    if (length(recommended_initial) != dimension) {
      stop(
        "`recommended_initial` must match the target dimension.",
        call. = FALSE
      )
    }
  }

  if (!is.null(recommended_proposal_covariance)) {
    recommended_proposal_covariance <- ensure_square_matrix(
      recommended_proposal_covariance,
      dimension,
      "recommended_proposal_covariance"
    )
  }

  if (!is.null(mode_centers)) {
    mode_centers <- as.matrix(mode_centers)

    if (!is.numeric(mode_centers) || ncol(mode_centers) != dimension) {
      stop(
        "`mode_centers` must have one column per target dimension.",
        call. = FALSE
      )
    }
  }

  structure(
    list(
      name = name,
      description = description,
      scenario = scenario,
      target_type = target_type,
      dimension = dimension,
      mean = mean,
      covariance = covariance,
      density = density,
      log_density = function(x) log(density(x)),
      recommended_initial = recommended_initial,
      recommended_proposal_covariance = recommended_proposal_covariance,
      recommended_adapt_start = recommended_adapt_start,
      recommended_delayed_rejection_scale = recommended_delayed_rejection_scale,
      plot_limits = plot_limits,
      mode_centers = mode_centers
    ),
    class = "adaptive_mcmc_target"
  )
}

#' Create a Gaussian Teaching Target
#'
#' This constructor creates a multivariate Gaussian target object. The object
#' stores a density closure and enough metadata to keep the teaching
#' examples reproducible.
#'
#' @param mean Numeric mean vector of the Gaussian target.
#' @param covariance Positive-definite covariance matrix of the Gaussian
#'   target.
#' @param name Short display name used in summaries and plots.
#' @param description Plain-language description of why this target is useful
#'   in the teaching examples.
#' @param scenario Short scenario label used by the teaching examples.
#' @param recommended_initial Optional starting point used in the teaching
#'   examples.
#' @param recommended_proposal_covariance Optional starting proposal covariance
#'   used in the teaching examples.
#' @param recommended_adapt_start Iteration at which adaptation is recommended
#'   to begin in the teaching examples.
#' @param recommended_delayed_rejection_scale Second-stage scale factor
#'   recommended for the DRAM example.
#' @param plot_limits Optional plot limits. For a two-dimensional target this
#'   should be a list with components `x` and `y`.
#' @param target_type Short internal label describing the target family.
#' @param mode_centers Optional matrix of mode centres, one row per mode.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' gaussian_target <- make_gaussian_target(
#'   mean = c(0, 0),
#'   covariance = matrix(c(1, 0.8, 0.8, 1), nrow = 2),
#'   name = "Example Gaussian",
#'   description = "A simple correlated Gaussian target."
#' )
#' gaussian_target$density(c(0, 0))
#'
#' @export
make_gaussian_target <- function(mean,
                                 covariance,
                                 name = "Gaussian target",
                                 description = "A Gaussian target for MCMC examples.",
                                 scenario = "gaussian",
                                 recommended_initial = NULL,
                                 recommended_proposal_covariance = NULL,
                                 recommended_adapt_start = 50L,
                                 recommended_delayed_rejection_scale = 0.5,
                                 plot_limits = NULL,
                                 target_type = "gaussian",
                                 mode_centers = NULL) {
  mean <- as.numeric(mean)

  if (length(mean) == 0L) {
    stop("`mean` must have at least one entry.", call. = FALSE)
  }

  covariance <- ensure_square_matrix(covariance, length(mean), "covariance")

  if (inherits(try(chol(covariance), silent = TRUE), "try-error")) {
    stop("`covariance` must be positive definite.", call. = FALSE)
  }

  precision <- solve(covariance)
  determinant_value <- as.numeric(determinant(covariance, logarithm = FALSE)$modulus)
  dimension <- length(mean)
  marginal_sd <- sqrt(diag(covariance))

  density <- function(x) {
    x <- as.numeric(x)

    if (length(x) != dimension) {
      stop(sprintf("The state must have length %d.", dimension), call. = FALSE)
    }

    centered <- x - mean
    quadratic_form <- sum(centered * drop(precision %*% centered))
    normalising_constant <- 1 / sqrt((2 * pi)^dimension * determinant_value)
    normalising_constant * exp(-0.5 * quadratic_form)
  }

  if (is.null(plot_limits)) {
    if (dimension == 1L) {
      plot_limits <- c(mean - 4 * marginal_sd, mean + 4 * marginal_sd)
    } else if (dimension == 2L) {
      plot_limits <- list(
        x = c(mean[1] - 4 * marginal_sd[1], mean[1] + 4 * marginal_sd[1]),
        y = c(mean[2] - 4 * marginal_sd[2], mean[2] + 4 * marginal_sd[2])
      )
    }
  }

  new_teaching_target(
    density = density,
    dimension = dimension,
    name = name,
    description = description,
    scenario = scenario,
    recommended_initial = recommended_initial,
    recommended_proposal_covariance = recommended_proposal_covariance,
    recommended_adapt_start = recommended_adapt_start,
    recommended_delayed_rejection_scale = recommended_delayed_rejection_scale,
    plot_limits = plot_limits,
    mean = mean,
    covariance = covariance,
    target_type = target_type,
    mode_centers = mode_centers
  )
}

#' Standard Normal Teaching Target
#'
#' This is the cleanest one-dimensional baseline in the package. It is useful
#' when you want to explain the mechanics of the samplers without also worrying
#' about correlation or multimodality.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' target <- make_standard_normal_target()
#' target$density(0)
#'
#' @export
make_standard_normal_target <- function() {
  make_gaussian_target(
    mean = 0,
    covariance = matrix(1, nrow = 1),
    name = "Standard normal target",
    description = paste(
      "A one-dimensional standard Gaussian used when the goal is to explain",
      "the proposal, acceptance ratio, and adaptation steps one line at a",
      "time without extra geometric complications."
    ),
    scenario = "standard_normal",
    recommended_initial = 3,
    recommended_proposal_covariance = matrix(1, nrow = 1),
    recommended_adapt_start = 30L,
    recommended_delayed_rejection_scale = 0.5,
    plot_limits = c(-4, 4)
  )
}

#' Correlated Gaussian Ridge Target
#'
#' This two-dimensional Gaussian target stays in the package because it is the
#' cleanest first example of why covariance adaptation matters. If you want a
#' non-Gaussian ridge, use `make_correlated_t_target()`.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' ridge_target <- make_correlated_gaussian_target()
#' ridge_target$density(c(0, 0))
#'
#' @export
make_correlated_gaussian_target <- function() {
  make_gaussian_target(
    mean = c(0, 0),
    covariance = matrix(c(1, 0.98, 0.98, 1), nrow = 2),
    name = "Correlated Gaussian ridge",
    description = paste(
      "A strongly correlated two-dimensional Gaussian.",
      "It is the simplest clean example of how a fixed isotropic proposal",
      "struggles to move along a narrow ridge."
    ),
    scenario = "correlated_gaussian",
    recommended_initial = c(-3, 3),
    recommended_proposal_covariance = diag(2),
    recommended_adapt_start = 60L,
    recommended_delayed_rejection_scale = 0.5,
    plot_limits = list(x = c(-5, 5), y = c(-5, 5)),
    mode_centers = matrix(c(0, 0), nrow = 1)
  )
}

#' Correlated Student-t Ridge Target
#'
#' This target answers the obvious teaching question: the target does not have
#' to be Gaussian. The geometry is still a correlated ridge, but the heavier
#' tails make the target visibly non-Gaussian while keeping the mathematics and
#' contour plots easy to explain.
#'
#' @param df Degrees of freedom of the multivariate Student-t target. Smaller
#'   values give heavier tails.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' target <- make_correlated_t_target()
#' target$density(c(0, 0))
#'
#' @export
make_correlated_t_target <- function(df = 4) {
  df <- ensure_positive_number(df, "df")
  scale_matrix <- matrix(c(1, 0.97, 0.97, 1), nrow = 2)
  precision <- solve(scale_matrix)

  density <- function(x) {
    x <- as.numeric(x)

    if (length(x) != 2L) {
      stop("The correlated Student-t target is two-dimensional.", call. = FALSE)
    }

    quadratic_form <- sum(x * drop(precision %*% x))
    (1 + quadratic_form / df)^(-0.5 * (df + 2))
  }

  new_teaching_target(
    density = density,
    dimension = 2L,
    name = sprintf("Correlated Student-t ridge (df = %s)", format(df)),
    description = paste(
      "A non-Gaussian correlated ridge with heavier tails than a Gaussian.",
      "It is useful when you want the adaptive Metropolis demonstrations to",
      "show that covariance learning is about geometry, not about assuming",
      "the target itself is Gaussian."
    ),
    scenario = "correlated_student_t",
    recommended_initial = c(-3, 3),
    recommended_proposal_covariance = 0.25 * diag(2),
    recommended_adapt_start = 20L,
    recommended_delayed_rejection_scale = 0.5,
    plot_limits = list(x = c(-6, 6), y = c(-6, 6)),
    covariance = scale_matrix * df / (df - 2),
    target_type = "student_t",
    mode_centers = matrix(c(0, 0), nrow = 1)
  )
}

#' Poor Proposal Scale Target
#'
#' This one-dimensional Gaussian target is paired with an intentionally tiny
#' proposal variance. The result is usually a very high acceptance rate but very
#' slow exploration, which makes the example useful for discussing
#' autocorrelation.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' small_scale_target <- make_small_scale_target()
#' small_scale_target$recommended_proposal_covariance
#'
#' @export
make_small_scale_target <- function() {
  make_gaussian_target(
    mean = 0,
    covariance = matrix(1, nrow = 1),
    name = "Poor proposal scale example",
    description = paste(
      "A one-dimensional standard Gaussian target paired with a very small",
      "random-walk proposal variance. The example highlights the difference",
      "between high acceptance and efficient exploration."
    ),
    scenario = "small_scale",
    recommended_initial = 4,
    recommended_proposal_covariance = matrix(0.0025, nrow = 1),
    recommended_adapt_start = 40L,
    recommended_delayed_rejection_scale = 0.5,
    plot_limits = c(-4, 4),
    mode_centers = matrix(0, ncol = 1)
  )
}

#' Multimodal Gaussian Mixture Target
#'
#' This target is a simple symmetric Gaussian mixture. It is intentionally used
#' to show a limitation as well as a strength: adaptive random-walk methods can
#' improve local scale and covariance choices, but they do not automatically
#' solve multimodal trapping.
#'
#' @param separation Distance between the two component means.
#' @param component_sd Standard deviation of each Gaussian mixture component.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' target <- make_multimodal_target()
#' target$density(0)
#'
#' @export
make_multimodal_target <- function(separation = 4,
                                   component_sd = 1) {
  separation <- ensure_positive_number(separation, "separation")
  component_sd <- ensure_positive_number(component_sd, "component_sd")

  mean_left <- -separation
  mean_right <- separation
  variance <- component_sd^2

  density <- function(x) {
    x <- as.numeric(x)

    if (length(x) != 1L) {
      stop("The multimodal teaching target is one-dimensional.", call. = FALSE)
    }

    0.5 * stats::dnorm(x, mean = mean_left, sd = component_sd) +
      0.5 * stats::dnorm(x, mean = mean_right, sd = component_sd)
  }

  new_teaching_target(
    density = density,
    dimension = 1L,
    name = "Multimodal Gaussian mixture",
    description = paste(
      "A symmetric two-mode Gaussian mixture. It is useful for showing that",
      "adaptive random-walk tuning can help local efficiency without solving",
      "the basic difficulty of moving between separated modes."
    ),
    scenario = "multimodal",
    recommended_initial = mean_left,
    recommended_proposal_covariance = matrix(1, nrow = 1),
    recommended_adapt_start = 40L,
    recommended_delayed_rejection_scale = 0.5,
    plot_limits = c(mean_left - 4 * component_sd, mean_right + 4 * component_sd),
    mean = 0,
    covariance = matrix(separation^2 + variance, nrow = 1),
    target_type = "gaussian_mixture",
    mode_centers = matrix(c(mean_left, mean_right), ncol = 1)
  )
}

#' Delayed Rejection Teaching Target
#'
#' This two-dimensional Gaussian target is paired with a deliberately ambitious
#' first-stage proposal covariance. The setup is designed so that DRAM has a
#' genuine opportunity to rescue rejected large moves with a smaller second
#' stage.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' dram_target <- make_dram_rescue_target()
#' dram_target$recommended_delayed_rejection_scale
#'
#' @export
make_dram_rescue_target <- function() {
  base_covariance <- matrix(c(1, 0.7, 0.7, 1), nrow = 2)

  make_gaussian_target(
    mean = c(0, 0),
    covariance = base_covariance,
    name = "Delayed rejection example",
    description = paste(
      "A moderately correlated Gaussian target used with an intentionally",
      "large first-stage proposal covariance so that DRAM can show the value",
      "of a smaller rescue proposal after rejection."
    ),
    scenario = "dram_rescue",
    recommended_initial = c(3, -3),
    recommended_proposal_covariance = 6 * base_covariance,
    recommended_adapt_start = 50L,
    recommended_delayed_rejection_scale = 0.35,
    plot_limits = list(x = c(-5, 5), y = c(-5, 5)),
    mode_centers = matrix(c(0, 0), nrow = 1)
  )
}

#' Higher-Dimensional Gaussian Teaching Target
#'
#' This constructor creates a correlated Gaussian target in dimension `d`. It is
#' meant for demonstrating why manual proposal tuning becomes awkward in higher
#' dimensions and why RAM's acceptance-rate targeting can be attractive there.
#'
#' @param dimension Target dimension.
#' @param rho Correlation parameter used in the AR(1)-style covariance matrix.
#'
#' @return An object of class `"adaptive_mcmc_target"`.
#'
#' @examples
#' target <- make_high_dimensional_gaussian_target(dimension = 5)
#' target$dimension
#'
#' @export
make_high_dimensional_gaussian_target <- function(dimension = 8L,
                                                  rho = 0.7) {
  dimension <- ensure_positive_integer(dimension, "dimension")
  rho <- as.numeric(rho)

  if (length(rho) != 1L || is.na(rho) || abs(rho) >= 1) {
    stop("`rho` must be a single number strictly between -1 and 1.", call. = FALSE)
  }

  covariance <- outer(
    seq_len(dimension),
    seq_len(dimension),
    function(i, j) rho^abs(i - j)
  )

  make_gaussian_target(
    mean = rep(0, dimension),
    covariance = covariance,
    name = sprintf("High-dimensional Gaussian (d = %d)", dimension),
    description = paste(
      "A correlated Gaussian target in moderate dimension.",
      "It is used to show that manual covariance tuning becomes increasingly",
      "awkward as dimension grows, and to contrast AM's empirical covariance",
      "learning with RAM's acceptance-rate targeting."
    ),
    scenario = "high_dimensional_gaussian",
    recommended_initial = rep(3, dimension),
    recommended_proposal_covariance = diag(dimension),
    recommended_adapt_start = 80L,
    recommended_delayed_rejection_scale = 0.5
  )
}
