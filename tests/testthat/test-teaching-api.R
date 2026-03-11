test_that("teaching wrappers and summaries are reproducible and usable", {
  target <- make_standard_normal_target()

  set.seed(123)
  wrapper_run <- run_metropolis(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 80,
    proposal_covariance = target$recommended_proposal_covariance
  )

  set.seed(123)
  direct_run <- basic_metropolis(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 80,
    proposal_covariance = target$recommended_proposal_covariance
  )

  expect_equal(wrapper_run$draws, direct_run$draws)
  expect_equal(wrapper_run$accepted, direct_run$accepted)

  set.seed(456)
  adaptive_one <- run_adaptive_metropolis(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 80,
    initial_proposal_covariance = target$recommended_proposal_covariance,
    adapt_start = 20
  )

  set.seed(456)
  adaptive_two <- run_adaptive_metropolis(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 80,
    initial_proposal_covariance = target$recommended_proposal_covariance,
    adapt_start = 20
  )

  expect_equal(adaptive_one$draws, adaptive_two$draws)
  expect_equal(adaptive_one$accepted, adaptive_two$accepted)
  expect_true(all(c("algorithm", "coordinate", "effective_sample_size") %in%
    names(effective_sample_size_summary(adaptive_one))))
})

test_that("each sampler accepts a plain R target density function", {
  standard_normal_density <- function(x) {
    x <- as.numeric(x)
    exp(-0.5 * sum(x^2))
  }

  basic_run <- basic_metropolis(
    target_density = standard_normal_density,
    initial_state = c(0, 0),
    n_iter = 40,
    proposal_covariance = diag(2)
  )

  adaptive_run <- adaptive_metropolis(
    target_density = standard_normal_density,
    initial_state = c(0, 0),
    n_iter = 40,
    initial_proposal_covariance = diag(2),
    adapt_start = 20
  )

  ram_run <- robust_adaptive_metropolis(
    target_density = standard_normal_density,
    initial_state = c(0, 0),
    n_iter = 40,
    initial_proposal_covariance = diag(2)
  )

  dram_run <- dram(
    target_density = standard_normal_density,
    initial_state = c(0, 0),
    n_iter = 40,
    initial_proposal_covariance = diag(2),
    adapt_start = 20
  )

  expect_s3_class(basic_run, "adaptive_mcmc_run")
  expect_s3_class(adaptive_run, "adaptive_mcmc_run")
  expect_s3_class(ram_run, "adaptive_mcmc_run")
  expect_s3_class(dram_run, "adaptive_mcmc_run")
})

test_that("additional targets and high-level comparison helpers return expected objects", {
  multimodal_target <- make_multimodal_target()
  high_dimensional_target <- make_high_dimensional_gaussian_target(dimension = 4L)

  expect_s3_class(multimodal_target, "adaptive_mcmc_target")
  expect_true(is.finite(multimodal_target$log_density(0)))
  expect_equal(high_dimensional_target$dimension, 4L)
  expect_true(is.finite(high_dimensional_target$log_density(rep(0, 4))))

  comparison_result <- compare_adaptive_methods(
    target = high_dimensional_target,
    n_iter = 100,
    burn_in = 20,
    seed = 1,
    benchmark = FALSE
  )

  expect_s3_class(comparison_result, "adaptive_method_comparison")
  expect_s3_class(comparison_result$comparison, "adaptive_mcmc_comparison")
  expect_true(all(c("basic", "adaptive", "ram", "dram") %in% names(comparison_result$comparison$runs)))
  expect_true(all(c("algorithm", "effective_sample_size") %in% names(comparison_result$summary)))

  ram_run <- run_ram(
    log_target = high_dimensional_target$log_density,
    initial_state = high_dimensional_target$recommended_initial,
    n_iter = 80,
    initial_proposal_covariance = high_dimensional_target$recommended_proposal_covariance
  )
  dram_run <- run_dram(
    log_target = high_dimensional_target$log_density,
    initial_state = high_dimensional_target$recommended_initial,
    n_iter = 80,
    initial_proposal_covariance = high_dimensional_target$recommended_proposal_covariance
  )

  expect_s3_class(ram_run, "adaptive_mcmc_run")
  expect_s3_class(dram_run, "adaptive_mcmc_run")

  validity_summary <- summarize_adaptive_validity(comparison_result$comparison)

  expect_equal(nrow(validity_summary), 3L)
  expect_true(all(c(
    "algorithm",
    "kernel_change_ratio",
    "max_proposal_condition_number"
  ) %in% names(validity_summary)))
})

test_that("samplers remain stable under early adaptation and pathological rejection periods", {
  rejection_target <- function(x) {
    if (all(abs(as.numeric(x)) < 1e-12)) {
      return(0)
    }

    -Inf
  }

  rejection_run <- run_metropolis(
    log_target = rejection_target,
    initial_state = 0,
    n_iter = 20,
    proposal_covariance = matrix(1, nrow = 1)
  )

  expect_false(any(rejection_run$accepted))
  expect_true(all(rejection_run$draws == 0))
  expect_true(all(rejection_run$acceptance_probability == 0))

  ridge_target <- make_correlated_gaussian_target()
  near_singular_covariance <- matrix(c(1, 0.999999, 0.999999, 1), nrow = 2)

  adaptive_run <- run_adaptive_metropolis(
    log_target = ridge_target$log_density,
    initial_state = ridge_target$recommended_initial,
    n_iter = 25,
    initial_proposal_covariance = near_singular_covariance,
    adapt_start = 2,
    regularization = 1e-6
  )

  smallest_eigenvalue <- vapply(
    seq_len(dim(adaptive_run$proposal_covariance_history)[3]),
    function(index) {
      min(eigen(adaptive_run$proposal_covariance_history[, , index, drop = TRUE],
        symmetric = TRUE,
        only.values = TRUE
      )$values)
    },
    numeric(1)
  )

  expect_true(all(is.finite(adaptive_run$adaptation_magnitude)))
  expect_true(all(is.finite(smallest_eigenvalue)))
  expect_true(all(smallest_eigenvalue > 0))
})

test_that("new single-run plotting helpers work on compact graphics devices", {
  multimodal_target <- make_multimodal_target()
  one_dimensional_target <- make_standard_normal_target()

  multimodal_run <- run_metropolis(
    log_target = multimodal_target$log_density,
    initial_state = multimodal_target$recommended_initial,
    n_iter = 120,
    proposal_covariance = multimodal_target$recommended_proposal_covariance
  )
  one_dimensional_adaptive_run <- run_adaptive_metropolis(
    log_target = one_dimensional_target$log_density,
    initial_state = one_dimensional_target$recommended_initial,
    n_iter = 120,
    initial_proposal_covariance = one_dimensional_target$recommended_proposal_covariance,
    adapt_start = one_dimensional_target$recommended_adapt_start
  )

  ridge_comparison <- compare_adaptive_methods(
    target = make_correlated_gaussian_target(),
    n_iter = 120,
    burn_in = 20,
    seed = 1,
    benchmark = FALSE
  )

  output_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(output_file, width = 5, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(plot_acf_chain(multimodal_run, dim = 1, lag_max = 15, burn_in = 20))
  expect_no_error(plot_acceptance(ridge_comparison$comparison$runs$ram))
  expect_no_error(plot_covariance_evolution(ridge_comparison$comparison$runs$adaptive))
  expect_no_error(plot_covariance_evolution(one_dimensional_adaptive_run))
  expect_no_error(plot_mode_switching_behavior(multimodal_run, target = multimodal_target, burn_in = 20))
})

test_that("multi-panel plotting helpers work on a compact graphics device", {
  comparison <- run_sampler_suite(make_correlated_gaussian_target(), n_iter = 120, seed = 1)

  output_file <- tempfile(fileext = ".pdf")

  grDevices::pdf(output_file, width = 5, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(plot_sampler_comparison(comparison, burn_in = 20))
  expect_no_error(plot_adaptation_diagnostic(comparison$runs$adaptive))
  expect_no_error(plot_containment_diagnostic(comparison$runs$adaptive))
  expect_no_error(plot_trace(comparison$runs$adaptive, burn_in = 20))
})
