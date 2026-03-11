test_that("target constructors return usable log densities", {
  ridge_target <- make_correlated_gaussian_target()
  student_t_target <- make_correlated_t_target()
  small_target <- make_small_scale_target()
  dram_target <- make_dram_rescue_target()

  expect_s3_class(ridge_target, "adaptive_mcmc_target")
  expect_true(is.finite(ridge_target$log_density(c(0, 0))))
  expect_true(is.finite(student_t_target$log_density(c(0, 0))))
  expect_true(is.finite(small_target$log_density(0)))
  expect_true(is.finite(dram_target$log_density(c(0, 0))))
})

test_that("suite summaries and diagnostics run on canonical examples", {
  comparison <- run_sampler_suite(make_correlated_gaussian_target(), n_iter = 120, seed = 1)

  summary_table <- summarize_sampler_suite(comparison, burn_in = 20)
  diminishing <- diminishing_adaptation_diagnostic(comparison$runs$adaptive)
  containment <- containment_diagnostic(comparison$runs$adaptive)
  benchmark <- benchmark_sampler_suite(make_small_scale_target(), n_iter = 60, repetitions = 2, seed = 1)

  expect_true(all(c("basic", "adaptive", "ram", "dram") %in% names(comparison$runs)))
  expect_true(all(c("algorithm", "effective_sample_size") %in% names(summary_table)))
  expect_equal(nrow(diminishing), comparison$runs$adaptive$n_iter)
  expect_equal(nrow(containment), comparison$runs$adaptive$n_iter)
  expect_true(all(c("raw", "summary") %in% names(benchmark)))
})

test_that("the non-Gaussian ridge target works with the comparison helpers", {
  comparison_result <- compare_adaptive_methods(
    make_correlated_t_target(),
    n_iter = 500,
    burn_in = 100,
    seed = 1,
    benchmark = FALSE
  )

  summary_table <- comparison_result$summary
  adaptive_summary <- summary_table[summary_table$algorithm == "Adaptive Metropolis", ]
  basic_summary <- summary_table[summary_table$algorithm == "Basic Metropolis", ]
  validity_summary <- summarize_adaptive_validity(comparison_result$comparison)

  expect_true(any(comparison_result$comparison$runs$adaptive$adaptation_magnitude > 0))
  expect_true(all(is.finite(adaptive_summary$effective_sample_size)))
  expect_gt(
    mean(adaptive_summary$effective_sample_size),
    mean(basic_summary$effective_sample_size)
  )
  expect_equal(nrow(validity_summary), 3L)
  expect_true(all(validity_summary$min_proposal_eigenvalue > 0))
})
