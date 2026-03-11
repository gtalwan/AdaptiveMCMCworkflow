test_that("basic and adaptive samplers return the expected shapes", {
  set.seed(1)
  target <- make_correlated_gaussian_target()

  basic_run <- basic_metropolis(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 80,
    proposal_covariance = target$recommended_proposal_covariance
  )

  adaptive_run <- adaptive_metropolis(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 80,
    initial_proposal_covariance = target$recommended_proposal_covariance,
    adapt_start = 20
  )

  expect_s3_class(basic_run, "adaptive_mcmc_run")
  expect_equal(dim(basic_run$draws), c(80, 2))
  expect_true(is.logical(basic_run$accepted))
  expect_true(all(basic_run$acceptance_probability >= 0 & basic_run$acceptance_probability <= 1))
  expect_equal(dim(basic_run$proposal_covariance_history), c(2, 2, 80))

  expect_equal(dim(adaptive_run$draws), c(80, 2))
  expect_true(any(adaptive_run$adaptation_magnitude > 0))
  expect_equal(dim(adaptive_run$proposal_covariance_history), c(2, 2, 80))
})

test_that("RAM keeps a positive proposal covariance history", {
  set.seed(1)
  target <- make_small_scale_target()

  ram_run <- robust_adaptive_metropolis(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 100,
    initial_proposal_covariance = target$recommended_proposal_covariance
  )

  smallest_eigenvalue <- vapply(
    seq_len(dim(ram_run$proposal_covariance_history)[3]),
    function(index) {
      min(eigen(ram_run$proposal_covariance_history[, , index, drop = TRUE],
        symmetric = TRUE,
        only.values = TRUE
      )$values)
    },
    numeric(1)
  )

  expect_true(all(smallest_eigenvalue > 0))
  expect_true(any(ram_run$adaptation_magnitude > 0))
})

test_that("DRAM stage-two bookkeeping is internally consistent", {
  set.seed(1)
  target <- make_dram_rescue_target()

  dram_run <- dram(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 120,
    initial_proposal_covariance = target$recommended_proposal_covariance,
    adapt_start = target$recommended_adapt_start,
    delayed_rejection_scale = target$recommended_delayed_rejection_scale
  )

  expect_true(any(dram_run$second_stage_attempted))
  expect_identical(dram_run$second_stage_attempted, !dram_run$stage_one_accepted)
  expect_identical(dram_run$accepted, dram_run$stage_one_accepted | dram_run$stage_two_accepted)
  expect_true(all(dram_run$stage_two_acceptance_probability >= 0 & dram_run$stage_two_acceptance_probability <= 1))
})

test_that("run summaries use post-burn-in rates and conditional DRAM rescue rates", {
  set.seed(1)
  target <- make_dram_rescue_target()

  dram_run <- dram(
    log_target = target$log_density,
    initial_state = target$recommended_initial,
    n_iter = 120,
    initial_proposal_covariance = target$recommended_proposal_covariance,
    adapt_start = target$recommended_adapt_start,
    delayed_rejection_scale = target$recommended_delayed_rejection_scale
  )

  burn_in <- 20L
  kept <- (burn_in + 1L):dram_run$n_iter
  attempted <- dram_run$second_stage_attempted[kept]
  summary_table <- summarize_sampler_run(dram_run, burn_in = burn_in)

  expect_equal(summary_table$acceptance_rate[1], mean(dram_run$accepted[kept]))
  expect_equal(summary_table$stage_one_accept_rate[1], mean(dram_run$stage_one_accepted[kept]))
  expect_equal(summary_table$second_stage_attempt_rate[1], mean(attempted))
  expect_equal(summary_table$second_stage_move_rate[1], mean(dram_run$stage_two_accepted[kept]))
  expect_equal(
    summary_table$second_stage_accept_rate[1],
    mean(dram_run$stage_two_accepted[kept][attempted])
  )
})

test_that("compiled sampler entry points run on a user-supplied density", {
  target_density <- function(theta) {
    theta <- as.numeric(theta)
    exp(-0.5 * sum(theta^2))
  }

  set.seed(1)
  basic_cpp_run <- AdaptiveMCMCWorkflow:::basic_metropolis_cpp(
    target_density = target_density,
    initial_state = c(0, 0),
    n_iter = 40L,
    proposal_covariance = diag(2)
  )

  set.seed(1)
  adaptive_cpp_run <- AdaptiveMCMCWorkflow:::adaptive_metropolis_cpp(
    target_density = target_density,
    initial_state = c(0, 0),
    n_iter = 40L,
    initial_proposal_covariance = diag(2),
    adapt_start = 20L,
    epsilon = 1e-8
  )

  set.seed(1)
  ram_cpp_run <- AdaptiveMCMCWorkflow:::robust_adaptive_metropolis_cpp(
    target_density = target_density,
    initial_state = c(0, 0),
    n_iter = 40L,
    initial_proposal_covariance = diag(2),
    target_acceptance = 0.234,
    adapt_exponent = 0.67
  )

  set.seed(1)
  dram_cpp_run <- AdaptiveMCMCWorkflow:::dram_cpp(
    target_density = target_density,
    initial_state = c(0, 0),
    n_iter = 40L,
    initial_proposal_covariance = diag(2),
    adapt_start = 20L,
    delayed_rejection_scale = 0.5,
    epsilon = 1e-8
  )

  expect_equal(dim(basic_cpp_run$draws), c(40, 2))
  expect_equal(dim(adaptive_cpp_run$draws), c(40, 2))
  expect_equal(dim(ram_cpp_run$draws), c(40, 2))
  expect_equal(dim(dram_cpp_run$draws), c(40, 2))
  expect_true(all(basic_cpp_run$acceptance_probability >= 0 & basic_cpp_run$acceptance_probability <= 1))
  expect_true(all(adaptive_cpp_run$acceptance_probability >= 0 & adaptive_cpp_run$acceptance_probability <= 1))
  expect_true(all(ram_cpp_run$acceptance_probability >= 0 & ram_cpp_run$acceptance_probability <= 1))
  expect_true(all(dram_cpp_run$stage_one_acceptance_probability >= 0 & dram_cpp_run$stage_one_acceptance_probability <= 1))
  expect_true(any(adaptive_cpp_run$adaptation_magnitude > 0))
  expect_true(any(ram_cpp_run$adaptation_magnitude > 0))
})
