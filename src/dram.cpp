// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using Rcpp::Function;
using Rcpp::List;
using Rcpp::Named;
using arma::accu;
using arma::cube;
using arma::dot;
using arma::eye;
using arma::ivec;
using arma::mat;
using arma::mean;
using arma::mvnrnd;
using arma::norm;
using arma::rowvec;
using arma::solve;
using arma::uword;
using arma::vec;
namespace fill = arma::fill;

// Delayed Rejection Adaptive Metropolis (DRAM)
//
// DRAM combines two ideas:
//
// 1. Adaptive Metropolis:
//    update the proposal covariance using the realised chain history.
//
// 2. Delayed Rejection:
//    if the first proposal is rejected, try a second smaller proposal
//    instead of staying put immediately.
//
// [[Rcpp::export]]
List dram_cpp(Function target_density,
              const vec& initial_state,
              const int n_iter,
              const mat& initial_proposal_covariance,
              const int adapt_start = 50,
              const double delayed_rejection_scale = 0.5,
              const double epsilon = 1e-8) {
  const uword dimension = initial_state.n_elem;
  const mat I_d = eye<mat>(dimension, dimension);
  
  // Storage for the chain and diagnostics.
  mat chain_draws(n_iter, dimension, fill::zeros);
  ivec accepted_move(n_iter, fill::zeros);
  ivec stage_one_accepted(n_iter, fill::zeros);
  ivec stage_two_accepted(n_iter, fill::zeros);
  ivec second_stage_attempted(n_iter, fill::zeros);
  vec stage_one_acceptance_probability(n_iter, fill::zeros);
  vec stage_two_acceptance_probability(n_iter, fill::zeros);
  vec target_density_values(n_iter, fill::zeros);
  cube proposal_covariance_history(dimension, dimension, n_iter, fill::zeros);
  vec adaptation_magnitude(n_iter, fill::zeros);
  vec jump_distance(n_iter, fill::zeros);
  
  // Start the chain at the initial state X_n.
  vec X_n = initial_state;
  double pi_X_n = Rcpp::as<double>(target_density(Rcpp::wrap(X_n)));
  
  // Initial proposal covariance for the first-stage Gaussian proposal.
  mat S_n = initial_proposal_covariance + epsilon * I_d;
  
  for (int iteration = 0; iteration < n_iter; ++iteration) {
    const vec X_previous = X_n;
    
    // Store the covariance used at this iteration.
    proposal_covariance_history.slice(iteration) = S_n;
    
    // ------------------------------------------------------------
    // Step 1: First-stage proposal
    //
    // Propose:
    //   Y_1 ~ N(X_n, S_n)
    // ------------------------------------------------------------
    const vec Y_1 = mvnrnd(X_n, S_n);
    const double pi_Y_1 = Rcpp::as<double>(target_density(Rcpp::wrap(Y_1)));
    
    // First-stage Metropolis acceptance probability:
    //   alpha_1 = min(1, pi(Y_1) / pi(X_n))
    double alpha_1 = 0.0;
    if (R_finite(pi_Y_1) && pi_Y_1 > 0.0 && pi_X_n > 0.0) {
      alpha_1 = std::min(1.0, pi_Y_1 / pi_X_n);
    }
    
    
    stage_one_acceptance_probability(iteration) = alpha_1;
    
    const int accepted_stage_one = (R::runif(0.0, 1.0) < alpha_1);
    stage_one_accepted(iteration) = accepted_stage_one;
  
    if (accepted_stage_one) {
      // If the first proposal is accepted, move to Y_1.
      X_n = Y_1;
      pi_X_n = pi_Y_1;
      accepted_move(iteration) = 1;
    } else {
      
      // ----------------------------------------------------------
      // Step 2: Second-stage proposal
      //
      // If Y_1 is rejected, try a smaller rescue proposal:
      //   Y_2 ~ N(X_n, gamma S_n)
      //
      // where gamma = delayed_rejection_scale and gamma < 1.
      // ----------------------------------------------------------
      second_stage_attempted(iteration) = 1;
      
      
      const mat stage_two_covariance = delayed_rejection_scale * S_n;
      const vec Y_2 = mvnrnd(X_n, stage_two_covariance);
      const double pi_Y_2 = Rcpp::as<double>(target_density(Rcpp::wrap(Y_2)));
      
      
      // ----------------------------------------------------------
      // Step 3: Second-stage acceptance probability
      //
      // This is NOT the ordinary Metropolis ratio.
      //
      // Why?
      // Because Y_2 is only proposed after Y_1 has already been rejected.
      // So we need a correction term to preserve the correct stationary
      // distribution.
      //
      // The delayed rejection formula uses:
      // 1. a reverse first-stage acceptance probability
      // 2. a first-stage proposal density ratio
      // ----------------------------------------------------------
    
      // Reverse first-stage acceptance probability:
      // this is alpha_1(Y_2, Y_1), meaning:
      // if the chain were at Y_2, how likely would Y_1 be accepted
      // as a first-stage proposal?
      double reverse_alpha_1 = 0.0;
      if (R_finite(pi_Y_1) && pi_Y_1 > 0.0 && pi_Y_2 > 0.0) {
        reverse_alpha_1 = std::min(1.0, pi_Y_1 / pi_Y_2);
      }
      
      // First-stage proposal density ratio:
      //
      //   q_1(Y_1 | Y_2) / q_1(Y_1 | X_n)
      //
      // This compares how likely the first proposal Y_1 would be:
      // - if the chain started at Y_2
      // - versus if the chain started at X_n
      //
      // Because q_1 is Gaussian in both cases, the normalising constants
      // cancel, so only the quadratic terms remain.
      const double q1_ratio = std::exp(
        -0.5 * dot(Y_1 - Y_2, solve(S_n, Y_1 - Y_2)) +
          0.5 * dot(Y_1 - X_n, solve(S_n, Y_1 - X_n))
      );
      
      // Delayed rejection acceptance probability:
      //
      //                pi(Y_2) q_1(Y_1 | Y_2) [1 - alpha_1(Y_2, Y_1)]
      // alpha_2 = ---------------------------------------------------------
      //                 pi(X_n) q_1(Y_1 | X_n) [1 - alpha_1(X_n, Y_1)]
      //
      // The q_1 ratio above is exactly:
      //   q_1(Y_1 | Y_2) / q_1(Y_1 | X_n)
      //
      // So we compute alpha_2 using that corrected ratio.
      double alpha_2 = 0.0;
      const double alpha_2_numerator =
        pi_Y_2 * q1_ratio * (1.0 - reverse_alpha_1);
      const double alpha_2_denominator =
        pi_X_n * (1.0 - alpha_1);
      
      
      if (R_finite(alpha_2_numerator) && alpha_2_numerator > 0.0 &&
          alpha_2_denominator > 0.0) {
        alpha_2 = std::min(1.0, alpha_2_numerator / alpha_2_denominator);
      }
      
      
      stage_two_acceptance_probability(iteration) = alpha_2;
      
      
      const int accepted_stage_two = (R::runif(0.0, 1.0) < alpha_2);
      stage_two_accepted(iteration) = accepted_stage_two;
      
      
      if (accepted_stage_two) {
        // If the second-stage proposal is accepted, move to Y_2.
        X_n = Y_2;
        pi_X_n = pi_Y_2;
        accepted_move(iteration) = 1;
      }
      
      // If Y_2 is also rejected, the chain stays at X_n.
    }
    
    // Store the realised chain value after delayed rejection is complete.
    chain_draws.row(iteration) = X_n.t();
    target_density_values(iteration) = pi_X_n;
    jump_distance(iteration) = norm(X_n - X_previous, 2);
  
    // ------------------------------------------------------------
    // Step 4: Adaptive Metropolis covariance update
    //
    // After enough iterations, update S_n using the empirical
    // covariance of the realised chain history.
    // ------------------------------------------------------------
    if ((iteration + 1) >= adapt_start && iteration >= 1) {
      const mat S_previous = S_n;
      const mat theta_history = chain_draws.rows(0, iteration);
      const rowvec theta_bar = mean(theta_history, 0);
      const mat centered_history = theta_history.each_row() - theta_bar;
      
      S_n =
        (centered_history.t() * centered_history) /
          static_cast<double>(iteration + 1) +
            epsilon * I_d;
      
      adaptation_magnitude(iteration) = norm(S_n - S_previous, "fro");
    }
  }
  
  return List::create(
    Named("algorithm") = "Delayed Rejection Adaptive Metropolis",
    Named("draws") = chain_draws,
    Named("accepted") = accepted_move,
    Named("acceptance_rate") = accu(accepted_move) / static_cast<double>(n_iter),
    Named("stage_one_accepted") = stage_one_accepted,
    Named("stage_two_accepted") = stage_two_accepted,
    Named("second_stage_attempted") = second_stage_attempted,
    Named("stage_one_acceptance_probability") = stage_one_acceptance_probability,
    Named("stage_two_acceptance_probability") = stage_two_acceptance_probability,
    Named("target_values") = target_density_values,
    Named("proposal_covariance_history") = proposal_covariance_history,
    Named("adaptation_magnitude") = adaptation_magnitude,
    Named("jump_distance") = jump_distance
  );
}

