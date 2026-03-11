// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using Rcpp::Function;
using Rcpp::List;
using Rcpp::Named;
using arma::accu;
using arma::chol;
using arma::cube;
using arma::dot;
using arma::eye;
using arma::ivec;
using arma::mat;
using arma::norm;
using arma::randn;
using arma::uword;
using arma::vec;
namespace fill = arma::fill;

// Robust Adaptive Metropolis (RAM)
//
// Proposal:
//   Y_n = X_{n-1} + S_{n-1} U_n
//   where U_n ~ N(0, I)
//
// Adaptation:
//   S_n S_n^T
//   = S_{n-1} [ I + eta_n (alpha_n - alpha^*) (U_n U_n^T / ||U_n||^2) ] S_{n-1}^T
//
// [[Rcpp::export]]
List robust_adaptive_metropolis_cpp(Function target_density,
                                    const vec& initial_state,
                                    const int n_iter,
                                    const mat& initial_proposal_covariance,
                                    const double target_acceptance = 0.234,
                                    const double adapt_exponent = 0.67) {
  const uword dimension = initial_state.n_elem;
  const mat I_d = eye<mat>(dimension, dimension);
  
  mat chain_draws(n_iter, dimension, fill::zeros);
  ivec accepted_move(n_iter, fill::zeros);
  vec acceptance_probability(n_iter, fill::zeros);
  vec target_density_values(n_iter, fill::zeros);
  cube proposal_covariance_history(dimension, dimension, n_iter, fill::zeros);
  vec adaptation_magnitude(n_iter, fill::zeros);
  vec jump_distance(n_iter, fill::zeros);
  
  
  // Start the chain at X_0.
  vec X_n = initial_state;
  double pi_X_n = Rcpp::as<double>(target_density(Rcpp::wrap(X_n)));
  
  
  // S_n is the proposal factor and Sigma_n = S_n S_n^T is the proposal covariance.
  mat S_n;
  if (!chol(S_n, initial_proposal_covariance, "lower")) {
    Rcpp::stop("initial_proposal_covariance must be positive definite.");
  }
  
  mat Sigma_n = initial_proposal_covariance;

  for (int iteration = 0; iteration < n_iter; ++iteration) {
    const vec X_n_minus_1 = X_n;
    const double pi_X_n_minus_1 = pi_X_n;
    const mat Sigma_n_minus_1 = Sigma_n;
    const mat S_n_minus_1 = S_n;
  
    // Store the covariance used at this iteration.
    proposal_covariance_history.slice(iteration) = Sigma_n_minus_1;
  
    // Step 1: Proposal
    // Draw U_n ~ N(0, I) and propose Y_n = X_{n-1} + S_{n-1} U_n.
    const vec U_n = randn<vec>(dimension);
    const vec Y_n = X_n_minus_1 + S_n_minus_1 * U_n;
    const double pi_Y_n = Rcpp::as<double>(target_density(Rcpp::wrap(Y_n)));
  
    // Step 2: Accept/Reject
    double alpha_n = 0.0;
    if (R_finite(pi_Y_n) && pi_Y_n > 0.0 && pi_X_n_minus_1 > 0.0) {
      alpha_n = std::min(1.0, pi_Y_n / pi_X_n_minus_1);
    }
  
    acceptance_probability(iteration) = alpha_n;
  
    const int accepted = (R::runif(0.0, 1.0) < alpha_n);
    accepted_move(iteration) = accepted;
  
    if (accepted) {
      X_n = Y_n;
      pi_X_n = pi_Y_n;
    } else {
      X_n = X_n_minus_1;
      pi_X_n = pi_X_n_minus_1;
    }
  
    // Store the realised state.
    chain_draws.row(iteration) = X_n.t();
    target_density_values(iteration) = pi_X_n;
    jump_distance(iteration) = norm(X_n - X_n_minus_1, 2);
  
    // Step 3: Adapt Proposal
    const double eta_n =
      1.0 / std::pow(static_cast<double>(iteration + 1), adapt_exponent);
  
    const double U_n_norm_squared = dot(U_n, U_n);
    const mat U_n_direction = (U_n * U_n.t()) / U_n_norm_squared;
    
    
    const mat update_matrix =
      I_d + eta_n * (alpha_n - target_acceptance) * U_n_direction;
    
    
    // Update the proposal covariance exactly as written on the slide:
    // Sigma_n = S_n S_n^T
    Sigma_n = S_n_minus_1 * update_matrix * S_n_minus_1.t();
    
    
    // Recover a factor S_n for the next proposal step.
    if (!chol(S_n, Sigma_n, "lower")) {
      Rcpp::stop("RAM produced a non-positive proposal covariance.");
    }
    
    
    adaptation_magnitude(iteration) = norm(Sigma_n - Sigma_n_minus_1, "fro");
  }
  
  
  return List::create(
    Named("algorithm") = "Robust Adaptive Metropolis",
    Named("draws") = chain_draws,
    Named("accepted") = accepted_move,
    Named("acceptance_probability") = acceptance_probability,
    Named("acceptance_rate") = accu(accepted_move) / static_cast<double>(n_iter),
    Named("target_values") = target_density_values,
    Named("proposal_covariance_history") = proposal_covariance_history,
    Named("adaptation_magnitude") = adaptation_magnitude,
    Named("jump_distance") = jump_distance
  );
}
