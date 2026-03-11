// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using Rcpp::Function;
using Rcpp::List;
using Rcpp::Named;
using arma::accu;
using arma::cube;
using arma::ivec;
using arma::mat;
using arma::mvnrnd;
using arma::norm;
using arma::uword;
using arma::vec;
namespace fill = arma::fill;

// Random-Walk Metropolis sampler with a fixed Gaussian proposal.
//
// Proposal:
//   theta* ~ N(theta_n, S)
//
// Because the proposal is symmetric, the Metropolis-Hastings acceptance
// probability simplifies to
//   alpha = min(1, pi(theta*) / pi(theta_n)).
//
// [[Rcpp::export]]
List basic_metropolis_cpp(Function target_density,
                          const vec& initial_state,
                          const int n_iter,
                          const mat& proposal_covariance) {
  const uword dimension = initial_state.n_elem;
  
  // Storage for the Markov chain and diagnostics.
  mat chain_draws(n_iter, dimension, fill::zeros);
  ivec accepted_move(n_iter, fill::zeros);
  vec acceptance_probability(n_iter, fill::zeros);
  vec target_density_values(n_iter, fill::zeros);
  cube proposal_covariance_history(dimension, dimension, n_iter, fill::zeros);
  vec adaptation_magnitude(n_iter, fill::zeros);
  vec jump_distance(n_iter, fill::zeros);

  // Start the chain at the user-supplied initial state.
  vec theta_n = initial_state;

  // Evaluate the target density at the initial state.
  // This code assumes target_density returns the density itself,
  // not the log-density.
  double pi_theta_n = Rcpp::as<double>(target_density(Rcpp::wrap(theta_n)));

  for (int iteration = 0; iteration < n_iter; ++iteration) {
    // Save the current state so we can measure the realised jump size.
    const vec theta_previous = theta_n;
    
    // Step 1: propose a new state using a Gaussian random walk- from armadillo
    const vec theta_star = mvnrnd(theta_n, proposal_covariance);
  
    // Evaluate the target density at the proposal.
    const double pi_theta_star =
      Rcpp::as<double>(target_density(Rcpp::wrap(theta_star)));
  
    // Step 2: compute the Metropolis acceptance probability.
    // If either density is invalid, force rejection.
    double alpha_n = 0.0;
    if (R_finite(pi_theta_star) && pi_theta_star > 0.0 &&
        R_finite(pi_theta_n) && pi_theta_n > 0.0) {
      alpha_n = std::min(1.0, pi_theta_star / pi_theta_n);
    }
  
    acceptance_probability(iteration) = alpha_n;
  
    // Step 3: accept or reject the proposal.
    const int accepted = (R::runif(0.0, 1.0) < alpha_n);
    accepted_move(iteration) = accepted;
  
    if (accepted) {
      theta_n = theta_star;
      pi_theta_n = pi_theta_star;
    }
  
    // Step 4: store the post-decision state and diagnostics.
    chain_draws.row(iteration) = theta_n.t();
    target_density_values(iteration) = pi_theta_n;
  
    // The proposal covariance is fixed for Basic Metropolis,
    // so the same matrix is stored at every iteration.
    proposal_covariance_history.slice(iteration) = proposal_covariance;
  
    // No adaptation occurs in this algorithm, so this stays zero.
    adaptation_magnitude(iteration) = 0.0;
  
    // Realised Euclidean jump distance after accept/reject.
    jump_distance(iteration) = norm(theta_n - theta_previous, 2);
  }

  return List::create(
    Named("algorithm") = "Basic Metropolis",
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