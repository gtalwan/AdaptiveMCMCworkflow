// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using Rcpp::Function;
using Rcpp::List;
using Rcpp::Named;
using arma::accu;
using arma::cube;
using arma::eye;
using arma::ivec;
using arma::mat;
using arma::mean;
using arma::mvnrnd;
using arma::norm;
using arma::rowvec;
using arma::uword;
using arma::vec;
namespace fill = arma::fill;

// Adaptive Metropolis sampler.
//
// Proposal at iteration n:
//   theta* ~ N(theta_n, S_n)
//
// The proposal covariance S_n is updated from the empirical covariance
// of the chain history, plus a small ridge term epsilon * I to stabilize
// the matrix.
//
// [[Rcpp::export]]
List adaptive_metropolis_cpp(Function target_density,
                             const vec& initial_state,
                             const int n_iter,
                             const mat& initial_proposal_covariance,
                             const int adapt_start = 50,
                             const double epsilon = 1e-8) {
  const uword dimension = initial_state.n_elem;
  const mat I_d = eye<mat>(dimension, dimension);
  
  // Storage for draws and diagnostics.
  mat chain_draws(n_iter, dimension, fill::zeros);
  ivec accepted_move(n_iter, fill::zeros);
  vec acceptance_probability(n_iter, fill::zeros);
  vec target_density_values(n_iter, fill::zeros);
  cube proposal_covariance_history(dimension, dimension, n_iter, fill::zeros);
  vec adaptation_magnitude(n_iter, fill::zeros);
  vec jump_distance(n_iter, fill::zeros);

  // Start the chain at the user-supplied initial value.
  vec theta_n = initial_state;

  // Evaluate the target density at the starting point.
  // This assumes target_density returns the density itself, not log-density.
  double pi_theta_n = Rcpp::as<double>(target_density(Rcpp::wrap(theta_n)));

  // Initial proposal covariance.
  mat S_n = initial_proposal_covariance + epsilon * I_d;

  for (int iteration = 0; iteration < n_iter; ++iteration) {
    // Save the current state so we can measure the realised jump size.
    const vec theta_previous = theta_n;
  
    // Store the proposal covariance used at this iteration.
    proposal_covariance_history.slice(iteration) = S_n;
  
    // Step 1: propose a new state from the current Gaussian random walk.
    const vec theta_star = mvnrnd(theta_n, S_n);
  
    // Evaluate the target density at the proposal.
    const double pi_theta_star =
      Rcpp::as<double>(target_density(Rcpp::wrap(theta_star)));
  
    // Step 2: Metropolis acceptance probability.
    double alpha_n = 0.0;
    if (R_finite(pi_theta_star) && pi_theta_star > 0.0 && pi_theta_n > 0.0) {
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
  
    // Store the realised state of the Markov chain.
    chain_draws.row(iteration) = theta_n.t();
    target_density_values(iteration) = pi_theta_n;
    jump_distance(iteration) = norm(theta_n - theta_previous, 2);
  
    // Step 4: adapt the proposal covariance after enough iterations - process called burn in
    if ((iteration + 1) >= adapt_start && iteration >= 1) {
      const mat S_previous = S_n;
    
      // Use the realised chain history, including repeated states from rejections.
      const mat theta_history = chain_draws.rows(0, iteration);
    
      // Empirical mean of the chain history.
      const rowvec theta_bar = mean(theta_history, 0);
    
      // Center each draw by subtracting the empirical mean.
      const mat centered_history = theta_history.each_row() - theta_bar;
    
      // Empirical covariance plus epsilon * I for numerical stability.
      S_n =
        (centered_history.t() * centered_history) /
          static_cast<double>(iteration + 1) +
            epsilon * I_d;
    
      // Frobenius norm of the covariance change measures adaptation size.
      adaptation_magnitude(iteration) = norm(S_n - S_previous, "fro");
    }
  }

  return List::create(
    Named("algorithm") = "Adaptive Metropolis",
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