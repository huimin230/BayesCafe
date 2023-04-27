#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

// static arma::uvec which(IntegerVector x, int c);
// static int num(IntegerVector x, int c);
// static arma::uvec organize_label(IntegerVector x);
// static arma::uvec remove(NumericVector x, int c);

// [[Rcpp::export]]
Rcpp::List mcmc(arma::mat Y, int K, NumericVector s, arma::mat P, double f, int iter, int burn) {
  // Read data information
  int n = Y.n_rows;
  int p = Y.n_cols;

  // Set the hyperparameters
  double pi = 0.5;
  double a_omega = 0.1;
  double b_omega = 1.9;
  double omega = 0.05;
  double a_phi = 0.001;
  double b_phi = 0.001;
  double a_0 = 0.001;
  double b_0 = 0.001;
  double a_mu = 0.001;
  double b_mu = 0.001;
  
  // Markov random field (MRF) prior setting 
  IntegerVector e(K);
  for (int k = 0; k < K; k++){
    e(k) = 1;
  }
  int n_neighbor = P.n_cols;
  int min_element = 10;
  
  // Set tuning parameters and algorithm settings
  double tau_log_mu = 0.1;
  double tau_log_mu0 = 0.1;
  double tau_log_phi = 1;  
  double tun_mu;
  bool findFeature = true;

  int W = p*0.2;

  arma::mat Mu(n, p);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < p; j++){
      Mu(i, j) = Y(i, j)/s(i);
    }
  }
  
  // Create the spaces to store the results 
  IntegerMatrix H(n, p);
  NumericMatrix H_ppi(n, p);
  IntegerVector H_sum(iter);
  IntegerMatrix z_store(iter, n);
  NumericVector map_z_store(iter);
  IntegerMatrix gamma_store(iter, p);
  NumericVector gamma_ppi(p);
  NumericVector gamma_BF(p);
  NumericMatrix logBF(iter, p);
  IntegerVector gamma_sum(iter);
  NumericVector map_gamma_store(iter);
  NumericMatrix phi_store(iter, p);
  NumericVector map_store(iter);   // to store full posterior 
  arma::cube u_store(K, p, iter);
  
  // Set temporary variables
  int t, i, j, w, k, gamma_new, count = 10, H_temp, gamma_sum_temp = 0, neighbor_index, n_neighbor_temp;
  double c_pi = 3.1415926, map_z_temp, map_gamma_temp;
  double hastings, phi_new, mu_new;
  double try_mu = 0.0, try_phi = 0.0, try_gamma = 0.0;
  double accept_mu = 0.0, accept_phi = 0.0, accept_gamma = 0.0;
  IntegerVector gamma_map, z_map, H_sum_temp(n, 0);
  NumericVector prob_temp(2);
  NumericVector prob(K);

  // Initialization
  NumericVector clusters(K);
  NumericVector prob_init(K);
  for (int k = 0; k < K; k++){
    clusters(k) = k;
    prob_init(k) = pow(K, -1);
  }
  IntegerVector z(n);
  for (int i = 0; i < n; i++){
    z(i) = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(clusters, 1, TRUE, prob_init));
  }
  NumericMatrix u(K, p);
  for (k = 0; k < K; k++){
    for (j = 0; j < p; j++){
      u(k, j) = 1;
      
    }
  }

  IntegerVector gamma(p);
  NumericVector phi(p);
  double phi_start = 10;
  for (j = 0; j < p; j++) { 
    gamma(j) = 0;
    // gamma(j) = rbinom(1, 1, 0.1)(0);
    phi(j) = phi_start;
    gamma_sum_temp = gamma_sum_temp + gamma(j);
    
    for (i = 0; i < n; i++) {
      H_ppi(i, j) = 0;
      if (Y(i, j) != 0) {
        H(i, j) = 0;
      } else {
        H(i, j) = 0;
        H_sum_temp(i) = H_sum_temp(i) + H(i, j);
      }
      }
    }

  
  // Start MCMC algorithms =====================================================
  for (t = 0; t < iter; t++) {
    // Gibbs Sampling for updating cluster allocation indicator z ==============
    for (i = 0; i < n; i++) {
      // Compute I(z_i = z_i')
      IntegerVector neighbor_temp(K);

      for (int ii = 0; ii < n_neighbor; ii++){
        if (P(i, ii) != 0) {
          neighbor_index = P(i, ii) - 1;
          neighbor_temp(z(neighbor_index)) = neighbor_temp(z(neighbor_index)) + 1;
          }
        }

      // Calculate posterior probability of z
      NumericVector loglklh(K);

      for (k = 0; k < K; k++) {
        for (j = 0; j < p; j++) {
          if (gamma(j) == 1 && H(i, j) == 0) {
            loglklh(k) = loglklh(k) + phi(j)*(-log(u(k, j)*s(i) + phi(j))) + Y(i, j)*(log(u(k, j)) - log(u(k, j)*s(i) + phi(j)));

          }
        }

        loglklh(k) = loglklh(k) + e(k) + f*(neighbor_temp(k));

      }

      for (k = 0; k < K; k++) {
        double diff = 0;
        for (int m = 0; m < K; m++) {
          if (m != k) {
            diff = diff + exp(loglklh(m) - loglklh(k));
          }
        }
        prob(k) = 1/(1 + diff);
      }

      // Generate new sample
      int z_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(clusters, 1, TRUE, prob));

      // Make sure each group has at least multiple observations
      k = z(i);

      if (sum(z == k) == min_element && z_new != k) {
        z(i) = k;
      }else {
        z(i) = z_new;
      }
    } // End of updating z

    
     // Gibbs Sampler for updating zero-inflation indicator H ==================
     for (i = 0; i < n; i++) {
       for (j = 0; j < p; j++) {

         if (Y(i, j) == 0) {
           int index = z(i);
           prob_temp(0) = phi(j)*(log(phi(j)) - log(s(i)*u(index, j) + phi(j))) + log(1 - pi);
           prob_temp(1) = log(pi);
           prob_temp(1) = 1/(1 + exp(prob_temp(0) - prob_temp(1)));

          H_temp = H(i, j);
          H(i, j) = rbinom(1, 1, prob_temp(1))(0);
          H_sum_temp(i) = H_sum_temp(i) - H_temp + H(i, j);

         }
       }
     }  // End of updating H


    // Random walk Metropolis-Hastings (RWMH) algorithm for updating dispersion phi
    for (j = 0; j < p; j++) {
      phi_new = exp(r_truncnorm(log(phi(j)), tau_log_phi, log(1), log(100)));

      hastings = 0;

      for (i = 0; i < n; i++) {
        if (H(i, j) == 0) {
          int index = z(i);
          hastings = hastings + lgamma(Y(i, j) + phi_new) - lgamma(phi_new) +
            phi_new*(log(phi_new) - log(u(index, j)*s(i) + phi_new)) + Y(i, j)*(- log(u(index, j)*s(i) + phi_new));
          hastings = hastings - (lgamma(Y(i, j) + phi(j)) - lgamma(phi(j)) +
            phi(j)*(log(phi(j)) - log(u(index, j)*s(i) + phi(j))) + Y(i, j)*(- log(u(index, j)*s(i) + phi(j))));
        }

      }

      hastings = hastings + (a_phi - 1)*log(phi_new) - b_phi*phi_new;
      hastings = hastings - ((a_phi - 1)*log(phi(j)) - b_phi*phi(j));

      if (t >= burn) {
        try_phi++;
      }

      if(hastings >= log(double(rand()%10001)/10000)) {
        phi(j) = phi_new;

        if (t >= burn) {
          accept_phi++;
        }
      }
    } // End of updating phi


    // Add-delete algorithm for updating feature selection indicator gamma =====
    if (findFeature) {
      for (w = 0; w < W; w++){
        j = rand()%p;
        gamma_new = 1 - gamma(j);
        tun_mu = log(mean(Mu.col(j)));
        
        if (gamma_new == 0) { // Delete step
          hastings = 0;
          mu_new = exp(rnorm(1, tun_mu, tau_log_mu0)(0));
          
          // Log-likelihood by feature
          for (i = 0; i < n; i++) {
            if (H(i, j) == 0){
              int index = z(i);
              hastings = hastings + phi(j)*(- log(mu_new*s(i) + phi(j))) + Y(i, j)*(log(mu_new*s(i)) - log(mu_new*s(i) + phi(j)));
              hastings = hastings - (phi(j)*(- log(u(index, j)*s(i) + phi(j))) + Y(i, j)*(log(u(index, j)*s(i)) - log(u(index, j)*s(i) + phi(j))));
            }
          }
          
          // Proposal distribution
          for (k = 0; k < K; k++){
            double mu = u(k, j);
            hastings = hastings + (-log(2*c_pi)/2 - log(10*tau_log_mu) - pow((log(mu) - tun_mu), 2)/2/pow(10*tau_log_mu, 2));
          }
          hastings = hastings - (-log(2*c_pi)/2 - log(tau_log_mu0) - pow((log(mu_new) - tun_mu),2)/2/pow(tau_log_mu0, 2));
          
          // Conditional prior distribution
          hastings = hastings + a_0*log(b_0) - lgamma(a_0) + (a_0 - 1)*log(mu_new) - b_0*mu_new;
          
          for (k = 0; k < K; k ++){
            double mu = u(k, j);
            hastings = hastings - (a_mu*log(b_mu) - lgamma(a_mu) + (a_mu - 1)*log(mu) - b_mu*mu);
          }
          
          // Prior
          hastings = hastings - log(omega) + log(1 - omega);
          
          if (t >= burn) {
            try_gamma++;
          }
          
          if(hastings >= log(double(rand()%10001)/10000)) {
            gamma(j) = gamma_new;
            gamma_sum_temp--;
            
            for (k = 0; k < K; k++) {
              u(k, j) = mu_new;
            }
            
            if (t >= burn) {
              accept_gamma++;
            }
          }
          
        } // End of delete step
        
        else { //Add step
          hastings = 0;
          double mu = u(0, j);
          NumericVector mu_new_vector(K);
          
          for (int k = 0; k < K; k++) {
            mu_new_vector(k) = exp(rnorm(1, tun_mu, 10*tau_log_mu0)(0));
          }
          
          // Log-likelihood by feature
          for (i = 0; i < n; i++) {
            if (H(i, j) == 0){
              int index = z(i);
              hastings = hastings + phi(j)*(- log(mu_new_vector(index)*s(i) + phi(j))) + Y(i, j)*(log(mu_new_vector(index)*s(i)) - log(mu_new_vector(index)*s(i) + phi(j)));
              hastings = hastings - (phi(j)*(- log(mu*s(i) + phi(j))) + Y(i, j)*(log(mu*s(i)) - log(mu*s(i) + phi(j))));
            }
          }
          
          // Proposal distribution
          hastings = hastings + (-log(2*c_pi)/2 - log(tau_log_mu0) - pow((log(mu) - tun_mu), 2)/2/pow(tau_log_mu0, 2));
          for (k = 0; k < K; k++){
            hastings = hastings - (-log(2*c_pi)/2 - log(10*tau_log_mu) - pow((log(mu_new_vector(k)) - tun_mu), 2)/2/pow(10*tau_log_mu, 2));
          }
          
          
          // Conditional prior distribution
          for (k = 0; k < K; k++) {
            hastings = hastings + a_mu*log(b_mu) - lgamma(a_mu) + (a_mu - 1)*log(mu_new_vector(k)) - b_mu*mu_new_vector(k);
          }
          
          hastings = hastings - (a_0*log(b_0) - lgamma(a_0) + (a_0 - 1)*log(mu) - b_0*mu);
          
          // Prior
          hastings = hastings + log(omega) - log(1 - omega);
          
          if (t >= burn) {
            try_gamma++;
          }
          
          if(hastings >= log(double(rand()%10001)/10000)) {
            gamma(j) = gamma_new;
            gamma_sum_temp++;
            
            for (k = 0; k < K; k++) {
              u(k, j) = mu_new_vector(k);
            }
            
            if (t >= burn) {
              accept_gamma++;
            }
          }
          
        } // End of add step
        
      } 
    } // End of updating gamma

    gamma_sum(t) = gamma_sum_temp;

    // Random walk Metropolis-Hastings (RWMH) algorithm for updating M =========
    for (j = 0; j < p; j++) {
      if (gamma(j) == 0) {
        mu_new = exp(rnorm(1, log(u(0, j)), tau_log_mu0)(0));
        hastings = 0;

        for (i = 0; i < n; i++) {
         if (H(i, j) == 0){
           hastings = hastings + phi(j)*(-log(mu_new*s(i) + phi(j))) + Y(i, j)*(log(mu_new) - log(mu_new*s(i) + phi(j)));
           hastings = hastings - (phi(j)*(-log(u(0, j)*s(i) + phi(j))) + Y(i, j)*(log(u(0, j)) - log(u(0, j)*s(i) + phi(j))));
         }

        }

        hastings = hastings + (a_0 - 1)*log(mu_new) - b_0*mu_new;
        hastings = hastings - ((a_0 - 1)*log(u(0, j)) - b_0*u(0, j));

        if (t >= burn) {
          try_mu++;
        }

        if(hastings >= log(double(rand()%10001)/10000)) {
          for (k = 0; k < K; k++) {
            u(k, j) = mu_new;
          }

          if (t >= burn) {
            accept_mu++;
          }
        }

      } // End of updating mu0

      else{
        for (int k = 0; k < K; k++) {
          mu_new = exp(rnorm(1, log(u(k, j)), tau_log_mu)(0));

          for (int i = 0; i < n; i++) {
            if (z(i) == k && H(i, j) == 0) {
              hastings = hastings + phi(j)*(-log(mu_new*s(i) + phi(j))) + Y(i, j)*(log(mu_new) - log(mu_new*s(i) + phi(j)));
              hastings = hastings - (phi(j)*(-log(u(k, j)*s(i) + phi(j))) + Y(i, j)*(log(u(k, j)) - log(u(k, j)*s(i) + phi(j))));
            }

          }

          hastings = hastings + (a_0 - 1)*log(mu_new) - b_0*mu_new;
          hastings = hastings - ((a_0 - 1)*log(u(k, j)) - b_0*u(k, j));

          if (t >= burn) {
            try_mu++;
          }

          if(hastings >= log(double(rand()%10001)/10000)) {
            u(k, j) = mu_new;

            if (t >= burn) {
              accept_mu++;
            }
          }
        }
      } // End of updating muk
    } // End of updating mu


    // Monitor the process =====================================================
    if((t*100/(iter-1)) == count) {
      Rcout<<count<< "% has been done\n";
      count = count + 10;

    }

    // Calculate the marginal posterior probability to obtain MAP ==============
    double map_z = 0;
    
    for (i = 0; i < n; i++) {
      n_neighbor_temp = 0;
      
      // Compute I(z_i = z_i')
      for (int ii = 0; ii < n_neighbor; ii++){
        if (P(i, ii) != 0) {
          if(z(i) == z(P(i, ii) - 1)){
            n_neighbor_temp = n_neighbor_temp + 1;
          }
        }
      }
      
      map_z = map_z + e(z(i)) + f*n_neighbor_temp;
      
      for (j = 0; j < p; j++){
        if (H(i, j) == 0) {
          int index = z(i);
          map_z = map_z + lgamma(Y(i, j) + phi(j)) - lgamma(phi(j)) - lgamma(Y(i, j) + 1) +
            phi(j)*(log(phi(j)) - log(u(index, j)*s(i) + phi(j))) + Y(i, j)*(log(u(index, j)*s(i)) - log(u(index, j)*s(i) + phi(j)));
        }
      }
    }
    
    map_z_store(t) = map_z;
    
    
    // Calculate the marginal posterior probability of gamma to obtain MAP =====
    double map_gamma = 0;
    
    for (j = 0; j < p; j++) {
      map_gamma = map_gamma + gamma(j)*(log(omega) - log(1 - omega));
      
      for (i = 0; i < n; i++) {
        if (H(i, j) == 0) {
          int index = z(i);
          map_gamma = map_gamma + lgamma(Y(i, j) + phi(j)) - lgamma(phi(j))  - lgamma(Y(i, j) + 1) +
            phi(j)*(log(phi(j)) - log(u(index, j)*s(i) + phi(j))) + Y(i, j)*(log(u(index, j)*s(i)) - log(u(index, j)*s(i) + phi(j)));
        }
      }
    }
    
    map_gamma_store(t) = map_gamma;
    

    // Obtain the MAP estimates of z and gamma =================================
    if (t == burn) {
      z_map = z;
      map_z_temp = map_z;

      gamma_map = gamma;
      map_gamma_temp = map_gamma;
      
    }else if (t > burn) {
      if (map_z > map_z_temp) {
        z_map = z;
        map_z_temp = map_z;
      }

      if (map_gamma > map_gamma_temp) {
        gamma_map = gamma;
        map_gamma_temp = map_gamma;
      }
    }

    // Obtain full posterior ===================================================
    double map = 0;
    
    for (i = 0; i < n; i++) {
      n_neighbor_temp = 0;
      
      // Compute I(z_i = z_i')
      for (int ii = 0; ii < n_neighbor; ii++){
        if (P(i, ii) != 0) {
          if(z(i) == z(P(i, ii) - 1)){
            n_neighbor_temp = n_neighbor_temp + 1;
          }
        }
      }
      
      map_z = map_z + e(z(i)) + f*n_neighbor_temp;
      
    }
    
    for (j = 0; j < p; j++) {
      map = map + gamma(j)*(log(omega) - log(1 - omega));
      map = map + (a_phi - 1)*log(phi(j)) - b_phi*phi(j);
      
      if (gamma(j) == 0) {
        map = map + (a_0 - 1)*log(u(0, j)) - b_0*u(0, j);
      }else{
        for (k = 0; k < K; k++) {
          map = map + (a_mu - 1)*log(u(k, j)) - b_mu*u(k, j);
          
        }
      }
      
      for (i = 0; i < n; i++) {
        map = map + H(i, j)*(log(pi) - log(1 - pi));
        int index = z(i);
        
        if (H(i, j) == 0) {
          map = map + lgamma(Y(i, j) + phi(j)) - lgamma(phi(j))  - lgamma(Y(i, j) + 1) +
            phi(j)*(log(phi(j)) - log(u(index, j)*s(i) + phi(j))) + Y(i, j)*(log(u(index, j)*s(i)) - log(u(index, j)*s(i) + phi(j)));
        }
        
      }
    }
    
    map_store(t) = map;
    
    
    // Store the results =======================================================
    for (i = 0; i < n; i++){
      z_store(t, i) = z(i);
    }
    u_store.slice(t) = as<arma::mat>(u);
    
    H_sum(t) = 0;
    for (i = 0; i < n; i++){
      H_sum(t) = H_sum(t) + H_sum_temp(i);
    }
    
    for (j = 0; j < p; j++) {
      phi_store(t, j) = phi(j);
      gamma_store(t, j) = gamma(j);
      
      if (t >= burn){
        gamma_ppi(j) = gamma_ppi(j) + gamma(j);
        for (i = 0; i < n; i++) {
          H_ppi(i, j) = H_ppi(i, j) + H(i, j);
      
        }
      }
    }
    
  } // End of iterations
  
  for (j = 0; j < p; j++) {
    gamma_ppi(j) = gamma_ppi(j)/(iter - burn);
    gamma_BF(j) = gamma_ppi(j)/(1 - gamma_ppi(j))*b_omega/a_omega;
    
    for (i = 0; i < n; i++) {
      H_ppi(i, j) = H_ppi(i, j)/(iter - burn);
    }
  }
  
  
  accept_phi = accept_phi/try_phi;
  accept_gamma = accept_gamma/try_gamma;
  accept_mu = accept_mu/try_mu;

  return Rcpp::List::create(Rcpp::Named("phi_store") = phi_store, 
                            Rcpp::Named("z_store") = z_store, Rcpp::Named("map_z_store") = map_z_store, Rcpp::Named("z_map") = z_map,
                            Rcpp::Named("H_sum") = H_sum, Rcpp::Named("H_ppi") = H_ppi,           
                            Rcpp::Named("u_store") = u_store, 
                            Rcpp::Named("gamma_store") = gamma_store, Rcpp::Named("gamma_sum") = gamma_sum, Rcpp::Named("gamma_ppi") = gamma_ppi,
                            Rcpp::Named("gamma_BF") = gamma_BF, Rcpp::Named("map_gamma_store") = map_gamma_store, Rcpp::Named("gamma_map") = gamma_map,
                            Rcpp::Named("accept_mu") = accept_mu, Rcpp::Named("accept_phi") = accept_phi, Rcpp::Named("accept_gamma") = accept_gamma,
                            Rcpp::Named("map_store") = map_store);
}



// int num(IntegerVector x, int c) {
//   int n = x.size();
//   int count = 0;
//   int i;
//   for (i = 0; i < n; i++) {
//     if (x(i) == c) {
//       count++;
//     }
//   }
//   return count;
// }

// arma::uvec which(IntegerVector x, int c) {
//   int n = x.size();
//   int count = 0;
//   int i;
//   int m = num(x, c);
//   arma::uvec index(m);
//   for (i = 0; i < n; i++) {
//     if (x(i) == c) {
//       index(count) = i;
//       count++;
//     }
//   }
//   return index;
// }


// arma::uvec remove(NumericVector x, int c) {
//   int n = x.size();
//   int count = 0;
//   arma::uvec y(n-1);
//   for (int i = 0; i < n; i++) {
//     if (i != c) {
//       y(count) = x(i);
//       count++;
//     }
//   }
//   return y;
// }


// Rcpp::Environment base("package:base");
// Function do_unique = base["unique"];
// 
// arma::uvec organize_label(IntegerVector x) {
//   int n = x.size();
//   arma::uvec x_new(n);
//   int count = 1;
//   int K = (Rcpp::unique(x)).size();
//   IntegerVector cluster = do_unique(x);
// 
//   for (int j = 0; j < K; j++) {
//     for (int i = 0; i < n; i++) {
//       if (x(i) == cluster[j]) {
//         x_new(i) = count;
//       }
//     }
//     count = count + 1;
//   }
//   return x_new;
// }
// 
