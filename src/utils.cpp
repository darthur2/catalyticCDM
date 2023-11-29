#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double compute_prior_rate(NumericMatrix sg,
                          NumericVector p,
                          double s_hat,
                          double g_hat,
                          NumericVector p_hat,
                          NumericMatrix R){
  
  int J = sg.nrow();
  int C = p_hat.size();
  double prior_rate = 0;
  
  NumericVector s1(J);
  NumericVector s2(J);
  NumericVector g1(J);
  NumericVector g2(J);
  
  for (int j = 0; j < J; j++){
    s1[j] = log(sg(j,0));
    s2[j] = log(1-sg(j,0));
    g1[j] = log(sg(j,1));
    g2[j] = log(1-sg(j,1));
  }
  
  for (int c = 0; c < C; c++){
    prior_rate += p_hat[c]*log(p[c]);
  }
  
  for (int c = 0; c < C; c++){
    for (int j = 0; j < J; j++){
      prior_rate += s_hat*p_hat[c]*R(c,j)*s1[j] + 
        (1-s_hat)*p_hat[c]*R(c,j)*s2[j] +
        g_hat*p_hat[c]*(1-R(c,j))*g1[j] +
        (1-g_hat)*p_hat[c]*(1-R(c,j))*g2[j];
    }
  }
  
  prior_rate = -1*prior_rate;
  
  return prior_rate;
}

// [[Rcpp::export]]
NumericMatrix compute_z_post_unnorm(NumericMatrix Y,
                                    NumericMatrix R,
                                    NumericVector p,
                                    NumericMatrix sg){
  int N = Y.nrow();
  int J = R.ncol();
  int C = R.nrow();
  
  NumericMatrix PZ(N, C);
  NumericVector s1(J);
  NumericVector s2(J);
  NumericVector g1(J);
  NumericVector g2(J);
  
  for (int j = 0; j < J; j++){
    s1[j] = log(sg(j,0));
    s2[j] = log(1-sg(j,0));
    g1[j] = log(sg(j,1));
    g2[j] = log(1-sg(j,1));
  }
  
  for (int i = 0; i < N; i++){
    for (int c = 0; c < C; c++){
      for (int j = 0; j < J; j++){
        PZ(i,c) += (1-Y(i,j))*R(c,j)*s1[j] + Y(i,j)*R(c,j)*s2[j] + 
          Y(i,j)*(1-R(c,j))*g1[j] + (1-Y(i,j))*(1-R(c,j))*g2[j];
      }
      PZ(i,c) *= exp(PZ(i,c))*p[c];
    }
  }
  
  return PZ;
}


// [[Rcpp::export]]
NumericVector row_sums(NumericMatrix X) {
  int nrow = X.nrow();
  int ncol = X.ncol();
  NumericVector out(nrow);
  
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += X(i, j);
    }
    out[i] = total;
  }
  return out;
}


// [[Rcpp::export]]
IntegerVector sample_one_multinomial(NumericVector p) {
  int k = p.size();
  IntegerVector multinom(k);
  rmultinom(1, p.begin(), k, multinom.begin());
  return multinom;
}


// [[Rcpp::export]]
IntegerMatrix sample_all_multinomial(NumericMatrix P){
  int nrow = P.nrow();
  int ncol = P.ncol();
  
  IntegerMatrix Z(nrow, ncol);
  
  for (int i = 0; i < nrow; i++){
    Z(i,_) = sample_one_multinomial(P(i,_));
  }
  return Z;
}


// [[Rcpp::export]]
NumericMatrix compute_dina_probs(NumericMatrix eta, NumericVector slip, NumericVector guess){
  int nrow = eta.nrow();
  int ncol = eta.ncol();
  NumericMatrix probs(nrow, ncol);
  
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      probs(i,j) = pow(1-slip[j], eta(i,j))*pow(guess[j], 1-eta(i,j));
    }
  }
  
  return probs;
}


// [[Rcpp::export]]
IntegerMatrix sample_dina_data(NumericMatrix probs){
  int nrow = probs.nrow();
  int ncol = probs.ncol();
  IntegerMatrix data(nrow, ncol);
  
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      data(i,j) = R::rbinom(1.0, probs(i,j));
    }
  } 
  
  return data;
}

// [[Rcpp::export]]
NumericMatrix compute_dina_probs_marginal(NumericMatrix Y, 
                                          NumericMatrix R, 
                                          NumericVector slip, 
                                          NumericVector guess, 
                                          NumericVector p) {
  int N = Y.nrow();
  int J = Y.ncol();
  int C = R.nrow();
  
  NumericMatrix probs(N, J);
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < J; j++) {
      double prob_ij = 0.0;
      for(int c = 0; c < C; c++) {
        prob_ij += p(c) * pow(slip(j), (1 - Y(i, j)) * R(c, j)) * pow(1 - slip(j), Y(i, j) * R(c, j)) *
          pow(guess(j), Y(i, j) * (1 - R(c, j))) * pow(1 - guess(j), (1 - Y(i, j)) * (1 - R(c, j)));
      }
      probs(i, j) = prob_ij;
    }
  }
  
  return probs;
}

double log_gamma_function(double x) {
  return R::lgammafn(x);
}

// [[Rcpp::export]]
double compute_tau_lp(double tau,
                      NumericVector p,
                      NumericVector p_hat,
                      NumericVector slip,
                      double s_hat, 
                      NumericVector guess,
                      double g_hat,
                      NumericVector pj_hat,
                      double tau_a, 
                      double tau_b) {
  
  int K = log2(p.size());
  int J = slip.size();
  
  double log_p = log_gamma_function(tau + pow(2, K));
  
  for (int c = 0; c < pow(2, K); c++) {
    log_p -= log_gamma_function(tau * p_hat[c] + 1);
    log_p += tau * p_hat[c] * log(p[c]);
  }
  
  double log_sg = 0.0;
  
  for (int j = 0; j < J; j++) {
    double a_s = tau * s_hat * pj_hat[j] + 1;
    double b_s = tau * (1 - s_hat) * pj_hat[j] + 1;
    double a_g = tau * g_hat * (1 - pj_hat[j]) + 1;
    double b_g = tau * (1 - g_hat) * (1 - pj_hat[j]) + 1;
    
    log_sg += log_gamma_function(a_s + b_s) - log_gamma_function(a_s) - log_gamma_function(b_s);
    log_sg += (a_s - 1) * log(slip[j]) + (b_s - 1) * log(1.0 - slip[j]);
    
    log_sg += log_gamma_function(a_g + b_g) - log_gamma_function(a_g) - log_gamma_function(b_g);
    log_sg += (a_g - 1) * log(guess[j]) + (b_g - 1) * log(1.0 - guess[j]);
  }
  
  double log_tau_prior = (tau_a - 1) * log(tau) - tau * tau_b;
  
  double log_tau_posterior = log_p + log_sg + log_tau_prior;
  
  return log_tau_posterior;
}


