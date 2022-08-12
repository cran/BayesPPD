
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#define NDEBUG
#include <RcppEigen.h>
#include <RcppNumerical.h>





using namespace Rcpp;
using namespace Numer;



// This class defined the integrand, later used by the numerical integration function. 

class Mintegrand: public Func {

private:
  const double b_c1_post;
  const double b_c2_post;
  const double b_t1_post;
  const double b_t2_post;
  const double delta;
  const std::string ns;
  const double upper_inf;
  const std::string dType;

public:
  Mintegrand(double b_c1_post_, double b_c2_post_, double b_t1_post_, double b_t2_post_, double delta_, std::string ns_, double upper_inf_, std::string dType_) :
  b_c1_post(b_c1_post_), b_c2_post(b_c2_post_), b_t1_post(b_t1_post_), b_t2_post(b_t2_post_), delta(delta_), ns(ns_), upper_inf(upper_inf_), dType(dType_){}



  double operator()(const double& x) const{

    double res= 0;
    if(ns == ">"){
      if(dType == "Bernoulli"){
        res = R::dbeta(x, b_c1_post, b_c2_post,FALSE) * R::pbeta(x+delta, b_t1_post, b_t2_post, TRUE, FALSE);
      }
      if(dType == "Poisson"){
        res = R::dgamma(x, b_c1_post, 1/b_c2_post,FALSE) * R::pgamma(x+delta, b_t1_post, 1/b_t2_post, TRUE, FALSE);
      }
      if(dType == "Exponential"){
        res = R::dgamma(x, b_c1_post, 1/b_c2_post,FALSE) * R::pgamma(x*delta, b_t1_post, 1/b_t2_post, TRUE, FALSE);
      }
    }else{
      if(dType == "Bernoulli"){
        res = R::dbeta(x, b_c1_post, b_c2_post,FALSE) * (1 - R::pbeta(x+delta, b_t1_post, b_t2_post, TRUE, FALSE));
      }
      if(dType == "Poisson"){
        res = R::dgamma(x, b_c1_post, 1/b_c2_post,FALSE) * (1 - R::pgamma(x+delta, b_t1_post, 1/b_t2_post, TRUE, FALSE));
      }
      if(dType == "Exponential"){
        res = R::dgamma(x, b_c1_post, 1/b_c2_post,FALSE) * (1 - R::pgamma(x*delta, b_t1_post, 1/b_t2_post, TRUE, FALSE));
      }
    }
    
    return(res);
  }

};


// This function performs numerical integration for Bernoulli, Poisson and Exponential responses. 
double num_integrate(std::string dType, double b_c1_post, double b_c2_post, double b_t1_post, double b_t2_post, double delta, std::string ns, double upper_inf){
  Mintegrand f(b_c1_post, b_c2_post, b_t1_post, b_t2_post, delta, ns, upper_inf, dType);
  double err_est;
  int err_code;
  double result = 0;
  if(dType == "Bernoulli"){
    result = integrate(f, 0.0, 1.0, err_est, err_code);
  }
  if(dType == "Poisson"){
    result = integrate(f, 0.0, upper_inf, err_est, err_code);
  }
  if(dType == "Exponential"){
    result = integrate(f, 0.0, upper_inf, err_est, err_code);
  }
  return(result);
}


// This function computes the parameters of the posterior distribution of mu_c
// for Bernoulli, Poisson and Exponential responses.
// [[Rcpp::export]]
arma::vec two_grp_fixed_a0(std::string dType, double & y_c, double & n_c,
                           arma::mat & historical, double b_01, double b_02){
  double b_c1_post = 0;
  double b_c2_post = 0;
  
  if(dType == "Bernoulli"){
    
    b_c1_post = (double) y_c + arma::dot(historical.col(0), historical.col(2)) + b_01;
    b_c2_post = (double) n_c - (double) y_c +  arma::dot((historical.col(1) - historical.col(0)), historical.col(2)) + b_02;
  }
  if(dType == "Poisson"){
    // y here is the sum of all the individual y's
    b_c1_post = (double) y_c + arma::dot(historical.col(0), historical.col(2)) + b_01;
    b_c2_post = (double) n_c + arma::dot(historical.col(1), historical.col(2)) + b_02;
  }
  if(dType == "Exponential"){
    b_c1_post = (double) n_c + arma::dot(historical.col(1), historical.col(2)) + b_01;
    b_c2_post = (double) y_c + arma::dot(historical.col(0), historical.col(2)) + b_02;
  }
  
  arma::vec result(2);
  result[0] = b_c1_post;
  result[1] = b_c2_post;
  
  return(result);
  
}


// This function performs sample size determination for Bernoulli, Poisson and Exponential responses.
// This function calls the num_integrate function for computing the posterior probability. 
// [[Rcpp::export]]
Rcpp::List power_two_grp_fixed_a0(std::string dType, double n_t, double n_c,
                                  arma::mat historical, std::string ns,
                                  NumericVector p_t_prior_samps,
                                  NumericVector p_c_prior_samps,
                                  double b_t1, double b_t2, double b_01,
                                  double b_02, double delta, double gamma,
                                  int N, double upper_inf) {

  RNGScope scope;



  NumericVector post_samples(N);

  for (int i=0;i<N;i++){

    // sample a p_t and p_c from the sampling distribution
    int ind_t = floor(R::runif(0,p_t_prior_samps.size()));
    int ind_c = floor(R::runif(0,p_c_prior_samps.size()));
    double p_t = p_t_prior_samps[ind_t];
    double p_c = p_c_prior_samps[ind_c];

    double b_t1_post = 0;
    double b_t2_post = 0;
    double b_c1_post = 0;
    double b_c2_post = 0;

    if(dType == "Bernoulli"){

      //generate data
      int y_t = R::rbinom(n_t, p_t);
      int y_c = R::rbinom(n_c, p_c);

      b_t1_post = (double) y_t + b_t1;
      b_t2_post = (double) n_t - (double) y_t + b_t2;
      b_c1_post = (double) y_c + arma::dot(historical.col(0), historical.col(2)) + b_01;
      b_c2_post = (double) n_c - (double) y_c +  arma::dot((historical.col(1) - historical.col(0)), historical.col(2)) + b_02;
    }
    if(dType == "Poisson"){
      int y_t = R::rpois(p_t*n_t); // y here is the sum of all the individual y's
      int y_c = R::rpois(p_c*n_c);

      b_t1_post = (double) y_t + b_t1;
      b_t2_post = (double) n_t + b_t2;
      b_c1_post = (double) y_c + arma::dot(historical.col(0), historical.col(2)) + b_01;
      b_c2_post = (double) n_c + arma::dot(historical.col(1), historical.col(2)) + b_02;
    }
    if(dType == "Exponential"){
      double y_t = R::rgamma(n_t, 1/p_t);
      double y_c = R::rgamma(n_c, 1/p_c);

      b_t1_post = (double) n_t + b_t1;
      b_t2_post = (double) y_t + b_t2;
      b_c1_post = (double) n_c + arma::dot(historical.col(1), historical.col(2)) + b_01;
      b_c2_post = (double) y_c + arma::dot(historical.col(0), historical.col(2)) + b_02;
    }

    double post_result = num_integrate(dType, b_c1_post, b_c2_post, b_t1_post, b_t2_post, delta, ns, upper_inf);


    post_samples[i] = post_result;

  }

  double beta = mean(post_samples >= gamma);

  return Rcpp::List::create(
    Rcpp::Named("power/type I error")    = beta);
}


  

// This function generates posterior samples of mu_c, tau, and tau_0k 
// using Gibbs sampling for normal responses.
// [[Rcpp::export]]
Rcpp::List two_grp_fixed_a0_normal(double & y_c, double & n_c, double & v, arma::mat & historical, int & nMC, int & nBI) {
  Rcpp::RNGScope scope;

  arma::vec gibbs_mu_c(nMC+nBI);
  gibbs_mu_c.zeros();
  arma::vec gibbs_tau(nMC+nBI);
  gibbs_tau.ones();
  arma::mat gibbs_tau0(nMC+nBI,historical.n_rows);
  gibbs_tau0.ones();



  for (int i=1;i<nMC+nBI;i++){

    //sampling mu
    double num = 0;
    double denom = 0;
    for(int j=0;j<historical.n_rows;j++){
      num = num + historical(j,0) * historical(j,3) * gibbs_tau0(i-1,j);
      denom = denom + historical(j,1) * historical(j,3) * gibbs_tau0(i-1,j);
    }
    double mu = (y_c*gibbs_tau[i-1] + num)/(n_c*gibbs_tau[i-1] + denom);
    double var = 1/(n_c*gibbs_tau[i-1] + denom);
    gibbs_mu_c[i] = R::rnorm(mu, sqrt(var));


    //sampling tau
    double y2sum = (n_c-1)*v + n_c*pow(y_c/n_c,2);
    double ss = y2sum - 2*gibbs_mu_c[i]*y_c + n_c*pow(gibbs_mu_c[i],2);
    gibbs_tau[i] = R::rgamma(n_c/2, 1/(ss/2));


    //sampling tau0
    if(historical(0,1)!=0){ // sample size is not zero
      for(int k=0;k<historical.n_rows;k++){
        double y_h = historical(k,0);
        double n_h = historical(k,1);
        double v_h = historical(k,2);
        double y2sum = (n_h-1)*v_h + n_h*pow(y_h/n_h,2);
        double ss_h = y2sum - 2*gibbs_mu_c[i]*y_h + n_h*pow(gibbs_mu_c[i],2);
        gibbs_tau0(i,k) = R::rgamma(historical(k,3)*n_h/2, 1/(historical(k,3)*ss_h/2));
      }
    }


  }
  arma::vec gibbs_mu_c_sub = gibbs_mu_c.subvec(nBI, nMC+nBI-1);
  arma::vec gibbs_tau_sub = gibbs_tau.subvec(nBI, nMC+nBI-1);


  arma::mat gibbs_tau0_sub = gibbs_tau0.rows(nBI, nMC+nBI-1);

  Rcpp::List result;
  // if no historical data is provided, no tau0 samples are provided
  if(historical(0,1)==0){
    result = Rcpp::List::create(
      Rcpp::Named("posterior samples of mu_c")    = gibbs_mu_c_sub,
      Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub);
  }else{
    result = Rcpp::List::create(
      Rcpp::Named("posterior samples of mu_c")    = gibbs_mu_c_sub,
      Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub,
      Rcpp::Named("posterior samples of tau_0")    = gibbs_tau0_sub);
  }




  return(result);
}



// This function performs sample size determination for normal responses.
// This function calls the two_grp_fixed_a0_normal() function to obtain posterior samples of mu_c. 
// [[Rcpp::export]]
Rcpp::List power_two_grp_fixed_a0_normal(double n_t, double n_c, arma::mat historical, std::string ns,
                                         NumericVector mu_t_prior_samps,
                                         NumericVector mu_c_prior_samps,
                                         NumericVector var_t_prior_samps,
                                         NumericVector var_c_prior_samps,
                                         double delta, double gamma, int nMC, int nBI,
                                         int N) {

  RNGScope scope;



  NumericVector post_samples(N);
  NumericVector mu_c_samples(N);
  NumericVector mu_t_samples(N);
  NumericVector tau_samples(N);
  arma::mat tau0_samples(N, historical.n_rows);

  for (int i=0;i<N;i++){

    int ind_t = floor(R::runif(0,mu_t_prior_samps.size()));
    int ind_c = floor(R::runif(0,mu_c_prior_samps.size()));
    int ind_vt = floor(R::runif(0,var_t_prior_samps.size()));
    int ind_vc = floor(R::runif(0,var_c_prior_samps.size()));
    double mu_t = mu_t_prior_samps[ind_t];
    double mu_c = mu_c_prior_samps[ind_c];
    double var_t = var_t_prior_samps[ind_vt];
    double var_c = var_c_prior_samps[ind_vc];


    NumericVector y_t = Rcpp::rnorm(n_t, mu_t, sqrt(var_t));
    arma::vec y_c = Rcpp::rnorm(n_c, mu_c, sqrt(var_c));

    double s_c = sum(y_c);
    double v_c = var(y_c);

    Rcpp::List lis = two_grp_fixed_a0_normal(s_c, n_c, v_c, historical, nMC, nBI);
    NumericVector mu_c_post = lis[0];
    NumericVector mu_t_post = Rcpp::rt(nMC, n_t-1) * sd(y_t)/sqrt(n_t) + mean(y_t);
    NumericVector tau_post = lis[1];
    arma::mat tau0_post;
    if(historical(0,1)!=0){
      arma::mat m = lis[2];
      tau0_post = m;
    }


    mu_c_samples[i] = mean(mu_c_post);
    mu_t_samples[i] = mean(mu_t_post);
    tau_samples[i] = mean(tau_post);


    // calculate posterior mean of tau0
    arma::vec tau0_vec(tau0_post.n_cols);
    if(historical(0,1)!=0){
      for(int k = 0; (unsigned)k < tau0_post.n_cols; k++){
        tau0_vec(k) = arma::mean(tau0_post.col(k));
      }
      tau0_samples.row(i) = tau0_vec.t();
    }


    if(ns == ">"){
      post_samples[i] = mean(mu_t_post - mu_c_post < delta);
    }else{
      post_samples[i] = mean(mu_t_post - mu_c_post > delta);
    }
    

  }

  double mean_mu_t = mean(mu_t_samples);
  double mean_mu_c = mean(mu_c_samples);
  double mean_tau = mean(tau_samples);

  // calculate average posterior means of tau0
  arma::vec final_tau0_vec(tau0_samples.n_cols);
  if(historical(0,1)!=0){
    for(int k = 0; (unsigned)k < tau0_samples.n_cols; k++){
      final_tau0_vec(k) = arma::mean(tau0_samples.col(k));
    }
  }

  double power = mean(post_samples >= gamma);

  if(historical(0,1)==0){
    return Rcpp::List::create(
      Rcpp::Named("average posterior mean of mu_t")  = mean_mu_t,
      Rcpp::Named("average posterior mean of mu_c")  = mean_mu_c,
      Rcpp::Named("average posterior mean of tau")  = mean_tau,
      Rcpp::Named("power/type I error")    = power);
  }else{
    return Rcpp::List::create(
      Rcpp::Named("average posterior mean of mu_t")  = mean_mu_t,
      Rcpp::Named("average posterior mean of mu_c")  = mean_mu_c,
      Rcpp::Named("average posterior mean of tau")  = mean_tau,
      Rcpp::Named("average posterior mean of tau0")  = final_tau0_vec,
      Rcpp::Named("power/type I error")    = power);
  }
}




