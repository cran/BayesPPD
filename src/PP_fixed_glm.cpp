// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>





using namespace Rcpp;

//------------------------------------ Class Specific Functions -------------------------------------------//
class glm{

public:

  // data
  arma::vec y;  // outcome data
  arma::vec n;  // number of trials (for binomial data)
  arma::mat x;  // covariates
  Rcpp::List historical; // historical data



  // model definition
  std::string dType;
  std::string dLink;
  bool dCurrent; // if false, then current likelihood contribution is zero

  int P;


  arma::vec                lower_limits;
  arma::vec                upper_limits;
  arma::vec                slice_widths;
  int m;


  // public member functions;
  glm(std::string & dType0, std::string & dLink0, arma::vec & y0, arma::vec & n0, arma::mat & x0, Rcpp::List & historical0,
           arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, bool & dCurrent);

  double logFC(const arma::vec & parm0, const int & p);

};

glm::glm(	std::string & dType0, std::string & dLink0, arma::vec & y0, arma::vec & n0, arma::mat & x0, Rcpp::List & historical0,
                    arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, bool & dCurrent0)
{


  dType = dType0;
  dLink = dLink0;

  y     = y0;
  if ((dType=="Bernoulli")) {n.resize(y.size()); n.ones();} else {n = n0;}
  x = x0;
  historical = historical0;


  dCurrent = dCurrent0;

  if (dCurrent0==TRUE){
    P = x.n_cols; // here an intercept has been added to x already

  }else{
    Rcpp::List dat = historical0[0];
    arma::mat x_h = dat["x0"];
    P = x_h.n_cols + 1; // add intercept
  }


  lower_limits = lower_limits0;
  upper_limits = upper_limits0;
  slice_widths = slice_widths0;

  m=10;
}

// Define the log likelihood of the posterior using the power prior
// for GLMs for Bernoulli, Poisson and Exponential responses.
double glm::logFC(const arma::vec & parm0, const int & p)
{
  // extract regression parameters;
  arma::vec beta0 = parm0;
  arma::vec mean;

  arma::vec beta_h = beta0;


  // compute likelihood contribution;
  double lp = 0;

  if (dCurrent==TRUE){

    arma::uvec ind;
    ind << 0;
    beta_h = join_cols(parm0.elem(ind), parm0.subvec(2,P-1));

    // compute mean parameter;
    mean = x*beta0;


    if (dLink=="Logistic") 		    { mean = exp(mean) / (1 + exp(mean)); 					}
    if (dLink=="Probit")
      { for (int k=0; (unsigned)k< x.n_rows;k++)   { mean[k] = R::pnorm(mean[k], 0.0, 1.0, true, false);  } }
    if (dLink=="Log") 	   		                     { mean = exp(mean); 				    				}
    if (dLink=="Identity-Positive")
      { 	for (int k=0; (unsigned)k< x.n_rows;k++) { mean[k] = std::max(mean[k],0.00001); }	}
    if (dLink=="Identity-Probability")
      { 	for (int k=0; (unsigned)k< x.n_rows;k++) { mean[k] = std::min(std::max(mean[k],0.00001),0.99999); }	}
    if (dLink=="Complementary Log-Log")    {  mean = 1 - exp(-exp(mean)); }

    
    // compute current data likelihood
    if ((dType=="Bernoulli") | (dType=="Binomial"))	{lp += sum( y % log(mean) + (n-y) % log(1-mean) );}
    if (dType=="Poisson")	{ lp += sum(  y % log(mean) - mean );}
    if (dType=="Exponential")	{ lp += sum(  log(mean) - y % mean );}

  }


  // compute historical data likelihood
  for(int i = 0; i < historical.size(); i++){
    Rcpp::List dat = historical[i];
    arma::vec y_h = dat["y0"];
    arma::mat x_h = dat["x0"];
    arma::vec v(x_h.n_rows);
    v.ones();
    x_h.insert_cols(0, v);
    double a0 = dat["a0"];
    arma::vec n_h;
    n_h.zeros();
    if (dType=="Bernoulli") {n_h.resize(y_h.size()); n_h.ones();} 
    if (dType=="Binomial") {n_h = Rcpp::as<arma::vec>(dat["n0"]);}

    arma:: vec mean_h = x_h*beta_h;


    if (dLink=="Logistic") 		    { mean_h = exp(mean_h) / (1 + exp(mean_h)); 					}
    if (dLink=="Probit")
      { for (int k=0; (unsigned)k< x_h.n_rows;k++)   { mean_h[k] = R::pnorm(mean_h[k], 0.0, 1.0, true, false);  } }
    if (dLink=="Log") 	   		    { mean_h = exp(mean_h); 				    				}
    if (dLink=="Identity-Positive")
      { 	for (int k=0; (unsigned)k< x_h.n_rows;k++) { mean_h[k] = std::max(mean_h[k],0.00001); }						}
    if (dLink=="Identity-Probability")
      { 	for (int k=0; (unsigned)k<x_h.n_rows;k++) { mean_h[k] = std::min(std::max(mean_h[k],0.00001),0.99999); }	}
    if (dLink=="Complementary Log-Log")    {  mean_h = 1 - exp(-exp(mean_h)); }


    if ((dType=="Bernoulli") | (dType=="Binomial")) {lp += a0 * sum( y_h % log(mean_h) + (n_h - y_h) % log(1-mean_h) );}
    if (dType=="Poisson")	{ lp += a0 * sum( y_h % log(mean_h) - mean_h );}
    if (dType=="Exponential")	{lp += a0 * sum( log(mean_h) - y_h % mean_h );}


  }

  return  lp;
}


// slice sampler
void slice( arma::vec & parms, glm & b)
{

  double b0, f0, f0_L, f0_R, f0_x1, h0, L, R, V, J, K,w,lower,upper;
  arma::vec parm0;

  for (int p = 0; p < b.P; p++)
  {
    // create vector of parameters to modify for slice sampling;
    parm0 = parms;

    // extract slice width and parameter bounds;
    w     = b.slice_widths[p];
    lower = b.lower_limits[p];
    upper = b.upper_limits[p];

    // skip over fixed parameter values;
    if (lower==upper){parms(p) = lower;}
    else
    {
      // current value of the parameter in question;
      b0 = parm0(p);

      // calculate current full conditional value;
      f0 = b.logFC(parm0,p);

      // calculate height of the horizontal slice;
      h0 = f0 - R::rexp(1.0);

      // Calculate initial horizontal interval;
      L = parm0(p) - R::runif(0.0,1.0)*w;
      R = L+w;

      // Truncate bounds to support of the parameter space;
      L = std::max(L,lower);
      R = std::min(R,upper);

      // Step out;
      V = R::runif(0.0,1.0);
      J = floor(b.m*V);
      K = (b.m-1)-J;

      // compute log of full conditional at current boundaries;
      parm0(p) = L; f0_L = b.logFC(parm0,p);
      parm0(p) = R; f0_R = b.logFC(parm0,p);

      while(J>0 and h0<f0_L and L>=lower)
      {
        L        = L-w; if (L<=lower) {L=lower;}
        J        = J-1;
        parm0(p) = L;
        f0_L     = b.logFC(parm0,p);
      }
      while(K>0 and h0<f0_R and R<=upper)
      {
        R        = R+w; if (R>=upper) {R=upper;}
        K        = K-1;
        parm0(p) = R;
        f0_R     = b.logFC(parm0,p);
      }

      // perform rejection sampling;
      int stop  = 0;
      while(stop == 0)
      {
        parm0(p)     = L + R::runif(0.0,1.0)*(R-L);
        f0_x1        = b.logFC(parm0,p);

        if      ( h0       <  f0_x1 ) { parms(p) = parm0(p); stop = 1;  }
        else if ( parm0(p) <  b0    ) { L = parm0(p);                     }
        else if ( parm0(p) >= b0    ) { R = parm0(p);                     }

        if (-0.0000000001 <= L-R and L-R <= 0.0000000001)
        {
          parms(p)= 0.5*(L+R);
          stop      = 1;
        }
      }
    }
  }
}


// This function generates posterior samples of beta using slice sampling 
// for GLMs with Bernoulli, Poisson and Exponential responses.
// [[Rcpp::export]]
arma::mat glm_fixed_a0(std::string & dType0, std::string & dLink0, arma::vec & y0, arma::vec & n0, arma::mat & x0, Rcpp::List & historical0,
                    arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                    int nMC, int nBI, bool & dCurrent0)
{

  Rcpp::RNGScope scope;

  arma::vec v(x0.n_rows);
  v.ones();
  x0.insert_cols(0, v);


  int P;

  if (dCurrent0==TRUE){
    P = x0.n_cols;
  }else{
    Rcpp::List dat = historical0[0];
    arma::mat x_h = dat["x0"];
    P = x_h.n_cols + 1; // add intercept
  }


  // declare object and set values;
  glm b(dType0,dLink0,y0,n0,x0,historical0,lower_limits0,upper_limits0,slice_widths0,dCurrent0);



  // Construct container for mcmc samples;
  arma::mat samples(nMC,P);

  // create parameter vector container and initial values;
  arma::vec parms(P);
  for (int p=0;p<P;p++)
  {

    parms[p]= R::runif(0,1);

  }

  for (int s=-nBI;s<nMC;s++)
  {
    slice(parms,b);

    if (s>=0){	samples.row(s) = parms.t();	}
  }

  return(samples);

}


// This function generates posterior samples of beta, tau and tau_0 using Gibbs sampling 
// for the normal linear model.
// The power prior with fixed a_0 is used. 
// [[Rcpp::export]]
Rcpp::List glm_fixed_a0_normal(arma::vec & y, arma::mat & x, Rcpp::List & historical, int & nMC, int & nBI) {
  Rcpp::RNGScope scope;


  arma::vec v(x.n_rows);
  v.ones();
  x.insert_cols(0, v);

  arma::mat gibbs_beta(nMC+nBI, x.n_cols);
  gibbs_beta.zeros();
  arma::vec gibbs_tau(nMC+nBI);
  gibbs_tau.ones();
  arma::mat gibbs_tau0(nMC+nBI,historical.size());
  gibbs_tau0.ones();



  for (int i=1;i<nMC+nBI;i++){

    //sampling beta
    arma::mat num(x.n_cols,x.n_cols);
    arma::vec denom(x.n_cols);
    num.zeros();
    denom.zeros();
    for(int j=0;j<historical.size();j++){
      Rcpp::List dat = historical[j];
      arma::vec y_h = dat["y0"];
      arma::mat x_h = dat["x0"];
      arma::vec v(x_h.n_rows);
      v.ones();
      x_h.insert_cols(0, v); //insert intercept
      arma::vec ins(x_h.n_rows);
      ins.zeros();
      x_h.insert_cols(1, ins); //insert treatment indicator
      double a0 = dat["a0"];
      num = num + a0 * gibbs_tau0(i-1,j) * x_h.t() * x_h;
      denom = denom + a0 * gibbs_tau0(i-1,j) * x_h.t() * y_h;
    }
    arma::mat A = gibbs_tau[i-1]*x.t()*x + num;
    arma::vec B = gibbs_tau[i-1]*x.t()*y + denom;
    arma::mat A_i = inv(A);
    arma::mat m = arma::randn(1, x.n_cols);

    gibbs_beta.row(i) = arma::repmat(A_i*B, 1, 1).t() + m * arma::chol(A_i);


    arma::mat var = (y-x*gibbs_beta.row(i).t()).t() * (y-x*gibbs_beta.row(i).t());

    //sampling tau
    gibbs_tau[i] = R::rgamma(y.size()/2, 1/(var(0,0)/2));


    //sampling tau0
    for(int k=0;k<historical.size();k++){
      Rcpp::List dat = historical[k];
      arma::vec y_h = dat["y0"];
      arma::mat x_h = dat["x0"];
      arma::vec v(x_h.n_rows);
      v.ones();
      x_h.insert_cols(0, v);
      arma::vec ins(x_h.n_rows);
      ins.zeros();
      x_h.insert_cols(1, ins);
      double a0 = dat["a0"];
      arma::mat var = (y_h-x_h*gibbs_beta.row(i).t()).t() * (y_h-x_h*gibbs_beta.row(i).t());
      gibbs_tau0(i,k) = R::rgamma(a0*y_h.size()/2, 1/(a0*var(0,0)/2));
    }

  }

  arma::mat gibbs_beta_sub = gibbs_beta.rows(nBI, nMC+nBI-1);
  arma::vec gibbs_tau_sub = gibbs_tau.subvec(nBI, nMC+nBI-1);
  arma::mat gibbs_tau0_sub = gibbs_tau0.rows(nBI, nMC+nBI-1);

  Rcpp::List result;
  if(historical.size()==0){
    result = Rcpp::List::create(
      Rcpp::Named("posterior samples of beta")    = gibbs_beta_sub,
      Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub);
  }else{
    result = Rcpp::List::create(
      Rcpp::Named("posterior samples of beta")    = gibbs_beta_sub,
      Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub,
      Rcpp::Named("posterior samples of tau_0")    = gibbs_tau0_sub);
  }


  return(result);
}


// This function performs sample size determination for all four data types using the power prior with fixed a_0.
// This function calls glm_fixed_a0() and glm_fixed_a0_normal() to obtain posterior samples of beta. 
// [[Rcpp::export]]
Rcpp::List power_glm_fixed_a0(std::string & dType0, std::string & dLink0, double & n_total, arma::vec & n0,
                              Rcpp::List & historical0, arma::mat & x_samps,
                              arma::mat & beta_c_prior_samps, arma::vec & var_prior_samps,
                          arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                          double & delta, double & gamma,
                          int nMC, int nBI, int N, bool & dCurrent0) {

  NumericVector power(N);


  arma::mat beta_mean(N, beta_c_prior_samps.n_cols);
  arma::vec tau_samples(N);
  arma::mat tau0_samples(N, historical0.size());

  for (int i=0;i<N;i++){

    // data simulation for the control group
    int beta_ind = floor(Rcpp::runif(1, 0, beta_c_prior_samps.n_rows)(0));
    arma::rowvec beta_s = beta_c_prior_samps.row(beta_ind);


    // bootstrap x from historial data
    arma::mat x_prior_samps;

    if(historical0.size()==0){
      x_prior_samps = x_samps;
    }

    for(int j = 0; j < historical0.size(); j++){
      Rcpp::List dat = historical0[j];
      arma::mat x_h = dat["x0"];
      x_prior_samps = join_cols(x_prior_samps, x_h);
    }

    NumericVector ind = floor(Rcpp::runif(n_total,0,x_prior_samps.n_rows));
    arma::uvec v = Rcpp::as<arma::uvec>(ind);
    arma::mat x_s = x_prior_samps.rows(v);
    arma::vec vect(x_s.n_rows);
    vect.ones();
    x_s.insert_cols(0, vect);
    arma::vec x_trt = Rcpp::rbinom(n_total, 1, 0.5);
    x_s.insert_cols(1, x_trt);



    arma::vec mean = x_s*beta_s.t();


    if (dLink0=="Logistic") 		    { mean = exp(mean) / (1 + exp(mean)); 					}
    if (dLink0=="Probit")
      { for (int k=0; (unsigned)k< x_s.n_rows;k++)   { mean[k] = R::pnorm(mean[k], 0.0, 1.0, true, false);  } }
    if (dLink0=="Log") 	   		    { mean = exp(mean); 				    				}
    if (dLink0=="Identity-Positive")
      { 	for (int k=0; (unsigned)k< x_s.n_rows;k++) { mean[k] = std::max(mean[k],0.00001); }						}
    if (dLink0=="Identity-Probability")
      { 	for (int k=0; (unsigned)k< x_s.n_rows;k++) { mean[k] = std::min(std::max(mean[k],0.00001),0.99999); }	}
    if (dLink0=="Complementary Log-Log")    {  mean = 1 - exp(-exp(mean)); }


    // simulate y's
    arma::vec y_s(n_total);
    arma::vec n0(n_total);

    if (dType0=="Bernoulli")	{ for(int k = 0; k < n_total; k++){ y_s[k] = R::rbinom(1, mean(k)); }}
    if (dType0=="Binomial")	{ for(int k = 0; k < n_total; k++){ y_s[k] = R::rbinom(n0(k), mean(k)); }}
    if (dType0=="Poisson")	{ for(int k = 0; k < n_total; k++){ y_s[k] = R::rpois(mean(k)); }}
    if (dType0=="Exponential")	{ for(int k = 0; k < n_total; k++){ y_s[k] = R::rexp(1/mean(k)); }}
    if (dType0=="Normal")	{
      int ind_v = floor(R::runif(0,var_prior_samps.size()));
      double var = var_prior_samps[ind_v];
      for(int k = 0; k < n_total; k++){ y_s[k] = R::rnorm(mean(k),sqrt(var)); }
    }


    // need to remove intercept column when plugging into following functions
    arma::mat x_sub = x_s.cols(1, x_s.n_cols-1);
  

    arma::mat beta_samps(nMC,x_s.n_cols);
    if(dType0=="Normal")	{
      Rcpp::List lis = glm_fixed_a0_normal(y_s, x_sub, historical0, nMC, nBI);
      arma::mat a = lis[0];
      beta_samps = a;
      arma::vec tau_post = lis[1];
      arma::mat tau0_post = lis[2];
      tau_samples[i] = arma::mean(tau_post);

      arma::vec tau0_vec(tau0_post.n_cols);
      for(int k = 0; (unsigned)k < tau0_post.n_cols; k++){
        tau0_vec(k) = arma::mean(tau0_post.col(k));
      }

      tau0_samples.row(i) = tau0_vec.t();

    }else{
      
      beta_samps = glm_fixed_a0(dType0, dLink0, y_s, n0, x_sub, historical0,
                             lower_limits0,upper_limits0,slice_widths0,
                             nMC, nBI, dCurrent0);
    }

    for(int k = 0; (unsigned)k < beta_samps.n_cols; k++){
      beta_mean(i,k) = arma::mean(beta_samps.col(k));
    }



    arma::vec beta1 = beta_samps.col(1);
    arma::vec beta_sub = beta1.elem( find( beta1 < delta ) );
    double r = beta_sub.size()/double(beta1.size());

    power[i] = r;
  }

  double res = mean(power >= gamma);

  arma::vec mean_beta_vec(beta_mean.n_cols);
  for(int k = 0; (unsigned)k < beta_mean.n_cols; k++){
    mean_beta_vec(k) = arma::mean(beta_mean.col(k));
  }


  if(dType0 == "Normal"){
    double mean_tau = mean(tau_samples);

    // compute average posterior means
    arma::vec final_tau0_vec(tau0_samples.n_cols);
    for(int k = 0; (unsigned)k < tau0_samples.n_cols; k++){
      final_tau0_vec(k) = arma::mean(tau0_samples.col(k));
    }
    return Rcpp::List::create(
      Rcpp::Named("average posterior mean of beta") = mean_beta_vec,
      Rcpp::Named("average posterior mean of tau")  = mean_tau,
      Rcpp::Named("average posterior mean of tau0")  = final_tau0_vec,
      Rcpp::Named("power/type I error")    = res);
  }else{
    return Rcpp::List::create(
      Rcpp::Named("average posterior mean of beta") = mean_beta_vec,
      Rcpp::Named("power/type I error")    = res);
    }
}


// This function computes the mean and covariance matrix of beta based on asymptotics for normal responses. 
// This function is later called by power_glm_fixed_a0_approx().
Rcpp::List glm_fixed_a0_normal_approx(arma::vec & y, arma::mat & x, Rcpp::List & historical, double & delta) {
  
  arma::vec v(x.n_rows);
  v.ones();
  x.insert_cols(0, v);
  
  // compute beta_hat
  arma::mat xx(x.n_cols,x.n_cols);
  arma::mat a2xx(x.n_cols,x.n_cols);
  arma::vec xy(x.n_cols);
  xx.zeros();
  xy.zeros();
  a2xx.zeros();
  for(int j=0;j<historical.size();j++){
    Rcpp::List dat = historical[j];
    arma::vec y_h = dat["y0"];
    arma::mat x_h = dat["x0"];
    arma::vec v(x_h.n_rows);
    v.ones();
    x_h.insert_cols(0, v); //insert intercept
    arma::vec ins(x_h.n_rows);
    ins.zeros();
    x_h.insert_cols(1, ins); //insert treatment indicator
    double a0 = dat["a0"];
    xx = xx + a0 * x_h.t() * x_h;
    xy = xy + a0 * x_h.t() * y_h;
    a2xx = a2xx + a0 * a0 * x_h.t() * x_h;
  }
  arma::mat A = x.t()*x + xx;
  arma::vec B = x.t()*y + xy;
  arma::mat A_i = inv(A);
  
  arma::vec beta_hat = A_i*B;
  
  // compute tau_hat
  double ayx = 0;
  double an = 0;
  for(int j=0;j<historical.size();j++){
    Rcpp::List dat = historical[j];
    arma::vec y_h = dat["y0"];
    arma::mat x_h = dat["x0"];
    arma::vec v(x_h.n_rows);
    v.ones();
    x_h.insert_cols(0, v); //insert intercept
    arma::vec ins(x_h.n_rows);
    ins.zeros();
    x_h.insert_cols(1, ins); //insert treatment indicator
    double a0 = dat["a0"];
    arma::mat m = a0*(y_h-x_h*beta_hat).t() * (y_h-x_h*beta_hat);
    ayx = ayx + m(0,0);
    an = an + y_h.size();
    an = an + a0 * (y_h.size()- x.n_cols);
  }
  
  arma::mat yx = (y-x*beta_hat).t() * (y-x*beta_hat);
  double tau_hat_inv =  (yx(0,0) + ayx) / (y.size() - x.n_cols + an);
  
  arma::mat beta_cov = tau_hat_inv * A_i * (x.t()*x + a2xx) * A_i.t();
  
  double prob = R::pnorm(delta, beta_hat(1,0), sqrt(beta_cov(1,1)), TRUE, FALSE);
  
  return Rcpp::List::create(
    Rcpp::Named("coefficient") = beta_hat,
    Rcpp::Named("covariance") = beta_cov,
    Rcpp::Named("prob")=prob);
}



// This function computes the mean and covariance matrix of beta based on asymptotics using the Newton-Raphson algorithm 
// for Bernoulli, Poisson and Exponential responses.
// This function is later called by power_glm_fixed_a0_approx().
Rcpp::List newton_raphson(std::string & dType, arma::vec & y, arma::mat & x, Rcpp::List & historical, double & delta, int & n, double & tol) {
  
  arma::vec v(x.n_rows);
  v.ones();
  x.insert_cols(0, v);
  
  // extract regression parameters;
  arma::vec mu, var, dmu;
  arma::vec beta(x.n_cols);
  beta.zeros();
  arma::vec u(x.n_cols);
  arma::mat H(x.n_cols, x.n_cols);
  arma::mat H_i(x.n_cols, x.n_cols);
  
  for(int k = 0; k < n; k++){
    
    u.zeros();
    H.zeros();
    
    arma::vec eta = x*beta;

    
    if (dType=="Bernoulli") { mu=exp(eta)/(1+exp(eta)); var=mu % (1-mu); dmu=exp(eta)/pow((1+exp(eta)),2); }
    if (dType=="Poisson") 	{ mu=exp(eta); var=mu; dmu=exp(eta); 	}
    if (dType=="Exponential")	{ mu=1/eta; var=mu%mu; dmu=-1/(eta%eta); }
    
    
    u += ((((y-mu)/var)%dmu).t() * x).t();
    
    
    arma::mat x1(x.n_rows, x.n_cols);
    for(int j = 0; j < x.n_cols; j++){
      x1.col(j) = x.col(j) % (dmu%dmu/var);
    }
    H += x1.t() * x;

    // compute historical data likelihood
    for(int i = 0; i < historical.size(); i++){
      Rcpp::List dat = historical[i];
      arma::vec y_h = dat["y0"];
      arma::mat x_h = dat["x0"];
      double a0 = dat["a0"];
      //add intercept
      arma::vec v(x_h.n_rows);
      v.ones();
      x_h.insert_cols(0, v);
      // add treatment indicator
      arma::vec v2(x_h.n_rows);
      v2.zeros();
      x_h.insert_cols(1, v2);
      
      arma:: vec eta_h = x_h*beta;
      
      if (dType=="Bernoulli") { mu=exp(eta_h)/(1+exp(eta_h)); var=mu%(1-mu); dmu=exp(eta_h)/pow((1+exp(eta_h)),2); }
      if (dType=="Poisson") 	{ mu=exp(eta_h); var=mu; dmu=exp(eta_h); 	}
      if (dType=="Exponential")	{ mu=1/eta_h; var=mu%mu; dmu=-1/(eta_h%eta_h); }
      
      u += a0 * ((((y_h-mu)/var)%dmu).t() * x_h).t();
      
      arma::mat x2(x_h.n_rows, x_h.n_cols);
      for(int j = 0; j < x_h.n_cols; j++){
        x2.col(j) = x_h.col(j) % (dmu%dmu/var);
      }
      
      
      H += a0 * x2.t() * x_h;

    }

    H_i = inv(H);
    arma::vec diff = H_i*u;

    if (max(abs(diff)) < tol) {
      k=n;
    }else{
      beta += diff;
    }
    
  }
  
  double prob = R::pnorm(delta, beta(1,0), sqrt(H_i(1,1)), TRUE, FALSE);
  
  return Rcpp::List::create(
    Rcpp::Named("coefficient") = beta.col(0),
    Rcpp::Named("covariance") = H_i,
    Rcpp::Named("prob")=prob);
}



// This function performs approximate sample size determination for all four data types using the power prior with fixed a_0.
// [[Rcpp::export]]
double power_glm_fixed_a0_approx(std::string & dType0, double & n_total,
                              Rcpp::List & historical0, arma::mat & x_samps,
                              arma::mat & beta_c_prior_samps, arma::vec & var_prior_samps,
                              double & delta, double & gamma,
                              int nNR, double tol, int N) {

  NumericVector power(N);


  for (int i=0;i<N;i++){
    
    // data simulation for the control group
    int beta_ind = floor(Rcpp::runif(1, 0, beta_c_prior_samps.n_rows)(0));
    arma::rowvec beta_s = beta_c_prior_samps.row(beta_ind);

    // bootstrap x from historical data
    arma::mat x_prior_samps;

    if(historical0.size()==0){
      x_prior_samps = x_samps;
    }

    for(int j = 0; j < historical0.size(); j++){
      Rcpp::List dat = historical0[j];
      arma::mat x_h = dat["x0"];
      x_prior_samps = join_cols(x_prior_samps, x_h);
    }

    NumericVector ind = floor(Rcpp::runif(n_total,0,x_prior_samps.n_rows));
    arma::uvec v = Rcpp::as<arma::uvec>(ind);
    arma::mat x_s = x_prior_samps.rows(v);
    arma::vec vect(x_s.n_rows);
    vect.ones();
    x_s.insert_cols(0, vect);
    arma::vec x_trt = Rcpp::rbinom(n_total, 1, 0.7);
    x_s.insert_cols(1, x_trt);

    arma::vec mean = x_s*beta_s.t();
    

    // simulate y's
    arma::vec y_s(n_total);

    if (dType0=="Bernoulli")	{
      mean = exp(mean) / (1 + exp(mean));
      for(int k = 0; k < n_total; k++){ y_s[k] = R::rbinom(1, mean(k)); }
    }
    if (dType0=="Poisson")	{
      mean = exp(mean);
      for(int k = 0; k < n_total; k++){ y_s[k] = R::rpois(mean(k)); }
    }
    if (dType0=="Exponential")	{
      for(int k = 0; k < n_total; k++){ y_s[k] = R::rexp(1/mean(k)); }
    }
    if (dType0=="Normal")	{
      int ind_v = floor(R::runif(0,var_prior_samps.size()));
      double var = var_prior_samps[ind_v];
      for(int k = 0; k < n_total; k++){ y_s[k] = R::rnorm(mean(k),sqrt(var)); }
    }

    // need to remove intercept column when plugging into following functions
    arma::mat x_sub = x_s.cols(1, x_s.n_cols-1);


    double prob = 0;
    if(dType0 == "Normal"){
      prob = glm_fixed_a0_normal_approx(y_s, x_sub, historical0, delta)[2];
    }else{
      prob = newton_raphson(dType0, y_s, x_sub, historical0, delta, nNR, tol)[2];
    }
    power[i] = prob;
  }
  
  double res = mean(power >= gamma);

  return res;

}

  

