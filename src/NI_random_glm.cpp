// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>





using namespace Rcpp;

//------------------------------------ Class Specific Functions -------------------------------------------//
class random_a0_glm{

public:

  // data
  arma::vec y;  // outcome data
  arma::vec n;  // number of trials (for binomial data)
  arma::mat x;  // covariates
  Rcpp::List historical; // historical data



  // model definition
  std::string dType;
  std::string dLink;

  int P;


  arma::vec coef;

  // beta prior for a_0 hyperparameters;
  double c_1;
  double c_2;


  arma::vec                lower_limits;
  arma::vec                upper_limits;
  arma::vec                slice_widths;
  int m;



  // public member functions;
  random_a0_glm(std::string & dType0, std::string & dLink0, arma::vec & y0, arma::vec & n0, arma::mat & x0, Rcpp::List & historical0,
           double & c_1,double & c_2,arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, arma::vec & coef0);

  double logFC(const arma::vec & parm0, const int & p);

};

random_a0_glm::random_a0_glm(	std::string & dType0, std::string & dLink0, arma::vec & y0, arma::vec & n0, arma::mat & x0, Rcpp::List & historical0,
                    double & c_10,double & c_20, arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, arma::vec & coef0)
{

  dType = dType0;
  dLink = dLink0;

  y = y0;
  //n = n0;
  if ((dType=="Bernoulli")) {n.resize(y.size()); n.ones();} else {n = n0;}
  x     = x0;
  historical = historical0;



  P = x.n_cols;

  coef = coef0;

  c_1 = c_10;
  c_2 = c_20;


  lower_limits = lower_limits0;
  upper_limits = upper_limits0;
  slice_widths = slice_widths0;


  m=10;
}

// beta has P elements; a0 has historical.size elements
double random_a0_glm::logFC(const arma::vec & parm0, const int & p)
{
  // extract regression parameters;
  arma::vec beta0 = parm0.subvec(0,P-1);
  arma::vec mean;
  arma::vec a0_vec = parm0.subvec(P, parm0.size()-1);

  //arma::vec beta_h = beta0;
  arma::uvec ind;
  ind << 0;
  arma::vec beta_h = join_cols(parm0.elem(ind), parm0.subvec(2,P-1));



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


  // compute likelihood contribution;
  double lp = 0;

  if ((dType=="Bernoulli") | (dType=="Binomial"))	{lp += sum( y % log(mean) + (n-y) % log(1-mean) );}
  if (dType=="Poisson")	{ lp += sum(  y % log(mean) - mean );}
  if (dType=="Exponential")	{ lp += sum(  log(mean) - y % mean );}



  for(int i = 0; i < historical.size(); i++){
      Rcpp::List dat = historical[i];
      arma::vec y_h = dat["y0"];
      arma::mat x_h = dat["x0"];
      arma::vec v(x_h.n_rows);
      v.ones();
      x_h.insert_cols(0, v);

      double a0 = a0_vec[i];
      arma::vec n_h;
      n_h.zeros();
      if (dType=="Bernoulli") {n_h.resize(y_h.size()); n_h.ones();} 
      if (dType=="Binomial") {n_h = Rcpp::as<arma::vec>(dat["n0"]);}


      arma:: vec mean_h = exp(x_h*beta_h)/(1 + exp(x_h*beta_h));

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


  // sampling a0
  if(p > P-1){

    double degree = (coef.size()-1) / historical.size();
    arma::vec vec(coef.size());

    vec[0] = 1;
    for(int k = 0;k < degree; k++){
      vec.subvec(k*historical.size()+1, historical.size()*(k+1)) = pow(a0_vec, k+1);
    }

    double f_a0 = sum(vec % coef);

    lp += (c_1-1)*log(parm0[p]) + (c_2-1)*log(1-parm0[p]) - f_a0;

  }

  return  lp;
}


void slice( arma::vec & parms, random_a0_glm & b)
{

  double b0, f0, f0_L, f0_R, f0_x1, h0, L, R, V, J, K,w,lower,upper;
  arma::vec parm0;

  for (int p = 0; p < b.P+b.historical.size(); p++)
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


// [[Rcpp::export]]
Rcpp::List glm_random_a0(std::string & dType0, std::string & dLink0, arma::vec & y0, arma::vec & n0, arma::mat & x0, Rcpp::List & historical0,
                   double & c_10,double & c_20, arma::vec & coef0, arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                   int nMC, int nBI)
{

  Rcpp::RNGScope scope;

  arma::vec v(x0.n_rows);
  v.ones();
  x0.insert_cols(0, v); //insert intercept


  int P = x0.n_cols;

  // declare regression object and set values;
  random_a0_glm b(dType0,dLink0,y0,n0,x0,historical0,c_10,c_20,lower_limits0,upper_limits0,slice_widths0,coef0);

  // sample betas and alphas;
  arma::mat samples(nMC, P+historical0.size());


  // create parameter vector container and initial values;
  arma::vec parms(P+historical0.size());
  for (int p=0;p<P+historical0.size();p++)
  {

    parms[p]= R::runif(0,1);

  }


  for (int s=-nBI;s<nMC;s++)
  {
    slice(parms,b);
    if (s>=0){	 samples.row(s) = parms.t(); }
  }

  return Rcpp::List::create(
    //Rcpp::Named("post_samples")    = post_samples,
    Rcpp::Named("posterior samples of beta")    = samples.cols(0,P-1),
    Rcpp::Named("posterior samples of a0")    = samples.cols(P,P+historical0.size()-1));

}





// [[Rcpp::export]]
Rcpp::List power_glm_random_a0(std::string & dType0, std::string & dLink0, double & n_total, arma::vec & n0, Rcpp::List & historical0,
                               arma::mat & beta_c_prior_samps, double & c_10, double & c_20, arma::vec & coef0,
                               arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                              double & delta, double & gamma,
                              int nMC, int nBI, int N) {

  NumericVector power(N);

  //arma::mat y_samp(N, n_total);
  //arma::mat x_samp;
  arma::mat beta_mean(N, beta_c_prior_samps.n_cols);
  arma::mat a0_mean(N, historical0.size());

  for (int i=0;i<N;i++){

    // data simulation for the control group
    int beta_ind = floor(Rcpp::runif(1, 0, beta_c_prior_samps.n_rows)(0));
    arma::rowvec beta_s = beta_c_prior_samps.row(beta_ind);


    // bootstrap x from historial data
    arma::mat x_prior_samps;
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



    //y_samp.row(i) = y_s.t();
    //x_samp = join_cols(x_samp, x_s);

    // need to remove intercept column when plugging into following functions
    arma::mat x_sub = x_s.cols(1, x_s.n_cols-1);

    arma::mat beta_samps(nMC,x_s.n_cols);

    Rcpp::List lis = glm_random_a0(dType0, dLink0, y_s, n0, x_sub, historical0, c_10, c_20, coef0,
                                     lower_limits0, upper_limits0, slice_widths0,
                                     nMC, nBI);

    arma::mat a = lis[0];
    beta_samps = a;
    arma::mat a0_mat = lis[1];



    for(int k = 0; (unsigned)k < beta_samps.n_cols; k++){
      beta_mean(i,k) = arma::mean(beta_samps.col(k));
    }

    arma::vec beta1 = beta_samps.col(1);
    arma::vec beta_sub = beta1.elem( find( beta1 < delta ) );
    double r = beta_sub.size()/double(beta1.size());

    power[i] = r;

    arma::vec mean_a0_vec(a0_mat.n_cols);

    for(int k = 0; (unsigned)k < a0_mat.n_cols; k++){
      mean_a0_vec(k) = arma::mean(a0_mat.col(k));
    }

    a0_mean.row(i) = mean_a0_vec.t();
  }

  double res = mean(power >= gamma);

  arma::vec mean_beta_vec(beta_mean.n_cols);
  for(int k = 0; (unsigned)k < beta_mean.n_cols; k++){
    mean_beta_vec(k) = arma::mean(beta_mean.col(k));
  }

  arma::vec final_a0_vec(a0_mean.n_cols);
  for(int k = 0; (unsigned)k < a0_mean.n_cols; k++){
    final_a0_vec(k) = arma::mean(a0_mean.col(k));
  }


  return Rcpp::List::create(
    //Rcpp::Named("post_samples")    = post_samples,
    Rcpp::Named("average posterior mean of beta") = mean_beta_vec,
    Rcpp::Named("average posterior means of a0")  = final_a0_vec,
    Rcpp::Named("power/type I error")    = res
  );
}

