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
  //n = n0;
  if ((dType=="Bernoulli")) {n.resize(y.size()); n.ones();} else {n = n0;}
  x = x0;
  //arma::vec v(x0.n_rows);
  //v.ones();
  //x.insert_cols(0, v);
  historical = historical0;


  dCurrent = dCurrent0;

  if (dCurrent0==TRUE){
    P = x.n_cols;

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


double glm::logFC(const arma::vec & parm0, const int & p)
{
  // extract regression parameters;
  arma::vec beta0 = parm0, mean;

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

    //if ((dType=="Bernoulli")) {n.resize(y.size()); n.ones();} else {n = n0;}


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
    if (dType=="Bernoulli") {n_h.resize(y_h.size()); n_h.ones();} else {n_h = Rcpp::as<arma::vec>(dat["n"]);}

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

  return  lp;
}

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
      //Rcpp::Named("post_samples")    = post_samples,
      Rcpp::Named("posterior samples of beta")    = gibbs_beta_sub,
      Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub);
  }else{
    result = Rcpp::List::create(
      //Rcpp::Named("post_samples")    = post_samples,
      Rcpp::Named("posterior samples of beta")    = gibbs_beta_sub,
      Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub,
      Rcpp::Named("posterior samples of tau_0")    = gibbs_tau0_sub);
  }


  return(result);
}





// [[Rcpp::export]]
Rcpp::List power_glm_fixed_a0(std::string & dType0, std::string & dLink0, int & n_total, arma::vec & n0,
                              Rcpp::List & historical0, arma::mat & x_samps,
                              arma::mat & beta_c_prior_samps, arma::vec & var_prior_samps,
                          arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                          double & delta, double & gamma,
                          int nMC, int nBI, int N, bool & dCurrent0) {

  NumericVector power(N);

  //arma::mat y_samp(N, n_total);
  //arma::mat x_samp;
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


    //y_samp.row(i) = y_s.t();
    //x_samp = join_cols(x_samp, x_s);

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
