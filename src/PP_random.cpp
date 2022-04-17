// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

///////////////////////////////////////////////////////////////////////////////
  ///////////////////// --------  two group case with random a0  --------- ////////////////
  class random_a0 {
    public:

      // data containers;
    arma::vec y_normal; // for normal data with covariates
    arma::mat x_normal; // for normal case with covariates
    double y_c;
    double n_c;
    double v_c;

    arma::mat historical;
    Rcpp::List historical_normal; // for normal data with covariates

    // model definition
    std::string dType;
    // parameter containers;
    int P;
    // intial prior hyperparameters (not needed for normal case);
    double b_01;
    double b_02;
    // beta prior for a_0 hyperparameters;
    arma::vec c_1;
    arma::vec c_2;

    // slice sampling containers;
    arma::vec lower_limits;
    arma::vec upper_limits;
    arma::vec slice_widths;
    int m;



    // public member functions;
    random_a0(std::string dType0, double y0, double n0, double v0, arma::vec y_normal0,
              arma::mat x_normal0, arma::mat historical0,
              Rcpp::List historical_normal0,
              double b_010, double b_020, arma::vec c_10, arma::vec c_20,
              arma::vec & lower_limits0, arma::vec & upper_limits0,
              arma::vec & slice_widths0);


    double logFC(const arma::vec & parm0, const int & p);
  };



// define constructor
random_a0::random_a0(std::string dType0, double y0, double n0, double v0,
                     arma::vec y_normal0, arma::mat x_normal0, arma::mat historical0,
                     Rcpp::List historical_normal0,
                     double b_010, double b_020, arma::vec c_10, arma::vec c_20,
                     arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0)
{

  y_c = y0;

  n_c = n0;
  v_c = v0;

  y_normal     = y_normal0;
  x_normal     = x_normal0;

  dType = dType0;



  if((dType=="Normal") & (x_normal0.n_rows!=0)){
    historical_normal = historical_normal0;
    P = historical_normal0.size();
  }else{
    historical = historical0;
    P = historical.n_rows;
  }

  b_01 = b_010;
  b_02 = b_020;
  c_1 = c_10;
  c_2 = c_20;


  lower_limits = lower_limits0;
  upper_limits = upper_limits0;
  slice_widths = slice_widths0;

  m=10;
}

// Define the log likelihood of the posterior using the normalized power prior
// for two group cases, and for the normal linear model.
double random_a0::logFC(const arma::vec & parm0, const int & p)
{
  double log_C_1 = 0;
  double log_C_0 = 0;
  if (dType=="Bernoulli")	{
    log_C_0 = R::lbeta(arma::dot(historical.col(0), parm0) + b_01,
                       arma::dot((historical.col(1) - historical.col(0)), parm0) + b_02);
    log_C_1 = R::lbeta(y_c + arma::dot(historical.col(0), parm0) + b_01,
                       n_c - y_c + arma::dot((historical.col(1) -
                                                historical.col(0)), parm0) + b_02);
  }
  if (dType=="Poisson")	{
    log_C_0 = lgamma(arma::dot(historical.col(0), parm0) + b_01) -
      (arma::dot(historical.col(0), parm0) + b_01) * log(arma::dot(historical.col(1), parm0) + b_02);


    log_C_1 = lgamma(y_c + arma::dot(historical.col(0), parm0) + b_01) -
      (y_c + arma::dot(historical.col(0), parm0) + b_01) * log(n_c + arma::dot(historical.col(1), parm0) + b_02);

  }
  if (dType=="Exponential"){
    log_C_0 = lgamma(arma::dot(historical.col(1), parm0) + b_01) -
      (arma::dot(historical.col(1), parm0) + b_01) * log(arma::dot(historical.col(0), parm0) + b_02);

    log_C_1 = lgamma(n_c + arma::dot(historical.col(1), parm0) + b_01) -
      (n_c + arma::dot(historical.col(1), parm0) + b_01) * log(y_c + arma::dot(historical.col(0), parm0) + b_02);

  }
  if (dType=="Normal"){
    if (x_normal.n_rows==0){ // if there are no covariates
      double a_n = 0;
      double a_y = 0;
      double a_y2 = 0;
      for(int j=0;(unsigned)j<historical.n_rows;j++){
        double y_h = historical(j,0);
        double n_h = historical(j,1);
        double v_h = historical(j,2);
        a_n = a_n + n_h*parm0[j];
        double y2sum_h = (n_h-1)*v_h + n_h*pow(y_h/n_h,2);
        a_y2 = a_y2 + parm0[j]*y2sum_h;
        a_y = a_y + parm0[j]*y_h;
      }

      double log_c_a0 = lgamma(0.5*(a_n-1))-0.5*log(a_n)+0.5*(a_n-1)*log(2.0)-
        0.5*(a_n-1)*log(a_y2-pow(a_y,2)/a_n);

      double y_h = historical(p,0);
      double n_h = historical(p,1);
      double v_h = historical(p,2);
      double y2sum_h = (n_h-1)*v_h + n_h*pow(y_h/n_h,2);
      double ss_h = y2sum_h - 2*parm0[P]*y_h + n_h*pow(parm0[P],2);

      log_C_0 = parm0[p] * parm0[P+1] / 2 * ss_h + log_c_a0;
      log_C_1 = parm0[p] * n_h / 2 * log(parm0[P+1]);

    }else{ // normal linear model


      double a_y2 = 0;
      double a_n = 0;
      arma::mat a_xy(x_normal.n_cols-1,1);
      a_xy.zeros();
      arma::mat a_x2(x_normal.n_cols-1,x_normal.n_cols-1);
      a_x2.zeros();
      for(int j=0;j<historical_normal.size();j++){
        Rcpp::List dat = historical_normal[j];
        arma::vec y_j = dat["y0"];
        arma::mat x_j = dat["x0"];
        arma::vec v(x_j.n_rows);
        v.ones();
        x_j.insert_cols(0, v);

        a_n = a_n + y_j.size()*parm0[j];
        a_y2 = a_y2 + parm0[j]*sum(y_j%y_j);
        a_xy = a_xy + parm0[j]*x_j.t()*y_j;
        a_x2 = a_x2 + parm0[j]*x_j.t()*x_j;
      }
      arma::mat amat = a_xy.t()*inv(a_x2)*a_xy;
      double val;
      double sign;

      log_det(val, sign, a_x2);
      double log_c_a0 = -0.5*log(exp(val)*sign) + lgamma(0.5*(a_n-(x_normal.n_cols-1)))-
        0.5*(a_n-(x_normal.n_cols-1))*log(0.5*(a_y2-amat(0,0)));

      Rcpp::List dat = historical_normal[p];
      arma::vec y_h = dat["y0"];
      arma::mat x_h = dat["x0"];
      arma::vec v(x_h.n_rows);
      v.ones();
      x_h.insert_cols(0, v);
      arma::vec ins(x_h.n_rows);
      ins.zeros();
      x_h.insert_cols(1, ins); // x_h does not have treatment indicator
      log_C_0 = parm0[p] * parm0[P+x_normal.n_cols] / 2 * sum((y_h - x_h*parm0.subvec(P,P+x_normal.n_cols-1))%(y_h - x_h*parm0.subvec(P,P+x_normal.n_cols-1))) + log_c_a0;
      log_C_1 = parm0[p] * y_h.size() / 2 * log(parm0[P+x_normal.n_cols]);
    }


  }

  return log_C_1 - log_C_0 + (c_1[p]-1)*log(parm0[p]) + (c_2[p]-1)*log(1-parm0[p]);

}

////////////////////////////////////////////////////////////////////////////////
  ////////////////// --------------  Slice Sampler  --------------- //////////////


  void slice( arma::vec & parms, random_a0 & b)
{

  double b0, f0, f0_L, f0_R, f0_x1, h0, L, R, V, J, K, w, lower, upper;
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




// This function generates the parameters of the posterior distribution of a_0 (using slice sampling) 
// and mu_c (the distribution of mu_c given a_0 has closed form solutions) for Bernoulli, Poisson and Exponential responses.
// [[Rcpp::export]]
Rcpp::List two_grp_random_a0(std::string & dType0, double & y0, double & n0, arma::mat & historical0,
                            double & b_010, double & b_020, arma::vec & c_10, arma::vec & c_20,
                            arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                            int nMC, int nBI)
{

  Rcpp::RNGScope scope;


  Rcpp::List historical_normal0;
  arma::vec y_normal0(1);
  arma::mat x_normal0(1,1);
  x_normal0.zeros();
  random_a0 b(dType0, y0, n0, 0.1, y_normal0, x_normal0, historical0, historical_normal0, b_010, b_020, c_10, c_20, lower_limits0, upper_limits0, slice_widths0);



  int P = historical0.n_rows;

  arma::mat a0_samples(nMC,P);

  arma::vec parms(P);
  for (int p=0;p < P;p++)
  {
    parms[p]= R::runif(0,1);
  }

  for (int s=-nBI;s<nMC;s++)
  {
    // call the slice sampler for samples of a_0
    slice(parms,b);

    if (s>=0){	a0_samples.row(s) = parms.t();	}
  }

  
  // sample mu_c given samples of a_0
  arma::vec mu_c_post(nMC);

  if(dType0=="Bernoulli"){
    for (int j=0;(unsigned)j < a0_samples.n_rows;j++){
      double alpha_c = (double) y0 + arma::dot(historical0.col(0), a0_samples.row(j)) + b_010;
      double beta_c = (double) n0 - (double) y0 +  arma::dot((historical0.col(1) - historical0.col(0)), a0_samples.row(j)) + b_020;
      double mu_c_samp = R::rbeta(alpha_c, beta_c);
      mu_c_post[j] = mu_c_samp;
    }
  }

  if(dType0=="Poisson"){
    for (int j=0;(unsigned)j < a0_samples.n_rows;j++){
      double alpha_c = (double) y0 + arma::dot(historical0.col(0), a0_samples.row(j)) + b_010;
      double beta_c = (double) n0 + arma::dot(historical0.col(1), a0_samples.row(j)) + b_020;
      double mu_c_samp = R::rgamma(alpha_c, 1/beta_c);
      mu_c_post[j] = mu_c_samp;
    }
  }

  if(dType0=="Exponential"){
    for (int j=0;(unsigned)j < a0_samples.n_rows;j++){
      double alpha_c = (double) n0 + arma::dot(historical0.col(1), a0_samples.row(j)) + b_010;
      double beta_c = (double) y0 + arma::dot(historical0.col(0), a0_samples.row(j)) + b_020;
      double mu_c_samp = R::rgamma(alpha_c, 1/beta_c);
      mu_c_post[j] = mu_c_samp;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("posterior samples of mu_c")    = mu_c_post,
    Rcpp::Named("posterior samples of a0")    = a0_samples);
}



// This function generates posterior samples of mu_c, tau, and a_0 
// using Gibbs and slice sampling for normal responses.
// [[Rcpp::export]]
Rcpp::List two_grp_random_a0_normal(double y0, double n0, double v0,
                                    arma::mat historical0,
                                    arma::vec & c_10, arma::vec & c_20,
                                    arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                                    int nMC, int nBI) {
  Rcpp::RNGScope scope;



  arma::vec gibbs_mu_c(nMC+nBI);
  gibbs_mu_c.zeros();
  arma::vec gibbs_tau(nMC+nBI);
  gibbs_tau.ones();
  arma::mat a0_samples(nMC+nBI,historical0.n_rows);


  for (int i=1;i<nMC+nBI;i++){


    double num = 0;
    double denom = 0;
    for(int j=0;(unsigned)j<historical0.n_rows;j++){
      double a0 = a0_samples(i-1,j);
      num = num + historical0(j,0) * a0;
      denom = denom + historical0(j,1) * a0;
    }

    //sampling mu

    double mu = (y0 + num)/(n0 + denom);
    double var = 1/(n0*gibbs_tau[i-1] + gibbs_tau[i-1]*denom);
    gibbs_mu_c[i] = R::rnorm(mu, sqrt(var));

    //sampling tau

    double a = 0;
    for(int k=0;(unsigned)k<historical0.n_rows;k++){
      double a0 = a0_samples(i-1,k);
      double y_h = historical0(k,0);
      double n_h = historical0(k,1);
      double v_h = historical0(k,2);
      double y2sum_h = (n_h-1)*v_h + n_h*pow(y_h/n_h,2);
      double ss_h = y2sum_h - 2*gibbs_mu_c[i]*y_h + n_h*pow(gibbs_mu_c[i],2);

      a = a + a0*ss_h;
    }
    double y2sum = (n0-1)*v0 + n0*pow(y0/n0,2);
    double ss = y2sum - 2*gibbs_mu_c[i]*y0 + n0*pow(gibbs_mu_c[i],2);

    gibbs_tau[i] = R::rgamma((n0+denom)/2, 1/((ss + a)/2));


    //sampling a0

    std::string dType0 = "Normal";
    Rcpp::List historical_normal0;
    arma::vec y_normal0(1);
    arma::mat x_normal0(0,1);
    x_normal0.zeros();
    random_a0 b(dType0, y0, n0, v0, y_normal0, x_normal0, historical0,
                historical_normal0, 0.1, 0.1, c_10, c_20,
                lower_limits0, upper_limits0, slice_widths0);

    int P = historical0.n_rows;

    arma::vec parms(P+2);
    for (int p=0;p < P;p++)
    {
      parms[p]= R::runif(0,1);
    }
    parms[P] = gibbs_mu_c[i];
    parms[P+1] = gibbs_tau[i];



    slice(parms,b);


    a0_samples.row(i) = parms.subvec(0,P-1).t();

  }
  arma::vec gibbs_mu_c_sub = gibbs_mu_c.subvec(nBI, nMC+nBI-1);
  arma::vec gibbs_tau_sub = gibbs_tau.subvec(nBI, nMC+nBI-1);
  arma::mat a0_samples_sub = a0_samples.rows(nBI, nMC+nBI-1);

  return Rcpp::List::create(
    Rcpp::Named("posterior samples of mu_c")    = gibbs_mu_c_sub,
    Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub,
    Rcpp::Named("posterior samples of a0")    = a0_samples_sub);
}



// This function generates posterior samples of beta, tau, and a_0 
// using Gibbs and slice sampling for the normal linear model.
// The normalized power prior is used. 
// [[Rcpp::export]]
Rcpp::List glm_random_a0_normal(arma::vec y_normal0, arma::mat x_normal0,
                                Rcpp::List historical_normal0,
                                arma::vec & c_10, arma::vec & c_20,
                                arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                                int nMC, int nBI) {
  Rcpp::RNGScope scope;

  arma::vec v(x_normal0.n_rows);
  v.ones();
  x_normal0.insert_cols(0, v);

  arma::mat gibbs_beta(nMC+nBI, x_normal0.n_cols);
  gibbs_beta.zeros();
  arma::vec gibbs_tau(nMC+nBI);
  gibbs_tau.ones();
  int P = historical_normal0.size();
  arma::mat a0_samples(nMC+nBI,P);


  for (int i=1;i<nMC+nBI;i++){

    //sampling beta
    arma::mat num(x_normal0.n_cols,x_normal0.n_cols);
    arma::vec denom(x_normal0.n_cols);
    num.zeros();
    denom.zeros();
    for(int j=0;j<historical_normal0.size();j++){
      Rcpp::List dat = historical_normal0[j];
      arma::vec y_h = dat["y0"];
      arma::mat x_h = dat["x0"];
      arma::vec v(x_h.n_rows);
      v.ones();
      x_h.insert_cols(0, v);
      arma::vec ins(x_h.n_rows);
      ins.zeros();
      x_h.insert_cols(1, ins);
      double a0 = a0_samples(i-1,j);
      num = num + a0 * x_h.t() * x_h;
      denom = denom + a0 * x_h.t() * y_h;
    }
    arma::mat A = gibbs_tau[i-1]*(x_normal0.t()*x_normal0 + num);
    arma::vec B = gibbs_tau[i-1]*(x_normal0.t()*y_normal0 + denom);
    arma::mat A_i = inv(A);
    arma::mat m = arma::randn(1, x_normal0.n_cols);

    gibbs_beta.row(i) = arma::repmat(A_i*B, 1, 1).t() + m * arma::chol(A_i);


    //sampling tau

    double num1 = 0;
    double a_n = 0;
    for(int k=0;k<historical_normal0.size();k++){
      Rcpp::List dat = historical_normal0[k];
      arma::vec y_h = dat["y0"];
      arma::mat x_h = dat["x0"];
      arma::vec v(x_h.n_rows);
      v.ones();
      x_h.insert_cols(0, v);
      arma::vec ins(x_h.n_rows);
      ins.zeros();
      x_h.insert_cols(1, ins);
      double a0 = a0_samples(i-1,k);
      arma::mat m2 = (y_h-x_h*gibbs_beta.row(i).t()).t() * (y_h-x_h*gibbs_beta.row(i).t());
      num1 = num1 + a0*m2(0,0);
      a_n = a_n + a0*y_h.size();
    }

    double alpha = (y_normal0.size()+a_n)/2;
    arma::mat m1 = (y_normal0-x_normal0*gibbs_beta.row(i).t()).t() * (y_normal0-x_normal0*gibbs_beta.row(i).t());
    double beta = (m1(0,0) + num1)/2;

    gibbs_tau[i] = R::rgamma(alpha, 1/beta);



    //sampling a0

    std::string dType0 = "Normal";
    arma::mat historical0(1,1);
    historical0.zeros();
    random_a0 b(dType0, 0.1, 1, 0.1, y_normal0, x_normal0, historical0,
                historical_normal0, 0.1, 0.1, c_10, c_20,
                lower_limits0, upper_limits0, slice_widths0);


    arma::vec parms(P+x_normal0.n_cols+1);
    for (int p=0;p < P;p++)
    {
      parms[p]= R::runif(0,1);
    }
    parms.subvec(P,P+x_normal0.n_cols-1) = gibbs_beta.row(i).t();
    parms[P+x_normal0.n_cols] = gibbs_tau[i];


    slice(parms,b);


    a0_samples.row(i) = parms.subvec(0,P-1).t();

  }
  arma::mat gibbs_beta_sub = gibbs_beta.rows(nBI, nMC+nBI-1);
  arma::vec gibbs_tau_sub = gibbs_tau.subvec(nBI, nMC+nBI-1);
  arma::mat a0_samples_sub = a0_samples.rows(nBI, nMC+nBI-1);

  return Rcpp::List::create(
    Rcpp::Named("posterior samples of beta")    = gibbs_beta_sub,
    Rcpp::Named("posterior samples of tau")    = gibbs_tau_sub,
    Rcpp::Named("posterior samples of a0")    = a0_samples_sub);
}





// This function performs sample size determination for all four data types using the normalized power prior.
// This function calls two_grp_random_a0() and two_grp_random_a0_normal() to obtain posterior samples of a_0 and mu_c. 
// [[Rcpp::export]]
Rcpp::List power_two_grp_random_a0(std::string & dType0, double & n_t, double & n0, arma::mat & historical0,
                                   arma::vec & mu_t_prior_samps, arma::vec & mu_c_prior_samps, arma::vec & var_t_prior_samps, arma::vec & var_c_prior_samps,
                                   double & b_t1, double & b_t2, double & b_010, double & b_020, arma::vec & c_10, arma::vec & c_20,
                                   arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                                   double & delta, double & gamma,
                                   int & nMC, int & nBI, int & N) {



  NumericVector post_probs(N);
  NumericVector mu_t_samples(N);
  NumericVector mu_c_samples(N);
  NumericVector tau_samples(N);
  arma::mat a0_mean(N, historical0.n_rows);

  for (int i=0;i<N;i++){

    // data simulation
    int ind_t = floor(R::runif(0,mu_t_prior_samps.size()));
    int ind_c = floor(R::runif(0,mu_c_prior_samps.size()));
    double mu_t = mu_t_prior_samps[ind_t];
    double mu_c = mu_c_prior_samps[ind_c];

    arma::mat a0_mat;

    NumericVector mu_t_post(nMC);
    NumericVector mu_c_post(nMC);

    if(dType0 == "Bernoulli"){
      int y_t = R::rbinom(n_t, mu_t);
      int y_c = R::rbinom(n0, mu_c);

      double y_cd = (double) y_c;

      Rcpp::List lis = two_grp_random_a0(dType0, y_cd, n0, historical0, b_010, b_020, c_10, c_20,
                                 lower_limits0, upper_limits0, slice_widths0,nMC, nBI);


      arma::mat a = lis[1];
      a0_mat = a;
      NumericVector b = lis[0];
      mu_c_post = b;


      double alpha_t = (double) y_t + b_t1;
      double beta_t = (double) n_t - (double) y_t + b_t2;
      mu_t_post = Rcpp::rbeta(a0_mat.n_rows, alpha_t, beta_t);
      post_probs[i] = mean(mu_t_post - mu_c_post < delta);

    }
    if(dType0 == "Poisson"){
      int y_t = R::rpois(mu_t*n_t);
      int y_c = R::rpois(mu_c*n0);

      double y_cd = (double) y_c;

      Rcpp::List lis = two_grp_random_a0(dType0, y_cd, n0, historical0, b_010, b_020, c_10, c_20,
                                 lower_limits0, upper_limits0, slice_widths0,nMC, nBI);

      arma::mat a = lis[1];
      a0_mat = a;
      NumericVector b = lis[0];
      mu_c_post = b;

      double alpha_t = (double) y_t + b_t1;
      double beta_t = (double) n_t + b_t2;
      mu_t_post = Rcpp::rgamma(a0_mat.n_rows, alpha_t, 1/beta_t);
      post_probs[i] = mean(mu_t_post - mu_c_post < delta);


    }

    if(dType0 == "Exponential"){
      double y_t = R::rgamma(n_t, 1/mu_t); //distribution of sum of y is Gamma, parametrized by scale parameter
      double y_c = R::rgamma(n0, 1/mu_c);

      double y_cd = (double) y_c;

      Rcpp::List lis = two_grp_random_a0(dType0, y_cd, n0, historical0, b_010, b_020, c_10, c_20,
                                 lower_limits0, upper_limits0, slice_widths0,nMC, nBI);

      arma::mat a = lis[1];
      a0_mat = a;
      NumericVector b = lis[0];
      mu_c_post = b;

      double alpha_t = (double) n_t + b_t1;
      double beta_t = (double) y_t + b_t2;
      mu_t_post = Rcpp::rgamma(a0_mat.n_rows, alpha_t, 1/beta_t);
      post_probs[i] = mean(mu_t_post/mu_c_post < delta);

    }
    if(dType0 == "Normal"){
      int ind_t = floor(R::runif(0,var_t_prior_samps.size()));
      int ind_c = floor(R::runif(0,var_c_prior_samps.size()));
      double var_t = var_t_prior_samps[ind_t];
      double var_c = var_c_prior_samps[ind_c];

      NumericVector y_t = Rcpp::rnorm(n_t, mu_t, sqrt(var_t));
      arma::vec y_c = Rcpp::rnorm(n0, mu_c, sqrt(var_c));

      double s_c = sum(y_c);
      double v_c = var(y_c);

      Rcpp::List lis = two_grp_random_a0_normal(s_c, n0, v_c, historical0, c_10, c_20,
                                                lower_limits0, upper_limits0, slice_widths0,nMC, nBI);

      arma::mat a = lis[2];
      a0_mat = a;
      NumericVector b = lis[0];
      mu_c_post = b;
      mu_t_post = Rcpp::rt(nMC, n_t-1) * sd(y_t)/sqrt(n_t) + mean(y_t);
      NumericVector tau_post = lis[1];
      tau_samples[i] = mean(tau_post);


      post_probs[i] = mean(mu_t_post - mu_c_post < delta);

    }


    mu_c_samples[i] = mean(mu_c_post);
    mu_t_samples[i] = mean(mu_t_post);

    arma::vec mean_a0_vec(a0_mat.n_cols);

    for(int k = 0; (unsigned)k < a0_mat.n_cols; k++){
      mean_a0_vec(k) = arma::mean(a0_mat.col(k));
    }

    a0_mean.row(i) = mean_a0_vec.t();


  }

  arma::vec final_a0_vec(a0_mean.n_cols);
  for(int k = 0; (unsigned)k < a0_mean.n_cols; k++){
    final_a0_vec(k) = arma::mean(a0_mean.col(k));
  }

  double res = mean(post_probs >= gamma);
  double mean_mu_t = mean(mu_t_samples);
  double mean_mu_c = mean(mu_c_samples);
  if(dType0 == "Normal"){
    double mean_tau = mean(tau_samples);
    return Rcpp::List::create(
      Rcpp::Named("average posterior mean of mu_t")  = mean_mu_t,
      Rcpp::Named("average posterior mean of mu_c")  = mean_mu_c,
      Rcpp::Named("average posterior mean of tau")  = mean_tau,
      Rcpp::Named("average posterior means of a0")  = final_a0_vec,
      Rcpp::Named("power/type I error")    = res);
  }else{
      return Rcpp::List::create(
        Rcpp::Named("average posterior mean of mu_t")  = mean_mu_t,
        Rcpp::Named("average posterior mean of mu_c")  = mean_mu_c,
        Rcpp::Named("average posterior means of a0")  = final_a0_vec,
        Rcpp::Named("power/type I error")    = res);
      }


}



// This function performs sample size determination for the normal linear model using the normalized power prior.
// This function calls glm_random_a0_normal() to obtain posterior samples of a_0, tau and beta. 
// [[Rcpp::export]]
Rcpp::List power_glm_random_a0_normal(double & n_total, Rcpp::List & historical0,
                                      arma::mat & beta_c_prior_samps, arma::vec & var_prior_samps,
                                      arma::vec & c_10, arma::vec & c_20,
                                      arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0,
                                      double & delta, double & gamma,
                                      int nMC, int nBI, int N) {

  NumericVector power(N);


  arma::mat beta_mean(N, beta_c_prior_samps.n_cols);
  NumericVector tau_samples(N);
  arma::mat a0_mean(N, historical0.size());

  for (int i=0;i<N;i++){

    // data simulation for the control group
    int beta_ind = floor(Rcpp::runif(1, 0, beta_c_prior_samps.n_rows)(0));
    arma::rowvec beta_s = beta_c_prior_samps.row(beta_ind);


    // bootstrap x from historical data
    arma::mat x_prior_samps;
    for(int j = 0; j < historical0.size(); j++){
      Rcpp::List dat = historical0[j];
      arma::mat x_h = dat["x0"];
      x_prior_samps = join_cols(x_prior_samps, x_h);
    }

    //Rcpp::IntegerVector pool = Rcpp::seq(0, x_prior_samps.n_rows-1);
    //std::random_shuffle(pool.begin(), pool.end()););
    //Rcpp::IntegerVector subpool = pool[Rcpp::seq(0, n_total-1)];
    //arma::uvec v = Rcpp::as<arma::uvec>(subpool);
    NumericVector ind = floor(Rcpp::runif(n_total,0,x_prior_samps.n_rows));
    arma::uvec v = Rcpp::as<arma::uvec>(ind);
    arma::mat x_s = x_prior_samps.rows(v);
    arma::vec vect(x_s.n_rows);
    vect.ones();
    x_s.insert_cols(0, vect);
    arma::vec x_trt = Rcpp::rbinom(n_total, 1, 0.5);
    x_s.insert_cols(1, x_trt);


    arma::vec mean = x_s*beta_s.t();


    // simulate y's
    arma::vec y_s(n_total);


    int index = floor(R::runif(0,var_prior_samps.size()));
    double var = var_prior_samps[index];
    for(int k = 0; k < n_total; k++){ y_s[k] = R::rnorm(mean(k),sqrt(var)); }


    arma::mat beta_samps(nMC,x_s.n_cols);

    // need to remove intercept column when plugging into following functions
    arma::mat x_sub = x_s.cols(1, x_s.n_cols-1);

    Rcpp::List lis = glm_random_a0_normal(y_s, x_sub,
                                         historical0,
                                         c_10, c_20,
                                         lower_limits0, upper_limits0, slice_widths0,
                                         nMC, nBI);
    arma::mat a = lis[0];
    beta_samps = a;
    arma::vec tau_post = lis[1];
    arma::mat a0_mat = lis[2];


    for(int k = 0; (unsigned)k < beta_samps.n_cols; k++){
      beta_mean(i,k) = arma::mean(beta_samps.col(k));
    }

    arma::vec beta1 = beta_samps.col(1);
    arma::vec beta_sub = beta1.elem( find( beta1 < delta ) );
    double r = beta_sub.size()/double(beta1.size());

    power[i] = r;

    tau_samples[i] = arma::mean(tau_post);

    arma::vec mean_a0_vec(a0_mat.n_cols);

    for(int k = 0; (unsigned)k < a0_mat.n_cols; k++){
      mean_a0_vec(k) = arma::mean(a0_mat.col(k));
    }

    a0_mean.row(i) = mean_a0_vec.t();
  }

  double res = mean(power >= gamma);

  // calculate average posterior means
  arma::vec mean_beta_vec(beta_mean.n_cols);
  for(int k = 0; (unsigned)k < beta_mean.n_cols; k++){
    mean_beta_vec(k) = arma::mean(beta_mean.col(k));
  }

  arma::vec final_a0_vec(a0_mean.n_cols);
  for(int k = 0; (unsigned)k < a0_mean.n_cols; k++){
    final_a0_vec(k) = arma::mean(a0_mean.col(k));
  }

  double mean_tau = mean(tau_samples);

  return Rcpp::List::create(
    Rcpp::Named("average posterior mean of beta") = mean_beta_vec,
    Rcpp::Named("average posterior mean of tau")  = mean_tau,
    Rcpp::Named("average posterior means of a0")  = final_a0_vec,
    Rcpp::Named("power/type I error")    = res
  );

}
