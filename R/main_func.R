

#' Model fitting for two groups (treatment and control group, no covariates) with fixed a0 when outcome follows Normal distribution
#'
#' @description Model fitting using power priors for two groups (treatment and control group, no covariates) with fixed \eqn{a_0} when outcome follows Normal distribution
#'
#' @param y.c Sum of responses (assumed to follow Normal distribution) for the control group.
#' @param n.c Sample size of the control group.
#' @param v.c Variance of responses for the control group.
#' @param historical (Optional) matrix of historical dataset(s) with four columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the variance of responses for the control group.
#' \item The fourth column contains the discounting parameter value \eqn{a_0} (between 0 and 1).
#' }
#' Each row represents a historical dataset.
#' @param nMC Number of iterations (excluding burn-in samples) for the Gibbs sampler. The default is 10,000.
#' @param nBI Number of burn-in samples for the Gibbs sampler. The default is 250.
#'
#' @details The power prior is applied on the data of the control group only.
#' Therefore, only summaries of the responses of the control group need to be entered.
#'
#' The responses are assumed to follow \eqn{N(\mu_c, \tau^{-1})} where \eqn{\mu_c} is the mean of responses for the control group
#' and \eqn{\tau} is the precision parameter. Each historical dataset \eqn{D_{0k}} is assumed to have a different precision parameter \eqn{\tau_k}.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}, and the initial prior for \eqn{\tau_k} is \eqn{\tau_k^{-1}}. The initial prior for the \eqn{\mu_c} is the uniform improper prior.
#' Posterior samples are obtained through Gibbs sampling.
#'
#' @return Posterior samples of \eqn{\mu_c}, \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are returned.
#' @references Chen, Ming-Hui, et al. "Bayesian design of noninferiority trials for medical devices using historical data." Biometrics 67.3 (2011): 1163-1170.
#' @seealso \code{\link{power.two.grp.fixed.a0}}
#' @examples
#' y.c <- 200 # The responses are assumed to follow normal distribution
#' n.c <- 100
#' v.c <- 2
#'
#' # Simulate three historical datasets
#' historical <- matrix(0, ncol=4, nrow=3)
#' historical[1,] <- c(200, 100, 2, 0.3)
#' historical[2,] <- c(300, 100, 2, 0.5)
#' historical[3,] <- c(400, 100, 2, 0.7)
#'
#' set.seed(1)
#' result <- two.grp.fixed.a0(y.c, n.c, v.c, historical, nMC=10000, nBI=250)

#' @export
two.grp.fixed.a0 <- function(y.c, n.c, v.c, historical=matrix(0,1,4), nMC=10000, nBI=250) {

  return(two_grp_fixed_a0_normal(y.c, n.c, v.c, historical, nMC, nBI))

}



#' Power/type I error calculation for data with two groups (treatment and control group, no covariates) with fixed a0
#'
#' @description Power/type I error calculation for data with two groups (treatment and control group, no covariates) with fixed \eqn{a_0} using power priors
#' @param data.type Character string specifying the type of response. The options are "Normal", "Bernoulli", "Poisson" and "Exponential".
#' @param n.t Sample size of the treatment group for the simulated datasets.
#' @param n.c Sample size of the control group for the simulated datasets.
#' @param historical (Optional) matrix of historical dataset(s). If \code{data.type} is "Normal", \code{historical} is a matrix with four columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the variance of responses for the control group.
#' \item The fourth column contains the discounting parameter value \eqn{a_0} (between 0 and 1).
#' }
#' For all other data types, \code{historical} is a matrix with three columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the discounting parameter value \eqn{a_0} (between 0 and 1).
#' }
#' Each row represents a historical dataset.
#' @param samp.prior.mu.t Vector of possible values of \eqn{\mu_t} to sample (with replacement) from. The vector contains realizations from the sampling prior (e.g. normal distribution) for \eqn{\mu_t}.
#' @param samp.prior.mu.c Vector of possible values of \eqn{\mu_c} to sample (with replacement) from. The vector contains realizations from the sampling prior (e.g. normal distribution) for \eqn{\mu_c}.
#' @param samp.prior.var.t Vector of possible values of \eqn{\sigma^2_t} to sample (with replacement) from. Only applies if \code{data.type} is "Normal". The vector contains realizations from the sampling prior (e.g. inverse-gamma distribution) for \eqn{\sigma^2_t}.
#' @param samp.prior.var.c Vector of possible values of \eqn{\sigma^2_c} to sample (with replacement) from. Only applies if \code{data.type} is "Normal". The vector contains realizations from the sampling prior (e.g. inverse-gamma distribution) for \eqn{\sigma^2_c}
#' @param prior.mu.t.shape1 First hyperparameter of the initial prior for \eqn{\mu_t}. The default is 1. Does not apply if \code{data.type} is "Normal".
#' @param prior.mu.t.shape2 Second hyperparameter of the initial prior for \eqn{\mu_t}. The default is 1. Does not apply if \code{data.type} is "Normal".
#' @param prior.mu.c.shape1 First hyperparameter of the initial prior for \eqn{\mu_c}. The default is 1. Does not apply if \code{data.type} is "Normal".
#' @param prior.mu.c.shape2 Second hyperparameter of the initial prior for \eqn{\mu_c}. The default is 1. Does not apply if \code{data.type} is "Normal".
#' @inheritParams power.glm.fixed.a0
#'
#'
#' @details If \code{data.type} is "Bernoulli", "Poisson" or "Exponential", a single response from the treatment group is assumed to follow Bern(\eqn{\mu_t}), Pois(\eqn{\mu_t}) or Exp(rate=\eqn{\mu_t}), respectively,
#' where \eqn{\mu_t} is the mean of responses for the treatment group. If \code{data.type} is "Normal", a single response from the treatment group is assumed to follow \eqn{N(\mu_t, \tau^{-1})}
#' where \eqn{\tau} is the precision parameter.
#' The distributional assumptions for the control group data are analogous.
#'
#' \code{samp.prior.mu.t} and \code{samp.prior.mu.c} can be generated using the sampling priors (see example).
#'
#' If \code{data.type} is "Bernoulli", the initial prior for \eqn{\mu_t} is
#' beta(\code{prior.mu.t.shape1}, \code{prior.mu.t.shape2}).
#' If \code{data.type} is "Poisson", the initial prior for \eqn{\mu_t} is
#' Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' If \code{data.type} is "Exponential", the initial prior for \eqn{\mu_t} is
#' Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' The initial priors used for the control group data are analogous.
#'
#' If \code{data.type} is "Normal", each historical dataset \eqn{D_{0k}} is assumed to have a different precision parameter \eqn{\tau_k}.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}, and the initial prior for \eqn{\tau_k} is \eqn{\tau_k^{-1}}.
#' The initial prior for the \eqn{\mu_c} is the uniform improper prior.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#'
#' If \code{data.type} is "Normal", Gibbs sampling is used for model fitting. For all other data types,
#' numerical integration is used for modeling fitting.
#'
#' @return Power or type I error is returned, depending on the sampling prior used. If \code{data.type} is "Normal", average posterior means of \eqn{\mu_c}, \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are also returned.
#' @references Yixuan Qiu, Sreekumar Balan, Matt Beall, Mark Sauder, Naoaki Okazaki and Thomas Hahn (2019). RcppNumerical: 'Rcpp' Integration for Numerical Computing Libraries. R package version 0.4-0. https://CRAN.R-project.org/package=RcppNumerical
#'
#' Chen, Ming-Hui, et al. "Bayesian design of noninferiority trials for medical devices using historical data." Biometrics 67.3 (2011): 1163-1170.
#' @seealso \code{\link{two.grp.fixed.a0}}
#' @examples
#' data.type <- "Bernoulli"
#' n.t <- 100
#' n.c <- 100
#'
#' # Simulate three historical datasets
#' historical <- matrix(0, ncol=3, nrow=3)
#' historical[1,] <- c(70, 100, 0.3)
#' historical[2,] <- c(60, 100, 0.5)
#' historical[3,] <- c(50, 100, 0.7)
#'
#' # Generate sampling priors
#' set.seed(1)
#' b_st1 <- b_st2 <- 1
#' b_sc1 <- b_sc2 <- 1
#' samp.prior.mu.t <- rbeta(50000, b_st1, b_st2)
#' samp.prior.mu.c <- rbeta(50000, b_st1, b_st2)
#' # The null hypothesis here is H0: mu_t - mu_c >= 0. To calculate power,
#' # we can provide samples of mu.t and mu.c such that the mass of mu_t - mu_c < 0.
#' # To calculate type I error, we can provide samples of mu.t and mu.c such that
#' # the mass of mu_t - mu_c >= 0.
#' sub_ind <- which(samp.prior.mu.t < samp.prior.mu.c)
#' # Here, mass is put on the alternative region, so power is calculated.
#' samp.prior.mu.t <- samp.prior.mu.t[sub_ind]
#' samp.prior.mu.c <- samp.prior.mu.c[sub_ind]
#'
#' N <- 1000 # N should be larger in practice
#' result <- power.two.grp.fixed.a0(data.type=data.type, n.t=n.t, n.c=n.t, historical=historical,
#'                                  samp.prior.mu.t=samp.prior.mu.t, samp.prior.mu.c=samp.prior.mu.c,
#'                                  delta=0, N=N)
#'
#' @export
power.two.grp.fixed.a0 <- function(data.type, n.t, n.c, historical=matrix(0,1,4), samp.prior.mu.t,
                                   samp.prior.mu.c, samp.prior.var.t, samp.prior.var.c, prior.mu.t.shape1=1,
                                   prior.mu.t.shape2=1, prior.mu.c.shape1=1, prior.mu.c.shape2=1, delta=0,
                                   gamma=0.95, nMC=10000, nBI=250, N=10000) {

  if(data.type == "Normal"){
    return(power_two_grp_fixed_a0_normal(n.t, n.c, historical, samp.prior.mu.t, samp.prior.mu.c,
                                         samp.prior.var.t, samp.prior.var.c, delta, gamma, nMC, nBI, N))
  }else{
    return(power_two_grp_fixed_a0(data.type, n.t, n.c, historical, samp.prior.mu.t, samp.prior.mu.c,
                                  prior.mu.t.shape1, prior.mu.t.shape2, prior.mu.c.shape1,
                                  prior.mu.c.shape2, delta, gamma, N, Inf))
  }
}








#' Model fitting for generalized linear models with fixed a0
#'
#' @description Model fitting using power priors for generalized linear models with fixed \eqn{a_0}
#'
#' @param y Vector of responses.
#' @param x Matrix of covariates. The first column should be the treatment indicator with 1 indicating treatment group. The number of rows should equal the length of the response vector \code{y}.
#'
#' @inheritParams power.glm.fixed.a0
#'
#' @details If \code{data.type} is "Normal", the response \eqn{y_i} is assumed to follow \eqn{N(x_i'\beta, \tau^{-1})} where \eqn{x_i} is the vector of covariates for subject \eqn{i}.
#' Each historical dataset \eqn{D_{0k}} is assumed to have a different precision parameter \eqn{\tau_k}.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}, and the initial prior for \eqn{\tau_k} is \eqn{\tau_k^{-1}}.
#' The initial prior for \eqn{\beta} is the uniform improper prior. Posterior samples are obtained through Gibbs sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling. The initial prior for \eqn{\beta} is the uniform improper prior.
#' The default lower limits for the parameters are -100. The default upper limits
#' for the parameters are 100. The default slice widths for the parameters are 1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' @return If \code{data.type} is "Normal", posterior samples of \eqn{\beta}, \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are returned.
#' For all other data types, a matrix of posterior samples of \eqn{\beta} is returned. The first column contains posterior samples of the intercept.
#' The second column contains posterior samples of \eqn{\beta_1}, the parameter for the treatment indicator.
#' @references Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#' @seealso \code{\link{power.glm.fixed.a0}}
#' @examples
#' data.type <- "Bernoulli"
#' data.link <- "Logistic"
#'
#' # Simulate current data
#' set.seed(1)
#' p <- 3
#' n_total <- 100
#' y <- rbinom(n_total,size=1,prob=0.6)
#' # The first column of x is the treatment indicator.
#' x <- cbind(rbinom(n_total,size=1,prob=0.5),
#'            matrix(rnorm(p*n_total),ncol=p,nrow=n_total))
#'
#' # Simulate two historical datasets
#' # Note that x0 does not have the treatment indicator
#' historical <- list(list(y0=rbinom(n_total,size=1,prob=0.2),
#'                         x0=matrix(rnorm(p*n_total),ncol=p,nrow=n_total), a0=0.2),
#'                    list(y0=rbinom(n_total, size=1, prob=0.5),
#'                         x0=matrix(rnorm(p*n_total),ncol=p,nrow=n_total), a0=0.3))
#'
#' # Set parameters of the slice sampler
#' lower.limits <- rep(-100, 5) # The dimension is the number of columns of x plus 1 (intercept)
#' upper.limits <- rep(100, 5)
#' slice.widths <- rep(1, 5)
#'
#' nMC <- 1000 # nMC should be larger in practice
#' nBI <- 250
#' result <- glm.fixed.a0(data.type=data.type, data.link=data.link, y=y, x=x, historical=historical,
#'                        lower.limits=lower.limits, upper.limits=upper.limits,
#'                        slice.widths=slice.widths, nMC=nMC, nBI=nBI)
#'
#' colMeans(result) # posterior mean of beta
#'
#' @export
glm.fixed.a0 <- function(data.type, data.link, y, n=1, x, historical=list(),
                         lower.limits=rep(-100, 50), upper.limits=rep(100, 50),
                         slice.widths=rep(1, 50), nMC=10000, nBI=250) {

  if(data.type == "Normal"){
    return(glm_fixed_a0_normal(y, x, historical, nMC, nBI))
  }else{
    return(glm_fixed_a0(data.type, data.link, y, n, x, historical, lower.limits, upper.limits, slice.widths, nMC, nBI, TRUE))
  }
}



#' Power/type I error calculation for generalized linear models with fixed a0
#'
#' @description Power/type I error calculation for generalized linear models with fixed \eqn{a_0} using power priors
#'
#' @param data.type Character string specifying the type of response. The options are "Normal", "Bernoulli", "Binomial", "Poisson" and "Exponential".
#' @param data.link Character string specifying the link function. The options are "Logistic", "Probit", "Log", "Identity-Positive", "Identity-Probability" and "Complementary Log-Log". Does not apply if \code{data.type} is "Normal".
#' @param n (For binomial data only) vector of integers specifying the number of subjects who have a particular value of the covariate vector. If the data is binary and all covariates are discrete, collapsing Bernoulli data into a binomial structure can make the slice sampler much faster.
#' @param historical (Optional) list of historical dataset(s). East historical dataset is stored in a list which constains three \emph{named} elements: \code{y0}, \code{x0} and \code{a0}.
#' \itemize{
#' \item \code{y0} is a vector of responses.
#' \item \code{x0} is a matrix of covariates. \code{x0} should NOT have the treatment indicator. Apart from missing the treatent/control indicator, \code{x0} should have the same set of covariates in the same order as \code{x}.
#' \item \code{a0} is a number between 0 and 1 indicating the discounting parameter value for that historical dataset.
#' }
#' @param x.samples Matrix of possible values of covariates from which covariate vectors are sampled with replacement. Only applies when there is no historical dataset. The matrix should not include the treatment indicator.
#' @param samp.prior.beta Matrix of possible values of \eqn{\beta} to sample (with replacement) from. Each row is a possible \eqn{\beta} vector (a realization from the sampling prior for \eqn{\beta}), where the first element is the coefficient for the intercept and the second element is the coefficient for the treatment indicator.
#' The length of the vector should be equal to the total number of parameters, i.e. P+2 where P is the number of columns of \code{x0} in \code{historical}.
#' @param samp.prior.var Vector of possible values of \eqn{\sigma^2} to sample (with replacement) from. Only applies if \code{data.type} is "Normal". The vector contains realizations from the sampling prior (e.g. inverse-gamma distribution) for \eqn{\sigma^2}.
#' @param data.size Sample size of the simulated datasets.
#' @param lower.limits Vector of lower limits for parameters to be used by the slice sampler. The length of the vector should be equal to the total number of parameters, i.e. P+1 where P is the number of covariates. The default is -100 for all parameters (may not be appropriate for all situations). Does not apply if \code{data.type} is "Normal".
#' @param upper.limits Vector of upper limits for parameters to be used by the slice sampler. The length of the vector should be equal to the total number of parameters, i.e. P+1 where P is the number of covariates. The default is 100 for all parameters (may not be appropriate for all situations). Does not apply if \code{data.type} is "Normal".
#' @param slice.widths Vector of initial slice widths for parameters to be used by the slice sampler. The length of the vector should be equal to the total number of parameters, i.e. P+1 where P is the number of covariates. The default is 1 for all parameter (may not be appropriate for all situations). Does not apply if \code{data.type} is "Normal".
#' @param nMC Number of iterations (excluding burn-in samples) for the slice sampler or Gibbs sampler. The default is 10,000.
#' @param nBI Number of burn-in samples for the slice sampler or Gibbs sampler. The default is 250.
#' @param delta Prespecified constant that defines the boundary of the null hypothesis. The default is zero.
#' @param gamma Posterior probability threshold for rejecting the null. The null hypothesis is rejected if posterior probability is greater \code{gamma}. The default is 0.95.
#' @param N Number of simulated datasets to generate. The default is 10,000.
#'
#' @details If historical datasets are provided, the algorithm samples with replacement from the historical covariates to construct the simulated datasets.
#' Otherwise, the algorithm samples with replacement from \code{x.samples}. One of the arguments \code{historical} and \code{x.samples} must be provided.
#'
#' \code{samp.prior.beta} can be generated using the sampling priors (see example).
#' \code{samp.prior.var} is necessary for generating normally distributed data.
#'
#' If \code{data.type} is "Normal", the response \eqn{y_i} is assumed to follow \eqn{N(x_i'\beta, \tau^{-1})} where \eqn{x_i} is the vector of covariates for subject \eqn{i}.
#' Each historical dataset \eqn{D_{0k}} is assumed to have a different precision parameter \eqn{\tau_k}.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}, and the initial prior for \eqn{\tau_k} is \eqn{\tau_k^{-1}}.
#' The initial prior for \eqn{\beta} is the uniform improper prior. Posterior samples are obtained through Gibbs sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling. The initial prior for \eqn{\beta} is the uniform improper prior.
#' The default lower limits for the parameters are -100. The default upper limits
#' for the parameters are 100. The default slice widths for the parameters are 1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#'
#' @return Power or type I error is returned, depending on the sampling prior used. If \code{data.type} is "Normal", average posterior means of \eqn{\beta}, \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are also returned.
#' For all other data types, the average posterior mean of \eqn{\beta} is also returned. The first column of \eqn{\beta} contains posterior samples of the intercept. The second column contains posterior samples of \eqn{\beta_1}, the parameter for the treatment indicator.
#' @references Chen, Ming-Hui, et al. "Bayesian design of noninferiority trials for medical devices using historical data." Biometrics 67.3 (2011): 1163-1170.
#'
#' Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#'
#' @seealso \code{\link{glm.fixed.a0}}
#' @examples
#'
#' data.type <- "Bernoulli"
#' data.link <- "Logistic"
#' data.size <- 100
#'
#' # Simulate two historical datasets
#' p <- 3
#' historical <- list(list(y0=rbinom(data.size,size=1,prob=0.2),
#'                         x0=matrix(rnorm(p*data.size),ncol=p,nrow=data.size), a0=0.2),
#'                    list(y0=rbinom(data.size, size=1, prob=0.5),
#'                         x0=matrix(rnorm(p*data.size),ncol=p,nrow=data.size), a0=0.3))
#'
#' # Generate sampling priors
#'
#' # The null hypothesis here is H0: beta_1 >= 0. To calculate power,
#' # we can provide samples of beta_1 such that the mass of beta_1 < 0.
#' # To calculate type I error, we can provide samples of beta_1 such that
#' # the mass of beta_1 >= 0.
#' samp.prior.beta1 <- rnorm(100, mean=-3, sd=1)
#' # Here, mass is put on the alternative region, so power is calculated.
#' samp.prior.beta <- cbind(rnorm(100), samp.prior.beta1, matrix(rnorm(100*p), 100, p))
#'
#' nMC <- 100 # nMC should be larger in practice
#' nBI <- 50
#' N <- 5 # N should be larger in practice
#' result <- power.glm.fixed.a0(data.type=data.type, data.link=data.link,
#'                              data.size=data.size, historical=historical,
#'                              samp.prior.beta=samp.prior.beta,
#'                              delta=0, nMC=nMC, nBI=nBI, N=N)
#'
#' @export
power.glm.fixed.a0 <- function(data.type, data.link="", data.size, n=1, historical=list(), x.samples=matrix(),
                               samp.prior.beta, samp.prior.var=0,
                               lower.limits=rep(-100, 50), upper.limits=rep(100, 50),
                               slice.widths=rep(1, 50),
                               delta=0, gamma=0.95, nMC=10000, nBI=250, N=10000) {

  return(power_glm_fixed_a0(data.type, data.link, data.size, n, historical, x.samples, samp.prior.beta, samp.prior.var, lower.limits, upper.limits, slice.widths, delta, gamma, nMC, nBI, N, TRUE))
}



#' Model fitting for two groups (treatment and control group, no covariates) with random a0
#'
#' @description Model fitting using normalized power priors for two groups (treatment and control group, no covariates) with random \eqn{a_0}
#'
#' @param y.c Sum of responses for the control group.
#' @param n.c Sample size of the control group.
#' @param v.c (For normal data only) variance of responses for the control group.
#' @param historical Matrix of historical dataset(s). If \code{data.type} is "Normal", \code{historical} is a matrix with three columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the variance of responses for the control group.
#' }
#' For all other data types, \code{historical} is a matrix with two columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' }
#' Each row represents a historical dataset.
#' @param lower.limits Vector of lower limits for parameters to be used by the slice sampler. The length of the vector should be equal to the number of historical datasets. The default is 0 for all parameters (may not be appropriate for all situations).
#' @param upper.limits Vector of upper limits for parameters to be used by the slice sampler. The length of the vector should be equal to the number of historical datasets. The default is 1 for all parameters (may not be appropriate for all situations).
#' @param slice.widths Vector of initial slice widths used by the slice sampler. The length of the vector should be equal to the number of historical datasets. The default is 0.1 for all parameter (may not be appropriate for all situations).
#' @param prior.a0.shape1 First shape parameter of the beta prior for \eqn{a_0}. The default is 1.
#' @param prior.a0.shape2 Second shape parameter of the beta prior for \eqn{a_0}. The default is 1.
#'
#' @inheritParams power.two.grp.fixed.a0
#' @inheritParams power.glm.fixed.a0
#'
#' @details If \code{data.type} is "Bernoulli", "Poisson" or "Exponential", a single response from the treatment group is assumed to follow Bern(\eqn{\mu_t}), Pois(\eqn{\mu_t}) or Exp(rate=\eqn{\mu_t}), respectively,
#' where \eqn{\mu_t} is the mean of responses for the treatment group. If \code{data.type} is "Normal", a single response from the treatment group is assumed to follow \eqn{N(\mu_t, \tau^{-1})}
#' where \eqn{\tau} is the precision parameter.
#' The distributional assumptions for the control group data are analogous.
#'
#' If \code{data.type} is "Bernoulli", the initial prior for \eqn{\mu_t} is beta(\code{prior.mu.t.shape1}, \code{prior.mu.t.shape2}).
#' If \code{data.type} is "Poisson", the initial prior for \eqn{\mu_t} is Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' If \code{data.type} is "Exponential", the initial prior for \eqn{\mu_t} is Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' The initial priors used for the control group data are analogous.
#'
#' If \code{data.type} is "Normal", historical datasets are assumed to have the same precision parameter \eqn{\tau} as the current dataset for computational simplicity.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}. The initial prior for the \eqn{\mu_c} is the uniform improper prior.
#' Posterior samples of \eqn{\mu_c} and \eqn{\tau} are obtained through Gibbs sampling.
#'
#' Posterior samples of \eqn{a_0} are obtained through slice sampling. The default lower limits for the parameters are 0. The default upper limits
#' for the parameters are 1.  The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' @return If \code{data.type} is "Normal", posterior samples of \eqn{\mu_c}, \eqn{\tau} and \eqn{a_0} are returned.
#' For all other data types, posterior samples of \eqn{\mu} and \eqn{a_0} are returned. If there are \eqn{K} historical datasets,
#' then \eqn{a_0 = (a_{01},\cdots,a_{0K})}.
#'
#' @references Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#' @seealso \code{\link{power.two.grp.random.a0}}
#' @examples
#'
#' data.type <- "Bernoulli"
#' y.c <- 70
#' n.c <- 100
#'
#' # Simulate three historical datasets
#' historical <- matrix(0, ncol=2, nrow=3)
#' historical[1,] <- c(70, 100)
#' historical[2,] <- c(60, 100)
#' historical[3,] <- c(50, 100)
#'
#' # Set parameters of the slice sampler
#' lower.limits <- rep(0, 3) # The dimension is the number of historical datasets
#' upper.limits <- rep(1, 3)
#' slice.widths <- rep(0.1, 3)
#'
#' result <- two.grp.random.a0(data.type=data.type, y.c=y.c, n.c=n.c, historical=historical,
#'                             lower.limits=lower.limits, upper.limits=upper.limits,
#'                             slice.widths=slice.widths, nMC=10000, nBI=250)
#' @export
two.grp.random.a0 <- function(data.type, y.c, n.c, v.c, historical,
                        prior.mu.c.shape1=1,prior.mu.c.shape2=1,
                        prior.a0.shape1=1,prior.a0.shape2=1,
                        lower.limits=rep(0, 10), upper.limits=rep(1, 10),
                        slice.widths=rep(0.1, 10), nMC=10000, nBI=250) {

  if(data.type == "Normal"){
    return(two_grp_random_a0_normal(y.c, n.c, v.c, historical, prior.a0.shape1, prior.a0.shape2, lower.limits, upper.limits, slice.widths, nMC, nBI))
  }else{
    return(two_grp_random_a0(data.type, y.c, n.c, historical, prior.mu.c.shape1, prior.mu.c.shape2,
                             prior.a0.shape1, prior.a0.shape2, lower.limits, upper.limits, slice.widths, nMC, nBI))
  }
}





#' Power/type I error calculation for two groups (treatment and control group, no covariates) with random a0
#'
#' @description Power/type I error calculation using normalized power priors for two groups (treatment and control group, no covariates) with random \eqn{a_0}
#'
#' @param n.t Sample size of the treatment group for the simulated datasets.
#' @param n.c Sample size of the control group for the simulated datasets.
#' @inheritParams two.grp.random.a0
#' @inheritParams power.two.grp.fixed.a0
#' @inheritParams power.glm.fixed.a0
#'
#' @details If \code{data.type} is "Bernoulli", "Poisson" or "Exponential", a single response from the treatment group is assumed to follow Bern(\eqn{\mu_t}), Pois(\eqn{\mu_t}) or Exp(rate=\eqn{\mu_t}), respectively,
#' where \eqn{\mu_t} is the mean of responses for the treatment group. If \code{data.type} is "Normal", a single response from the treatment group is assumed to follow \eqn{N(\mu_t, \tau^{-1})}
#' where \eqn{\tau} is the precision parameter.
#' The distributional assumptions for the control group data are analogous.
#'
#' \code{samp.prior.mu.t} and \code{samp.prior.mu.c} can be generated using the sampling priors (see example).
#'
#' If \code{data.type} is "Bernoulli", the initial prior for \eqn{\mu_t} is beta(\code{prior.mu.t.shape1}, \code{prior.mu.t.shape2}).
#' If \code{data.type} is "Poisson", the initial prior for \eqn{\mu_t} is Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' If \code{data.type} is "Exponential", the initial prior for \eqn{\mu_t} is Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' The initial priors used for the control group data are analogous.
#'
#' If \code{data.type} is "Normal", historical datasets are assumed to have the same precision parameter as the current dataset for computational simplicity.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}. The initial prior for the \eqn{\mu_c} is the uniform improper prior.
#' Posterior samples of \eqn{\mu_c} and \eqn{\tau} are obtained through Gibbs sampling.
#'
#' Posterior samples of \eqn{a_0} are obtained through slice sampling. The default lower limits for the parameters are 0. The default upper limits
#' for the parameters are 1.  The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#'
#' @return Power or type I error is returned, depending on the sampling prior used. If \code{data.type} is "Normal", average posterior means of \eqn{\mu_t}, \eqn{\mu_c}, \eqn{\tau} and \eqn{a_0} are also returned.
#' For all other data types, average posterior means of \eqn{\mu_t}, \eqn{\mu_c} and \eqn{a_0} are also returned.
#' @references Chen, Ming-Hui, et al. "Bayesian design of noninferiority trials for medical devices using historical data." Biometrics 67.3 (2011): 1163-1170.
#'
#' Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#'
#' @seealso \code{\link{two.grp.random.a0}}
#' @examples
#'
#' data.type <- "Bernoulli"
#' n.t <- 100
#' n.c <- 100
#'
#' # Simulate three historical datasets
#' historical <- matrix(0, ncol=2, nrow=3)
#' historical[1,] <- c(70, 100)
#' historical[2,] <- c(60, 100)
#' historical[3,] <- c(50, 100)
#'
#' # Generate sampling priors
#' set.seed(1)
#' b_st1 <- b_st2 <- 1
#' b_sc1 <- b_sc2 <- 1
#' samp.prior.mu.t <- rbeta(50000, b_st1, b_st2)
#' samp.prior.mu.c <- rbeta(50000, b_st1, b_st2)
#' # The null hypothesis here is H0: mu_t - mu_c >= 0. To calculate power,
#' # we can provide samples of mu.t and mu.c such that the mass of mu_t - mu_c < 0.
#' # To calculate type I error, we can provide samples of mu.t and mu.c such that
#' # the mass of mu_t - mu_c >= 0.
#' sub_ind <- which(samp.prior.mu.t < samp.prior.mu.c)
#' # Here, mass is put on the alternative region, so power is calculated.
#' samp.prior.mu.t <- samp.prior.mu.t[sub_ind]
#' samp.prior.mu.c <- samp.prior.mu.c[sub_ind]
#'
#' N <- 10 # N should be larger in practice
#' result <- power.two.grp.random.a0(data.type=data.type, n.t=n.t, n.c=n.c, historical=historical,
#'                                   samp.prior.mu.t=samp.prior.mu.t, samp.prior.mu.c=samp.prior.mu.c,
#'                                   delta=0, nMC=10000, nBI=250, N=N)
#' @export
power.two.grp.random.a0 <- function(data.type, n.t, n.c, historical,
                              samp.prior.mu.t, samp.prior.mu.c,
                              samp.prior.var.t=0, samp.prior.var.c=0,
                              prior.mu.t.shape1=1,prior.mu.t.shape2=1,
                              prior.mu.c.shape1=1,prior.mu.c.shape2=1,
                              prior.a0.shape1=1,prior.a0.shape2=1,
                              lower.limits=rep(0, 10), upper.limits=rep(1, 10),
                              slice.widths=rep(0.1, 10), delta=0, gamma=0.95,
                              nMC=10000, nBI=250, N=10000) {

  return(power_two_grp_random_a0(data.type, n.t, n.c, historical,
                                 samp.prior.mu.t, samp.prior.mu.c,
                                 samp.prior.var.t, samp.prior.var.c,
                                 prior.mu.t.shape1, prior.mu.t.shape2,
                                 prior.mu.c.shape1, prior.mu.c.shape2,
                                 prior.a0.shape1, prior.a0.shape2,
                                 lower.limits, upper.limits, slice.widths,
                                 delta, gamma, nMC, nBI, N))

}



#' Model fitting for generalized linear models with random a0
#'
#' @description Model fitting using normalized power priors for generalized linear models with random \eqn{a_0}
#'
#' @param historical List of historical dataset(s). East historical dataset is stored in a list which constains two \emph{named} elements: \code{y0} and \code{x0}.
#' \itemize{
#' \item \code{y0} is a vector of responses.
#' \item \code{x0} is a matrix of covariates. \code{x0} should NOT have the treatment indicator. Apart from missing the treatent/control indicator, \code{x0} should have the same set of covariates in the same order as \code{x}.
#' }
#' @param lower.limits Vector of lower limits for parameters to be used by the slice sampler. If \code{data.type} is "Normal", slice sampling is used for \eqn{a_0}, and the length of the vector should be equal to the number of historical datasets.
#' For all other data types, slice sampling is used for \eqn{\beta} and \eqn{a_0}. The first P+1 elements apply to the sampling of \eqn{\beta} and the rest apply to the sampling of \eqn{a_0}.
#' The length of the vector should be equal to the sum of the total number of parameters (i.e. P+1 where P is the number of covariates) and the number of historical datasets.
#'  The default is -100 for all parameters (may not be appropriate for all situations).
#' @param upper.limits Vector of upper limits for parameters to be used by the slice sampler. If \code{data.type} is "Normal", slice sampling is used for \eqn{a_0}, and the length of the vector should be equal to the number of historical datasets.
#' For all other data types, slice sampling is used for \eqn{\beta} and \eqn{a_0}. The first P+1 elements apply to the sampling of \eqn{\beta} and the rest apply to the sampling of \eqn{a_0}.
#' The length of the vector should be equal to the sum of the total number of parameters (i.e. P+1 where P is the number of covariates) and the number of historical datasets.
#'  The default is 100 for all parameters (may not be appropriate for all situations).
#' @param slice.widths Vector of initial slice widths used by the slice sampler. If \code{data.type} is "Normal", slice sampling is used for \eqn{a_0}, and the length of the vector should be equal to the number of historical datasets.
#' For all other data types, slice sampling is used for \eqn{\beta} and \eqn{a_0}. The first P+1 elements apply to the sampling of \eqn{\beta} and the rest apply to the sampling of \eqn{a_0}.
#' The length of the vector should be equal to the sum of the total number of parameters (i.e. P+1 where P is the number of covariates) and the number of historical datasets.
#' The default is 0.1 for all parameter (may not be appropriate for all situations).
#' @param a0.coefficients Vector of coefficients for \eqn{a_0} returned by the function \code{\link{normalizing.constant}}. This is necessary for estimating the normalizing constant for the normalized power prior. Does not apply if \code{data.type} is "Normal".
#'
#' @inheritParams power.glm.fixed.a0
#' @inheritParams glm.fixed.a0
#' @inheritParams two.grp.random.a0
#'
#' @details
#' The user should use the function \code{\link{normalizing.constant}} to obtain \code{a0.coefficients} (does not apply if \code{data.type} is "Normal").
#'
#' If \code{data.type} is "Normal", the response \eqn{y_i} is assumed to follow \eqn{N(x_i'\beta, \tau^{-1})} where \eqn{x_i} is the vector of covariates for subject \eqn{i}.
#' Historical datasets are assumed to have the same precision parameter as the current dataset for computational simplicity.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}.
#' The initial prior for \eqn{\beta} is the uniform improper prior. Posterior samples for \eqn{\beta} and \eqn{\tau} are obtained through Gibbs sampling.
#' Posterior samples for \eqn{a_0} are obtained through slice sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling. The initial prior for \eqn{\beta} is the uniform improper prior.
#' The default lower limits for the parameters are -100. The default upper limits
#' for the parameters are 100. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#'
#' @return If \code{data.type} is "Normal", posterior samples of \eqn{\beta}, \eqn{\tau} and \eqn{a_0} are returned.
#' For all other data types, posterior samples of \eqn{\beta} and \eqn{a_0} are returned.
#' @references Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#' @seealso \code{\link{normalizing.constant}} and \code{\link{power.glm.random.a0}}
#' @examples
#'
#' data.type <- "Bernoulli"
#' data.link <- "Logistic"
#'
#' # Simulate current data
#' set.seed(1)
#' p <- 3
#' n_total <- 100
#' y <- rbinom(n_total,size=1,prob=0.6)
#' # The first column of x is the treatment indicator.
#' x <- cbind(rbinom(n_total,size=1,prob=0.5),
#'            matrix(rnorm(p*n_total),ncol=p,nrow=n_total))
#'
#' # Simulate two historical datasets
#' # Note that x0 does not have the treatment indicator
#' historical <- list(list(y0=rbinom(n_total,size=1,prob=0.2),
#'                         x0=matrix(rnorm(p*n_total),ncol=p,nrow=n_total)),
#'                    list(y0=rbinom(n_total, size=1, prob=0.5),
#'                         x0=matrix(rnorm(p*n_total),ncol=p,nrow=n_total)))
#'
#' # Please see function "normalizing.constant" for how to obtain a0.coefficients
#' a0.coefficients <- c(1, 0.5)
#'
#' # Set parameters of the slice sampler
#' lower.limits <- rep(-100, 5) # The dimension is the number of columns of x plus 1 (intercept)
#' upper.limits <- rep(100, 5)
#' slice.widths <- rep(0.1, 5)
#'
#' nMC <- 500 # nMC should be larger in practice
#' nBI <- 100
#' result <- glm.random.a0(data.type=data.type, data.link=data.link, y=y, x=x,
#'                         historical=historical, a0.coefficients=a0.coefficients,
#'                         lower.limits=lower.limits, upper.limits=upper.limits,
#'                         slice.widths=slice.widths, nMC=nMC, nBI=nBI)
#'
#' @export
glm.random.a0 <- function(data.type, data.link, y, n=1, x, historical,
                         prior.a0.shape1=1, prior.a0.shape2=1, a0.coefficients,
                         lower.limits=rep(-100, 50), upper.limits=rep(100, 50),
                         slice.widths=rep(0.1, 50), nMC=10000, nBI=250) {

  if(data.type == "Normal"){
    return(glm_random_a0_normal(y, x, historical, prior.a0.shape1, prior.a0.shape2, lower.limits, upper.limits, slice.widths, nMC, nBI))
  }else{
    return(glm_random_a0(data.type, data.link, y, n, x, historical,
                         prior.a0.shape1, prior.a0.shape2, a0.coefficients,
                         lower.limits, upper.limits,
                         slice.widths, nMC, nBI))
  }
}


#' Power/type I error calculation for generalized linear models with random a0
#'
#' @description Power/type I error calculation using normalized power priors for generalized linear models with random \eqn{a_0}
#'

#' @inheritParams glm.random.a0
#' @inheritParams power.glm.fixed.a0
#' @inheritParams two.grp.random.a0
#' @details  The user should use the function \code{\link{normalizing.constant}} to obtain \code{a0.coefficients} (does not apply if \code{data.type} is "Normal").
#'
#' \code{samp.prior.beta} can be generated using the sampling priors (see example).
#' \code{samp.prior.var} is necessary for generating normally distributed data.
#'
#' If \code{data.type} is "Normal", the response \eqn{y_i} is assumed to follow \eqn{N(x_i'\beta, \tau^{-1})} where \eqn{x_i} is the vector of covariates for subject \eqn{i}.
#' Historical datasets are assumed to have the same precision parameter as the current dataset for computational simplicity.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}.
#' The initial prior for \eqn{\beta} is the uniform improper prior. Posterior samples for \eqn{\beta} and \eqn{\tau} are obtained through Gibbs sampling.
#' Posterior samples for \eqn{a_0} are obtained through slice sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling. The initial prior for \eqn{\beta} is the uniform improper prior.
#' The default lower limits for the parameters are -100. The default upper limits
#' for the parameters are 100. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#'
#' @return Power or type I error is returned, depending on the sampling prior used. If \code{data.type} is "Normal", average posterior means of \eqn{\beta}, \eqn{\tau} and \eqn{a_0} are also returned.
#' For all other data types, average posterior means of \eqn{\beta} and \eqn{a_0} are also returned.
#' @references Chen, Ming-Hui, et al. "Bayesian design of noninferiority trials for medical devices using historical data." Biometrics 67.3 (2011): 1163-1170.
#'
#' Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#' @seealso \code{\link{normalizing.constant}} and \code{\link{glm.random.a0}}
#' @examples
#'
#' data.type <- "Bernoulli"
#' data.link <- "Logistic"
#' data.size <- 100
#'
#' # Simulate two historical datasets
#' p <- 3
#' historical <- list(list(y0=rbinom(data.size,size=1,prob=0.2),
#'                         x0=matrix(rnorm(p*data.size),ncol=p,nrow=data.size)),
#'                    list(y0=rbinom(data.size, size=1, prob=0.5),
#'                         x0=matrix(rnorm(p*data.size),ncol=p,nrow=data.size)))
#'
#' # Generate sampling priors
#'
#' # The null hypothesis here is H0: beta_1 >= 0. To calculate power,
#' # we can provide samples of beta_1 such that the mass of beta_1 < 0.
#' # To calculate type I error, we can provide samples of beta_1 such that
#' # the mass of beta_1 >= 0.
#' samp.prior.beta1 <- rnorm(100, mean=-3, sd=1)
#' # Here, mass is put on the alternative region, so power is calculated.
#' samp.prior.beta <- cbind(rnorm(100), samp.prior.beta1, matrix(rnorm(100*p), 100, p))
#'
#' # Please see function "normalizing.constant" for how to obstain a0.coefficients
#' a0.coefficients <- c(1, 0.5)
#'
#' nMC <- 100 # nMC should be larger in practice
#' nBI <- 50
#' N <- 3 # N should be larger in practice
#' result <- power.glm.random.a0(data.type=data.type, data.link=data.link,
#'                               data.size=data.size, historical=historical,
#'                               samp.prior.beta=samp.prior.beta, a0.coefficients=a0.coefficients,
#'                               delta=0, nMC=nMC, nBI=nBI, N=N)
#'
#' @export
power.glm.random.a0 <- function(data.type, data.link, data.size, n=1, historical,
                                samp.prior.beta, samp.prior.var,
                                prior.a0.shape1=1, prior.a0.shape2=1, a0.coefficients,
                                lower.limits=rep(-100, 50), upper.limits=rep(100, 50),slice.widths=rep(0.1, 50),
                                delta=0, gamma=0.95, nMC=10000, nBI=250, N=10000) {

  if(data.type == "Normal"){
    return(power_glm_random_a0_normal(data.size, historical,
                                      samp.prior.beta, samp.prior.var,
                                      prior.a0.shape1, prior.a0.shape2,
                                      lower.limits, upper.limits, slice.widths,
                                      delta, gamma, nMC, nBI, N))
  }else{
    return(power_glm_random_a0(data.type, data.link, data.size, n, historical,
                               samp.prior.beta, prior.a0.shape1, prior.a0.shape2,
                               a0.coefficients, lower.limits, upper.limits, slice.widths,
                               delta, gamma, nMC, nBI, N))
  }

}



