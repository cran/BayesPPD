### two group, random a0 #####




#' Model fitting for two groups (treatment and control group, no covariates) with random a0
#'
#' @description Model fitting using normalized power priors for two groups (treatment and control group, no covariates) with random \eqn{a_0}
#'
#' @param y.c Sum of responses for the control group.
#' @param n.c Sample size of the control group.
#' @param v.c (For normal data only) sample variance of responses for the control group.
#' @param historical Matrix of historical dataset(s). If \code{data.type} is "Normal", \code{historical} is a matrix with three columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the sample variance of responses for the control group.
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
#' @param prior.a0.shape1 Vector of the first shape parameters of the independent beta priors for \eqn{a_0}. The length of the vector should be equal to the number of historical datasets. The default is a vector of one's.
#' @param prior.a0.shape2 Vector of the second shape parameters of the independent beta priors for \eqn{a_0}. The length of the vector should be equal to the number of historical datasets. The default is a vector of one's.
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
#' Independent beta(\code{prior.a0.shape1},\code{prior.a0.shape1}) priors are used for \eqn{a_0}. Posterior samples of \eqn{a_0} are obtained through slice sampling. The default lower limits for the parameters are 0. The default upper limits
#' for the parameters are 1.  The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' @return The function returns a S3 object with a \code{summary} method. If \code{data.type} is "Normal", posterior samples of \eqn{\mu_c}, \eqn{\tau} and \eqn{a_0} are returned.
#' For all other data types, posterior samples of \eqn{\mu_c} and \eqn{a_0} are returned. If there are \eqn{K} historical datasets,
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
#' set.seed(1)
#' result <- two.grp.random.a0(data.type=data.type, y.c=y.c, n.c=n.c, historical=historical,
#'                             lower.limits=lower.limits, upper.limits=upper.limits,
#'                             slice.widths=slice.widths, nMC=10000, nBI=250)
#' summary(result)
#' @export
two.grp.random.a0 <- function(data.type, y.c, n.c, v.c, historical,
                              prior.mu.c.shape1=1,prior.mu.c.shape2=1,
                              prior.a0.shape1=rep(1,10),prior.a0.shape2=rep(1,10),
                              lower.limits=rep(0, 10), upper.limits=rep(1, 10),
                              slice.widths=rep(0.1, 10), nMC=10000, nBI=250) {
  
  if(data.type == "Normal"){
    result <- two_grp_random_a0_normal(y.c, n.c, v.c, historical, prior.a0.shape1, prior.a0.shape2, lower.limits, upper.limits, slice.widths, nMC, nBI)
  }else{
    result <- two_grp_random_a0(data.type, y.c, n.c, historical, prior.mu.c.shape1, prior.mu.c.shape2,
                             prior.a0.shape1, prior.a0.shape2, lower.limits, upper.limits, slice.widths, nMC, nBI)
  }
  out <- list(posterior.samples=result, data.type=data.type)
  structure(out, class=c("tgrandom"))
}

#' @importFrom stats quantile
#' @export
summary.tgrandom <- function(object, ...) {
  
  r <- object$posterior.samples
  
  # mu_c
  mu_c <- r$`posterior samples of mu_c`
  output1 <- cbind("mean"=mean(mu_c),"sd"=sd(mu_c),"2.5%"=quantile(mu_c,probs=0.025),"97.5%"=quantile(mu_c,probs=0.975))
  rownames(output1) <- "mu_c"
  output1 <- round(output1, digits=3)
  
  if(object$data.type=="Normal"){
    # tau
    taus <- r$`posterior samples of tau`
    output2 <- cbind("mean"=mean(taus),"sd"=sd(taus),"2.5%"=quantile(taus,probs=0.025),"97.5%"=quantile(taus,probs=0.975))
    rownames(output2) <- "tau"
    output2 <- round(output2, digits=3)
  }else{
    output2 <- NULL
  }
  
  # a0
  a0s <- r$`posterior samples of a0`
  m <- apply(a0s, 2, mean)
  std <- apply(a0s, 2, sd)
  q1 <- apply(a0s, 2, quantile, probs=0.025)
  q2 <- apply(a0s, 2, quantile, probs=0.975)
  output3 <- cbind("mean"=m,"sd"=std,"2.5%"=q1,"97.5%"=q2)
  rownames(output3) <- paste0("a0_", 1:nrow(output3))
  output3 <- round(output3, digits=3)
  
  out <- rbind(output1, output2, output3)
  out

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
#' Independent beta(\code{prior.a0.shape1},\code{prior.a0.shape1}) priors are used for \eqn{a_0}. Posterior samples of \eqn{a_0} are obtained through slice sampling. The default lower limits for the parameters are 0. The default upper limits
#' for the parameters are 1.  The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#'
#' @return The function returns a S3 object with a \code{summary} method. Power or type I error is returned, depending on the sampling prior used.
#' The posterior probabilities of the alternative hypothesis are returned.  
#' Average posterior means of \eqn{\mu_t} and \eqn{\mu_c} and their corresponding biases are returned. 
#' The average posterior mean of \eqn{a_0} is returned. 
#' If \code{data.type} is "Normal", the average posterior mean of \eqn{\tau} is also returned.
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
#' summary(result)
#' @export
power.two.grp.random.a0 <- function(data.type, n.t, n.c, historical,nullspace.ineq=">",
                                    samp.prior.mu.t, samp.prior.mu.c,
                                    samp.prior.var.t=0, samp.prior.var.c=0,
                                    prior.mu.t.shape1=1,prior.mu.t.shape2=1,
                                    prior.mu.c.shape1=1,prior.mu.c.shape2=1,
                                    prior.a0.shape1=rep(1,10),prior.a0.shape2=rep(1,10),
                                    lower.limits=rep(0, 10), upper.limits=rep(1, 10),
                                    slice.widths=rep(0.1, 10), delta=0, gamma=0.95,
                                    nMC=10000, nBI=250, N=10000) {
  
  out <- power_two_grp_random_a0(data.type, n.t, n.c, historical,nullspace.ineq,
                                 samp.prior.mu.t, samp.prior.mu.c,
                                 samp.prior.var.t, samp.prior.var.c,
                                 prior.mu.t.shape1, prior.mu.t.shape2,
                                 prior.mu.c.shape1, prior.mu.c.shape2,
                                 prior.a0.shape1, prior.a0.shape2,
                                 lower.limits, upper.limits, slice.widths,
                                 delta, gamma, nMC, nBI, N)
  
  structure(out, class=c("powertg"))
}

