### two group, fixed a0 #####



#' Model fitting for two groups (treatment and control group, no covariates) with fixed a0
#'
#' @description Model fitting using power priors for two groups (treatment and control group, no covariates) with fixed \eqn{a_0}
#'
#' @param historical (Optional) matrix of historical dataset(s). If \code{data.type} is "Normal", \code{historical} is a matrix with four columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the sample variance of responses for the control group.
#' \item The fourth column contains the discounting parameter value \eqn{a_0} (between 0 and 1).
#' }
#' For all other data types, \code{historical} is a matrix with three columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the discounting parameter value \eqn{a_0} (between 0 and 1).
#' }
#' Each row represents a historical dataset.
#' @param nMC (For normal data only) number of iterations (excluding burn-in samples) for the Gibbs sampler. The default is 10,000.
#' @param nBI (For normal data only) number of burn-in samples for the Gibbs sampler. The default is 250.

#' @inheritParams two.grp.random.a0
#' @inheritParams power.two.grp.fixed.a0
#'
#' @details The power prior is applied on the data of the control group only.
#' Therefore, only summaries of the responses of the control group need to be entered.
#'
#' If \code{data.type} is "Bernoulli", "Poisson" or "Exponential", a single response from the treatment group is assumed to follow Bern(\eqn{\mu_t}), Pois(\eqn{\mu_t}) or Exp(rate=\eqn{\mu_t}), respectively,
#' where \eqn{\mu_t} is the mean of responses for the treatment group. The distributional assumptions for the control group data are analogous.
#'
#' If \code{data.type} is "Bernoulli", the initial prior for \eqn{\mu_t} is beta(\code{prior.mu.t.shape1}, \code{prior.mu.t.shape2}).
#' If \code{data.type} is "Poisson", the initial prior for \eqn{\mu_t} is Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' If \code{data.type} is "Exponential", the initial prior for \eqn{\mu_t} is Gamma(\code{prior.mu.t.shape1}, rate=\code{prior.mu.t.shape2}).
#' The initial priors used for the control group data are analogous.
#'
#' If \code{data.type} is "Normal", the responses are assumed to follow \eqn{N(\mu_c, \tau^{-1})} where \eqn{\mu_c} is the mean of responses for the control group
#' and \eqn{\tau} is the precision parameter. Each historical dataset \eqn{D_{0k}} is assumed to have a different precision parameter \eqn{\tau_k}.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}, and the initial prior for \eqn{\tau_k} is \eqn{\tau_k^{-1}}. The initial prior for the \eqn{\mu_c} is the uniform improper prior.
#' Posterior samples are obtained through Gibbs sampling.
#'
#' @return The function returns a S3 object with a \code{summary} method. If \code{data.type} is "Normal", posterior samples of \eqn{\mu_c}, \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are returned
#' in the list item named \code{posterior.params}.
#' For all other data types, two scalars, \eqn{c_1} and \eqn{c_2}, are returned in the list item named \code{posterior.params}, representing the two parameters of the posterior distribution of \eqn{\mu_c}.
#' For Bernoulli responses, the posterior distribution of \eqn{\mu_c} is beta(\eqn{c_1}, \eqn{c_2}).
#' For Poisson responses, the posterior distribution of \eqn{\mu_c} is Gamma(\eqn{c_1}, \eqn{c_2}) where \eqn{c_2} is the rate parameter.
#' For exponential responses, the posterior distribution of \eqn{\mu_c} is Gamma(\eqn{c_1}, \eqn{c_2}) where \eqn{c_2} is the rate parameter.
#' @references Chen, Ming-Hui, et al. "Bayesian design of noninferiority trials for medical devices using historical data." Biometrics 67.3 (2011): 1163-1170.
#' @seealso \code{\link{power.two.grp.fixed.a0}}
#' @examples
#' data.type <- "Bernoulli"
#' y.c <- 70
#' n.c <- 100
#'
#' # Simulate three historical datasets
#' historical <- matrix(0, ncol=3, nrow=3)
#' historical[1,] <- c(70, 100, 0.3)
#' historical[2,] <- c(60, 100, 0.5)
#' historical[3,] <- c(50, 100, 0.7)
#'
#' set.seed(1)
#' result <- two.grp.fixed.a0(data.type=data.type, y.c=y.c, n.c=n.c, historical=historical)
#' summary(result)
#' @export
two.grp.fixed.a0 <- function(data.type, y.c, n.c, v.c, historical=matrix(0,1,4), prior.mu.c.shape1=1, prior.mu.c.shape2=1, nMC=10000, nBI=250) {
  if(data.type == "Normal"){
    result <- two_grp_fixed_a0_normal(y.c, n.c, v.c, historical, nMC, nBI)
  }else{
    result <- two_grp_fixed_a0(data.type, y.c, n.c,
                            historical, prior.mu.c.shape1, prior.mu.c.shape2)
  }
  out <- list(posterior.params=result, data.type=data.type)
  
  structure(out, class=c("tgfixed"))
  
}



#' @importFrom stats quantile
#' @export
summary.tgfixed <- function(object, ...) {
  
  r <- object$posterior.params
  
  if(object$data.type=="Normal"){
    # mu_c
    mu_c <- r$`posterior samples of mu_c`
    output1 <- cbind("mean"=mean(mu_c),"sd"=sd(mu_c),"2.5%"=quantile(mu_c,probs=0.025),"97.5%"=quantile(mu_c,probs=0.975))
    rownames(output1) <- "mu_c"
    output1 <- round(output1, digits=3)
    
    # tau
    taus <- r$`posterior samples of tau`
    output2 <- cbind("mean"=mean(taus),"sd"=sd(taus),"2.5%"=quantile(taus,probs=0.025),"97.5%"=quantile(taus,probs=0.975))
    rownames(output2) <- "tau"
    output2 <- round(output2, digits=3)
    # tau0
    if(length(r)==2){
      output3 <- NULL
    }else{
      tau0s <- r$`posterior samples of tau_0`
      m <- apply(tau0s, 2, mean)
      std <- apply(tau0s, 2, sd)
      q1 <- apply(tau0s, 2, quantile, probs=0.025)
      q2 <- apply(tau0s, 2, quantile, probs=0.975)
      output3 <- cbind("mean"=m,"sd"=std,"2.5%"=q1,"97.5%"=q2)
      rownames(output3) <- paste0("tau0_", 1:nrow(output3))
      output3 <- round(output3, digits=3)
    }
    
    out <- rbind(output1, output2, output3)
    out
    
  }else{
    c1 <- round(r[1],3)
    c2 <- round(r[2],3)
    if(object$data.type=="Bernoulli"){ cat("The posterior distribution for mu_c is beta(",c1,", ",c2,").")}
    if(object$data.type=="Poisson"){ cat("The posterior distribution for mu_c is gamma(",c1,", rate=",c2,").")}
    if(object$data.type=="Exponential"){ cat("The posterior distribution for mu_c is gamma(",c1,", rate=",c2,").")}

  }
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
#' \item The third column contains the sample variance of responses for the control group.
#' \item The fourth column contains the discounting parameter value \eqn{a_0} (between 0 and 1).
#' }
#' For all other data types, \code{historical} is a matrix with three columns:
#' \itemize{
#' \item The first column contains the sum of responses for the control group.
#' \item The second column contains the sample size of the control group.
#' \item The third column contains the discounting parameter value \eqn{a_0} (between 0 and 1).
#' }
#' Each row represents a historical dataset.
#' @param nullspace.ineq Character string specifying the inequality of the null hypothesis. The options are ">" and "<". If ">" is specified, the null hypothesis (for non-exponential data) is \eqn{H_0}: \eqn{\mu_t} - \eqn{\mu_c} \eqn{\ge} \eqn{\delta}. If "<" is specified, the null hypothesis is \eqn{H_0}: \eqn{\mu_t} - \eqn{\mu_c} \eqn{\le} \eqn{\delta}. The default choice is ">".
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
#' @return The function returns a S3 object with a \code{summary} method. Power or type I error is returned, depending on the sampling prior used.
#' The posterior probabilities of the alternative hypothesis are returned.  
#' Average posterior means of \eqn{\mu_t} and \eqn{\mu_c} and their corresponding biases are returned. 
#' If \code{data.type} is "Normal", average posterior means of \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are also returned.
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
#' summary(result)
#' @export
power.two.grp.fixed.a0 <- function(data.type, n.t, n.c, historical=matrix(0,1,4), nullspace.ineq=">", samp.prior.mu.t,
                                   samp.prior.mu.c, samp.prior.var.t, samp.prior.var.c, prior.mu.t.shape1=1,
                                   prior.mu.t.shape2=1, prior.mu.c.shape1=1, prior.mu.c.shape2=1, delta=0, 
                                   gamma=0.95, nMC=10000, nBI=250, N=10000) {
  
  if(data.type == "Normal"){
    out <- power_two_grp_fixed_a0_normal(n.t, n.c, historical, nullspace.ineq, samp.prior.mu.t, samp.prior.mu.c,
                                         samp.prior.var.t, samp.prior.var.c, delta, gamma, nMC, nBI, N)
  }else{
    out <- power_two_grp_fixed_a0(data.type, n.t, n.c, historical, nullspace.ineq, samp.prior.mu.t, samp.prior.mu.c,
                                  prior.mu.t.shape1, prior.mu.t.shape2, prior.mu.c.shape1,
                                  prior.mu.c.shape2, delta, gamma, N, Inf)
  }
  structure(out, class=c("powertg"))

}


#' @export
summary.powertg <- function(object, ...) {
  
  r <- round(object$`power/type I error`, digits=3)
  # beta
  mus <- c(object$`average posterior mean of mu_t`, object$`average posterior mean of mu_c`)
  bias <- c(object$`bias of the average posterior mean of mu_t`, object$`bias of the average posterior mean of mu_c`)
  postprob <- mean(object[[2]])
  output1 <- cbind(mus,bias)
  colnames(output1) <- c("average posterior mean","bias")
  rownames(output1) <- c("mu_t","mu_c")
  output1 <- round(output1, digits=3)
  print(output1)
  cat("The power/type I error rate is ",r,".\n")
  cat("The average of the", names(object[2]), "is",round(mean(object[[2]]),3),".")
  
}




