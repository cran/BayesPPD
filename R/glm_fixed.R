### glm, fixed a0 #####




#' Model fitting for generalized linear models with fixed a0
#'
#' @description Model fitting using power priors for generalized linear models with fixed \eqn{a_0}
#'
#' @param y Vector of responses.
#' @param x Matrix of covariates. The first column should be the treatment indicator with 1 indicating treatment group. The number of rows should equal the length of the response vector \code{y}.
#' @param n (For binomial data only) vector of integers specifying the number of subjects who have a particular value of the covariate vector. If the data is binary and all covariates are discrete, collapsing Bernoulli data into a binomial structure can make the slice sampler much faster.
#' The length of \code{n} should be equal to the number of rows of \code{x}. 
#' @param historical (Optional) list of historical dataset(s). East historical dataset is stored in a list which contains three \emph{named} elements: \code{y0}, \code{x0} and \code{a0}.
#' \itemize{
#' \item \code{y0} is a vector of responses.
#' \item \code{x0} is a matrix of covariates. If \code{borrow.treat} is FALSE (the default), \code{x0} should NOT have the treatment indicator. Apart from missing the treatment indicator, \code{x0} should have the same set of covariates in the same order as \code{x}.
#' If \code{borrow.treat} is TRUE, \code{x0} should have the same set of covariates in the same order as \code{x}, where the first column of \code{x0} must be the treatment indicator.
#' \item \code{a0} is a number between 0 and 1 indicating the discounting parameter value for that historical dataset.
#' }
#' For binomial data, an additional element \code{n0} is required.
#' \itemize{
#' \item \code{n0} is vector of integers specifying the number of subjects who have a particular value of the covariate vector. 
#' The length of \code{n0} should be equal to the number of rows of \code{x0}. 
#' }
#' @param current.data Logical value indicating whether current data is included. The default is TRUE. If FALSE, only historical data is included in the analysis,
#' and the posterior samples can be used as a discrete approximation to the sampling prior in \code{\link{power.glm.fixed.a0}}.
#' @param prior.beta.var Only applies if current.data = FALSE. If no current data is provided, the initial priors used for \eqn{\beta} are i.i.d. normal distributions with mean zero and variance equal to \code{prior.beta.var}. 
#' The length of the vector should be equal to the length of \eqn{\beta}. The default variance is 10. 
#' 
#' @inheritParams power.glm.fixed.a0
#'
#' @details If \code{data.type} is "Normal", the response \eqn{y_i} is assumed to follow \eqn{N(x_i'\beta, \tau^{-1})} where \eqn{x_i} is the vector of covariates for subject \eqn{i}.
#' Each historical dataset \eqn{D_{0k}} is assumed to have a different precision parameter \eqn{\tau_k}.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}, and the initial prior for \eqn{\tau_k} is \eqn{\tau_k^{-1}}.
#' The initial prior for \eqn{\beta} is the uniform improper prior. Posterior samples are obtained through Gibbs sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling. 
#' The default lower limits for the parameters are -100. The default upper limits
#' for the parameters are 100. The default slice widths for the parameters are 1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#' 
#' When \code{current.data} is set to FALSE, only historical data is included in the analysis,
#' and the posterior samples can be used as a discrete approximation to the sampling prior in \code{\link{power.glm.fixed.a0}}.
#'
#' @return The function returns a S3 object with a \code{summary} method. If \code{data.type} is "Normal", posterior samples of \eqn{\beta}, \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are returned.
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
#' summary(result)
#'
#' @export
glm.fixed.a0 <- function(data.type, data.link, y=0, x=matrix(), n=1, borrow.treat=FALSE, 
                         historical=list(), 
                         lower.limits=rep(-100, 50), upper.limits=rep(100, 50),
                         slice.widths=rep(1, 50), nMC=10000, nBI=250, current.data=TRUE, 
                         prior.beta.var = rep(10, 50)) {
  
  if(is.na(x[1,1])){
    x <- historical[[1]]$x0
  }
  if(is.null(colnames(x))){
    colnames(x) <- paste0("X", 1:ncol(x))
  }
  
  if(data.type == "Normal"){
    result <- glm_fixed_a0_normal(y, x, borrow.treat, historical, nMC, nBI)
    colnames(result$`posterior samples of beta`) <- c("intercept",colnames(x))
  }else{
    result <- glm_fixed_a0(data.type, data.link, y, n, x, borrow.treat, historical, prior.beta.var, lower.limits, upper.limits, slice.widths, nMC, nBI, current.data)
    colnames(result) <- c("intercept",colnames(x))
  }
  
  out <- list(posterior.samples=result, data.type=data.type)
  structure(out, class=c("glmfixed"))
}

#' @importFrom stats quantile
#' @export
summary.glmfixed <- function(object, ...) {
  
  r <- object$posterior.samples 
  
  if(object$data.type=="Normal"){
    
    # beta
    betas <- r$`posterior samples of beta`
    m <- apply(betas, 2, mean)
    std <- apply(betas, 2, sd)
    q1 <- apply(betas, 2, quantile, probs=0.025)
    q2 <- apply(betas, 2, quantile, probs=0.975)
    output1 <- cbind("mean"=m,"sd"=std,"2.5%"=q1,"97.5%"=q2)
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
    
    m <- apply(r, 2, mean)
    std <- apply(r, 2, sd)
    q1 <- apply(r, 2, quantile, probs=0.025)
    q2 <- apply(r, 2, quantile, probs=0.975)
    output <- cbind("mean"=m,"sd"=std,"2.5%"=q1,"97.5%"=q2)
    round(output, digits=3)
  }
  
}



#' Power/type I error calculation for generalized linear models with fixed a0
#'
#' @description Power/type I error calculation for generalized linear models with fixed \eqn{a_0} using power priors
#'
#' @param data.type Character string specifying the type of response. The options are "Normal", "Bernoulli", "Binomial", "Poisson" and "Exponential".
#' @param data.link Character string specifying the link function. The options are "Logistic", "Probit", "Log", "Identity-Positive", "Identity-Probability" and "Complementary Log-Log". Does not apply if \code{data.type} is "Normal".
#' @param n (For binomial data only) vector of integers specifying the number of subjects who have a particular value of the covariate vector. If the data is binary and all covariates are discrete, collapsing Bernoulli data into a binomial structure can make the slice sampler much faster.
#' The sum of \code{n} should be equal to \code{data.size}. The length of \code{n} should be equal to the number of rows of \code{x0}. 
#' @param borrow.treat Logical value indicating whether the historical information is used to inform the treatment effect parameter. The default value is FALSE. If TRUE, the first column of the historical covariate matrix must be the treatment indicator. 
#' If FALSE, the historical covariate matrix must NOT have the treatment indicator, since the historical data is assumed to be from the control group only.  
#' @param treat.assign.prob Probability of being assigned to the treatment group. The default value is 0.5. Only applies if \code{borrow.treat=FALSE}. 
#' @param historical (Optional) list of historical dataset(s). East historical dataset is stored in a list which contains three \emph{named} elements: \code{y0}, \code{x0} and \code{a0}.
#' \itemize{
#' \item \code{y0} is a vector of responses.
#' \item \code{x0} is a matrix of covariates. If \code{borrow.treat} is FALSE (the default), \code{x0} should NOT have the treatment indicator. 
#' If \code{borrow.treat} is TRUE, the first column of \code{x0} must be the treatment indicator.
#' \item \code{a0} is a number between 0 and 1 indicating the discounting parameter value for that historical dataset.
#' }
#' For binomial data, an additional element \code{n0} is required.
#' \itemize{
#' \item \code{n0} is vector of integers specifying the number of subjects who have a particular value of the covariate vector. 
#' The length of \code{n0} should be equal to the number of rows of \code{x0}. 
#' }
#' @param nullspace.ineq Character string specifying the inequality of the null hypothesis. The options are ">" and "<". If ">" is specified, the null hypothesis is \eqn{H_0}: \eqn{\beta_1} \eqn{\ge} \eqn{\delta}. If "<" is specified, the null hypothesis is \eqn{H_0}: \eqn{\beta_1} \eqn{\le} \eqn{\delta}. The default choice is ">".
#' @param x.samples (Only applies when there is no historical dataset) matrix of possible values of covariates from which covariate vectors are sampled with replacement. 
#' @param samp.prior.beta Matrix of possible values of \eqn{\beta} to sample (with replacement) from. Each row is a possible \eqn{\beta} vector (a realization from the sampling prior for \eqn{\beta}), where the first element is the coefficient for the intercept and the second element is the coefficient for the treatment indicator. 
#' The length of the vector should be equal to the total number of parameters. If P is the number of columns of \code{x0} in \code{historical}, the total number of parameters is P+2 if \code{borrow.treat=FALSE}, and is P+1 if \code{borrow.treat=TRUE}. 
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
#' @param approximate Logical value indicating whether the approximation method based on asymptotic theory is used. The default is FALSE. If TRUE, an approximation method based on the Newton-Raphson algorithm (assuming canonical links) is used.
#' This feature helps users quickly obtain a rough estimate of the sample size required for the desired level of power or type I error rate.
#' @param nNR (Only applies if \code{approximate=TRUE}) number of iterations of the Newton-Raphson algorithm. The default value is 10,000.
#' @param tol (Only applies if \code{approximate=TRUE}) absolute tolerance of the Newton-Raphson algorithm. The default value is 0.00001.
#'
#' @details If historical datasets are provided, the algorithm samples with replacement from the historical covariates to construct the simulated datasets.
#' Otherwise, the algorithm samples with replacement from \code{x.samples}. One of the arguments \code{historical} and \code{x.samples} must be provided.
#'
#' The sampling prior for the treatment parameter can be generated from a normal distribution (see examples). 
#' For example, suppose one wants to compute the power for the hypotheses \eqn{H_0: \beta_1 \ge 0} and \eqn{H_1: \beta_1 < 0.} 
#' To approximate the sampling prior for \eqn{\beta_1}, one can simply sample from a normal distribution with negative mean, 
#' so that the mass of the prior falls in the alternative space. Conversely, to compute the type I error rate, one can 
#' sample from a normal distribution with positive mean, so that the mass of the prior falls in the null space. 
#' The sampling prior for the other parameters can be generated by using the \code{glm.fixed.a0} function with \code{current.data} set to FALSE. 
#' The posterior samples based on only historical data can be used as a discrete approximation to the sampling prior. 
#' 
#' \code{samp.prior.var} is necessary for generating normally distributed data.
#'
#' If \code{data.type} is "Normal", the response \eqn{y_i} is assumed to follow \eqn{N(x_i'\beta, \tau^{-1})} where \eqn{x_i} is the vector of covariates for subject \eqn{i}.
#' Each historical dataset \eqn{D_{0k}} is assumed to have a different precision parameter \eqn{\tau_k}.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}, and the initial prior for \eqn{\tau_k} is \eqn{\tau_k^{-1}}.
#' The initial prior for \eqn{\beta} is the uniform improper prior. Posterior samples are obtained through Gibbs sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling. 
#' The default lower limits for the parameters are -100. The default upper limits
#' for the parameters are 100. The default slice widths for the parameters are 1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#'
#' Because running \code{power.glm.fixed.a0()} and \code{power.glm.random.a0()} is potentially time-consuming,
#' an approximation method based on asymptotic theory has been implemented for the model with fixed \eqn{a_0}.
#' In order to attain the exact sample size needed for the desired power, the user can start with the approximation
#' to get a rough estimate of the sample size required, using \code{power.glm.fixed.a0()} with \code{approximate=TRUE}.
#'
#' @return The function returns a S3 object with a \code{summary} method. Power or type I error is returned, depending on the sampling prior used.
#' The posterior probabilities of the alternative hypothesis are returned.  
#' The average posterior mean of \eqn{\beta} and its corresponding bias are returned. 
#' If \code{data.type} is "Normal", average posterior means of \eqn{\tau} and \eqn{\tau_k}'s (if historical data is given) are also returned. 
#' The first column of \eqn{\beta} contains posterior samples of the intercept. The second column contains posterior samples of \eqn{\beta_1}, the parameter for the treatment indicator.
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
#' summary(result)
#'
#' @export
power.glm.fixed.a0 <- function(data.type, data.link="", data.size, n=1, borrow.treat=FALSE, treat.assign.prob=0.5,  
                               historical=list(), nullspace.ineq=">", 
                               x.samples=matrix(), samp.prior.beta, samp.prior.var=0,
                               lower.limits=rep(-100, 50), upper.limits=rep(100, 50),
                               slice.widths=rep(1, 50),
                               delta=0, gamma=0.95, nMC=10000, nBI=250, N=10000, approximate=FALSE, nNR=10000, tol=0.00001) {
  
  if(approximate==TRUE){
    return(power_glm_fixed_a0_approx(data.type, data.size, treat.assign.prob, borrow.treat, 
                                     historical, nullspace.ineq, x.samples,
                                     samp.prior.beta, samp.prior.var,
                                     delta, gamma,
                                     nNR, tol, N))
  }else{
    if(length(historical)!=0){
      x <- historical[[1]]$x0
      if(is.null(colnames(x))){
        colnames(x) <- paste0("X", 1:ncol(x))
      }
    }else{
      if(is.null(colnames(x.samples))){
        colnames(x) <- paste0("X", 1:ncol(x.samples))
      }
    }
    
    out <- power_glm_fixed_a0(data.type, data.link, data.size, n, treat.assign.prob, borrow.treat, historical, nullspace.ineq, x.samples,
                              samp.prior.beta, samp.prior.var, lower.limits, upper.limits, slice.widths, delta, gamma, nMC, nBI, N, TRUE)
    if(borrow.treat==FALSE){
      rownames(out$`average posterior mean of beta`) <- c("intercept","treatment (beta_1)",colnames(x))
    }else{
      rownames(out$`average posterior mean of beta`) <- c("intercept", paste(colnames(x)[1],"(beta_1)"),colnames(x)[-1])
    }
    
    structure(out, class=c("powerglm"))
  }

}

#' @export
summary.powerglm <- function(object, ...) {
  
  r <- round(object$`power/type I error`, digits=3)
  # beta
  betas <- object$`average posterior mean of beta`
  bias <- object$`bias of the average posterior mean of beta`
  postprob <- mean(object[[2]])
  output1 <- cbind(betas,bias)
  colnames(output1) <- c("average posterior mean","bias")
  output1 <- round(output1, digits=3)
  print(output1)
  cat("The power/type I error rate is ",r,".\n")
  cat("The average of the", names(object[2]), "is",round(mean(object[[2]]),3),".")
  
}


