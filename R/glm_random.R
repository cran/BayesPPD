### glm, random a0 #####



#' Model fitting for generalized linear models with random a0
#'
#' @description Model fitting using normalized power priors for generalized linear models with random \eqn{a_0}
#'
#' @param historical List of historical dataset(s). East historical dataset is stored in a list which contains two \emph{named} elements: \code{y0} and \code{x0}.
#' \itemize{
#' \item \code{y0} is a vector of responses.
#' \item \code{x0} is a matrix of covariates. If \code{borrow.treat} is FALSE (the default), \code{x0} should NOT have the treatment indicator. Apart from missing the treatment indicator, \code{x0} should have the same set of covariates in the same order as \code{x}.
#' If \code{borrow.treat} is TRUE, \code{x0} should have the same set of covariates in the same order as \code{x}, where the first column of \code{x0} must be the treatment indicator.
#' }
#' For binomial data, an additional element \code{n0} is required.
#' \itemize{
#' \item \code{n0} is vector of integers specifying the number of subjects who have a particular value of the covariate vector. 
#' The length of \code{n0} should be equal to the number of rows of \code{x0}. 
#' }
#' @param prior.beta.var Vector of variances of the independent normal initial priors on \eqn{\beta} with mean zero. The length of the vector should be equal to the length of \eqn{\beta}. The default variance is 10.   
#' @param lower.limits Vector of lower limits for parameters to be used by the slice sampler. If \code{data.type} is "Normal", slice sampling is used for \eqn{a_0}, and the length of the vector should be equal to the number of historical datasets.
#' For all other data types, slice sampling is used for \eqn{\beta} and \eqn{a_0}. The first P+1 elements apply to the sampling of \eqn{\beta} and the rest apply to the sampling of \eqn{a_0}.
#' The length of the vector should be equal to the sum of the total number of parameters (i.e. P+1 where P is the number of covariates) and the number of historical datasets.
#'  The default is -100 for \eqn{\beta} and 0 for \eqn{a_0} (may not be appropriate for all situations).
#' @param upper.limits Vector of upper limits for parameters to be used by the slice sampler. If \code{data.type} is "Normal", slice sampling is used for \eqn{a_0}, and the length of the vector should be equal to the number of historical datasets.
#' For all other data types, slice sampling is used for \eqn{\beta} and \eqn{a_0}. The first P+1 elements apply to the sampling of \eqn{\beta} and the rest apply to the sampling of \eqn{a_0}.
#' The length of the vector should be equal to the sum of the total number of parameters (i.e. P+1 where P is the number of covariates) and the number of historical datasets.
#'  The default is 100 for \eqn{\beta} and 1 for \eqn{a_0}  (may not be appropriate for all situations).
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
#' Independent normal priors with mean zero and variance \code{prior.beta.var} are used for \eqn{\beta} to ensure the propriety of the normalized power prior. Posterior samples for \eqn{\beta} and \eqn{\tau} are obtained through Gibbs sampling.
#' Independent beta(\code{prior.a0.shape1}, \code{prior.a0.shape1}) priors are used for \eqn{a_0}. Posterior samples for \eqn{a_0} are obtained through slice sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling.
#' The default lower limits are -100 for \eqn{\beta} and 0 for \eqn{a_0}. The default upper limits
#' for the parameters are 100 for \eqn{\beta} and 1 for \eqn{a_0}. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#'
#' @return The function returns a S3 object with a \code{summary} method. If \code{data.type} is "Normal", posterior samples of \eqn{\beta}, \eqn{\tau} and \eqn{a_0} are returned.
#' For all other data types, posterior samples of \eqn{\beta} and \eqn{a_0} are returned.
#' The first column of the matrix of posterior samples of \eqn{\beta} contains posterior samples of the intercept.
#' The second column contains posterior samples of \eqn{\beta_1}, the parameter for the treatment indicator.
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
#' # Here, suppose one-degree polynomial regression is chosen by the "normalizing.constant"
#' # function. The coefficients are obtained for the intercept, a0_1 and a0_2.
#' a0.coefficients <- c(1, 0.5, -1)
#'
#' # Set parameters of the slice sampler
#' # The dimension is the number of columns of x plus 1 (intercept)
#' # plus the number of historical datasets
#' lower.limits <- c(rep(-100, 5), rep(0, 2))
#' upper.limits <- c(rep(100, 5), rep(1, 2))
#' slice.widths <- rep(0.1, 7)
#'
#' nMC <- 500 # nMC should be larger in practice
#' nBI <- 100
#' result <- glm.random.a0(data.type=data.type, data.link=data.link, y=y, x=x,
#'                         historical=historical, a0.coefficients=a0.coefficients,
#'                         lower.limits=lower.limits, upper.limits=upper.limits,
#'                         slice.widths=slice.widths, nMC=nMC, nBI=nBI)
#' summary(result)
#'
#' @export
glm.random.a0 <- function(data.type, data.link, y, x, n=1, borrow.treat=FALSE, historical, prior.beta.var=rep(10,50),
                          prior.a0.shape1=rep(1,10), prior.a0.shape2=rep(1,10), a0.coefficients,
                          lower.limits=NULL, upper.limits=NULL,
                          slice.widths=rep(0.1, 50), nMC=10000, nBI=250) {
  

  if(is.null(colnames(x))){
    colnames(x) <- paste0("X", 1:ncol(x))
  }
  if(data.type == "Normal"){
    if(is.null(lower.limits)){
      lower.limits = rep(0,length(historical))
    }
    if(is.null(upper.limits)){
      upper.limits = rep(1,length(historical))
    }
    result <- glm_random_a0_normal(y, x, borrow.treat, historical, prior.a0.shape1, prior.a0.shape2, lower.limits, upper.limits, slice.widths, nMC, nBI)
  }else{
    if(is.null(lower.limits)){
      lower.limits = c(rep(-100, ncol(x)+1), rep(0,length(historical)))
    }
    if(is.null(upper.limits)){
      upper.limits = c(rep(100, ncol(x)+1), rep(1,length(historical)))
    }
    result <- glm_random_a0(data.type, data.link, y, n, x, borrow.treat, historical, prior.beta.var, 
                         prior.a0.shape1, prior.a0.shape2, a0.coefficients,
                         lower.limits, upper.limits,
                         slice.widths, nMC, nBI)
  }
  colnames(result$`posterior samples of beta`) <- c("intercept",colnames(x))
  out <- list(posterior.samples=result, data.type=data.type, x=x)
  structure(out, class=c("glmrandom"))
}


#' @importFrom stats quantile
#' @export
summary.glmrandom <- function(object, ...) {
  
  r <- object$posterior.samples 
    
  # beta
  betas <- r$`posterior samples of beta`
  m <- apply(betas, 2, mean)
  std <- apply(betas, 2, sd)
  q1 <- apply(betas, 2, quantile, probs=0.025)
  q2 <- apply(betas, 2, quantile, probs=0.975)
  output1 <- cbind("mean"=m,"sd"=std,"2.5%"=q1,"97.5%"=q2)
  if(is.null(colnames(object$x))){
    colnames(object$x) <- paste0("V", 1:ncol(object$x))
  }
  rownames(output1) <- c("intercept",colnames(object$x))
  output1 <- round(output1, digits=3)
  
  # tau
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






#' Power/type I error calculation for generalized linear models with random a0
#'
#' @description Power/type I error calculation using normalized power priors for generalized linear models with random \eqn{a_0}
#'
#' @param historical List of historical dataset(s). East historical dataset is stored in a list which contains two \emph{named} elements: \code{y0} and \code{x0}.
#' \itemize{
#' \item \code{y0} is a vector of responses.
#' \item \code{x0} is a matrix of covariates. If \code{borrow.treat} is FALSE (the default), \code{x0} should NOT have the treatment indicator. 
#' If \code{borrow.treat} is TRUE, the first column of \code{x0} must be the treatment indicator.
#' }
#' For binomial data, an additional element \code{n0} is required.
#' \itemize{
#' \item \code{n0} is vector of integers specifying the number of subjects who have a particular value of the covariate vector. 
#' The length of \code{n0} should be equal to the number of rows of \code{x0}. 
#' }
#' @inheritParams glm.random.a0
#' @inheritParams power.glm.fixed.a0
#' @inheritParams two.grp.random.a0
#' @details  The user should use the function \code{\link{normalizing.constant}} to obtain \code{a0.coefficients} (does not apply if \code{data.type} is "Normal").
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
#' Historical datasets are assumed to have the same precision parameter as the current dataset for computational simplicity.
#' The initial prior for \eqn{\tau} is the Jeffery's prior, \eqn{\tau^{-1}}.
#' Independent normal priors with mean zero and variance \code{prior.beta.var} are used for \eqn{\beta} to ensure the propriety of the normalized power prior. Posterior samples for \eqn{\beta} and \eqn{\tau} are obtained through Gibbs sampling.
#' Independent beta(\code{prior.a0.shape1}, \code{prior.a0.shape1}) priors are used for \eqn{a_0}. Posterior samples for \eqn{a_0} are obtained through slice sampling.
#'
#' For all other data types, posterior samples are obtained through slice sampling. 
#' The default lower limits are -100 for \eqn{\beta} and 0 for \eqn{a_0}. The default upper limits
#' for the parameters are 100 for \eqn{\beta} and 1 for \eqn{a_0}. The default slice widths for the parameters are 0.1.
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
#' The average posterior mean of \eqn{a_0} is returned.
#' If \code{data.type} is "Normal", the average posterior mean of \eqn{\tau} is also returned. 
#' The first element of the average posterior means of \eqn{\beta} is the average posterior mean of the intercept.
#' The second element is the average posterior mean of \eqn{\beta_1}, the parameter for the treatment indicator.
#' 
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
#' # Please see function "normalizing.constant" for how to obtain a0.coefficients
#' # Here, suppose one-degree polynomial regression is chosen by the "normalizing.constant"
#' # function. The coefficients are obtained for the intercept, a0_1 and a0_2.
#' a0.coefficients <- c(1, 0.5, -1)
#'
#' nMC <- 100 # nMC should be larger in practice
#' nBI <- 50
#' N <- 3 # N should be larger in practice
#' result <- power.glm.random.a0(data.type=data.type, data.link=data.link,
#'                               data.size=data.size, historical=historical,
#'                               samp.prior.beta=samp.prior.beta, a0.coefficients=a0.coefficients,
#'                               delta=0, nMC=nMC, nBI=nBI, N=N)
#' summary(result)
#'
#' @export
power.glm.random.a0 <- function(data.type, data.link, data.size, n=1, treat.assign.prob=0.5, borrow.treat=FALSE, historical,nullspace.ineq=">",
                                samp.prior.beta, samp.prior.var, prior.beta.var=rep(10,50),
                                prior.a0.shape1=rep(1,10), prior.a0.shape2=rep(1,10), a0.coefficients,
                                lower.limits=NULL, upper.limits=NULL,slice.widths=rep(0.1, 50),
                                delta=0, gamma=0.95, nMC=10000, nBI=250, N=10000) {
  

  x <- historical[[1]]$x0
  if(is.null(colnames(x))){
    colnames(x) <- paste0("X", 1:ncol(x))
  }
  
  if(data.type == "Normal"){
    if(is.null(lower.limits)){
      lower.limits = rep(0,length(historical))
    }
    if(is.null(upper.limits)){
      upper.limits = rep(1,length(historical))
    }
    out <- power_glm_random_a0_normal(data.size, treat.assign.prob, borrow.treat, historical,nullspace.ineq,
                                      samp.prior.beta, samp.prior.var,
                                      prior.a0.shape1, prior.a0.shape2,
                                      lower.limits, upper.limits, slice.widths,
                                      delta, gamma, nMC, nBI, N)
  }else{
    if(is.null(lower.limits)){
      lower.limits = c(rep(-100, ncol(samp.prior.beta)), rep(0,length(historical)))
    }
    if(is.null(upper.limits)){
      upper.limits = c(rep(100, ncol(samp.prior.beta)), rep(1,length(historical)))
    }
    out <- power_glm_random_a0(data.type, data.link, data.size, n, treat.assign.prob, borrow.treat, historical,nullspace.ineq,
                               samp.prior.beta, prior.beta.var, prior.a0.shape1, prior.a0.shape2,
                               a0.coefficients, lower.limits, upper.limits, slice.widths,
                               delta, gamma, nMC, nBI, N)
  }
  if(borrow.treat==FALSE){
    rownames(out$`average posterior mean of beta`) <- c("intercept","treatment (beta_1)",colnames(x))
  }else{
    rownames(out$`average posterior mean of beta`) <- c("intercept", paste(colnames(x)[1],"(beta_1)"),colnames(x)[-1])
  }

  structure(out, class=c("powerglm"))
  
}



