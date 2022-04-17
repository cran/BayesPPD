
# Some of the following comments appeared in the source code for 
# Wang YB, Chen MH, Kuo L, Lewis PO (2018). “A New Monte Carlo Method for Estimating Marginal Likelihoods.” Bayesian Analysis, 13(2), 311–333.





calc_marg_l <- function(a0, mcmc, historical, data.type, data.link,init_var){

  sds <- apply(mcmc, 2, stats::sd)
  means <- apply(mcmc, 2, mean)

  # Standardize MCMC samples
  mcmc_scaled <- scale(mcmc, center=TRUE, scale = TRUE)

  # Determine r ("rcover") and K ("ncutslice") for PWK
  rcover <- sqrt(qchisq(0.95,df=ncol(mcmc)))
  ncutslice <- 100


  c_est_pp <- 0

  #Forming the partition subsets for the PWK estimation of c_0
  rings_pp <- LOR_partition_pp(r=rcover, nslice=ncutslice, mcmc=mcmc, a0=a0, historical, data.type, data.link,init_var)
  #Estimate c_0
  t1 <- dim(mcmc_scaled)[1]
  val <- rep(NA, t1)
  r_ <- sqrt(apply(mcmc_scaled^2, 1, sum )) # norm of mcmc samples
  kreprp <- rep(NA, t1)
  for (j in 1:t1 ){
    ring_position  <- which(r_[j]>=rings_pp[,1] & r_[j]<rings_pp[,2])
    if (sum(ring_position)>0){ # phi_t is in A_k

      partjo1 <- log(prod(sds))# log of jacobian

      kreprp[j] <- 0
      kreprp[j] <- kreprp[j] + logpowerprior(mcmc[j,], a0, historical, data.type, data.link,init_var) + partjo1
      val[j] <- rings_pp[ring_position, 3] - kreprp[j]

    }
  }
  lcons <- max(rings_pp[, 4])
  tcubevol <- log(sum(exp(rings_pp[,4] - lcons) ) ) + lcons 

  den <- denominator_control(val)
  c_est_pp <- NA
  c_est_pp <- tcubevol - den
  return(c_est_pp) # log(c0)

}


# This function computes the normalizing constant for a given a_0 value. 
calc_a0_func <- function(a0, historical, data.type, data.link,
                         init_var,lower_limits,upper_limits,slice_widths,nMC,nBI){
  dCurrent <- FALSE
  historical2 <- list()
  for(i in 1:length(historical)){
    l <- historical[[i]]
    l[["a0"]] = a0[i]
    historical2[[i]] = l
  }
  y <- rep(0,1)
  x <- matrix(0, nrow=1, ncol=1)
  n <- rep(0, 1)
  samples <- glm_fixed_a0(data.type,data.link,y,n,x,historical2,init_var,lower_limits,upper_limits,slice_widths,nMC,nBI,dCurrent)

  marg_l <- calc_marg_l(a0, samples, historical, data.type, data.link,init_var)


  return(marg_l)
}



#' Function for approximating the normalizing constant for generalized linear models with random a0
#'
#'
#'
#'
#' @description This function returns a vector of coefficients that defines a function \eqn{f(a_0)} that approximates the normalizing constant for generalized linear models with random \eqn{a_0}.
#' The user should input the values returned to \code{\link{glm.random.a0}} or \code{\link{power.glm.random.a0}}.
#'
#' @param grid Matrix of potential values for \eqn{a_0}, where the number of columns should equal the number of historial datasets. Note that the algorithm may fail if some grid values are close to zero. See \emph{Details} below.
#' @param historical List of historical dataset(s). East historical dataset is stored in a list which constains two \emph{named} elements: \code{y0} and \code{x0}.
#' \itemize{
#' \item \code{y0} is a vector of responses.
#' \item \code{x0} is a matrix of covariates. \code{x0} should NOT have the treatment/control group indicator. Apart from missing the treatent/control indicator, \code{x0} should have the same set of covariates in the same order as \code{x}.
#' }
#' For binomial data, an additional element \code{n0} is required. 
#' \itemize{
#' \item \code{n0} is vector of integers specifying the number of subjects who have a particular value of the covariate vector.
#' }
#' @param data.type Character string specifying the type of response. The options are "Bernoulli", "Binomial", "Poisson" and "Exponential".
#' 
#' @inheritParams power.glm.fixed.a0
#' @inheritParams glm.random.a0
#'
#' @details
#'
#' This function performs the following steps:
#'
#' \enumerate{
#'
#' \item	Suppose there are K historical datasets. The user inputs a grid of M rows and K columns of potential values for \eqn{a_0}. For example, one can choose the vector \code{v = c(0.1, 0.25, 0.5, 0.75, 1)}
#' and use \code{expand.grid(a0_1=v, a0_2=v, a0_3=v)} when \eqn{K=3} to get a grid with \eqn{M=5^3=125} rows and 3 columns. If there are more than three historical datasets, the dimension of \code{v} can be reduced
#' to limit the size of the grid. A large grid will increase runtime.
#' \item	For each row of \eqn{a_0} values in the grid, obtain \eqn{M} samples for \eqn{\beta} from the power prior associated with the current values of \eqn{a_0} using the slice sampler.
#' \item	For each of the M sets of posterior samples, execute the PWK algorithm (Wang et al., 2018) to estimate the log of normalizing constant \eqn{d_1,...,d_M} for the normalized power prior.
#' \item	At this point, one has a dataset with outcomes \eqn{d_1,...,d_M} and predictors corresponding to the rows of the \eqn{a_0} grid matrix. A polynomial regression is applied to estimate a function \eqn{d=f(a0)}.
#' The degree of the polynomial regression is determined by the algorithm to ensure \eqn{R^2 > 0.99}.
#' \item	The vector of coefficients from the polynomial regression model is returned by the function, which the user must input into \code{\link{glm.random.a0}} or \code{\link{power.glm.random.a0}}.
#'
#' }
#'
#' When a row of the \code{grid} contains elements that are close to zero, the resulting power prior will be flat and estimates of normalizing constants may be inaccurate.
#' Therefore, it is recommended that \code{grid} values should be at least 0.05.
#'
#' If one encounters the error message "some coefficients not defined because of singularities",
#' it could be due to the following factors: number of \code{grid} rows too large or too small, insufficient sample size of the historical data, insufficient number of iterations for the slice sampler,
#' or near-zero \code{grid} values.
#'
#'
#' @return Vector of coefficients for \eqn{a_0} that defines a function \eqn{f(a_0)} that approximates the normalizing constant, necessary for functions \code{\link{glm.random.a0}} and \code{\link{power.glm.random.a0}}.
#' The length of the vector is equal to 1+K*L where K is the number of historical datasets and L is the degree of the polynomial regression determined by the algorithm. 
#' @references Wang, Yu-Bo; Chen, Ming-Hui; Kuo, Lynn; Lewis, Paul O. A New Monte Carlo Method for Estimating Marginal Likelihoods. Bayesian Anal. 13 (2018), no. 2, 311--333.
#' @seealso \code{\link{glm.random.a0}} and \code{\link{power.glm.random.a0}}
#' @examples
#'
#' data.type <- "Bernoulli"
#' data.link <- "Logistic"
#' data.size <- 50
#'
#' # Simulate two historical datasets
#' p <- 1
#' set.seed(111)
#' x1 <- matrix(rnorm(p*data.size),ncol=p,nrow=data.size)
#' set.seed(222)
#' x2 <- matrix(rnorm(p*data.size),ncol=p,nrow=data.size)
#' beta <- c(1,2)
#' mean1 <- exp(x1*beta)/(1+exp(x1*beta))
#' mean2 <- exp(x2*beta)/(1+exp(x2*beta))
#' historical <- list(list(y0=rbinom(data.size,size=1,prob=mean1),x0=x1),
#'                    list(y0=rbinom(data.size, size=1, prob=mean2),x0=x2))
#'
#' # Create grid of possible values of a0 with two columns corresponding to a0_1 and a0_2
#' g <- c(0.1, 0.25, 0.5, 0.75, 1)
#' grid <- expand.grid(a0_1=g, a0_2=g)
#'
#' nMC <- 100 # nMC should be larger in practice
#' nBI <- 50
#' result <- normalizing.constant(grid=grid, historical=historical,
#'                                data.type=data.type, data.link=data.link,
#'                                nMC=nMC, nBI=nBI)
#' @importFrom stats sd qchisq pnorm dnorm lm
#' @export
normalizing.constant <- function(grid, historical, data.type, data.link,
                                 prior.beta.var=rep(10,50), lower.limits=rep(-100, 50), upper.limits=rep(100, 50), slice.widths=rep(1, 50),
                                 nMC=10000, nBI=250){


  d <- apply(grid, 1, calc_a0_func, historical, data.type, data.link,
             prior.beta.var,lower.limits,upper.limits,slice.widths,nMC,nBI) # not a number when a0 too low

  m <- as.matrix(grid)
  r2 <- 0
  degree <- 0
  while(r2 < 0.99){
    degree = degree + 1
    mat <- matrix(0, nrow=nrow(m), ncol=ncol(m)*degree)

    # create m degree polynomials
    for(i in 1:degree){
      mat[,((i-1)*ncol(m)+1):(ncol(m)*i)] <- m^i
    }
    fit <- lm(d ~ mat)
    r2 <- summary(fit)$r.squared

    if(NA %in% fit$coefficients){
      stop("some coefficients not defined because of singularities. Please adjust the grid.")
    }
  }
  
  result <- fit$coefficients

  #return(list("coef"=result, "logc"=d))
  return(result)
}

