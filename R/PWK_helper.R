
# The following comments appeared in the source code for
# Wang YB, Chen MH, Kuo L, Lewis PO (2018). “A New Monte Carlo Method for Estimating Marginal Likelihoods.” Bayesian Analysis, 13(2), 311–333.

library(stats)
logpowerprior <- function(mcmc, a0, historical, data.type, data.link,init_var){
  beta = mcmc
  lp = 0;

  # add independent normal priors for beta
  for(i in 1:length(beta)){
    lp = lp + dnorm(beta[i], mean=0, sd=sqrt(init_var[i]), log=TRUE)
  }

  for(i in 1:length(historical)){
      dat = historical[[i]]
      y_h = dat[["y0"]]
      x_h = dat[["x0"]]
      x_h = cbind(1,x_h)

      if (data.type=="Bernoulli") {n_h = rep(1,length(y_h))}
      if (data.type=="Binomial") {n_h = dat[["n"]]}

      a0_i = a0[i]

      mean = x_h%*%beta


      if (data.link=="Logistic") 		    { mean = exp(mean) / (1 + exp(mean)) 	}
      if (data.link=="Probit") 	   		  {  mean = pnorm(mean, 0.0, 1.0)  }
      if (data.link=="Log") 	   		    { mean = exp(mean)				    				}
      if (data.link=="Identity-Positive")
      { 	for (j in 1:length(mean)) { mean[j] = max(mean[j],0.00001)}						}
      if (data.link=="Identity-Probability")
      { 	for (j in 1:length(mean)) { mean[j] = min(max(mean[j],0.00001),0.99999) }	}
      if (data.link=="Complementary Log-Log") {  mean = 1 - exp(-exp(mean)) }



      #mean = min(max(mean, 10^(-4)), 0.9999)
      if (data.type=="Bernoulli"|data.type=="Binomial"){
        lp = lp + a0_i * sum(y_h * log(mean) + (n_h - y_h) * log1p(-mean) )


      }
      if (data.type=="Poisson")	{
        lp = lp + a0_i * sum(y_h * log(mean) - mean)
      }
      if (data.type=="Exponential")	{
        lp = lp + a0_i * sum(log(mean) - y_h * mean)
      }

  }


  return(lp)
}


#Purpose: floating control for summing all "ratios"
#Input: "ratios" are a vector saving all log(q(\theta^*_k)/q(\theta_t) ) for PWK
#Output: sum of "ratios" divided by an MC/MCMC sample size in log scale
denominator_control <- function(ratios)
{
  tot <- length(ratios)
  b <- ratios[!is.na(ratios)]

  est_d_r <- 0
  b.max <- max(b)
  est_d_r <- log(sum(exp(b-b.max)))+b.max - log(tot)
  return(est_d_r)
}



#Purpose: forming the rings for the estimation of c_0
#Input: "r" denotes the maximum radius to form the working parameter space, and "nslice" is the number of partition subsets
#Ouput: a matrix recording the interval of each ring,
#the posterior kernel (log-scale) of the representative point in each ring,
#and the volume of each ring

LOR_partition_pp <- function(r,nslice,mcmc,a0,historical, data.type, data.link,init_var){
  interval <- seq(0, r, length=(nslice+1) )
  rings <- cbind(interval[-(nslice+1)],interval[-1])
  P <- ncol(mcmc)
  reprp <- apply(rings, 1, mean)/sqrt(P) # square root of number of parameters
  sds <- apply(mcmc, 2, sd)
  means <- apply(mcmc, 2, mean)

  partjo1 <- log(prod(sds))
  kreprp <- rep(NA, nslice)
  for (i in 1:nslice ){
    rpp <- means+sds*reprp[i]
    kreprp[i] <- 0
    kreprp[i] <- kreprp[i] + logpowerprior(rpp, a0, historical, data.type, data.link,init_var) # not plugging in actual mcmc

    kreprp[i] <- kreprp[i] + partjo1
  }

  rings <- cbind(rings, kreprp)
  rarea <- pi^(P/2)*interval^P/gamma(P/2+1)  # pi^(p/2)/gamma(p/2+1) 0.4
  rvol <- log(rarea[-1]-rarea[-(nslice+1)]) + kreprp

  rings <- cbind(rings, rvol)
  return(rings)
}



