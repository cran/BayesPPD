
## replicate the non-inferiority design application of Chen et al. (2011) 
## test that the power.two.grp.fixed.a0() function produces similar power and type I error
## calculations (absolute difference < 0.01) compared to Table 2 in Chen et al. (2011). 
## The outcome is binary. 



library(testthat)

set.seed(1)

eps <- 1e-2
powers_chen <- c(0.84, 0.856, 0.884)
TIEs_chen <- c(0.03, 0.027, 0.028)

n.ts <- c(750, 810, 900)
n.cs <- c(250, 270, 300)
N <- 10000
historical <- matrix(0, ncol=3, nrow=2)
historical[1,] <- c(44, 535, 0.3)
historical[2,] <- c(33, 304, 0.3)

delta <- 0.041
powers_ppd <- NULL
TIEs_ppd <- NULL

test_tg_power <- function(){
  for(i in 1:length(powers_chen)){
    
    
    power <- power.two.grp.fixed.a0(data.type="Bernoulli",
                                    n.t=n.ts[i], n.c=n.cs[i], historical=historical,
                                    samp.prior.mu.t=0.092, samp.prior.mu.c=0.092,
                                    prior.mu.t.shape1=0.0001, prior.mu.t.shape2=0.0001,
                                    prior.mu.c.shape1=0.0001,prior.mu.c.shape2=0.0001,
                                    delta=delta, N=N)
    powers_ppd <- c(powers_ppd, power$`power/type I error`)
    TIE <- power.two.grp.fixed.a0(data.type="Bernoulli",
                                  n.t=n.ts[i], n.c=n.cs[i], historical=historical,
                                  samp.prior.mu.t=0.092+delta, samp.prior.mu.c=0.092,
                                  prior.mu.t.shape1=0.0001, prior.mu.t.shape2=0.0001,
                                  prior.mu.c.shape1=0.0001,prior.mu.c.shape2=0.0001,
                                  delta=delta, N=N)
    TIEs_ppd <- c(TIEs_ppd, TIE$`power/type I error`)
    
    
  }
  
  power_diff <- abs(powers_ppd - powers_chen)
  TIE_diff <- abs(TIEs_ppd - TIEs_chen)
  expect_true(all(power_diff < eps))
  expect_true(all(TIE_diff < eps))
  
}


test_that("power / type I error calculation evaluates correctly for two group cases with binary outcomes and fixed a0", test_tg_power())




## test that the power.two.grp.random.a0() function produces similar power and type I error
## calculations (absolute difference < 0.02) with a beta(1000, 1000) prior on a0 compared to 
## the power.two.grp.fixed.a0() function when a0 = 0.5.
## The outcome is normally distributed. 


set.seed(1)
N <- 10000
eps <- 0.02
historical <- matrix(0, ncol=4, nrow=1)
historical[1,] <- c(100, 100, 1, 0.5)
historical_random <- matrix(0, ncol=3, nrow=1)
historical_random[1,] <- c(100, 100, 1)
n.t <- n.c <- 100
delta <- 0

test_tg_power_random <- function(){

    skip("This test is skipped due to its slow speed.")
    
    power_ran <- power.two.grp.random.a0(data.type="Normal",
                                    n.t=n.t, n.c=n.c, historical=historical_random,
                                    samp.prior.mu.t=1, samp.prior.mu.c=1.5,
                                    samp.prior.var.t = 1, samp.prior.var.c = 1,
                                    prior.a0.shape1 = 1000, prior.a0.shape2 = 1000,
                                    delta=delta, N=N)

    power_fixed <- power.two.grp.fixed.a0(data.type="Normal",
                                  n.t=n.t, n.c=n.c, historical=historical,
                                  samp.prior.mu.t=1, samp.prior.mu.c=1.5,
                                  samp.prior.var.t = 1, samp.prior.var.c = 1,
                                  delta=delta, N=N)
    power_diff <- abs(power_ran$`power/type I error` - power_fixed$`power/type I error`)
     
    expect_true(power_diff < eps)
    
    
}
  
test_that("power / type I error calculation evaluates correctly for two group cases with normal outcomes and random a0", test_tg_power_random())






