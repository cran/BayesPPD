
## test that the glm.fixed.a0() function produces similar posterior mean estimates of beta  
## (all absolute differences < 0.02) compared to the results produced by PROC MCMC on SAS. Two 
## historical datasets are used. 
## The outcome is binary. 


library(testthat)

## generate data 

eps <- 1e-2
# Simulate current data
set.seed(1)
data.type <- "Bernoulli"
data.link <- "Logistic"

p <- 3
n_total <- 100
y <- rbinom(n_total,size=1,prob=0.5)
# The first column of x is the treatment indicator.
x <- cbind(rbinom(n_total,size=1,prob=0.5),
           matrix(rnorm(p*n_total),ncol=p,nrow=n_total))
colnames(x) <- c("x1","x2","x3","x4")
df_new <- cbind(y,x,a0=1)


# Simulate two historical datasets
y01 <- rbinom(n_total,size=1,prob=0.5)
x01 <- cbind(0,matrix(rnorm(p*n_total),ncol=p,nrow=n_total))
colnames(x01) <- c("x1","x2","x3","x4")
hist1 <- cbind(y01,x01,a0=0.3)
y02 <- rbinom(n_total,size=1,prob=0.5)
x02 <- cbind(0,matrix(rnorm(p*n_total),ncol=p,nrow=n_total))
colnames(x02) <- c("x1","x2","x3","x4")
hist2 <- cbind(y02,x02,a0=0.6)

df_glm_bern <- rbind(df_new, hist1, hist2)
#write.csv(df_glm_bern, file="df_glm_bern.csv")


# SAS posterior mean estimates: 0.0284, -0.1993, -0.201, -0.0474, -0.0445

sas_results <- c(0.0284, -0.1993, -0.201, -0.0474, -0.0445)
historical <- list(list(y0=y01, x0=x01[,-1], a0=0.3),
                   list(y0=y02, x0=x02[,-1], a0=0.6))
nMC <- 100000 
nBI <- 250

test_glm_mcmc <- function(){
  result <- glm.fixed.a0(data.type=data.type, data.link=data.link, y=y, x=x, historical=historical,
                       nMC=nMC, nBI=nBI)
  expect_true(all(abs(colMeans(result$posterior.samples)-sas_results) < eps))
}

test_that("the glm.fixed.a0() function evaluates correctly for GLMs with binary outcomes and fixed a0 compared to 
          the results produced by PROC MCMC", test_glm_mcmc())



## test that the glm.random.a0() function produces similar posterior mean estimates of beta
## (all differences < 0.02) with a beta(1000, 1000) prior on a0 compared to 
## the glm.fixed.a0() function when a0 = 0.5.
## The outcome is binary. 

## generate data 

set.seed(1)
eps <- 1e-2
data.type <- "Bernoulli"
data.link <- "Logistic"

p <- 3
n_total <- 100
y <- rbinom(n_total,size=1,prob=0.5)
# The first column of x is the treatment indicator.
x <- cbind(rbinom(n_total,size=1,prob=0.5),
           matrix(rnorm(p*n_total),ncol=p,nrow=n_total))
# Simulate two historical datasets
y01 <- rbinom(n_total,size=1,prob=0.5)
x01 <- matrix(rnorm(p*n_total),ncol=p,nrow=n_total)

historical <- list(list(y0=y01, x0=x01))
historical_fixed <- list(list(y0=y01, x0=x01, a0=0.5))

nMC <- 10000
nBI <- 250

test_glm_mcmc_random <- function(){
  
  grid <- as.matrix(seq(0.1,1,by=0.1))
  
  nMC <- 10000 # nMC should be larger in practice
  nBI <- 50
  a0.coefficients <- normalizing.constant(grid=grid, historical=historical,
                                 data.type=data.type, data.link=data.link,
                                 nMC=nMC, nBI=nBI)
  
  result_ran <- glm.random.a0(data.type=data.type, data.link=data.link, y=y, x=x, historical=historical,
                          prior.a0.shape1 = 1000, prior.a0.shape2 = 1000,a0.coefficients = a0.coefficients,
                         nMC=nMC, nBI=nBI)
  beta_ran <- colMeans(result_ran$posterior.samples[[1]])
  result_fix <- glm.fixed.a0(data.type=data.type, data.link=data.link, y=y, x=x, historical=historical_fixed,
                              nMC=nMC, nBI=nBI)
  beta_fixed <- colMeans(result_fix$posterior.samples)
  expect_true(all(abs(beta_ran-beta_fixed) < eps))
}

test_that("the glm.random.a0() function evaluates correctly for GLMs with binary outcomes and random a0", test_glm_mcmc_random())






