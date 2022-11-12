## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE--------------------------------------------------------------
library(kableExtra)
df <- data.frame(Cat = c("Historical Trial 1", "Historical Trial 2"), 
                 Sad = c("8.2% (44/535)", "10.9% (33/304)"))
kable(df, col.names = c("", "% TLF (# of failure/sample size)"), escape = F, caption = "Summary of the historical data.") %>%
  kable_styling(latex_options = "hold_position")


## ---- eval=TRUE---------------------------------------------------------------

historical <- matrix(0, ncol=3, nrow=2)
historical[1,] <- c(44, 535, 0.3)
historical[2,] <- c(33, 304, 0.3)


## ---- eval=TRUE---------------------------------------------------------------

library(BayesPPD)

set.seed(1)
n.t_vals <- seq(from=600, to=1000, by=50)
powers <- NULL

for(i in 1:length(n.t_vals)){
  n.t <- n.t_vals[i]
  results <- power.two.grp.fixed.a0(data.type="Bernoulli", 
      n.t=n.t, n.c=round(n.t/3), historical=historical,
      samp.prior.mu.t=0.092, samp.prior.mu.c=0.092,
      prior.mu.t.shape1=0.0001, prior.mu.t.shape2=0.0001, 
      prior.mu.c.shape1=0.0001,prior.mu.c.shape2=0.0001,
      delta=0.041, N=10000)
  power <- results$`power/type I error`
  powers <- c(powers, power)
}

powers




## ---- eval=TRUE---------------------------------------------------------------
library(ggplot2)

df <- data.frame(sample_size=n.t_vals, power=powers)
ggplot(data=df, aes(x=sample_size, y=powers)) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
  geom_point() +
  xlab("Sample Size") +
  ylab("Power")


## ---- eval=TRUE---------------------------------------------------------------
TIEs <- NULL

for(i in 1:length(n.t_vals)){
  n.t <- n.t_vals[i]
  results <- power.two.grp.fixed.a0(data.type="Bernoulli", 
      n.t=n.t, n.c=round(n.t/3), historical=historical,
      samp.prior.mu.t=0.092+0.041, samp.prior.mu.c=0.092,
      prior.mu.t.shape1=0.0001, prior.mu.t.shape2=0.0001, 
      prior.mu.c.shape1=0.0001,prior.mu.c.shape2=0.0001,
      delta=0.041, N=10000)
  TIE <- results$`power/type I error`
  TIEs <- c(TIEs, TIE)
}

TIEs


## ---- eval=FALSE--------------------------------------------------------------
#  historical <- matrix(0, ncol=2, nrow=2)
#  historical[1,] <- c(44, 535)
#  historical[2,] <- c(33, 304)

## ---- eval=FALSE--------------------------------------------------------------
#  n.t <- 750
#  results <- power.two.grp.random.a0(data.type="Bernoulli",
#        n.t=n.t, n.c=round(n.t/3),historical=historical,
#        samp.prior.mu.t=0.092, samp.prior.mu.c=0.092,
#        prior.mu.t.shape1=0.0001, prior.mu.t.shape2=0.0001,
#        prior.mu.c.shape1=0.0001,prior.mu.c.shape2=0.0001,
#        prior.a0.shape1=1,prior.a0.shape2=1,
#        delta=0.041, gamma=0.95,
#        nMC=10000, nBI=250, N=10000)
#  summary(results)

## ---- eval=TRUE---------------------------------------------------------------
data.type <- "Normal"
n.t <- 100
n.c <- 100

# Simulate three historical datasets
K <- 3
historical <- matrix(0, ncol=4, nrow=K)
# The columns are the sum of the responses, the sample size, the sample variance and a_0
historical[1,] <- c(50, 50, 1, 0.3)
historical[2,] <- c(30, 50, 1, 0.5)
historical[3,] <- c(20, 50, 1, 0.7)

## ---- eval=TRUE---------------------------------------------------------------
# Generate sampling priors
set.seed(1)
samp.prior.mu.t <- rnorm(50000)
samp.prior.mu.c <- rnorm(50000)

sub_ind <- which(samp.prior.mu.t < samp.prior.mu.c) 
# Here, mass is put on the alternative region, so power is calculated. 
samp.prior.mu.t <- samp.prior.mu.t[sub_ind]
samp.prior.mu.c <- samp.prior.mu.c[sub_ind]
samp.prior.var.t <- rgamma(100, 1, 1)
samp.prior.var.c <- rgamma(100, 1, 1)


## ---- eval=TRUE---------------------------------------------------------------
set.seed(1)
results <- power.two.grp.fixed.a0(data.type=data.type, n.t=n.t, n.c=n.t, historical=historical,  
           samp.prior.mu.t=samp.prior.mu.t, samp.prior.mu.c=samp.prior.mu.c, 
           samp.prior.var.t=samp.prior.var.t, samp.prior.var.c=samp.prior.var.c,
           delta=0, nMC=10000, nBI=250, N=100) 

summary(results)

## ---- eval=TRUE---------------------------------------------------------------
# Generate sampling priors
set.seed(1)
samp.prior.mu.t <- rnorm(50000)
samp.prior.mu.c <- rnorm(50000)

sub_ind <- which(samp.prior.mu.t >= samp.prior.mu.c) 
# Here, mass is put on the null region, so type I error rate is calculated. 
samp.prior.mu.t <- samp.prior.mu.t[sub_ind]
samp.prior.mu.c <- samp.prior.mu.c[sub_ind]

set.seed(1)
results <- power.two.grp.fixed.a0(data.type=data.type, n.t=n.t, n.c=n.t, historical=historical,  
           samp.prior.mu.t=samp.prior.mu.t, samp.prior.mu.c=samp.prior.mu.c, 
           samp.prior.var.t=samp.prior.var.t, samp.prior.var.c=samp.prior.var.c,
           delta=0, nMC=10000, nBI=250, N=100) 

summary(results)

## ---- eval=TRUE---------------------------------------------------------------
data.type <- "Normal"
n.t <- 100
n.c <- 100

# Simulate three historical datasets
K <- 3
historical <- matrix(0, ncol=3, nrow=K)
# The columns are the sum of the responses, the sample size, and the sample variance 
historical[1,] <- c(50, 50, 1)
historical[2,] <- c(30, 50, 1)
historical[3,] <- c(20, 50, 1)

## ---- eval=TRUE---------------------------------------------------------------
# Generate sampling priors
set.seed(1)
samp.prior.mu.t <- rnorm(50000)
samp.prior.mu.c <- rnorm(50000)

sub_ind <- which(samp.prior.mu.t < samp.prior.mu.c) 
# Here, mass is put on the alternative region, so power is calculated. 
samp.prior.mu.t <- samp.prior.mu.t[sub_ind]
samp.prior.mu.c <- samp.prior.mu.c[sub_ind]
samp.prior.var.t <- rgamma(100, 1, 1)
samp.prior.var.c <- rgamma(100, 1, 1)


## ---- eval=TRUE---------------------------------------------------------------
set.seed(1)
results <- power.two.grp.random.a0(data.type=data.type, n.t=n.t, n.c=n.t, historical=historical,  
           samp.prior.mu.t=samp.prior.mu.t, samp.prior.mu.c=samp.prior.mu.c, 
           samp.prior.var.t=samp.prior.var.t, samp.prior.var.c=samp.prior.var.c,
           delta=0, nMC=10000, nBI=250, N=100) 

summary(results)
results$`average posterior means of a0`
results$`average posterior mean of tau`

## ---- eval=TRUE---------------------------------------------------------------
# Generate sampling priors
set.seed(1)
samp.prior.mu.t <- rnorm(50000)
samp.prior.mu.c <- rnorm(50000)

sub_ind <- which(samp.prior.mu.t >= samp.prior.mu.c) 
# Here, mass is put on the null region, so type I error rate is calculated. 
samp.prior.mu.t <- samp.prior.mu.t[sub_ind]
samp.prior.mu.c <- samp.prior.mu.c[sub_ind]

set.seed(1)
results <- power.two.grp.random.a0(data.type=data.type, n.t=n.t, n.c=n.t, historical=historical,  
           samp.prior.mu.t=samp.prior.mu.t, samp.prior.mu.c=samp.prior.mu.c, 
           samp.prior.var.t=samp.prior.var.t, samp.prior.var.c=samp.prior.var.c,
           delta=0, nMC=10000, nBI=250, N=100) 

summary(results)

