
# BayesPPD 1.1.2 - Nov 25, 2023

## Bug Fixes

* Fix bug for glm.random.a0() and power.glm.random.a0() for binomial data.
* Update default values for lower.limits and upper.limits for glm.random.a0() and power.glm.random.a0().



# BayesPPD 1.1.1 - May 8, 2023

## Enhancements

* Update link to the package review article published in R Journal in DESCRIPTION.
* Update syntax that is deprecated by RcppArmadillo.



# BayesPPD 1.1.0 - Nov 11, 2022

## Enhancements

* Allow borrowing for treatment effect parameter for GLMs.
* Return values are S3 objects with `summary` methods implemented. 
* Add vignette for two group cases for binary and normal outcomes.
* Add automated testing. 
* Expand the return values of the power calculation functions to include the posterior probabilities of the alternative hypothesis and the biases of the average posterior mean estimates. 



# BayesPPD 1.0.7 - Oct 12, 2022

## Enhancements
* Allow the user to specify the direction of the null hypothesis. 
