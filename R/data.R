
#' AIDS Clinical Trial ACTG019 (1990).
#'
#' A dataset containing the ACTG019 clinical trial placebo group data (1990) in adults with asymptomatic HIV.
#'
#' @format A data frame with 404 rows and 4 variables:
#' \describe{
#'   \item{outcome}{binary variable with 1 indicating death, development of AIDS or ARC and 0 otherwise}
#'   \item{age}{patient age in years}
#'   \item{race}{binary variable with 1 indicating white and 0 otherwise}
#'   \item{T4count}{CD4 cell count (cell count per cubicmillimetre of serum)}
#' }
#' @source Chen, Ming-Hui, et al. "Prior Elicitation, Variable Selection and Bayesian Computation for Logistic Regression Models." Journal of the Royal Statistical Society. Series B, vol. 61, no. 1, 1999, pp. 223-242.
"actg019"


#' AIDS Clinical Trial ACTG036 (1991).
#'
#' A dataset containing the ACTG036 clinical trial data (1991) comparing zidovudine (AZT) with a placebo in asymptomatic patients with hereditary coagulation disorders and HIV infection.
#' The ACTG036 trial had the same response variable and covariates as the ACTG019 study. The ATCG019 data can be used as a historical dataset.
#'
#' 
#'
#' @format A data frame with 183 rows and 5 variables:
#' \describe{
#'   \item{outcome}{binary variable with 1 indicating death, development of AIDS or ARC and 0 otherwise}
#'   \item{treat}{binary variable with 1 indicating Zidovudine (AZT) treatment and 0 indicating placebo}
#'   \item{age}{patient age in years}
#'   \item{race}{binary variable with 1 indicating white and 0 otherwise}
#'   \item{T4count}{CD4 cell count (cell count per cubicmillimetre of serum)}
#' }
#' @source Chen, Ming-Hui, et al. "Prior Elicitation, Variable Selection and Bayesian Computation for Logistic Regression Models." Journal of the Royal Statistical Society. Series B, vol. 61, no. 1, 1999, pp. 223-242.
"actg036"