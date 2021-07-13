
#' AIDS Clinical Trial ACTG019 (1990).
#'
#' A dataset containing the ACTG019 clinical trial data (1990) comparing zidovudine (AZT) with a placebo in adults with asymptomatic HIV.
#'
#' @format A data frame with 822 rows and 8 variables:
#' \describe{
#'   \item{outcome}{binary variable with 1 indicating death, development of AIDS or ARC and 0 otherwise}
#'   \item{treat}{binary variable with 1 indicating Zidovudine (AZT) treatment and 0 indicating placebo}
#'   \item{age}{patient age in years}
#'   \item{race}{binary variable with 1 indicating white and 0 otherwise}
#'   \item{T4count}{CD4 cell count (cell count per cubicmillimetre of serum)}
#'   \item{sex}{binary variable with 1 indicating male and 0 indicating female}
#'   \item{sex_pref}{sexual orientation. Binary variable with 1 indicating homosexual/bisexual and 0 otherwise}
#'   \item{ivdu}{IV drug use. Binary variable with 1 indicating current or past user and 0 indicating never used}
#' }
#' @source Chen, Ming-Hui, et al. "Prior Elicitation, Variable Selection and Bayesian Computation for Logistic Regression Models." Journal of the Royal Statistical Society. Series B, vol. 61, no. 1, 1999, pp. 223-242.
"actg019"


#' AIDS Clinical Trial ACTG036 (1991).
#'
#' A dataset containing the ACTG036 clinical trial data (1991) comparing zidovudine (AZT) with a placebo in asymptomatic patients with hereditary coagulation disorders and HIV infection.
#' The ACTG036 trial had the same response variable and had many covariates in common with the ACTG019 study. The ATCG019 data can be used as a historical dataset.
#'
#' 
#'
#' @format A data frame with 183 rows and 7 variables:
#' \describe{
#'   \item{outcome}{binary variable with 1 indicating death, development of AIDS or ARC and 0 otherwise}
#'   \item{treat}{binary variable with 1 indicating Zidovudine (AZT) treatment and 0 indicating placebo}
#'   \item{age}{patient age in years}
#'   \item{race}{binary variable with 1 indicating white and 0 otherwise}
#'   \item{T4count}{CD4 cell count (cell count per cubicmillimetre of serum)}
#'   \item{hemtype}{haemophilia factor type. Binary variable with 1 indicating the haemophilia factor type was factor VIII and 0 otherwise}
#'   \item{mfc}{monoclonal factor concentrate use. Binary variable with 1 indicating the patient used a monoclonal factor concentrate and 0 otherwise}
#' }
#' @source Chen, Ming-Hui, et al. "Prior Elicitation, Variable Selection and Bayesian Computation for Logistic Regression Models." Journal of the Royal Statistical Society. Series B, vol. 61, no. 1, 1999, pp. 223-242.
"actg036"