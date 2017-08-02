#' Probability of a relative being affected with familial disease
#'
#' Given a series of pedigrees of different sizes, each with at least r affected relatives of the proband,
#' probAffectedRelative estimates the probability of a relative being affected with familial disease,
#' where this probability is an average over all observed relationships to probands.
#' This probability can be used at the pf parameter in the probFamilial function.
#'
#' Parameter r accounts for ascertainment of pedigrees that are highly likely to segregate familial disease.
#'
#' Each pedigree contributes to the likelihood the binomial probability of m successes in k trials,
#' conditional on at least r successes.  Estimation is then by maximum likelihood.
#'
#' @examples
#' # Familial nonmedullary thyroid cancer, table 2 in Charkes 2006.
#' m=c(2,2,2,4,2,2,2)
#' k=c(9,9,10,7,9,14,7)
#' probAffectedRelative(m,k,2)

#' @references Dudbridge F, Brown SJ, Ward L, Wilson SG, Walsh JP.
#' How many cases of disease in a pedigree imply familial disease?
#' Submitted.
#' @references Charkes ND (2006)
#' On the Prevalence of Familial Nonmedullary Thyroid Cancer in Multiply Affected Kindreds.
#' Thryoid 16:181-186
#'
#' @param m Vector in which each element corresponds to a pedigree and is the number of affected relatives of the proband.
#' @param k Vector in which element corresponds to a pedigree and is the number of relatives of the proband with known affection status.
#' @param r Minimum number of affected relatives in each pedigree.
#'
#' @return Probability of a relative being affected with familial disease.
#' @export
probAffectedRelative <- function(m,k,r) {
  llhd=function(p) {-sum(log(dbinom(m,k,p)/(1-pbinom(r-1,k,p))))}
  pf=optimise(llhd,c(0,1))$min
  pf
}
