#' Probability that a pedigree is segregating familial disease
#'
#' Given k relatives of the proband, among which m are affected, probFamilial calculates the probability that the pedigree is segregating familial disease as opposed to simply having m sporadic cases.
#'
#' pf can be understood as a combination of penetrance and relatedness.  For example, considering 1st degree relatives and a fully penetrant dominant mutation,
#' pf=0.5 since a full sibling of the proband has 0.5 prob of inheriting the same mutation, a child of the proband has
#' 0.5 prob of receiving the mutation, and a parent of the proband has 0.5 prob of being the one who transmitted the mutation.
#' For general pedigrees containing a mix of relationships, pf can be understood as the average over all observed relatives to the proband.
#'
#'
#' For a rare disease, pr is approximately the population lifetime risk.
#'
#' If parameter r equals 0, the returned value allFamilial is the probability that all m cases have the familial form of disease.
#' Otherwise allFamilial is the probability that at least (m-r) of the m cases have the familial form.
#'
#' @examples
#' # Familial nonmedullary thyroid cancer, figure 1 in Dudbridge et al.
#' pf <- 0.0981
#' pr <- 0.0031
#' priorF <- 0.013
#' for(k in 2:8)
#'   print(probFamilial(1,k,pf,pr,priorF)$probFamilial)
#' for(k in 2:8)
#'   print(probFamilial(2,k,pf,pr,priorF)$probFamilial)
#'
#' # Colorectal cancer, probabilities that all cases are familial
#' pf <- 0.2
#' pr <- 0.05
#' priorF <- 0.0865
#' probFamilial(2,8,pf,pr,priorF)$allFamilial
#' # [1] 0.6944444
#' probFamilial(8,8,pf,pr,priorF)$allFamilial
#' # [1] 0.232568
#' # Allowing up to one sporadic case in the pedigree
#' probFamilial(2,8,pf,pr,priorF,r=1)$allFamilial
#' # [1] 0.9722222
#' probFamilial(8,8,pf,pr,priorF,r=1)$allFamilial
#' # [1] 0.6046769
#'
#' @references Dudbridge F, Brown SJ, Ward L, Wilson SG, Walsh JP.
#' How many cases of disease in a pedigree imply familial disease?
#' Submitted.
#'
#' @param m Number of affected relatives of the proband.
#' @param k Number of relatives of the proband with known affection status.
#' @param pf Probability of relative being affected with familial disease.
#' @param pr Probability of relative being affected with sporadic disease.
#' @param priorF Prior probability that a pedigree is segregating familial disease.
#' @param r Number of sporadic cases allowed when calculating probability that all cases are familial.
#'
#' @return mprob Probability of m affected relatives in a pedigree segregating familial disease (equation 2 in Dudbridge et al).
#' @return rprob Probability of m affected relatives in a pedigree with only sporadic disease (calculated as a binomial probability).
#' @return probFamilial Posterior probability that the pedigree is segregating familial disease.
#' @return allFamilial Posterior probability that all the cases in a pedigree are familial, except for up to r sporadic cases.
#' @export
probFamilial <- function(m,k,pf,pr,priorF,r=0) {

# Prob of m affected out of k 1st degree relatives of the proband in a F-family
# Equation 2 in Dudbridge et al
  mprob <- sum(dbinom((0:m),k,pf)*dbinom(m-(0:m),k-(0:m),pr))

# Prob of m affected out of k 1st degree relatives of the proband in a R-family
  rprob <- dbinom(m,k,pr)

# Prob that a family is segregating familial disease when m out of k relatives are affected
# Equation 1 in Dudbridge et al
  famProb <- mprob*priorF/(mprob*priorF+dbinom(m,k,pr)*(1-priorF))

# Prob that all m cases are familial cases
  allFamilial <- sum(dbinom(((m-r):m),k,pf)*dbinom(m-((m-r):m),k-((m-r):m),pr))/mprob

  return(list(mprob=mprob,rprob=rprob,probFamilial=famProb,allFamilial=allFamilial))
}
