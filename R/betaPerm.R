# ------------------------------------------------------------
# Beta approximation of permutation p-values
#
# Under the null hypothesis, the smallest p-value of a gene across its
# tested markers approximately follows a Beta(a, b) distribution, with
# a close to 1 and b close to the effective number of independent tests
# (Ongen et al. 2016, FastQTL). Fitting this distribution to the smallest
# p-values from the permutations gives a permutation p-value,
# pbeta(observed smallest p-value, a, b), that is not limited by the
# number of permutations.
# ------------------------------------------------------------

# Permutations in which no marker gave a p-value below 1 are recorded
# as -1 and left out, matching the unpermuted data, where such a gene
# and model are treated as not tested; the fit is then the null
# distribution of the smallest p-value given that a test was possible.

# Fit Beta(a, b) by maximum likelihood to the permuted smallest
# p-values x. Returns c(a, b), or c(NA, NA) if fewer than min.n values
# are usable or the fit fails.
fitBetaPerm <- function(x, min.n=10)
{
  xo = x[is.finite(x) & x > 0 & x < 1]
  n  = length(xo)

  if(n < min.n){ return(c(NA, NA)) }

  S1 = sum(log(xo))
  S2 = sum(log1p(-xo))

  # starting values by the method of moments
  m = mean(xo)
  v = var(xo)
  k = if(is.finite(v) && v > 0) m*(1 - m)/v - 1 else NA
  if(is.finite(k) && k > 0){
    start = log(c(m*k, (1 - m)*k))
  }else{
    start = c(0, 0)
  }

  nll <- function(lp){
    a = exp(lp[1])
    b = exp(lp[2])
    -((a - 1)*S1 + (b - 1)*S2 - n*lbeta(a, b))
  }

  fit = tryCatch(optim(start, nll, method="BFGS",
                       control=list(maxit=1000, reltol=1e-12)),
                 error=function(e) NULL)

  if(is.null(fit) || fit$convergence != 0 || !is.finite(fit$value)){
    return(c(NA, NA))
  }

  exp(fit$par)
}

# Beta-approximated permutation p-values for each gene and model.
# minP: list of length nY*4 with the permuted smallest p-values of gene i
# and model k at position i + (k - 1)*nY; pval: nY x 4 matrix of observed
# smallest p-values; tested: nY x 4 logical matrix. Returns a list of
# nY x 4 matrices a, b and p.
betaPermP <- function(minP, pval, tested)
{
  nY = nrow(pval)
  a  = b = p = matrix(NA_real_, nY, 4)

  for(k in 1:4){
    for(i in which(tested[, k])){
      ab = fitBetaPerm(minP[[i + (k - 1)*nY]])
      if(!is.na(ab[1])){
        a[i, k] = ab[1]
        b[i, k] = ab[2]
        p[i, k] = if(pval[i, k] >= 1) 1 else pbeta(pval[i, k], ab[1], ab[2])
      }
    }
  }

  list(a=a, b=b, p=p)
}
