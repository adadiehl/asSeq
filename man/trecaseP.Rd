\name{trecaseP}
\alias{trecaseP}
\title{
  Permutation test for TReCASE model
}
\description{
  eQTL mapping using both Total Read Count (TReC) 
  and Allele Specific Expression (ASE).
}
\usage{
trecaseP(Y, Y1, Y2, X, Z, offset=NULL, 
      min.AS.reads=5, min.AS.sample=5, min.n.het=5, local.only=TRUE, 
      local.distance=2e+05, eChr=NULL, ePos=NULL, mChr=NULL, 
      mPos=NULL, converge=5e-5, convergeGLM=1e-8, scoreTestP=0.05, 
      transTestP=0.05, np.max=5000, np=c(20, 100, 500, 1000, 2500),
      aim.p=c(0.5, 0.2, 0.1, 0.05, 0.02), 
      confidence.p=0.01, trace=0, maxit=100)
}
\arguments{
  \item{Y}{ matrix of gene expression (total read count per gene per sample). 
    Each row corresponds to one sample and each column corresponds to one gene. 
  }
  
  \item{Y1}{ matrix of gene expression from one allele (haplotype). 
    Each row corresponds to one sample and each column corresponds to one gene. 
  }
  
  \item{Y2}{ matrix of gene expression from the other allele (haplotype). 
    Each row corresponds to one sample and each column corresponds to one gene. 
    Must be the same dimension as Y1.
  }

  \item{X}{ matrix of cofounding covariates. 
    Each row corresponds to one sample and each column corresponds to one 
    gene. Usually this include log(total reads per sample) 
    and PCs obtained from log transformed expression matrix.
  }
  
\item{Z}{ matrix of phased genotype data, Each row corresponds to one sample 
    and each column corresponds to one SNP. Must take value of 0 (AA), 1 (AB), 
    3 (BA) or 4(BB), where A and B are two alleles of the SNP. AB and BA are 
    different since AB means the 1st/2nd haplotype harbors the A/B allele, 
    respectively, while BA means the 1st/2nd haplotype harbors the B/A allele, 
    respectively. 
  }
    
  \item{offset}{
    A vector of offset of each sample in the model log(mu) = offset + Xb.
  }
  
  \item{min.AS.reads, min.AS.sample}{ one gene will be considered for 
  allele-specific association analysis only if at least min.AS.sample samples
  have at least min.AS.reads allele-specific reads.
  }
  
  \item{min.n.het}{ an allele-specific eQTL will be evaluated only if the 
  at least min.n.het samples have heterozygous genotype at the query SNP. 
  }
  
  \item{local.only}{ whether examine local eQTL only. }
  
  \item{local.distance}{ if the distance between a gene and a marker
   is smaller than local.distance, we call it a local eQTL. }
   
  \item{eChr}{ chromosomes of gene expression traits, must be
   a vector of integers, which means sex chromosomes (X, Y) need
   to be coded as numbers }
   
  \item{ePos}{ positions of gene expression traits, must be
   a vector of integers}
   
  \item{mChr}{ chromosomes of markers, must be a vector of integers,
   which means sex chromosomes (X, Y) need to be coded as numbers }
   
  \item{mPos}{ positions of markers, must be a vector of integers}
  
  \item{converge}{ convergence tolerance. The algorithm converges if 
  changes of all the parameters are smaller than value. }
   
  \item{convergeGLM}{ convergence tolerance to be passed to function glmNB.}

  \item{scoreTestP}{ p-value cutoff for score test of the overdispersion in 
  Poisson model. For each association regression, we first fit a Poisson 
  model, and move to negative binomial if the score test p-value smaller 
  than scoreTestP}

  \item{transTestP}{ p-value cutoff for the hypothesis testing of trans effect. 
    If p-value is smaller than transTestP, the eqTL is treated as trans-eQTL, 
    otherwise it is treated as cis-eQTL}

  \item{np.max}{ the maximum number of permutations }
   
  \item{np}{ number of permutations, a vector in ascending order, 
    see the explanation for \code{aim.p} }
   
  \item{aim.p}{ the permutation p-value threshold, corresponding
    to the number of permutations in \code{np}. We stop the 
    permutations if we can confident to say the permutation p-value 
    is larger than aim.p[i] after np[i] permutations. }
   
  \item{confidence.p}{ A binomial p-value determining whether
    the permutation p-value estimated from the np[i] permutations
    is larger than aim.p[i]. }

  \item{trace}{ If trace=0, do not output tracing information.  
  if trace = 1, output brief tracing information, and 
  if trace > 1, output detailed tracing information. }

  \item{maxit}{ maximum number of iterations of iterative 
  parameter estimate. It is also passed to glmNB to control 
  maximum iteration of IRLS. }
}
\details{
  Different numbers of permutations are performed for the each gene. 
  A larger number of permutations is performed if needed. 
}
\value{
  A data frame with columns gene ID, marker ID,
  nominal p-value, permutation p-value, number of permutations (nuse) 
  and number of permutations in which a test was possible (ntest), 
  for each of the TReC, ASE, TReCASE and TReC/TReCASE (trans/cis) 
  results. Permutations in which no marker could be tested (e.g. for 
  ASE, when too few permuted heterozygous samples have allele-specific 
  reads) carry no information and are left out: the permutation p-value 
  is the fraction of the ntest permutations whose smallest p-value is at 
  most the observed one, and NA if ntest is 0. 
  
  The permutation p-value cannot be smaller than one over the number 
  of permutations. The data frame therefore also has the columns 
  betaA, betaB and betaP for each result: the smallest p-values from 
  the permutations (excluding permutations in which no marker could be 
  tested) are fit by a Beta(betaA, betaB) distribution by maximum 
  likelihood, and betaP is the probability under this distribution of 
  a smallest p-value at most as large as the observed one (Ongen et al. 
  2016, Bioinformatics 32:1479). betaP is NA if fewer than 10 
  permutations were usable or the fit failed, and at least 
  .Machine$double.xmin (a smaller value, e.g. from an observed p-value 
  that underflows to 0, is reported as that). 
 }
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }
