\name{DoubleExpSeq-package}
\alias{DoubleExpSeq-package}
\alias{DoubleExpSeq}
\docType{package}
\title{
DoubleExpSeq is a package with application to RNA-Seq experiments that tests for differential exon usage.
}
\description{
The functions take inclusion and total counts. Inclusion counts are counts that express the exon. Exclusion counts are counts that explicitly skip the exon. The Total count is the sum of inclusion and exclusion.
The package provides 2 methods for the analysis of differential exon usage in RNA-Seq data: 1) DEB-Seq and 2) WEB-Seq. Each of these assume a double-binomial distribution with the normalization constant equal to 1. A simple empirical bayes strategy is used to shrink the dispersion parameter toward a common consensus from all exons.
DEB-Seq implements a 2-parameter empirical bayes strategy to estimate shrunken dispersion estimates. WEB-Seq implements a weighted likelihood approach and estimates the weight parameter via the empirical bayes strategy after reparameterizing the posterior predictive distribution in terms of only the weight parameter.
Once the dispersion estimates are found, the count data is fit with a double binomial GLM (ignoring the normalization constant) and the dispersion estimates used within the test statistic calculation, but have no bearing on the estimates of the proportion.
The main function returns p-values, adjusted p-values using the Benjamini-Hochberg procedure for multiple testing, and proportion estimates along with other relevant information for each exon.

}
\details{
\tabular{ll}{
Package: \tab DoubleExpSeq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-13\cr
License: \tab GPL-3 \cr
}
}
\author{
Sean Ruddy
Maintainer: Sean Ruddy <s.ruddy@yahoo.com>
}
\keyword{ package }
\examples{
## Toy exon data: "counts"=inclusion counts, "offsets"=total counts, "groups"=experiment design
  data(exon)
  
  ## Default will compare "G1" & "G2" using the WEB-Seq method
  ## and uses ALL groups to estimate dispersion
    results.G1G2.WEB <- DBGLM1( counts, offsets, groups)
    
    # Compare G1 & G3
    results.G1G3.WEB <- DBGLM1( counts, offsets, groups, contrast=c(1,3))
	
    # Compare G1 & G3. Does not use G2 for dispersion estimation.
    results.G1G3.noG2.WEB <- DBGLM1(counts, offsets, groups, contrast=c(1,3), use.all.groups=FALSE)
	    
    # Global minimum check of the weight parameter estimate in the WEB-Seq method.
    optimPlot( counts, offsets, groups, contrast=c(1,3), use.all.groups=FALSE)
		
## The DEB-Seq method. Less conservative, more powerful. Very similar ranks to WEB-Seq.
  results.G1G2.DEB <- DBGLM1( counts, offsets, groups, shrink.method="DEB")
		  
## M-A Plot
  WEB.sig <- rownames(results.G1G2.WEB$Sig)
  DB.MAPlot( counts, offsets, groups, de.tags=WEB.sig, main="WEB-Seq")
}
