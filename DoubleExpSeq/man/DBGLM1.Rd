\name{DBGLM1}
\alias{DBGLM1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Double Binomial Generalized Linear Model with Shrinkage of the Dispersion Parameter
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fits a double binomial GLM with the normalization constant set to 1 and uses shrinkage to obtain estimates of dispersion used for p-value calculation.
}
\usage{
DBGLM1( y, m, groups, shrink.method=c("WEB","DEB"),
        contrast=c(1,2), fdr.level=0.05, use.all.groups=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
numeric matrix of inclusion counts.
%%     ~~Describe \code{y} here~~
}
  \item{m}{
numeric matrix of total counts: inclusion + exclusion.
%%     ~~Describe \code{m} here~~
}
  \item{groups}{
vector or factor giving the experimental group/condition for each sample/library.
%%     ~~Describe \code{groups} here~~
}
  \item{shrink.method}{
for shinkage estimation of the dispersion parameter. "WEB" implements the WEB-Seq method. "DEB" implements the DEB-Seq method. Default is "WEB".
%%     ~~Describe \code{shrink.method} here~~
}
  \item{contrast}{
numeric vector of length 2 specifying which levels of the "groups" factor should be compared.
%%     ~~Describe \code{contrast} here~~
}
  \item{fdr.level}{
a numeric constant. The FDR level to determine the list of significant events. Default is 0.05.
%%     ~~Describe \code{fdr.level} here~~
}
  \item{use.all.groups}{
logical. If TRUE, all data in "y" is used to estimate dispersions. If FALSE, only the 2 groups given in "contrasts" are used to estimate dispersions. Only makes a difference if "y" contains more than 2 groups. Default is TRUE.
%%     ~~Describe \code{use.all.groups} here~~
}
}
\details{
This function tests for group differences for a two group comparison via a double binomial GLM with the normalization constant set to 1, and utilizes shrinkage estimates of the dispersion parameter for p-value calcuation which is done using a likelihood ratio test. The shrinkage estimates of the dispersion are found according to the selection of "shrink.method". "DEB" implements the DEB-Seq method which uses an empirical bayes strategy to obtain shrunken estimates of the dipersion parameter. "WEB" implements the WEB-Seq method which reparameterizes the empirical bayes strategy in terms of the weight parameter in the weighted liklelihood formulation. An emprical bayes estimate of the weight parameter is found and plugged into the weighted likelihood which is then maximized to obtain shrunken estimates of the dipsersion parameter. DEB-Seq has shown to be more powerful than WEB-Seq; however, WEB-Seq is more conservative thus being more robust against departures from assumptions and therefore maintains the required FDR better in moderate to larger sample sizes.

In the case when "groups" specifies more than two groups the default procedure is to use all groups to calculate the shrunken dispersion estimates. The argument "constrast" is used to specify a particular comparison of two of the groups. If "use.all.groups" is FALSE, only the data for the groups specified in "contrasts" are used to estimate the dispersions.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    \item{Sig}{a matrix consisting of the significant events at the specified FDR level. The matrix contains the proportion estimates, unadjusted and adjusted p-values, the effective sample size, mean total count and dispersion estimates. }
    \item{All}{the same matrix  as above but consisting of all events.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Sean Ruddy
%%  ~~who you are~~
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

