\name{EstimateDEBDisp}
\alias{EstimateDEBDisp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DEB-Seq: Empirical Bayes Estimates of Dispersion for a Double Binomial Distribution
%  ~~function to do ... ~~
}
\description{
Calculation of shrunken dispersion estimates via a 2-parameter empirical bayes method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EstimateDEBDisp(y,m,groups=NULL,neff=NULL,S=NULL,optim.method=c("BFGS","Nelder-Mead"))
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
vector or factor giving the experimental group/condition for each sample/library. Default NULL.
%%     ~~Describe \code{groups} here~~
}
  \item{neff}{
numeric vector of length equal to the number of rows of "y" where each value is the effective sample size for the event. Default is NULL in which case the effective sample size is calculated within the function.
%%     ~~Describe \code{neff} here~~
}
  \item{S}{
numeric vector of length equal to the number of rows of "y" where each value is the random variable for each event whose distribution across exons is gamma. Default is NULL in which case the vector is calculated internally.
%%     ~~Describe \code{S} here~~
}
  \item{optim.method}{
character string determining which optimization routine to use for estimating the parameters of the prior distribution. Default is "BFGS".
%%     ~~Describe \code{optim.method} here~~
}
}
\details{
The shrunken dispersion estimates are a function of 2 parameters of the generalized beta prime distribution which are estimated via maximum likelihood resulting in empricaly bayes shinkage of the dispersion parameter.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
vector of length equal to the number of rows of "y" where each value is the estimate of dispersion.
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
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
data(exon)
dispersions <- EstimateDEBDisp( counts, offsets, groups)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
