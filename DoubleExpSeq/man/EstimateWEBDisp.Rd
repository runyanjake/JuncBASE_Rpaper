\name{EstimateWEBDisp}
\alias{EstimateWEBDisp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WEB-Seq: Weighted Likelihood Empirical Bayes Estimates of Dispersion for a Double Binomial Distribution
%%  ~~function to do ... ~~
}
\description{
Calculation of shrunken dispersion estimates via weighted likelihood where the weight parameter is estimated using an empirical bayes strategy.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EstimateWEBDisp(y, m, groups, neff = NULL, S = NULL)
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
  \item{neff}{
numeric vector of length equal to the number of rows of "y" where each value is the effective sample size for the event. Default is NULL in which case the effective sample size is calculated within the function.
%%     ~~Describe \code{neff} here~~
}
  \item{S}{
numeric vector of length equal to the number of rows of "y" where each value is the random variable for each event whose distribution across exons is gamma. Default is NULL in which case the vector is calculated internally.
%%     ~~Describe \code{S} here~~
}
}
\details{
Shrunken dispersion estimates are obtained by maximizing the weighted sum of the likelihood for a given event and the sum of likelihoods for all events, the common likelihood. The weight given to the common likelihood is estimated via empirical bayes.
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
\examples{
data(exon)
dispersions <- EstimateWEBDisp( counts, offsets, groups)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
