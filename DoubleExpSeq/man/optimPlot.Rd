\name{optimPlot}
\alias{optimPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the WEB-Seq Maximum Likelihood Solution for the Weight Parameter in the Weighted Likelihood
%%  ~~function to do ... ~~
}
\description{
The MLE solution is signified on a plot as the minimum of the negative log likelihood of the generalized beta prime distribution, parameterized in terms of the weight parameter.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optimPlot(y, m, groups, contrast=c(1,2), use.all.groups=TRUE,...)
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
  \item{contrast}{
numeric vector of length 2 specifying which levels of the "groups" factor should be compared. This is only relevant if "use.all.groups" is FALSE.
%%     ~~Describe \code{contrast} here~~
}
  \item{use.all.groups}{
logical. If TRUE, all data in "y" is used to estimate dispersions. If FALSE, only the 2 groups given in "contrasts" are used to estimate dispersions. Only makes a difference if "y" contains more than 2 groups. Default is TRUE.
%%     ~~Describe \code{use.all.groups} here~~
}
  \item{\dots}{
further arguments passes to plot()
}
}
\details{
The MLE estimate of the (transformed) weight parameter in the WEB-Seq method is checked to be a true global minimum of the negative log likelihood of the generalized beta prime distribution. The weight parameter is transformed from an infinite range to the (0,1) range before optimization and this is the range on which the estimate is checked.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot to the current device
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
# If all groups were used to estimate dispersions
  optimPlot(counts, offsets, groups) 
# If only the 2 groups being compared were used to estimate dispersions
  optimPlot(counts, offsets, groups, contrast=c(1,3), use.all.groups=FALSE) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
