\name{DB.MAPlot}
\alias{DB.MAPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots Log-Fold Change versus Log-Concentration for Inclusion/Exclusion Data
%%  ~~function to do ... ~~
}
\description{
M-A Plot
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DB.MAPlot( y, m, groups, contrast=c(1,2), de.tags=NULL,
 col="lightgrey", deCol="red", deCex=0.2,
 xlab="Average Over Groups of log2 Mean Total Count",
 ylab="logFC of Odds Ratio", pch=19, cex=0.2,
 panel.last=grid(col = "red", lwd = 0.2) , ylim = c(-15, 15), ...)

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
numeric vector of length 2 specifying which levels of the "groups" factor should be compared.
%%     ~~Describe \code{contrast} here~~
}
  \item{de.tags}{
rownames for events identified as being differentially expressed.
%%     ~~Describe \code{de.tags} here~~
}
  \item{col}{
color given to the points.
%%     ~~Describe \code{col} here~~
}
  \item{deCol}{
color for the events given in "de.tags".
%%     ~~Describe \code{deCol} here~~
}
  \item{deCex}{
cex for the events given in "de.tags".
%%     ~~Describe \code{deCex} here~~
}
  \item{xlab}{
x-label of plot
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
y-label of plot
%%     ~~Describe \code{ylab} here~~
}
  \item{pch}{
pch given to the points.
%%     ~~Describe \code{pch} here~~
}
  \item{cex}{
cex given to the points.
%%     ~~Describe \code{cex} here~~
}
  \item{panel.last}{
an expression to be evaluated after plotting; the default grid() draws a background grid to aid interpretation of the plot.
%%     ~~Describe \code{panel.last} here~~
}
  \item{ylim}{
y-limits for the plot
%%     ~~Describe \code{ylim} here~~
}
  \item{\dots}{
further arguments passed to plot().
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The total counts are used to determine A, and the log-fold change of the odds ratio is used to determine M. In the case where a group has proportions all 1 or all 0, resulting in an infinite value for M, these points are plotted in orange and away from the main plot. Significant calls made for such events are still colored in red.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
data(exon)
results.G1G2 <- DBGLM1( counts, offsets, groups)
de.tags.G1G2 <- rownames(results.G1G2$Sig)
DB.MAPlot(counts, offsets, groups, contrast=c(1,2), de.tags=de.tags.G1G2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
