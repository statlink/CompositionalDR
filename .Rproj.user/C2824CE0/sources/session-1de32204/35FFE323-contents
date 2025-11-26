\name{cikmeans}
\alias{cikmeans}
\title{
The \eqn{K}--means algorithm with cluster indices computed
}

\description{
The \eqn{K}--means algorithm with cluster indices computed.
}

\usage{
cikmeans(y, ncl = 10, trim = 0, max.iters = 50, nstart = 10, all = FALSE)
}

\arguments{
\item{y}{
A matrix with numerical data.
}
\item{ncl}{
The maximum number of clusters to try. The minimum number of clusters is 2.
}
\item{trim}{
A number in [0, 1). If trim = 0, then the classical \eqn{K}--means algorithm is performed.
If you chose a number higher than 0 then the trimmed \eqn{K}--means of Garcia-Escudero et al. (2008) is performed.
}
\item{max.iters}{
The maximum number of iterations allowed during the \eqn{K}--means algortihm.
}
\item{nstart}{
How many random starts to perform?
}
\item{all}{
If this is TRUE, then the clustering indices of each observation for each number of clusters will be returned.
}
}

\details{
The \eqn{K}--means algorithm is performed and a series of cluster validity indices are computed.
}

\value{
A list inclusing:
\item{min_crit}{
A matrix with 9 columns and at least one row, where each column contains the value of a cluster validity index, whose minimal vlaue is preferred. Each row corresponds to a specific number of clusters, starting from 2 up to ncl.
}
\item{best_min}{
The number of clusters selected based upon the minimal valued cluster validity indices.
}
\item{max_crit}{
A matrix with 24 columns and at least one row, where each column contains the value of a cluster validity index, whose minimal vlaue is preferred. Each row corresponds to a specific number of clusters, starting from 2 up to ncl.
}
\item{best_max}{
The number of clusters selected based upon the maximal valued cluster validity indices.
}
\item{cluster}{
If the argument "all" is TRUE, then the clustering indices of each observation for each number of clusters will be returned in a matrix, where each column corresponds to the clustering of each number of clusters.
}
}

\references{
Garcia-Escudero Luis A., Gordaliza Alfonso, Matran Carlos, Mayo-Iscar Agustin. (2008). A general
trimming approach to robust cluster analysis. Annals of Statistics 36(3): 1324--1345.
}

\author{
Michail Tsagris and Nikolaos Kontemeniotis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Nikolaos Kontemeniotis \email{kontemeniotisn@gmail.com}.
.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{index_min}, \link{index_max}, \link{alfa.cikmeans}
}
}

\examples{
y <- as.matrix(iris[, 1:4])
y <- y / rowSums(y)
mod <- cikmeans(y, ncl = 5)
}

