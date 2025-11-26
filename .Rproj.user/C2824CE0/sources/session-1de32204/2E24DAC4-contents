\name{alfa.kmeans}
\alias{alfa.kmeans}
\title{
The \eqn{K}--means algorithm for compositional data using the \eqn{\alpha}--transformation
}

\description{
The \eqn{K}--means algorithm for compositional data using the \eqn{\alpha}--transformation.
}

\usage{
alfa.kmeans(x, ncl = 10, trim = 0, a = seq(-1, 1, by = 0.1), max.iters = 50,
nstart = 10)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{ncl}{
The maximum number of clusters to try. The minimum number of clusters is 2.
}
\item{trim}{
A number in [0, 1). If trim = 0, then the classical \eqn{K}--means algorithm is performed.
If you chose a number higher than 0 then the trimmed \eqn{K}--means of Garcia-Escudero et al. (2008) is performed.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be
between -1 and 1. If zero values are present it has to be greater than 0.
If a=0, the isometric log-ratio transformation is applied.
}
\item{max.iters}{
The maximum number of iterations allowed during the \eqn{K}--means algortihm.
}
\item{nstart}{
How many random starts to perform?
}
}

\details{
The \eqn{\alpha}--transformation is applied to the compositional data, and then the \eqn{K}--means algorithm is performed.
}

\value{
A list with the results of the kmeans function for each value of \eqn{\alpha}.
}

\references{
Tsagris M. and Kontemeniotis N. (2025). Simplicial clustering using the \eqn{\alpha}--transformation. https://arxiv.org/pdf/2509.05945.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation
for compositional data. In Proceedings of the 4th Compositional Data Analysis
Workshop, Girona, Spain. https://arxiv.org/pdf/1106.1451.pdf

Garcia-Escudero Luis A., Gordaliza Alfonso, Matran Carlos, Mayo-Iscar Agustin. (2008). A general
trimming approach to robust cluster analysis. Annals of Statistics 36(3): 1324--1345.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cikmeans}
}
}

\examples{
y <- as.matrix(iris[, 1:4])
y <- y / rowSums(y)
mod <- alfa.kmeans( y, ncl = 5, a = c(0, 0.5, 1) )
}

