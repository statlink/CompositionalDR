\name{index_min}
\alias{index_min}
\title{
Cluster indices (minimial valued) for the \eqn{K}--means algorithm
}

\description{
Cluster indices (minimal valued) for the \eqn{K}--means algorithm.
}

\usage{
index_min(y, mod)
}

\arguments{
\item{y}{
A matrix with numerical data.
}
\item{mod}{
An object with the result of the \code{\link{kmeans}} function.
}
}

\details{
A series of cluster validity indices (minimal valued) are computed.
}

\value{
A vector with 9 cluster validity indices.
}

%\references{
%Garcia-Escudero Luis A., Gordaliza Alfonso, Matran Carlos, Mayo-Iscar Agustin. (2008). A general
%trimming approach to robust cluster analysis. Annals of Statistics 36(3): 1324--1345.
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{index_max}, \link{cikmeans}, \link{alfa.cikmeans}
}
}

\examples{
y <- as.matrix(iris[, 1:4])
mod <- kmeans(y, 3)
mod <- index_min(y, mod)
}

