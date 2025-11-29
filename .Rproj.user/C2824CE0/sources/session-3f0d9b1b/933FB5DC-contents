\name{saa.kl}
\alias{saa.kl}
\title{
Simplicial archetypal analysis minimizing the Kullback-Leibler divergence
}

\description{
Simplicial archetypal analysis minimizing the Kullback-Leibler divergence.
}

\usage{
saa.kl(x, k, lr_w = 0.1, lr_h = 0.1, maxiter = 1000, tol = 1e-7, clip_exp = 50)
}

\arguments{
\item{x}{
An \eqn{n \ times D} matrix with compositional data. Zero values are allowed.
}
\item{k}{
The number of lower dimensions. It must be less than the dimensionality of the data, at most \eqn{D-1}.
}
\item{lr_w}{
The learning rate to multiply the gradient for the \eqn{W} matrix.
}
\item{lr_h}{
The learning rate to multiply the gradient for the \eqn{H} matrix.
}
\item{maxiter}{
The maximum number of iterations allowed.
}
\item{tol}{
The tolerance value to terminate the exponentiated gradient descent algorithm.
}
\item{clip_exp}{
A value to clip the gradient in case it becomes too large.
}
}

\details{
Simplicial archetypal analysis is performed. The typical archetypal analysis with simplicial constraints on the
matrices \eqn{W} and \eqn{H}. The objective function to be minimized is the Kullback-Leibler divergence, \eqn{KL(x||WH)}.
}

\value{
\item{W}{
The \eqn{W} matrix, an \eqn{n \times k} matrix with the mapped compositional data.
}
\item{H}{
The \eqn{H} matrix, an \eqn{k \times D} matrix.
}
\item{Z}{
The reconstructed data, \eqn{Z = WH}.
}
\item{KL}{
The Kullback-Leibler divergence \eqn{ KL(x||WH) }.
}
\item{iters}{
The number of iterations performed.
}
}

\references{
Cutler A. and Breiman L. (1994). Archetypal analysis. Technometrics, 36(4): 338--347.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{saa}, \link{dca}, link{ternary}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
mod <- saa.kl(x, 3)
group <- as.numeric(iris[, 5])
ternary(mod$W, colour = group)
}

