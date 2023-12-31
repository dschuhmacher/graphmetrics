% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{rperturb}
\alias{rperturb}
\title{Randomly perturb a spatial graph by scattering points and flipping edges}
\usage{
rperturb(g, scatter = 0, flip = 0.1)
}
\arguments{
\item{g}{an object of class \code{igraph} with edge attributes \code{x}and \code{y}.}

\item{scatter}{a non-negative number. The standard deviation of the perturbation.}

\item{flip}{the probability of flipping an edge.}
}
\value{
the perturbed graph, as an object of class \code{"igraph"}.
}
\description{
Points are scattered according to a centered bivariate normal distribution with
covariance matrix that is a multiple of the identity matrix.
}
\examples{
g1 <- rspatER(n=8,p=.3)
plot(g1, rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1),
     vertex.color="orange2", vertex.size=7, edge.color="orange2", edge.width=2)
box()
g2 <- rperturb(g1, scatter=0.02, flip=0.1)
plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue", vertex.size=5,
     edge.color="lightblue", edge.width=2)
}
