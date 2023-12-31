% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distcomp_iGraphMatch.R
\name{faq_match}
\alias{faq_match}
\alias{FW_match}
\title{Graph matching method using the Frank-Wolfe algorithm}
\usage{
faq_match(
  A,
  B,
  seeds = NULL,
  similarity = NULL,
  start = "bari",
  max_iter = 20,
  lap_method = NULL,
  dmat,
  type = "OSPA1",
  m = nrow(A[[1]]),
  n = nrow(B[[1]]),
  CV = 1,
  CE = 1,
  vpen = 1,
  verbose = verbose
)

FW_match(
  A,
  B,
  seeds = NULL,
  similarity = NULL,
  start = "bari",
  max_iter = 20,
  lap_method = NULL,
  dmat,
  type = "OSPA1",
  m = nrow(A[[1]]),
  n = nrow(B[[1]]),
  CV = 1,
  CE = 1,
  vpen = 1,
  verbose = verbose
)
}
\arguments{
\item{A}{A matrix, igraph object, or list of either.}

\item{B}{A matrix, igraph object, or list of either.}

\item{seeds}{A vector of integers or logicals, a matrix or a data frame. If
the seed pairs have the same indices in both graphs then seeds can be a
vector. If not, seeds must be a matrix or a data frame, with the first
column being the indices of \eqn{G_1} and the second column being the
corresponding indices of \eqn{G_2}.}

\item{similarity}{A matrix. An \code{n-by-n} matrix containing vertex
similarities. Mandatory for the "IsoRank" method.}

\item{start}{A matrix or a character. Any \code{nns-by-nns} matrix or
character value like "bari", "rds" or "convex" to initialize the starting matrix.}

\item{max_iter}{A number. Maximum number of replacing matches.}

\item{lap_method}{Choice for lap method. One of "lapjv", "lapmod", or "clue".}

\item{dmat}{A matrix. An n x n matrix  that contains vertex differences.}

\item{type}{One of \code{"OSPA1"}, \code{"OSPA2"} and \code{"TT"} determining the variant of the distance
used. See details.}

\item{m}{The initial size of the graph \code{g1}. Default is \eqn{m = nrow(A[[1]])}.}

\item{n}{The initial size of the graph \code{g2}. Default is \eqn{n = nrow(B[[1]])}.}

\item{CV}{An upper bound for the vertex distance in the OSPA distances. Default is \code{CV = 1}.}

\item{CE}{An upper bound for the edge distance in the OSPA distances. Default is \code{CE = 1}.}

\item{vpen}{A penalty for adding auxiliary vertices in the TT distance. Default is \code{vpen = 1}.}

\item{verbose}{Numerical. How much information on intermediate results should be printed.
Larger means more. Defaults to \code{verbose = 0} (print nothing).}
}
\value{
A  \code{\link[iGraphMatch]{graphMatch}} object.
}
\description{
Find the optimal matching w.r.t. a assignment based metric of two given graphs (of the same size)
using the Frank-Wolfe algorithm. \code{FW_match} is an alias for \code{faq_match} that is deprecated
and will be removed in future versions.
}
\details{
The function \code{faq_match} finds an optimal matching for two graphs of the same size by calling
the \code{\link[iGraphMatch]{graph_match_indefinite}} method which uses the Frank-Wolfe algorithm to solve the graph matching problem.
The vertex distances \code{dmat} are transformed w.r.t. the given \code{type} and passed to the \code{\link[iGraphMatch]{graph_match_indefinite}} method
by the \code{similarity} argument.

The graphs passed by the \code{\link[iGraphMatch]{gm}} function are already extended to graphs of the same size, see also \code{\link{gdist}}.
To find the optimal matching w.r.t. the given \code{type} the initial graph sizes need to be passed via \code{m} and \code{n}.
}
\seealso{
\code{\link{compdist}}
}
