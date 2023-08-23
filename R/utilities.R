#' Create a Spatial Graph (igraph with coordinates)
#'
#' @param x vector of x-coordinates of the vertices.
#' @param y vector of y-coordinates of the vertices.
#' @param from vector of indices from which to draw edges.
#' @param to vector of indices to which to draw edges.
#' @param sigma_adj the adjaceny matrix 
#'
#' @return an object of class \code{"igraph"}.
#' @export
#'
#' @examples
#' g1 <- spatgraph(c(-1,1,0), c(0,0,sqrt(3)), from=c(1,1), to=c(2,3))
#' g2 <- spatgraph(c(-0.5,0.5,0), c(0.5/sqrt(3),0.5/sqrt(3),2/sqrt(3)), from=c(1,2), to=c(2,3))
#' plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
#' plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue")
spatgraph <- function(x, y, from, to, sigma_adj=NULL) {
  if (!is.null(sigma_adj)) {
    sigma_adj[lower.tri(sigma_adj)] <- 0
    ind <- which(sigma_adj > 0, arr.ind=TRUE)
    from <- ind[,1]
    to <- ind[,2]
  }
  n <- length(x)
  stopifnot(length(y) == n)
  if (n > 0) {
    name <- 1:n
  } else {
    name <- integer(0)
  }
  xi_df <- data.frame(name=name, x=x, y=y)
  sigma_df <- data.frame(from=from, to=to)
  g <- igraph::graph_from_data_frame(d=sigma_df, directed=FALSE, vertices=xi_df)
  #  class(g) <- c("spatgraph", class(g))    # let's not do this for now, not sure how igraphs respond to that
}

# Convenience functions since DS never remembers the igraph commands involved :-)
get_coords <- function(g) {
  if (igraph::vcount(g) == 0) {
    return(matrix(0,0,2))
  }
  
  vattr <- igraph::vertex_attr(g)
  k <- which(names(vattr) == "name")
  if (length(k)) vattr <- vattr[-k]
  if(length(vattr) == 0){ # nothing left
    return(matrix(0,igraph::vcount(g), 2))
  }
  mat <- do.call(cbind, vattr)  # works for 1D as well
  
  return(mat)
}



# Extract Edges of a Spatial Graph
#
# @param g an \code{igraph} constructed with \code{\link{spatgraph}}.
# @param adj.mat logical. Whether to return the adjacency matrix rather than a from-to matrix
# @param edge_weighted logical. Whether to return edge weights in the adjacency matrix.
#
# @return A two-column matrix specifying in each row from which vertex to which other vertex
#   an edge exists. If \code{adj.mat = TRUE} the adjacency matrix of the graph is returned.
#   If \code{edge_weighted = TRUE} the entries of the adjacency matrix are given by the edge weights. 
#
#
get_edges <- function(g, adj.mat=FALSE, edge_weighted = FALSE) {
  if (adj.mat) {
    if(edge_weighted){
      return( as.matrix(igraph::as_adjacency_matrix(g, attr = "weight")) )
    } else {
      return( as.matrix(igraph::as_adjacency_matrix(g)) )
    }
  } else {
    return( igraph::as_edgelist(g, names=FALSE) )
  }
}

# could also add a function for geom graph,
# could also add a function for more general underlying point process
#' Sample a spatial Erdős--Rényi graph with underlying Poisson or Binomial point process
#'
#' @param lambda a non-negative number. The intensity of the underlying Poisson process. 
#' @param n a non-negative integer. The number of points in the underlying Binomial process.
#'   If \code{n} is specified \code{lambda} is ignored.
#' @param p a number \eqn{\in [0,1]}. The edge probability. 
#'
#' @return the sampled spatial ER graph, as an object of class \code{"igraph"}.
#' @export
#'
#' @examples
#' g <- rspatER(n=8,p=.3)
#' igraph::plot.igraph(g, rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), vertex.size=8)
#' box()
#
rspatER <- function(lambda, n=NULL, p=0.5) {
  if (is.null(n)) {
    n <- rpois(1,lambda)
  }
  x <- runif(n)
  y <- runif(n)
  if (n >= 2) {
    combos <- t(combn(n,2))
    wh <- as.logical(rbinom(choose(n,2),1,p))
    edgemat <- combos[wh, , drop=FALSE]
  } else {
    edgemat <- matrix(0,0,2)
  }
  spatgraph(x, y, edgemat[,1], edgemat[,2])
}


#' Randomly perturb a spatial graph by scattering points and flipping edges
#' 
#' Points are scattered according to a centered bivariate normal distribution with
#' covariance matrix that is a multiple of the identity matrix.
#'
#' @param g an object of class `igraph` with edge attributes `x`and `y`.
#' @param scatter a non-negative number. The standard deviation of the perturbation.
#' @param flip the probability of flipping an edge.
#'
#' @return the perturbed graph, as an object of class \code{"igraph"}.
#' @export
#'
#' @examples
#' g1 <- rspatER(n=8,p=.3)
#' plot(g1, rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1),
#'      vertex.color="orange2", vertex.size=7, edge.color="orange2", edge.width=2)
#' box()
#' g2 <- rperturb(g1, scatter=0.02, flip=0.1)
#' plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue", vertex.size=5,
#'      edge.color="lightblue", edge.width=2)
#
rperturb <- function(g, scatter=0, flip=0.1) {
  n <- igraph::vcount(g)
  coords <- get_coords(g)
  adjmat <- get_edges(g, adj.mat = TRUE)
  if (n >= 1) {
    x <- coords[,1] + rnorm(n, 0, scatter)
    y <- coords[,2] + rnorm(n, 0, scatter)
  } else {
    x <- numeric()
    y <- numeric()
  }
  if (n >= 2) {
    combos <- t(combn(n,2))
    wh1 <- as.logical(adjmat[lower.tri(adjmat)])
    wh2 <- xor(wh1, as.logical(rbinom(choose(n,2),1,flip)))
    edgemat2 <- combos[wh2, , drop=FALSE]  # combos is by row; wh1, wh2 are lower.tri by col = upper.tri by row
  } else {
    edgemat2 <- matrix(0,0,2) 
  }
  spatgraph(x, y, edgemat2[,1], edgemat2[,2])
}


# Compute Matrix of Distances Between Two Sequences of Points
#
# @param x,y two-column matrices specifying the 2-d coordinates of the points.
#
# @return a matrix of size \code{dim(x)[1]} times \code{dim(y)[1]}
# @export
#
# @examples
# x <- matrix(runif(8), 4, 2)
# y <- matrix(runif(10), 5, 2)
# crossdistR(x, y)
#
# I keep the legacy R function for now (the Rcpp function is faster by a factor of
# about 10; not that it really matters as we only call this once for each distance comparison)
crossdistR <- function(x,y) {
  if (length(x) == 2 && length(y) == 2) {
    return(sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2))
  } else {  # probably faster with dist and subsetting
    z1 <- outer(x[,1], y[,1], "-")
    z2 <- outer(x[,2], y[,2], "-")
    return(sqrt(z1^2 + z2^2))
  }
}

# Compute Matrix of Distances Between Two Sequences of Points in 3d
#
# @param x,y three-column matrices specifying the 3-d coordinates of the points.
#
# @return a matrix of size \code{dim(x)[1]} times \code{dim(y)[1]}
# @export
#
crossdistR_3d <- function(x,y) {
  z1 <- outer(x[,1], y[,1], "-")
  z2 <- outer(x[,2], y[,2], "-")
  z3 <- outer(x[,3], y[,3], "-")
  return(sqrt(z1^2 + z2^2 + z3^2))
}


# edgedist1 computes distance between individual edges, edgedist between vectors or
# matrices of edges (result is vector!)
edgedist1 <- function(e,f, maxedist=1, bound = TRUE, edge_weighted = FALSE) {
  if(!bound){
    # GTT case where we might have unbounded edgedistance
    # works the same for weighted edges
    return(abs(e-f))
  }
  else if (bound){
    if(edge_weighted){
      # for edge-weighted distance (or in general for OSPA2, TT)
      return(min(abs(e-f), maxedist))
    }
    else{
      # for 0,0.5,1 edges, i.e. for OSPA1 (and for unweighted OSPA2, TT)
      return(maxedist * abs(e-f))
    }
  }
}  
#
edgedistR <- Vectorize(edgedist1, vectorize.args = c("e", "f"))
# matrix(edgedist(sigma,tau), nn, nn)
# aber normalerweise nicht nötig, dass wir wieder eine Matrix drausmachen und die Diagonalen summieren wir mit
