#' Compute OSPA- and TT-like Metrics Between Two Graphs
#'
#' The graphs are assumed to be simple (no loops, undirected, no multiple edges) and without
#' weights, but may have different number of vertices.
#'
#' @param xi,eta two-column matrices. The vertex sets of the two graphs. 
#' @param sigma,tau square matrices. The adjacency matrices for the graphs, of size 
#'          n x n if the graph has n vertices.
#' @param g1,g2 objects of class \code{igraph} generated with \code{\link{spatgraph}}.
#'          Alternative format for passing the graphs (if \code{g1} is specified, \code{xi} and \code{sigma}
#'          are ignored; if \code{g2} is specified, \code{eta} and \code{tau} are ignored).
#' @param permlist a list of permutations. The search for the best match to compute the distance
#'          is restricted to this set. \code{permlist=NULL} means no restriction.
#' @param type one of "OSPA1", "OSPA2" and "TT".
#' @param edge_weighted Logical. Whether the graphs have edge_weights. Default is \code{edge_weighted = FALSE}.
#' @param CV,CE cut-off for distances between vertices and edges, respectively. Only relevant for OSPA types. 
#' @param vpen penalty for additional vertices (plus their incident vertices) in the TT metric. 
#'
#' @details This function computes metrics between two graphs that are based on a best match of their vertices. 
#'   The type argument specifies what we mean by a best match. ...
#'   
#'   All of these optimization problems may be seen as binary quadratic programming problems, which are hard to
#'   solve. The present function performs an exhaustive search, which is very slow. It may serve for very small
#'   graphs (\eqn{\leq} 10 points) and to verify the results of heuristic algorithms.
#'   
#'   The parameter \code{permlist} allows to restrict the search to a subset of permutations. For two graphs with
#'   m and n nodes, each permutation has to be specified as a vector of length \code{max(m,n)} of numbers from
#'   1 to \code{max(m,n)} for the OSPA-type metrics and a vector of length \code{m+n} of numbers from
#'   1 to \code{m+n} for the TT-type metric. These numbers give the images of the vertex numbers of the first graph
#'   in terms of vertex numbers of the second graph, where both graphs have been filled up with dummy points
#'   to create "virtual matches" for unmatched points. 
#'   
#'   \code{permlist} may be just a single one of these vectors, which is interpreted as a list of length 1.
#'   This is currently the main use of \code{permlist}: computing the metric if the optimal matching is 
#'   already (heuristically) known, e.g. from \code{\link[iGraphMatch]{gm}}. Future uses might be e.g.
#'   algorithms of the form: first find best match of vertices ignoring edges and then check locally
#'   for improvements, restricting e.g. the number of transpositions applied from the original matching)
#'   
#' @return A list with components \code{dist} and \code{perm} specifying the desired distance and the 
#'   (optimal) permutation on which it is based.
#' @export
#'
#' @examples
#' g1 <- rspatER(n=5, p=0.4)
#' plot(g1, rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), 
#'      vertex.color="orange2", vertex.size=7, edge.color="orange2", edge.width=2)
#' g2 <- rperturb(g1, scatter=0.05, flip=0.2)
#' plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue", vertex.size=5,
#'      edge.color="lightblue", edge.width=2)
#' compdist(g1=g1, g2=g2, type="OSPA1", CV=1, CE=1)
compdist <- function(xi, sigma, eta, tau, g1=NULL, g2=NULL, permlist = NULL,
                     type = c("OSPA1", "OSPA2", "TT"), edge_weighted = FALSE, CV = 1, CE = 1, vpen = 1) {

  type <- match.arg(type)   # picks one of "OSPA1", "OSPA2", "TT" based on the argument type
                            # or throws an error if the user's intention is not clear
  
  if (is(g1, "igraph")) {
    xi <- get_coords(g1)
    sigma <- get_edges(g1, adj.mat=TRUE)
  }
  if (is(g2, "igraph")) {
    eta <- get_coords(g2)
    tau <- get_edges(g2, adj.mat=TRUE)
  }
  
  m <- dim(xi)[1]
  n <- dim(eta)[1]
  stopifnot(dim(sigma)[1] == dim(sigma)[2] && dim(sigma)[1] == m)
  stopifnot(dim(tau)[1] == dim(tau)[2] && dim(tau)[1] == n)
  stopifnot(all(diag(sigma) == rep(0,m)) && all(diag(tau) == rep(0,n)))  # important check since we sum over the diagonals
  stopifnot("Edge weights are only implemented for \"OSPA2\" and \"TT\" " = !(type == "OSPA1" && edge_weighted == TRUE))
  
  if (m == 0 && n == 0) return(list(dist=0, perm=integer(0)))
  
  nn <- ifelse(type=="TT", m+n, max(m,n))
  
  if (is.null(permlist)) {
    permlist <- combinat::permn(nn)
      # if using the compdist function remains interesting: we can cut down on the run time quite a bit
      # if we only generate 1 permutation of each fixed assignment between non-dummy point
      # I guess this is one of the main reasons why cplex now does much better
    nperm <- factorial(nn)
  } else {
    if (!is.list(permlist)) {
      permlist <- list(permlist)  # assuming user has passed a single permutation vector
    }
    nperm <- length(permlist)
    stopifnot(nperm > 0 && nperm <= factorial(nn))  # minimial sanity check
  }
  
  
  if(type == "OSPA1"){
    
    ###### Fill-up OSPA 1 ######

    C <- ifelse(max(n,m) == 1, CV + 1/2 * CE, CV +  1/4*CE + min(m,n)/(max(m,n)-1) * CE/4 )
    # Fill up with vertex cost CV + half of all auxiliary edge cost
    # Edges are filled-up using auxiliary edges "0.5" -> yield constant cost (but only with factor "1/2" , rest in aux vertex cost) 
    dmat <- pmin(crossdist(xi, eta), CV)
    if (m < n) {
      dmat <- rbind(dmat, matrix(C, n-m, n)) 
      sigma <- cbind(sigma, matrix(0.5, m, n-m)) 
      sigma <- rbind(sigma, matrix(0.5, n-m, n))
      diag(sigma) <- rep(0, n) 
    } 
    else if (m > n) {
      dmat <- cbind(dmat, matrix(C, m, m-n)) 
      tau <- rbind(tau, matrix(0.5, m-n, n)) 
      tau <- cbind(tau, matrix(0.5, m, m-n))
      diag(tau) <- rep(0, m)
    }
    
    ###### Calculation OSPA1 ######
    
    cost <- rep(0,nperm)
    for (L in 1:nperm) {
      curpi <- permlist[[L]]
      if (nn == 1) {
        cost[L] <- sum(dmat[cbind(1:nn,curpi)]) 
      } else {
        cost[L] <- (1/nn) * sum(dmat[cbind(1:nn,curpi)]) + (1/(nn*(nn-1))) * (1/2) * sum(edgedist(sigma,tau[curpi,curpi], maxedist = CE, edge_weighted = edge_weighted))
      }
    }
  }
  
  
  else if (type == "OSPA2"){
    
    ###### Fill-up OSPA 2 #######
    
    C <- CV + CE  # penalty CV for additional points plus CE for its facility to generate edges
    
    # Fill-up using auxiliary elements in OSPA 2-like sense (see proof of metric)
    # Additional vertices have cost C (including cost for facility to generate edges)
    # Edges are filled up using 0's 
    
    dmat <- pmin(crossdist(xi, eta), CV)
    if (m < n) {
      dmat <- rbind(dmat, matrix(C, n-m, n)) 
      sigma <- cbind(sigma, matrix(0, m, n-m))
      sigma <- rbind(sigma, matrix(0, n-m, n))
      diag(sigma) <- rep(0, n) 
    } 
    else if (m > n) {
      dmat <- cbind(dmat, matrix(C, m, m-n)) 
      tau <- rbind(tau, matrix(0, m-n, n)) 
      tau <- cbind(tau, matrix(0, m, m-n))
      diag(tau) <- rep(0, m)
    }
    
    ###### Calculation OSPA2 ######
    
    cost <- rep(0,nperm)
    for (L in 1:nperm) {
      curpi <- permlist[[L]]
      if (nn == 1) {
        cost[L] <- sum(dmat[cbind(1:nn,curpi)])
      } else {
        cost[L] <- (1/nn) * sum(dmat[cbind(1:nn,curpi)]) + (1/(nn*(nn-1))) * (1/2) * sum(edgedist(sigma,tau[curpi,curpi], maxedist = CE, edge_weighted = edge_weighted))
      }
    }
  }
  
  
  else if (type == "TT"){
   
    ###### Fill-up TT ######
    
    C <- vpen # penalty for additional points (cost for facility to generate edges is non-constant)
    
    # Fill-up to size n+m using auxiliary elements in TT sense (see proof of metric)
    # Additional vertices have cost C = vpen 
    # Edges are filled-up using 0's 
    
    dmat <- crossdist(xi, eta)           
    dmat <- rbind(cbind(dmat,matrix(C, m, m)),cbind(matrix(C, n, n), matrix(0,n,m)))
    sigma <- rbind(cbind(sigma,matrix(0, m, n)),matrix(0, n, nn))
    tau <- rbind(cbind(tau,matrix(0, n, m)),matrix(0, m, nn))
    
    stopifnot(all(dim(dmat) == c(nn, nn)))
    
    ###### Calculation TT ######
    
    cost <- rep(0,nperm)
    for (L in 1:nperm) {
      curpi <- permlist[[L]]
      cost[L] <- sum(dmat[cbind(1:(nn),curpi)]) +  (1/2) * sum(edgedist(sigma,tau[curpi,curpi], bound = FALSE, edge_weighted = edge_weighted))
    }
  }
  
  
  else if (!(type %in% c("OSPA1","OSPA2","TT"))){
    stop("No valid type. Choose between OSPA1, OSPA2, TT.")
  }
  
  l <- which.min(cost)
  
  #largeones <- which(cost < cost[l] + sqrt(.Machine$double.eps))
  #if (length(largeones) > 1) { 
  #  warning("Maximum might not be unique")
  #  return(list(cost=cost[largeones], largepi=permlist[largeones]))
  #}
  return(list(dist=cost[l], perm=permlist[[l]]))
}

