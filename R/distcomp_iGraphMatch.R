#' Compute distance between two graphs 
#' 
#' Run an exact or heuristic algorithm to compute one of several assignment-based graph metrics
#' (a.k.a. GTT and GOSPA metrics) and its underlying vertex assignment. Based on [Schuhmacher and Wirth (2023)](https://arxiv.org/abs/2308.12165).\cr 
#' `gmspat` is an alias for `gdist`.
#' 
#' @param g1 An `igraph` object. 
#' @param g2 An `igraph` object.
#' @param method A function used for matching vertices. This function determines which (exact 
#' or heuristic algorithm) is used for computing the distance.  One of  `naive_match`,
#' `cplex_match`, `faq_match` and `auction_match`. See details.
#' @param type One of `"OSPA1"`, `"OSPA2"` and `"TT"` determining the variant of the distance
#' used. See details.  
#' @param CV An upper bound for the vertex distance in the OSPA distances. Default is `CV = 1`.
#' @param CE An upper bound for the edge distance in the OSPA distances. Default is `CE = 1`.
#' @param vpen A penalty for adding auxiliary vertices in the TT distance. Default is `vpen = 1`.
#' @param start_help Value in `[0,1]`. Determines to what extent the solution of the linear vertex matching
#' problem (ignoring the edges) influences the starting solution. Only relevant for method `faq_match`.
#' @param compensate Logical. Whether the computation of the value of a new assignments
#' for `auction_match` considers the influence of this new assignment to other 
#' bidders and objects. Only relevant for method `auction_match`.
#' @param stop_at Maximal number of full assignments to be reached. One of several stopping criteria
#' for `auction_match`.
#' @param maxiter Maximal number of iterations to be reached. One of several stopping criteria
#' for `auction_match`.
#' @param lang One of `"Cpp"` and `"R"`. Specifies whether the C++ or R implementation is called.
#' The C++-implementation is usually much faster. Only relevant for method `auction_match`. 
#' @param eps Small constant that is added to differences of personal values in the auction algorithm.
#' Only relevant for method `auction_match`.
#' @param verbose Numerical. How much information on intermediate results should be printed. 
#' Larger means more. Defaults to `verbose = 0` (print nothing). 
#'  
#' 
#' @details `gdist` expands the two graphs `g1` and `g2` of size `m` and `n` to graphs of the
#' same size in a way that depends on `type`:\cr
#' If `type` is `"TT"`, the graphs are filled up to a size `m+n` using auxiliary vertices, 
#' that are not connected to any other vertex. Adding an auxiliary vertex incurs cost `vpen`.\cr
#' If `type` is `"OSPA1"`, the graphs are filled up to a size \eqn{\max(m,n)}
#' using auxiliary vertices that are connected by auxiliary edges.\cr
#' If `type` is `"OSPA2"`, the graphs are filled up to a size \eqn{\max(m,n)}
#' using auxiliary vertices that are not connected to any other vertex.
#'  
#' Under the hood, the information about the vertex and edge structure of the extended graphs is treated separately: 
#' Euclidean distances between the vertex attributes of the two `igraph`s are passed as argument `dmat` to the `method`
#' function. The `igraph`s themselves are then stripped of their vertex attributes and passed to [iGraphMatch::gm] as
#' its main arguments together with the function specified by `method`. 
#' 
#' Available functions are [naive_match], [cplex_match], [faq_match] and [auction_match]. The first two solve the
#' optimization problem underlying the chosen graph metric `type` exactly. If the CPLEX solver is available on the
#' user's system and the R package `ROI.plugin.cplex` has been installed, it is advisable to use [cplex_match], as it
#' is usually much faster. It is overall still very slow and, depending on the example, may take a long time to return
#' even for relatively small graphs of only about 15 vertices. [faq_match] and [auction_match] provide faster heuristics.
#' We recommend [auction_match] for medium sized graphs (10-30 vertices) and [faq_match] for larger graphs.
#'  
#' The starting solution for `faq_match` in the set of doubly stochastic matrices is computed as the convex combination of 
#' the optimal solution `P_V` of the linear vertex matching problem (ignoring the edges) and the barycenter matrix,
#' which has each entry set to \eqn{1/n} (assuming we have \eqn{n \times n} matrices), corresponding to `start = "bari"` in [faq_match] and
#' [iGraphMatch::graph_match_indefinite] (which is what `faq_match` uses underneath).
#' More precisely, the starting solution is `start_help * P_V + (1- start_help) * bari`.
#' 
#' The function [iGraphMatch::gm] returns an [iGraphMatch::graphMatch] object that is transformed 
#' into a permutation. Additionally the distance/cost of this permutation is calculated using the
#' function [compdist].
#' 
#' @references Dominic Schuhmacher and Leoni Carla Wirth (2023).
#' Assignment Based Metrics for Attributed Graphs.
#' _ArXiv preprint_, \doi{10.48550/arXiv.2308.12165}
#' 
#' @return A list with components `dist` and `perm`, where `dist` is the graph distance corresponding
#' to the determined (near) optimal permutation `perm`. 
#' 
#' @export
#'  
#' @examples
#' set.seed(230820)
#' g1 <- rspatER(n=5, p=0.4)
#' plot(g1, rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), 
#'      vertex.color="orange2", vertex.size=7, edge.color="orange2", edge.width=2)
#' g2 <- rspatER(n=4, p=0.5)
#' plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue", vertex.size=5,
#'      edge.color="lightblue", edge.width=2)
#' gdist(g1, g2, method = naive_match, type="OSPA2")
#' if (requireNamespace("ROI.plugin.cplex", quietly = TRUE)) {
#'   gdist(g1, g2, method = cplex_match, type="OSPA2")
#' }
#' gdist(g1, g2, method = faq_match, type="OSPA2")
#' gdist(g1, g2, method = auction_match, type="OSPA2")
#'
gdist <- function(g1, g2, method = naive_match, type = c("OSPA1", "OSPA2", "TT"),
                   CV=1, CE=1, vpen=1, start_help = 0, compensate=TRUE,
                   stop_at=1, maxiter=100, lang = c("Cpp", "R"), eps = 0.01, verbose=0) {
  #edge_weighted ONLY for OSPA2 and TT! Not possible for OSPA1 (we are not able to get constant edges).
  #edge_weighted not for faq_match (not possible with calling min function)
  #start_help only used for faq_match
  #compensate, stop_at, maxiter, verbose only used for graph auction 
  type <- match.arg(type)
  lang <- match.arg(lang)
  
  xi <- get_coords(g1)
  eta <- get_coords(g2)

  edge_weighted <- igraph::is.weighted(g1) | igraph::is.weighted(g2)
  sigma <- get_edges(g1, adj.mat=TRUE, edge_weighted = edge_weighted)
  tau <- get_edges(g2, adj.mat=TRUE, edge_weighted = edge_weighted)

  m <- dim(xi)[1]
  n <- dim(eta)[1]

  stopifnot(dim(sigma)[1] == dim(sigma)[2] && dim(sigma)[1] == m)
  stopifnot(dim(tau)[1] == dim(tau)[2] && dim(tau)[1] == n)
  stopifnot(all(diag(sigma) == rep(0,m)) && all(diag(tau) == rep(0,n)))  
  stopifnot("Edge weights are only implemented for \"OSPA2\" and \"TT\" " = !(type == "OSPA1" && edge_weighted == TRUE))
  stopifnot("Edge weights not implemented for \"faq_match\"." = !(identical(method,faq_match) && edge_weighted == TRUE))
  
  if (m == 0 && n == 0) return(list(dist=0, perm=integer(0)))
  
  dmat <- crossdist(xi,eta)

  # Fill up the graphs to the same size corresponding to the chosen graph type. 
  if (type == "OSPA1"){
    # For OSPA1 we fill up the smaller graph to the size of the larger graph 
    # by inserting auxiliary vertices and connecting each auxiliary vertex using auxiliary edges.
    # Thus, every auxiliary vertex induces constant vertex cost CV and constant 
    # edge cost of 1/2*CE + min(m,n)/(max(m,n)-1) * CE/2. 
    # As we have "0" and "1" edges, we denote the auxiliary edge with "0.5". 
    # This yields constant edgecost (i.e. it does not matter whether we match the auxiliary edge with a "0"-edge or a "1" edge).
    # However, this yields a factor "1/2" for the constant edge cost (overall 1/4*CE + min(m,n)/(max(m,n)-1) * CE/4 ) when
    # calculating the edge cost via the edgedist function. 
    # The remaining constant edge cost (i.e. 1/4*CE + min(m,n)/(max(m,n)-1) * CE/4) is added to the cost of an auxiliary vertex.
    # Then the cost of one auxiliary vertex is given by C = CV +  1/4*CE + min(m,n)/(max(m,n)-1) * CE/4.
    # We need this new approach for the CPLEX and FW algorithms as they can not interpret "-9" correctly.
    # For the old approach ("-9" edges) see the comments in compdist (in file distcomp_independent.R). 
    
    dmat <- pmin(dmat, CV)
    C <- ifelse(max(n,m) == 1, CV + 1/2 * CE, CV +  1/4*CE + min(m,n)/(max(m,n)-1) * CE/4 )

    if (m < n) {
      dmat <- rbind(dmat, matrix(C, n-m, n)) 
      sigma <- rbind(cbind(sigma,matrix(0.5, m, n-m)),matrix(0.5, n-m, n))
      diag(sigma) <- rep(0, n) 
    } 
    else if (m > n) {
      dmat <- cbind(dmat, matrix(C, m, m-n)) 
      tau <- cbind(rbind(tau, matrix(0.5, m-n, n)),matrix(0.5, m, m-n)) 
      diag(tau) <- rep(0, m)
    }
  }
  
  else if (type == "OSPA2"){
    # For OSPA2 we fill up the smaller graph to the size of the larger graph 
    # by inserting auxiliary vertices that have no edges (i.e. "0" edges). 
    # The edge "fill-up" is done automatically when calling the gm-function. 
    # Hence we only have to fill up the dmat matrix using the cost for an 
    # auxiliary vertex that are given by C = CV + CE.
    
    dmat <- pmin(dmat, CV)
    C <- CV + CE
    
    if (m < n) {
      dmat <- rbind(dmat, matrix(C, n-m, n)) 
    } 
    else if (m > n) {
      dmat <- cbind(dmat, matrix(C, m, m-n)) 
    }
  }
  
  else if (type == "TT"){
    # For TT we fill up both graphs to a size n+m by inserting auxiliary vertices 
    # that have no edges (i.e. "0" edges). In this way we guarantee that each vertex
    # can be matched with an auxiliary vertex. 
    # We have to fill up the dmat matrix using the cost for an 
    # auxiliary vertex that is given by C = vpen.
    # The adjacency matrices are filled up with "0" edges. 
    
    C <- vpen
    
    dmat <- rbind(cbind(dmat,matrix(C, m, m)),cbind(matrix(C, n , n), matrix(0,n,m)))
    sigma <- rbind(cbind(sigma,matrix(0, m, n)),matrix(0, n , n+m))
    tau <- rbind(cbind(tau,matrix(0, n, m)),matrix(0, m , n+m))
  }
  
  if(edge_weighted == FALSE){
    stopifnot(all(sigma %in% c(0,0.5,1)))
    stopifnot(all(tau %in% c(0,0.5,1)))
  }
  
  g1_filled <- igraph::graph_from_adjacency_matrix(sigma, mode = "undirected", weighted = TRUE)
  g2_filled <- igraph::graph_from_adjacency_matrix(tau, mode = "undirected", weighted = TRUE)
  
  if (is.function(method)) {
    if(identical(method,faq_match) ){
      if(start_help>0){
        nn <- ifelse(type=="TT", m+n, max(m,n))
        g1_empty <- igraph::graph_from_data_frame(data.frame(from = "", to = "")[-1,], directed = FALSE, vertices = igraph::get.vertex.attribute(g1))
        g2_empty <- igraph::graph_from_data_frame(data.frame(from = "", to = "")[-1,], directed = FALSE, vertices = igraph::get.vertex.attribute(g2))
        opt_vertex <- gdist(g1_empty,g2_empty, method = faq_match, type = type,CV=CV, CE=CE, vpen=vpen)$perm
        opt_vert_mat <- diag(1,nn)[opt_vertex,]
        start <- start_help * opt_vert_mat + (1-start_help) * matrix(1/nn, nrow = nn, ncol = nn)
      }
      else{
        start = "bari"
      }
      res <- iGraphMatch::gm(g1_filled, g2_filled, seeds = NULL, similarity = NULL, method = method, start = start, dmat = dmat, type = type,
                             m = m, n = n, CV=CV, CE=CE, vpen=vpen, verbose = verbose)
    }
    else if (identical(method, auction_match)){
      res <- iGraphMatch::gm(g1_filled, g2_filled, seeds = NULL, similarity = NULL, method = method, dmat = dmat, type = type, compensate = compensate,
                             m = m, n = n, edge_weighted = edge_weighted, CV=CV, CE=CE, vpen=vpen, stop_at = stop_at, maxiter = maxiter, lang=lang, eps = eps, verbose = verbose)
    }
    else{
      res <- iGraphMatch::gm(g1_filled, g2_filled, seeds = NULL, similarity = NULL, method = method, dmat = dmat, type = type,
                             m = m, n = n, edge_weighted = edge_weighted, CV=CV, CE=CE, vpen=vpen, verbose = verbose)
    }
  } 

  corr <- res@corr
  if (!isTRUE(all.equal(sort(corr[,1]), corr[,1])) || !(dim(corr)[1] == max(corr[,1]))) {
    warning("Non-standard or unordered vertex names ", corr[,1])
  }
  perm <- corr[,2]
  
  dist <- compdist(g1 = g1, g2 = g2 , permlist = perm, type=type, edge_weighted = edge_weighted, CV = CV, CE = CE, vpen = vpen)$dist
  if (verbose > 0) {
    cat("cost", dist, "perm", perm, "\n")
  }
  return(list(dist=dist, perm=perm))
}



#' Graph matching method by exhaustive search
#' 
#' Find the optimal matching w.r.t. a assignment based metric of two given graphs (of the same size) by an exhaustive search. 
#' @inheritParams iGraphMatch::gm 
#' @param dmat A matrix. An n x n matrix  that contains vertex differences. 
#' @param m The initial size of the graph \code{g1}. Default is \eqn{m = nrow(A[[1]])}. 
#' @param n The initial size of the graph \code{g2}. Default is \eqn{n = nrow(B[[1]])}.
#' @param edge_weighted Whether the graphs have edge_weights. Defaults to `FALSE`.
#' @inheritParams gdist 
#' 
#' @details The function \code{naive_match} finds an optimal matching for two graphs of the same size by performing
#' an exhaustive search. It is a version of the \code{\link{compdist}} function but callable as \code{method} in the functions
#' \code{\link[iGraphMatch]{gm}} and \code{\link{gdist}}.
#' 
#' The graphs passed by the \code{\link[iGraphMatch]{gm}} function are already extended to graphs of the same size, see also \code{\link{gdist}}.
#' To find the optimal matching w.r.t. the given \code{type} the initial graph sizes need to be passed via \code{m} and \code{n}.
#' 
#' @seealso \code{\link{compdist}}
#' @return A  \code{\link[iGraphMatch]{graphMatch}} object. 
#' 
#' @export
#'
naive_match <- function(A, B, seeds = NULL, similarity = NULL, dmat, type = "OSPA1", m = nrow(A[[1]]), n = nrow(B[[1]]), edge_weighted = FALSE, CV=1, CE=1, vpen=1, verbose=verbose){ 
  # possible types: OSPA1, OSPA2, TT
  # m,n give the original size of the graph (gm function passes filled-up graphs)
  # CV is cut-off for distance (d_V) between points (default = 1)
  # CE is cut-off for distance (d_E) between edges (default = 1)
  # vpen is penalty in TT distance (default = 1)

  #if (!(n <= 10 && m <= 10)) stop("You do not want to wait for the result.")
  
  if (m == 0 && n == 0) {
    return(iGraphMatch::graphMatch(corr = data.frame(corr_A = integer(0), corr_B = integer(0)),
                                   nnodes=c(0L,0L),
                                   detail=list(dmat)))
  }

  Amat <- A[[1]]  # actually the method functions are supposed to work with lists of graphs
  Bmat <- B[[1]]  # as long as we only pass lists of length 1 all is fine
  
  if(type == "OSPA1" || type == "OSPA2"){
    nn <- max(n,m)
    stopifnot(all(dim(dmat) == c(nn,nn)))
    
    ###### Calculation OSPA1 and OSPA2 ######
    allpi <- combinat::permn(nn)
    nfac <- factorial(nn)
    cost <- rep(0, nfac)
    
    if (nn == 1){ 
      for (L in 1:nfac) {
        curpi <- allpi[[L]]
        cost[L] <- sum(dmat[cbind(1:nn,curpi)])  
      }
      l <- which.min(cost)
      
      if (verbose > 0) {
        cat("Minimal cost:", cost[l], "\n")
      }
      # to return graph match object
      corr <- data.frame(corr_A = 1:nn, corr_B = allpi[[l]])
      return(iGraphMatch::graphMatch(corr=corr, nnodes=c(nn,nn), detail=list(dmat)))
    }
    
    for (L in 1:nfac) {
      curpi <- allpi[[L]]
      Bmatnew <- Bmat[curpi, curpi] 
      # same as P <- Matrix::Diagonal(n)[curpi, ]; P %*% B %*% t(P) and not clear if this would be faster
      cost[L] <- (1/nn) * sum(dmat[cbind(1:nn,curpi)]) + 0.5 * (1/nn) * (1/(nn-1)) * sum(edgedist(as.vector(Amat),as.vector(Bmatnew), maxedist = CE, edge_weighted = edge_weighted)) 
      # *this* 0.5 is because every edgedist appears twice 
    }
  }
  
  else if (type == "TT"){
    nn <- n + m
    stopifnot(all(dim(dmat) == c(nn,nn)))
    
    ###### Calculation TT ######
    allpi <- combinat::permn(nn)
    nfac <- factorial(nn)
    cost <- rep(0,nfac)
    for (L in 1:nfac) {
      curpi <- allpi[[L]]
      Bmatnew <- Bmat[curpi, curpi] 
      cost[L] <- sum(dmat[cbind(1:(nn),curpi)]) +  (1/2) * sum(edgedist(as.vector(Amat),as.vector(Bmatnew), bound = FALSE, edge_weighted = edge_weighted)) 
    }
  }
  
  l <- which.min(cost)
  if (verbose > 0) {
    cat("Minimal cost:", cost[l], "\n")
  }
  
  # we want to return a graph Match object
  corr <- data.frame(corr_A = 1:nn, corr_B = allpi[[l]])
  iGraphMatch::graphMatch(corr=corr, nnodes=c(nn,nn), detail=list(dmat))
}  


# to not create confusion, we pass the distance matrix in an extra argument dmat
# (rather than using the standard argument similarity)
# unfortunately gm does not pass the vertex attributes of the igraphs to the method-function 
# the advantage is (I think) that we can pass g1, g2 also in other formats

#' Graph matching method using CPLEX
#' 
#' Find the optimal matching w.r.t. a assignment based metric of two given graphs (of the same size) using the \code{CPLEX} solver. 
#' @inheritParams iGraphMatch::gm 
#' @param dmat A matrix. An n x n matrix  that contains vertex differences. 
#' @param m The initial size of the graph \code{g1}. Default is \eqn{m = nrow(A[[1]])}. 
#' @param n The initial size of the graph \code{g2}. Default is \eqn{n = nrow(B[[1]])}.
#' @param edge_weighted Whether the graphs have edge_weights. Defaults to `FALSE`.
#' @inheritParams gdist 
#' 
#' @details The function \code{CPLEX_match} transforms the graph matching problem into a binary quadratic programming problem using the 
#' optimization problem constructor \code{\link[ROI]{OP}}. The optimization problem is then solved using the \code{\link[ROI]{ROI_solve}} function
#' with the \code{CPLEX} solver. 
#' 
#' The graphs passed by the \code{\link[iGraphMatch]{gm}} function are already extended to graphs of the same size, see also \code{\link{gdist}}.
#' To find the optimal matching w.r.t. the given \code{type} the initial graph sizes need to be passed via \code{m} and \code{n}. 
#' 
#' @return A  \code{\link[iGraphMatch]{graphMatch}} object. 
#' @references IBM ILOG CPLEX Optimization Studio documentation
#' 
#' @export
#'
cplex_match <- function(A, B, seeds = NULL, similarity = NULL, dmat,  type = "OSPA1",  m = nrow(A[[1]]) , n = nrow(B[[1]]), edge_weighted = FALSE, CV = 1, CE = 1, vpen = 1, verbose=verbose){ 
  # possible types: OSPA1, OSPA2, TT
  # m,n give the *original* size of the graph (gm function passes filled-up graphs, so nrow(A[[1]]), nrow(B[[1]])
  #                                            are usually not correct, but need to be passed correctly)
  # CV is cut-off for distance (d_V) between points
  # CE is cut-off for distance (d_E) between edges
  # vpen is penalty in TT distance

  if (!requireNamespace("ROI.plugin.cplex", quietly = TRUE)) {
    stop("package 'ROI.plugin.cplex' is not available. This method requires an installation of 
          IBM ILOG CPLEX Optimizer, as well as the R packages 'Rcplex' and 'ROI.plugin.cplex'.")
  }
  
  #if (!((n <= 12) && (m<=12))) stop("You do not want to wait for the result.")
  
  if (m == 0 && n == 0) {
    return(iGraphMatch::graphMatch(corr = data.frame(corr_A = integer(0), corr_B = integer(0)),
                                   nnodes=c(0L,0L),
                                   detail=list(dmat)))
  }
  
  Amat <- A[[1]]
  Bmat <- B[[1]]

  ###### Solving the quadratic problem ######
  sigma <- as.matrix(Amat)
  tau <- as.matrix(Bmat)
  if(type == "OSPA1" || type == "OSPA2"){
    nn <- max(n,m)
    temp <- outer(sigma, tau, edgedist, maxedist = CE, edge_weighted = edge_weighted)     
  }
  else if (type == "TT"){
    nn <- n+m
    temp <- outer(sigma, tau, edgedist, bound = FALSE, edge_weighted = edge_weighted)
  }
  #temp <- outer(sigma, tau, \(x,y) {as.numeric(xor(x,y))})
  dearray <- aperm(temp, c(1,3,2,4))   # (i,j,k,l) -> (i,k,j,l)
  #dearray
  # this has the advantage that we really stack the columns as allowed when
  # linearizing the tilde{pi} matrix; it has the disadvantage that it means
  # we keep the each target point (second pattern) fix in turn and vary
  # the source point, which seems less natural...
  # alternative would be c(3,1,4,2). Then advantage and disadvantage exactly the other way round
  dedge=matrix(dearray, nn^2, nn^2)
  dvert <- as.numeric(dmat)  

  # since R is col-major, entries of pi fix each destination and run through source
  # so marginal1 has the structure that usually marginal2 would have and vice versa
  marginal1 <- rep(c(1,rep(0,nn)), each=nn, times=nn+1)
  marginal1 <- matrix(utils::head(marginal1, nn^2*nn), nn, nn*nn, byrow=TRUE)
  marginal2 <- rep(diag(1,nn), nn)
  marginal2 <- matrix(marginal2, nn, nn*nn)
  A <- rbind(marginal1, marginal2)  # matrix for constraints
  b <- rep(1, 2*nn)

  if (type %in% c("OSPA1","OSPA2")){
    if (nn == 1){
      bqp <- ROI::OP(objective = ROI::Q_objective(Q =as.matrix(0,1,1), L = dvert),   
                     constraints = ROI::L_constraint(A, dir = rep("==", 2*nn), rhs = b),  # constraints are still missing
                     types = rep("B", nn^2)) 
      # the annoying CPLEX environment opened / Closed CPLEX environment
      # could be suppressed by saying
      # invisible(capture.output(resc <- ROI::ROI_solve(bqp, solver = "cplex"), type = "message"))
      # instead (not sure if we want that)
      resc <- ROI::ROI_solve(bqp, solver = "cplex")
      if (verbose > 0) {
        cat("Minimal cost:", resc$objval/(nn), "\n")
      }
    }
    else{
      bqp <- ROI::OP(objective = ROI::Q_objective(Q = dedge, L = (nn-1)*dvert),   
                     constraints = ROI::L_constraint(A, dir = rep("==", 2*nn), rhs = b),  # constraints are still missing
                     types = rep("B", nn^2)) 
      resc <- ROI::ROI_solve(bqp, solver = "cplex")
      if (verbose > 0) {
        cat("Minimal cost:", resc$objval/(nn*(nn-1)), "\n")
      }
    }
  }
  
  else if (type == "TT"){
    bqp <- ROI::OP(objective = ROI::Q_objective(Q = dedge, L =dvert),   
                   constraints = ROI::L_constraint(A, dir = rep("==", 2*nn), rhs = b),  # constraints are still missing
                   types = rep("B", nn^2))
    resc <- ROI::ROI_solve(bqp, solver = "cplex", control=list(epgap=1e-5)) 
      #  control=list(epagap=1e-15)   sets the absolute optimality gap tolerance (default 1e-6) 
      #                               --> nothing happens in the close cases if we decrease it !???
      #  control=list(epgap=1e-6)     sets the relative optimality gap tolerance (default 1e-4) 
      #                               --> solution correct in the in the close cases if we set to 1e-5
      #                               (we might want to turn this into a parameter passed to cplex_match at some point)
    #cat("Minimal cost:", resc$objval, "\n")
  }

  
  match <- which(matrix(as.logical(resc$solution), nn, nn), arr.ind=TRUE)

  if(nn!=1){
    match <- match[order(match[,1]),]
    if (!all(match[,1] == 1:nn)) warning("result does not seem to be a permutation")
  }
  
  corr <- data.frame(corr_A = 1:nn, corr_B = match[,2])
  iGraphMatch::graphMatch(corr=corr, nnodes=c(nn,nn), detail=list(dmat))
}  


# approximate OSPA- and TT-like metric between two graphs using Frank-Wolfe algorithm
# OSPA-like: smaller graph is filled up to larger, no tunnelling
# TT-like: both graphs are filled up to size n+m to allow for tunneling
# allows for special cases, such as one point pattern empty or max cardinality = 1,
# but these are not (well) tested
#' Graph matching method using the Frank-Wolfe algorithm
#' 
#' Find the optimal matching w.r.t. a assignment based metric of two given graphs (of the same size)
#' using the Frank-Wolfe algorithm. `FW_match` is an alias for `faq_match` that is deprecated
#' and will be removed in future versions.
#' @inheritParams iGraphMatch::gm 
#' @inheritParams iGraphMatch::graph_match_indefinite
#' @param dmat A matrix. An n x n matrix  that contains vertex differences. 
#' @param m The initial size of the graph \code{g1}. Default is \eqn{m = nrow(A[[1]])}. 
#' @param n The initial size of the graph \code{g2}. Default is \eqn{n = nrow(B[[1]])}.
#' @inheritParams gdist 
#' 
#' @details The function \code{faq_match} finds an optimal matching for two graphs of the same size by calling
#' the \code{\link[iGraphMatch]{graph_match_indefinite}} method which uses the Frank-Wolfe algorithm to solve the graph matching problem. 
#' The vertex distances \code{dmat} are transformed w.r.t. the given \code{type} and passed to the \code{\link[iGraphMatch]{graph_match_indefinite}} method
#' by the \code{similarity} argument. 
#' 
#' The graphs passed by the \code{\link[iGraphMatch]{gm}} function are already extended to graphs of the same size, see also \code{\link{gdist}}.
#' To find the optimal matching w.r.t. the given \code{type} the initial graph sizes need to be passed via \code{m} and \code{n}.
#' 
#' @seealso \code{\link{compdist}}
#' @return A  \code{\link[iGraphMatch]{graphMatch}} object. 
#' 
#' @export
#'
faq_match <- function(A, B, seeds = NULL, similarity = NULL, start = "bari", max_iter = 20, 
                                  lap_method = NULL, dmat, type = "OSPA1", m = nrow(A[[1]]) , n = nrow(B[[1]]) , CV = 1, CE = 1, vpen = 1, verbose = verbose){
  # possible types: OSPA1, OSPA2, TT
  # m,n give the original size of the graph (gm function passes filled-up graphs)
  # CV is cut-off for distance (d_V) between points
  # CE is cut-off for distance (d_E) between edges
  # this function applies gm with indefinite method
  # does not work with edge_weighted (not possible to take min)

  Amat <- A[[1]]  # actually the method functions are supposed to work with lists of graphs
  Bmat <- B[[1]]  # as long as we only pass lists of length 1 all is fine

  if(type == "OSPA1" || type == "OSPA2"){
    
    nn <- max(n,m)  
    stopifnot(all(dim(dmat) == c(nn,nn)))
    
    ###### Calculation OSPA1 and OSPA2 ######
    if(CE == 0){
      return(iGraphMatch::gm(matrix(0,nn,nn), matrix(0,nn,nn), seeds = NULL, similarity = -(nn-1)*dmat, start = start, max_iter = max_iter, lap_method = lap_method, method = "indefinite"))
    }
    else{
      #return(graph_match_indefinite(A, B, seeds = NULL, similarity = -1/CE*(nn-1)*dmat, start = start, max_iter = max_iter, lap_method = lap_method))
      return(iGraphMatch::gm(Amat, Bmat, seeds = NULL, similarity = -1/CE*(nn-1)*dmat, start = start, max_iter = max_iter, lap_method = lap_method, method = "indefinite"))
    }
  }
  
  else if(type == "TT"){
    
    nn <- n+m
    stopifnot(all(dim(dmat) == c(nn,nn)))
    
    ###### Calculation TT ######
    
    return(iGraphMatch::gm(Amat, Bmat, seeds = NULL, similarity = -dmat, start = start, max_iter = max_iter, lap_method = lap_method, method = "indefinite"))
    
  }
}


#' @rdname faq_match
#' @export
FW_match <- faq_match


#' @rdname gdist
#' @export
gmspat <- gdist

