#' Graph matching method by auction algorithm
#' 
#' Find the optimal matching w.r.t. a assignment based metric of two given graphs (of the same size) 
#' sing an auction algorithm. `graph_auction` is an alias for `auction_match` that is deprecated
#' and will be removed in future versions.
#' 
#' @inheritParams iGraphMatch::gm 
#' @param dmat A matrix. An n x n matrix that contains vertex differences. 
#' @param compensate Logical. Whether or not the computation of the value of a new assignment in the \code{auction_match} algorithm also 
#'considers influences of this new assignment to other bidders and objects. Default is \code{compensate = TRUE}.
#' @param edge_weighted Whether the graphs have edge_weights. Defaults to `FALSE`.
#' @param m The initial size of the graph \code{g1}. Default is \eqn{m = nrow(A[[1]])}. 
#' @param n The initial size of the graph \code{g2}. Default is \eqn{n = nrow(B[[1]])}.
#' @inheritParams gmspat 
#' @param stop_at Maximal number of full assignments that can be reached before the auction algorithm stops. Default is \eqn{stop_at = 1}.
#' @param maxiter Maximal number of iterations that can be reached before the auction algorithm stops. Default is \eqn{maxiter = 100}.
#' @param lang One of "Cpp" and "R". Specifies whether the C++ or R implementation is called. The C++-implementation is usually much faster.   
#' @param eps Small constant that is added to differences of personal values in the auction algorithm. Only relevant for method \code{"auction_match"}. Default is \code{eps = 0.01}.
#' @param verbose Numerical. How much information on intermediate results should be printed during the algorithm. The larger verbose the more information is printed. Default is \code{verbose = 0}. 
#' 
#' @details The function \code{auction_match} approximates an optimal matching for two graphs of the same size by performing
#' an heuristic auction algorithm. 
#' 
#' The graphs passed by the \code{\link[iGraphMatch]{gm}} function are already extended to graphs of the same size, see also \code{\link{gmspat}}.
#' To find the optimal matching w.r.t. the given \code{type} the initial graph sizes need to be passed via \code{m} and \code{n}.
#' 
#' In each iteration the personal object values of a bidder are calculated. These values depend on the chosen \code{type}
#' and if \code{compensate = TRUE} also the influence of a possible new assignment of the current bidder to the current object to other bidders
#' and objects (for example due to un-assignment) is considered. 
#'
#' An iteration is done if every bidder got the chance to bid. A full assignment is reached if every bidder is assigned to an object.
#' If the number of iterations is \eqn{\geq maxiter} or the number of full assignments is \eqn{\geq stop_at} the algorithm stops 
#' and returns the best full assignment. 
#' 
#' @return A  \code{\link[iGraphMatch]{graphMatch}} object. 
#' 
#' @export
#' 
auction_match <- function(A, B, seeds = NULL, similarity = NULL, dmat, compensate=TRUE, type="OSPA1", edge_weighted = FALSE,
                          CV = 1, CE = 1, vpen = 1, m = nrow(A[[1]]) , n = nrow(B[[1]]), stop_at=1, maxiter=100, lang=c("Cpp","R"), eps = 0.01, verbose=0) { 
  #edge_weighted only for R Version
  #eps <- 0.01 # normally 1/(n+1), here not clear; .Machine$double.eps^(1/3) will make things very slow
  
  lang <- match.arg(lang)
  
  sigma <- as.matrix(A[[1]])
  tau <- as.matrix(B[[1]])
  
  if (m == 0 && n == 0) return(list(dist=0, perm=integer(0)))

  if (type == "OSPA1") {
    nn <- max(m,n)
    maxcost <- CV + CE/2
    # needed in application 
    # eps = eps * (CV + CE/2) *2
    # print(eps)
  }
  else if (type == "OSPA2") {
    nn <- max(m,n)
    maxcost <- CV + CE
    # needed in application
    # eps = eps * (CV + CE) *1/2
    # print(eps)
  }
  else if (type == "TT") {
    nn <- m+n
    maxcost <- vpen 
    #eps = eps * (max(m,n) * (max(m,n) -1)) *2
    #print(eps)
  }
  
  if (lang == "Cpp") {
    best_assignment <- rep(NA, nn)
    best_cost <- Inf
    numtype <- switch(type, TT=0, OSPA1=1, OSPA2=2)
    res <- auction_match_Cpp(sigma, tau, dmat, nn, maxcost, compensate, numtype, CE, stop_at, maxiter, eps, verbose)
    res$best_assignment = res$best_assignment + 1
  } else {
    res <- auction_match_R(sigma, tau, dmat, nn, maxcost, compensate, type, edge_weighted, CE, stop_at, maxiter, eps, verbose)
  }
  
  # from a practical point of view it may be desirable in larger examples to keep track
  # of near-matches up to 2 or more points (instead of only one) and do assign the remaining
  # vertices by the auction algorithm for vertices only (or restart graph auction on the remaining set.
  # if it is reasonably small)
  # note: near-match_num includes full matches
  if (res$nearmatch_num == 0) stop("Maximal number of iterations reached without a near-match.")
  if (res$fullmatch_num==0) warning("No full matched reached, returning best of ", res$nearmatch_num, " near-match(es).")
  if (verbose > 0) {
    print(res$best_assignment)
    message("Internally computed best_cost was ", res$best_cost)
  }
  
  # we want to return a graph Match object
  corr <- data.frame(corr_A = 1:nn, corr_B = res$best_assignment)
  iGraphMatch::graphMatch(corr=corr, nnodes=c(nn,nn), detail=list(dmat))
  #return(list(dist = best_cost, perm = best_assignment))
}


auction_match_R <- function(sigma, tau, dmat, nn, maxcost, compensate, type, edge_weighted, CE, stop_at, maxiter, eps, verbose) {
  pers_to_obj <- rep(NA, nn)
  obj_to_pers <- rep(NA, nn)
  prices <- rep(0, nn)
  pers_values <- matrix(0, nn, nn)
  iter <- 1
  fullmatch_num <- 0
  fullmatch <- FALSE 
  assignments <- list()
  best_assignment <- numeric()
  costs <- numeric()
  best_cost <- Inf
  change <- TRUE
  while (change && iter <= maxiter && fullmatch_num < stop_at) {
    if (verbose > 0) cat("\n----------------------------------\n", iter, "\n", sep="")
    # if (is.na(pers_to_obj[bidder])) {   # bidder is not assigned
    # bidding
    change <- FALSE
    for (bidder in 1:nn) {
      if (compensate) {
        my_values <- pers_values[bidder, ] <- compute_pers_values2(bidder, nn, pers_to_obj, obj_to_pers, maxcost, dmat, sigma, tau, prices, type, edge_weighted, CE)
      } else {
        my_values <- pers_values[bidder, ] <- compute_pers_values(bidder, nn, pers_to_obj, obj_to_pers, maxcost, dmat, sigma, tau, prices, type, edge_weighted, CE)
      }
      my_values_endsort <- sort(my_values, partial=nn-1) 
      # guarantees that for the resulting vector v (my_values_decr) that v[k] is smaller than v[partial] if k < partial
      # and larger if k > partial --> the two largest value will be at the end of the vector 
      # for about 10 or more points typically faster than sorting the whole vector
      bidfor <- which.max(my_values)
      oldbidfor <- pers_to_obj[bidder]
      if (is.na(oldbidfor) || oldbidfor != bidfor) { ## non-unique max *should* be no problem(?)
        change <- TRUE
        # if (!is.na(oldbidfor))  # not needed (at best saves a little time)
        obj_to_pers[oldbidfor] <- NA
        oldbidder <- obj_to_pers[bidfor] 
        # if (!is.na(oldbidder))  # not needed (at best saves a little time)
        pers_to_obj[oldbidder] <- NA
        pers_to_obj[bidder] <- bidfor 
        obj_to_pers[bidfor] <- bidder
        bidincrement <- (my_values_endsort[nn]-my_values_endsort[nn-1]) + eps
        if (verbose > 0) {
          cat("eps", eps, "bidincrement", my_values_endsort[nn]-my_values_endsort[nn-1], "\n")
        }
        #cat("pers_val_dist", bidincrement - eps, "eps", eps, "\n")
        prices[bidfor] <- prices[bidfor] + bidincrement
        
        unassigned <- sum(is.na(pers_to_obj))
        if (unassigned <= 1) {
          curpi <- pers_to_obj
          if (unassigned == 1) {
            curpi[is.na(curpi)] <- setdiff(1:nn, curpi[!is.na(curpi)]) 
            fullmatch <- FALSE
          } else {  # unassigned == 0  (we count only real full matches as full matches)
            if (!fullmatch) {
              fullmatch <- TRUE
              fullmatch_num <- fullmatch_num + 1
            }
          }
          if (type == "OSPA1" || type == "OSPA2") { 
            # @Leoni: it seems aux-edge-info is in sigma and tau and correctly interpreted by edgedist --> formulae for OSPA1 and OSPA2 exactly the same!?
            # yes. We fill-up OSPA1 using 0.5 edges. Thus, the edge matching cost is constant (i.e. the same for matching 0.5 with 0,1 edge). 
            # This yields, that we can calculate the edge distance as in the OSPA2 case (where d_E(0.5,y) = 1/2 * CE for all possible edges y). 
            # Clearly, this reduces the cost for an auxiliary edge matching (by factor 1/2), this is why we add the remaining constant cost via dmat (when filling up with aux vertices). 
            tempcost <- (1/nn) * sum(dmat[cbind(1:nn, curpi)]) + (1/(nn*(nn-1))) * (1/2) * sum(edgedist(sigma, tau[curpi,curpi], maxedist = CE, edge_weighted = edge_weighted))
          } else if (type == "TT") {
            tempcost <- sum(dmat[cbind(1:(nn), curpi)]) +  (1/2) * sum(edgedist(sigma, tau[curpi,curpi], bound = FALSE, edge_weighted = edge_weighted)) 
          }
          # print(curpi)
          # print(dmat[cbind(1:(nn), curpi)])
          # print(matrix(edgedist(sigma, tau[curpi,curpi], maxedist = CE), nn, nn))
          # cat(1/nn, "  ", 1/(nn*(nn-1)))
          # message(tempcost, " = ", sum(dmat[cbind(1:nn, curpi)]), " + ", sum(edgedist(sigma, tau[curpi,curpi], maxedist = CE)))
          
          assignments <- c(assignments, list(curpi))
          costs <- c(costs, tempcost)
          if (tempcost < best_cost) {
            best_cost <- tempcost
            best_assignment <- curpi
          }
        } else {
          fullmatch <- FALSE
        }
        
        if (verbose > 0) {
          print(paste(bidder, "for", bidfor))
          if (verbose > 1) print(pers_values)
          print(pers_to_obj)
          print(prices)
        }
      }
    }
    iter <- iter + 1
  }
  
  if (verbose > 0) {
    message("Done at change = ", change, ", iter = ", iter, ", full matches = ", fullmatch_num)
    print(assignments)
    print(costs)
  }
  nearmatch_num <- length(costs)
  
  return(list(best_assignment=best_assignment, best_cost=best_cost,
              change=change, iter=iter, fullmatch_num=fullmatch_num, nearmatch_num=nearmatch_num))
}



# NO compensation and optimistic towards edges to unassigned nodes
compute_pers_values <- function(i, nn, pers_to_obj, obj_to_pers, maxcost, dmat, sigma, tau, prices, type, edge_weighted, CE) {
  edgediff <- rep(0, nn)
  # neutralize *CE if edge_weighted is TRUE
  edge_weight_CE <- ifelse(edge_weighted, 1/CE, 1)
  # print(pers_to_obj)
  # print(obj_to_pers)
  for (j in 1:nn) {  # investigate values for matching person i to obj j
    pers_to_obj_hypo <- pers_to_obj
    obj_to_pers_hypo <- obj_to_pers
    oldi <- obj_to_pers_hypo[j]
    if (!is.na(oldi)) {  # if obj j is already matched
      pers_to_obj_hypo[oldi] <- NA
    }
    pers_to_obj_hypo[i] <- j
    obj_to_pers_hypo[j] <- i
    
    #print(pers_to_obj_hypo)
    #print(obj_to_pers_hypo)
    
    # for assigned vertices k != i we match i--k and j--pi(k) edge by edge;
    # among the unassigned vertices we match freely in the best possible way
    unassigned <- is.na(pers_to_obj_hypo)  
    assigned <- !unassigned
    aobj <- pers_to_obj_hypo[assigned]
    if(type == "OSPA1"|| type == "OSPA2"){
      edgediff[j] <- sum(unlist(edgedist(sigma[i,assigned], tau[j,aobj], maxedist = CE, edge_weighted = edge_weighted)))#sum(abs(sigma[i,assigned] - tau[j,aobj]))
      # here we always add a zero from the diagonal (i is always assigned to j)
      edgediff[j] <- edgediff[j] + abs( sum(sigma[i,unassigned]) - sum(tau[j, setdiff(1:nn, aobj)]) )*CE*edge_weight_CE
    }
    else if (type == "TT"){
      edgediff[j] <- sum(unlist(edgedist(sigma[i,assigned], tau[j,aobj], bound = FALSE, edge_weighted = edge_weighted)))
      # here we always add a zero from the diagonal (i is always assigned to j)
      edgediff[j] <- edgediff[j] + abs( sum(sigma[i,unassigned]) - sum(tau[j, setdiff(1:nn, aobj)]) )
    }
  }
  if (type == "OSPA1" || type == "OSPA2") { 
    persval <- rep(maxcost, nn) - (1/nn) * dmat[i,] - 0.5*(1/(nn*(nn-1))) * edgediff - prices # - indemnity
  }
  else if (type == "TT") {
    persval <- rep(maxcost, nn) -  dmat[i,] - 0.5* edgediff - prices # - indemnity
  }
  #print(persval)
  #if (i== 7) stop("safety stop");
  return(persval)
}


# compensation (called transferpay below) + optimistic towards edges to unassigned nodes
compute_pers_values2 <- function(i, nn, pers_to_obj, obj_to_pers, maxcost, dmat, sigma, tau, prices, type, edge_weighted, CE) {
  edgediff <- rep(0, nn)
  transferpay <- rep(NA, nn)
  # neutralize *CE if edge_weighted is TRUE
  edge_weight_CE <- ifelse(edge_weighted, 1/CE, 1)
  for (j in 1:nn) {  # investigate values for matching person i to obj j
    pers_to_obj_hypo <- pers_to_obj
    obj_to_pers_hypo <- obj_to_pers
    oldi <- obj_to_pers_hypo[j]
    if (!is.na(oldi)) {  # if obj j is already matched
      pers_to_obj_hypo[oldi] <- NA
    }
    pers_to_obj_hypo[i] <- j
    obj_to_pers_hypo[j] <- i
    
    # for assigned vertices k != i we match i--k and j--pi(k) edge by edge;
    # among the unassigned vertices we match freely in the best possible way
    unassigned <- is.na(pers_to_obj_hypo)  
    assigned <- !unassigned
    aobj <- pers_to_obj_hypo[assigned]
    if(type == "OSPA1"|| type == "OSPA2"){
      edgediff[j] <- sum(unlist(edgedist(sigma[i,assigned], tau[j,aobj], maxedist = CE)))#sum(abs(sigma[i,assigned] - tau[j,aobj]))
      # here we always add a zero from the diagonal (i is always assigned to j)
      edgediff[j] <- edgediff[j] + abs( sum(sigma[i,unassigned]) - sum(tau[j, setdiff(1:nn, aobj)]) )*CE*edge_weight_CE # optimistic
    }
    else if (type == "TT"){
      edgediff[j] <- sum(unlist(edgedist(sigma[i,assigned], tau[j,aobj], bound = FALSE, edge_weighted = edge_weighted)))
      # here we always add a zero from the diagonal (i is always assigned to j)
      edgediff[j] <- edgediff[j] + abs( sum(sigma[i,unassigned]) - sum(tau[j, setdiff(1:nn, aobj)]) ) # optimistic
    }
    #print(edgediff[j])
    theothers <- setdiff((1:nn)[assigned], i)  # Knoten, die vor der (Neu-)Zuordnung i --> j, zugeordnet waren und auch gleich zugeordnet bleiben
    oldj <- pers_to_obj[i]
    
    # externalities for switching object oldj to j (for bidder i)
    if (is.na(oldj)) {  # i did not have an active bid
      externality_oldj <- 0 # Can we do something more appropriate here?
      #externality_unass_j <- sum(tau[pers_to_obj[theothers],setdiff(1:nn, pers_to_obj[theothers])])
     } else {
      if(type == "OSPA1"|| type == "OSPA2"){
        distbefore <- sum(unlist(edgedist(sigma[theothers, i], tau[pers_to_obj[theothers], oldj], maxedist = CE, edge_weighted = edge_weighted))) #sum(abs(sigma[theothers, i] - tau[pers_to_obj[theothers], oldj])) 
        distafter <- sum(unlist(edgedist(sigma[theothers, i], tau[pers_to_obj_hypo[theothers], j], maxedist = CE, edge_weighted = edge_weighted))) #sum(abs(sigma[theothers, i] - tau[pers_to_obj_hypo[theothers], j])) 
      }
      else if (type == "TT"){
        distbefore <- sum(unlist(edgedist(sigma[theothers, i], tau[pers_to_obj[theothers], oldj], bound = FALSE, edge_weighted = edge_weighted))) 
        distafter <- sum(unlist(edgedist(sigma[theothers, i], tau[pers_to_obj_hypo[theothers], j], bound = FALSE, edge_weighted = edge_weighted)))  
      }
        externality_oldj <- distafter - distbefore  
      #externality_unass_j <- sum(tau[c(pers_to_obj[theothers],oldj),setdiff(1:nn, c(pers_to_obj[theothers],oldj))])
     }
    #print(paste(externality_oldj, "oldj"))
    # externalities for switching bidder oldi to i (for object j)
    if (is.na(oldi)) {  # j was not currently bid for 
      externality_oldi <- 0 # Can we do something more appropriate here?
      #externality_unass_i <- sum(sigma[theothers,setdiff(1:nn, theothers)]) 
    } else {
      if(type == "OSPA1"|| type == "OSPA2"){
        distbefore <- sum(unlist(edgedist(sigma[theothers, oldi], tau[pers_to_obj[theothers], j], maxedist = CE, edge_weighted = edge_weighted))) #sum(abs(sigma[theothers, oldi] - tau[pers_to_obj[theothers], j]))
        distafter <- sum(unlist(edgedist(sigma[theothers, i], tau[pers_to_obj_hypo[theothers], j], maxedist = CE, edge_weighted = edge_weighted))) #sum(abs(sigma[theothers, i] - tau[pers_to_obj_hypo[theothers], j]))  
      } 
      else if (type == "TT"){
        distbefore <- sum(unlist(edgedist(sigma[theothers, oldi], tau[pers_to_obj[theothers], j], bound = FALSE, edge_weighted = edge_weighted))) 
        distafter <- sum(unlist(edgedist(sigma[theothers, i], tau[pers_to_obj_hypo[theothers], j], bound = FALSE, edge_weighted = edge_weighted)))   
      }
      externality_oldi <- distafter - distbefore
      #externality_unass_i <- sum(sigma[c(theothers,oldi),setdiff(1:nn, c(theothers,oldi))]) 
    }
    #print(paste(externality_oldi, "oldi", i, j))
    #externality_unass_edgecost_new <- abs(externality_unass_i-externality_unass_j)
    #externality_unass_edgecost_old <- abs(sum(sigma[theothers,setdiff(1:nn, c(theothers,i))]) - sum(tau[pers_to_obj[theothers],setdiff(1:nn, c(pers_to_obj[theothers],j))]))
    #externality_unass <- externality_unass_edgecost_new - externality_unass_edgecost_old
    transferpay[j] <- externality_oldi + externality_oldj# + externality_unass
  }
  if (type == "OSPA1" || type == "OSPA2") { 
    persval <- rep(maxcost, nn) - (1/nn) * dmat[i,] - 0.5*(1/(nn*(nn-1))) * edgediff - prices - 0.5*(1/(nn*(nn-1))) * transferpay
  }
  else if (type == "TT") {
    persval <- rep(maxcost, nn) -  dmat[i,] - 0.5* edgediff - prices - 0.5 * transferpay
  }
  return(persval)
}

# # with compensation (called transferpay below) and maximal cost for edges with unassigned nodes (pessimistic)
# however "pessimistic" seems to suppress competition between the nodes to much (one always gets a bonus for)
# bidding for an unassigned compared to bidding for an assigned one, which creates an additional unassigned one afterwards)
# compute_pers_values2b <- function(i, nn, pers_to_obj, obj_to_pers, maxcost, dmat, sigma, tau, prices) {
#   edgediff <- rep(0, nn)
#   transferpay <- rep(NA, nn)
#   for (j in 1:nn) {  # investigate values for matching person i to obj j
#     pers_to_obj_hypo <- pers_to_obj
#     obj_to_pers_hypo <- obj_to_pers
#     oldi <- obj_to_pers_hypo[j]
#     if (!is.na(oldi)) {  # if obj j is already matched
#       pers_to_obj_hypo[oldi] <- NA
#     }
#     pers_to_obj_hypo[i] <- j
#     obj_to_pers_hypo[j] <- i
#     
#     # for assigned vertices k != i we match i--k and j--pi(k) edge by edge;
#     # among the unassigned vertices we match freely in the best possible way
#     unassigned <- is.na(pers_to_obj_hypo)  
#     assigned <- !unassigned
#     aobj <- pers_to_obj_hypo[assigned]
#     edgediff[j] <- sum(abs(sigma[i,assigned] - tau[j,aobj]))
#     # we have no benefit from all the potential edges of the unassigned --> add to cost
#     edgediff[j] <- edgediff[j] + sum(unassigned)*(nn-1)
#     
#     theothers <- setdiff((1:nn)[assigned], i)  # Knoten, die vor der (Neu-)Zuordnung i --> j, zugeordnet waren und auch gleich zugeordnet bleiben
#     oldj <- pers_to_obj[i]
#     
#     # in any case bidder i gets a bonus for every matching edge it creates for others bidders
#       transferpay[j] <- -sum(1 - abs(sigma[theothers, i] - tau[pers_to_obj_hypo[theothers],j]))
#     if (!is.na(oldj)) {
#       # if bidder i was assigned to an object before it has to pay for each matching edge other bidders lose due to the switch to a new object
#       transferpay[j] <- transferpay[j] + sum(1-abs(sigma[theothers, i] - tau[pers_to_obj[theothers],oldj]))  # pers_to_obj[theothers] = pers_to_obj_hypo[theothers]
#     }
#       
#     if (!is.na(oldi)) {
#       # if somebody else (oldi) bid for j before, bidder i also has to pay for each matching edge other bidders lose due to unmatching oldi 
#       # wenn ein anderes Biety sich 
#       transferpay[j] <- transferpay[j] + sum(1 - abs(sigma[theothers, oldi] - tau[pers_to_obj[theothers],j]))  # pers_to_obj[theothers] = pers_to_obj_hypo[theothers]
#     }
#   }
#   
#   if (i == 4) browser()
#   persval <- rep(maxcost, nn) - (1/nn) * dmat[i,] - 0.5*(1/(nn*(nn-1))) * edgediff - prices - 0.5*(1/(nn*(nn-1))) * transferpay
#   return(persval)
#   # in later versions we should have some retrospect adaptation of prices if the assignment changes (I guess) 
# }


#' @rdname auction_match
#' @export
graph_auction <- auction_match

