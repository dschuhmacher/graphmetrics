# This file contains code by Raoul Müller (with some substantial changes by the package authors)
# Wir sollten Raoul fragen bevor wir Package MIT diesem Code veröffentlichen (Raoul als ctb)

#'
#' Anderson's PERMANOVA for objects in a general metric space
#'
#' Perform a permutation test to decide whether \eqn{k \geq 2} groups of objects (\eqn{n} in total)
#' all come from the same distribution. The test is mainly sensitive to differing
#' locations of elements across the groups.
#' 
#' @param distmat the full distance matrix between all data objects, arranged according to
#'                group membership (first all objects of group 1, then all objects of group 2, etc.)
#' @param groupsizes the vector of group sizes (length \eqn{k}, sums to \eqn{n}).
#' @param nperm the number of permutations to use for the test. The p-value is honest but random. 
#'              A higher `nperm` reduces its variance.
#' @param alpha the (significance) level of the test.
#'
#' @details See Müller, Schuhmacher and Mateu (2023+). ANOVA for Data in Metric Spaces, with Applications to Spatial Point Patterns.
#'          _ArXiv preprint_ \doi{10.48550/arXiv.2201.08664}.
#'
#' @seealso [`msm_levene`] for a test that is sensitive to differing scatter across groups.
#'
#' @return A list with components
#'  * `Fstat` the value of the test statistic
#'  * `pval` the Monte Carlo p-value based on random permutations 
#'  * `reject` logical. The rejection decision.
#' @export
#'
#' @examples
#' k <- 2
#' n <- 10
#' x <- rnorm(n)
#' y <- rnorm(n, 1)
#' distmat  <- as.matrix(dist(c(x,y)))
#' groupsizes <- rep(n, k)
#' anderson_anova(distmat, groupsizes, nperm=9999) 
#' # because we have univariate data from shifted normal distributions
#' # this is up to the random permutations the same as a t-test:
#' res <- t.test(x, y, var.equal=TRUE) 
#' res
#' res$statistic^2
anderson_anova <- function(distmat, groupsizes, nperm=999, alpha=0.05){

  k <- length(groupsizes)
  stopifnot(k >= 2)
  n <- sum(groupsizes)
  stopifnot(all(dim(distmat) == c(n,n)))
  
  groupvec <- matrix(0,n,k)
  for (i in 1:k) {
    u <- i-1
    index1 <- sum(groupsizes[0:u]) + 1
    index2 <- sum(groupsizes[0:i])
    groupvec[index1:index2,i] <- 1
  }
  
  A <- -0.5*distmat^2
  G <- A - matrix(rowSums(A)/n,n,n) - matrix(colSums(A)/n,n,n,byrow = TRUE) + mean(A)
  X <- groupvec
  invmat <- diag(1/groupsizes)
  H <- X%*%invmat%*%t(X)
  
  SSt <- sum(diag(G))
  SSa <- sum(diag(H%*%G))
  SSr <- SSt - SSa
  # message("TSS ", SSt, "; MSS ", SSa, "; RSS ", SSr)
  
  Fstat <- SSa/SSr * (n-k)/(k-1)
  
  Fvec <- rep(0, nperm)
  for (j in seq_len(nperm)) {
    perm <- sample(n,n)
    newX <- X[perm,]
    newH <- newX%*%invmat%*%t(newX)
    newSSa <- sum(diag(newH%*%G))
    newSSr <- SSt - newSSa
    newF <- newSSa/newSSr * (n-k)/(k-1)
    Fvec[j] <- newF
  }
  
  pval <- (sum(Fvec >= Fstat)+1)/(nperm+1)
  reject <- (pval <= alpha)
  
  return(list(Fstat=Fstat, pval=pval, reject=reject))
}


# 1-factor Levene, distance based for k \geq 2 groups, based on the distance matrix of the elements
# based on MSM2022
# Not sure what is going on here results are strange
msm_levene_R <- function(distmat, groupsizes, nperm=999, alpha=0.05){

  k <- length(groupsizes)
  sqsizes <- groupsizes * (groupsizes-1) / 2
  n <- sum(groupsizes)
  
  groupvec <- matrix(FALSE,n,k)
  for (i in 1:k) {
    u <- i-1
    index1 <- sum(groupsizes[0:u]) + 1
    index2 <- sum(groupsizes[0:i])
    groupvec[index1:index2,i] <- TRUE
  }
  
  Vs <- rep(0,k)  # s as in plural
  for (i in 1:k) {
    index <- groupvec[,i]
    Vs[i] <- 1/2*sum(distmat[index,index])
  }
  Vp <- sum(Vs)/sum(sqsizes) #Vs not scaled yet
  Vs <- Vs / sqsizes
  
  
  enum <- 0
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      enum <- enum + groupsizes[i]*groupsizes[j] * (Vs[i]-Vs[j])^2
    }
  }
  enum <- enum/n
  
  ms <- rep(0,k)
  for (i in 1:k) {
    index <- groupvec[,i]
    m <- (distmat[index,index]-Vs[i])^2
    ms[i] <- sum(m[upper.tri(m)])
  }
  denom <- sum(ms)
  Lstat <- enum/denom #statistic 3 = DISTANCE BASED LEVENE
  
  if(nperm > 0){
    Tmat <- matrix(0,nperm,1)
    for (l in 1:nperm) {
      perm <- sample(n,n)
      pgroupvec <- groupvec[perm,]
      
      pVs <- rep(0,k)
      for (i in 1:k) {
        index <- pgroupvec[,i]
        pVs[i] <- 1/2*sum(distmat[index,index])
      }
      pVp <- sum(pVs)/sum(sqsizes) #Vs not scaled yet
      pVs <- pVs / sqsizes
      
      penum <- 0
      for (i in 1:(k-1)) {
        for (j in (i+1):k) {
          penum <- penum + groupsizes[i]*groupsizes[j] * (pVs[i]-pVs[j])^2
        }
      }
      penum <- penum/n
      pms <- rep(0,k)
      for (i in 1:k) {
        index <- pgroupvec[,i]
        m <- (distmat[index,index]-pVs[i])^2
        pms[i] <- sum(m[upper.tri(m)])
      }
      denom <- sum(pms)
      
      Tmat[l,1] <- enum/denom
    }
  
    pval <- (sum(Tmat[,1] >= Lstat) + 1)/(nperm + 1)
  } else {
    factor <- (sum(sqsizes)-k)/(k-1)
    pval <- pchisq(factor*Lstat, k-1, lower.tail = FALSE)
  }
  
  reject <- (pval <= alpha)
  return(list(Lstat=Lstat, pval=pval, reject=reject))
}


#' The Levene's test from Müller, Schuhmacher and Mateu (2022) for objects in a general metric space
#' 
#' Perform a permutation test to decide whether two groups of objects (\eqn{n} in total)
#' both come from the same distribution. The test is sensitive to differing scatter of elements
#' across the groups, but not at all to differing locations.
#' 
#' @param distmat the full distance matrix between all data objects, arranged according to
#'                group membership (first all objects of group 1, then all objects of group 2, etc.).
#' @param n1,n2 group sizes (\eqn{n1+n2=n}).
#' @param type either "balanced" or "unbalanced". Determines which statistic to use for the test.
#'             If `n1` = `n2` the statistics are equivalent. If n1`, `n2` differ a lot 
#'             the unbalanced type is recommended.
#' @param nperm the number of permutations to use for the test. The p-value is honest but random. 
#'              A higher `nperm` reduces its variance.
#' @param alpha the (significance) level of the test.
#'
#' @details See Müller, Schuhmacher and Mateu (2023+). ANOVA for Data in Metric Spaces, with Applications to Spatial Point Patterns.
#'          _ArXiv preprint_ \doi{10.48550/arXiv.2201.08664}.
#'
#' @seealso [`anderson_anova`] for a test that is sensitive to differing locations across groups.
#'
#' @return A list with components
#'  * `Fstat` the value of the test statistic
#'  * `pval` the Monte Carlo p-value based on random permutations 
#'  * `reject` logical. The rejection decision.
#' @export
#'
#' @examples
#' k <- 2
#' n <- 10
#' x <- rnorm(n)
#' y <- rnorm(n, 0, 3)
#' distmat  <- as.matrix(dist(c(x,y)))
#' groupsizes <- rep(n, k)
#' msm_levene(distmat/2, n, n, "balanced", nperm=9999, alpha=0.05) 
#' msm_levene(distmat/2, n, n, "unbalanced", nperm=9999, alpha=0.05)
#' # results are the same because groups have the same size
msm_levene <- function(distmat, n1, n2, type=c("balanced", "unbalanced"), nperm=999, alpha=0.05) {
  type <- match.arg(type)
  res <- msm_levene_Cpp(distmat, n1, n2, (type=="balanced"), nperm, alpha)

  if (nperm > 0) {
    Lstat <- res$Lstat
    pval <- res$pval  
  } else {
    Lstat <- res$Lstat
    pval <- pchisq(Lstat, 1, lower.tail = FALSE)
  }
  
  reject <- (pval <= alpha)
  return(list(Lstat=Lstat, pval=pval, reject=reject))
}
