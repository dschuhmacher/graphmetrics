###### anderson_anova ######
set.seed(230429)
k <- 2
n <- 10
x <- rnorm(n)
y <- rnorm(n, 1)
distmat  <- as.matrix(dist(c(x,y)))
groupsizes <- c(n, n)

test_that("anderson_anova, balanced data", {
  Fstat <- unname(t.test(x,y,var.equal=TRUE)$statistic^2)
  set.seed(111)
  res <- anderson_anova(distmat, c(n,n), nperm=9999)
  expect_equal(res$Fstat, Fstat)
  expect_equal(res$pval, 0.1897)
  expect_false(res$reject)
})

test_that("msm_levene, balanced data", {
  set.seed(112)
  resb <- msm_levene(distmat, n, n, "balanced", nperm=9999)
  set.seed(112)
  resub <- msm_levene(distmat, n, n, "unbalanced", nperm=9999)
  expect_equal(resb$Lstat, 0.0128512335)  # checked by hand
  expect_equal(resub$Lstat, 0.0128512335)
  expect_equal(resb$pval, 0.9336)
  expect_equal(resub$pval, 0.9336)
  expect_false(resb$reject)
  expect_false(resub$reject)
})


###### msm_levene ######
set.seed(230430)
k <- 2
n1 <- 9
n2 <- 10
x <- rnorm(n1)
y <- rnorm(n2, 0, 3)
distmat  <- as.matrix(dist(c(x,y)))
groupsizes <- c(n1, n2)

test_that("msm_levene, unbalanced data", {
  set.seed(222)
  resb <- msm_levene(distmat, n1, n2, "balanced", nperm=9999)
  set.seed(222)
  resub <- msm_levene(distmat, n1, n2, "unbalanced", nperm=9999)
  expect_equal(resb$Lstat, 3.883511176)  # checked by hand
  expect_equal(resub$Lstat, 3.268201613) # not checked but plausible
  expect_equal(resb$pval, 0.0049)  
  expect_equal(resub$pval, 0.0049)  # luck I suppose
  expect_true(resb$reject)
  expect_true(resub$reject)
})
