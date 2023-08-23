###### crossdist ######
test_that("crossdist, 2D", {
  set.seed(23082311) 
  x <- matrix(rnorm(20), 10, 2)
  y <- matrix(rnorm(18,5), 9, 2)
  z <- rbind(x,y)
  dfull <- as.matrix(dist(z))
  d_expected <- dfull[1:10, 11:19]
  attr(d_expected, "dimnames") <- NULL
  expect_equal(crossdist(x,y), d_expected)
})

test_that("crossdist, 3D", {
  set.seed(23082312) 
  x <- matrix(rnorm(15), 5, 3)
  y <- matrix(rnorm(18,0,10), 6, 3)
  z <- rbind(x,y)
  dfull <- as.matrix(dist(z))
  d_expected <- dfull[1:5, 6:11]
  attr(d_expected, "dimnames") <- NULL
  expect_equal(crossdist(x,y), d_expected)
})



###### Graph auction, Cpp versus R ######
set.seed(230424)
g1 <- rspatER(n=6, p=0.2)
g2 <- rspatER(n=5, p=0.3)

test_that("OSPA1, Cpp == R", {
  res_comp_Cpp <- gmspat(g1, g2, method = auction_match, compensate=TRUE, type="OSPA1", 
                         CV = 1, CE = 1, vpen = 1, stop_at=3, maxiter=100, lang="Cpp", verbose=0)
  res_comp_R <- gmspat(g1, g2, method = auction_match, compensate=TRUE, type="OSPA1",
                         CV = 1, CE = 1, vpen = 1, stop_at=3, maxiter=100, lang="R", verbose=0)
  expect_true(res_comp_Cpp$dist == res_comp_R$dist)
  all.equal(res_comp_Cpp, res_comp_R)
})

test_that("OSPA2, Cpp == R", {
  res_comp_Cpp <- gmspat(g1, g2, method = auction_match, compensate=TRUE, type="OSPA2", 
                         CV = 1, CE = 1, vpen = 1, stop_at=3, maxiter=100, lang="Cpp", verbose=0)
  res_comp_R <- gmspat(g1, g2, method = auction_match, compensate=TRUE, type="OSPA2",
                       CV = 1, CE = 1, vpen = 1, stop_at=3, maxiter=100, lang="R", verbose=0)
  expect_true(res_comp_Cpp$dist == res_comp_R$dist)
  all.equal(res_comp_Cpp, res_comp_R)
})

test_that("TT, Cpp == R", {
  res_comp_Cpp <- gmspat(g1, g2, method = auction_match, compensate=TRUE, type="TT", 
                         CV = 1, CE = 1, vpen = 1, stop_at=3, maxiter=500, lang="Cpp", verbose=0)
  res_comp_R <- gmspat(g1, g2, method = auction_match, compensate=TRUE, type="TT",
                       CV = 1, CE = 1, vpen = 1, stop_at=3, maxiter=500, lang="R", verbose=0)
  expect_true(res_comp_Cpp$dist == res_comp_R$dist)
  all.equal(res_comp_Cpp, res_comp_R)
})

