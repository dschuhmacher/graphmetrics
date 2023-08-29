# random "medium-sized" examples for heuristic algorithms, 
# verified based on closeness (non-neg. difference) to cplex solution

set.seed(230822)

g1 <- rspatER(n=12, p=0.4)
g2 <- rperturb(g1, scatter=0.1, flip=0.3)
h1 <- rspatER(n=11, p=0.25)
h2 <- rspatER(n=8, p=0.5)
# plot(g1, rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), 
#      vertex.color="orange2", vertex.size=7, edge.color="orange2", edge.width=2)
# plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue", vertex.size=5,
#      edge.color="lightblue", edge.width=2)
# res_exact <- list(perturbed=list(), independent=list())
# res_exact$perturbed$gtt <- gmspat(g1, g2, method = cplex_match, type="TT") # takes more than a minute
# res_exact$perturbed$gospa1 <- gmspat(g1, g2, method = cplex_match, type="OSPA1", CV=0.3)  
# res_exact$perturbed$gospa2 <- gmspat(g1, g2, method = cplex_match, type="OSPA2", CV=0.3)
# plot(h1, rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), 
#      vertex.color="orange2", vertex.size=7, edge.color="orange2", edge.width=2)
# plot(h2, add=TRUE, rescale=FALSE, vertex.color="lightblue", vertex.size=5,
#      edge.color="lightblue", edge.width=2)
# res_exact$independent$gtt <- gmspat(h1, h2, method = cplex_match, type="TT") 
# res_exact$independent$gospa1 <- gmspat(h1, h2, method = cplex_match, type="OSPA1", CV=0.3) # takes more than a minute
# res_exact$independent$gospa2 <- gmspat(h1, h2, method = cplex_match, type="OSPA2", CV=0.3)
# saveRDS(res_exact, file="gmspat_res_exact.rds")
# res_exact <- readRDS(file.path("_testdata/gmspat_res_exact.rds"))  # for checking plausibility

test_that("cross-section test FAQ", {
  res_faq <- list(perturbed=list(), independent=list())
  suppressWarnings(res_faq$perturbed$gtt <- gmspat(g1, g2, method = faq_match, type="TT"))  # reaches max number of iterations (which we cannot adapt yet)
  res_faq$perturbed$gospa1 <- gmspat(g1, g2, method = faq_match, type="OSPA1", CV=0.3)  
  res_faq$perturbed$gospa2 <- gmspat(g1, g2, method = faq_match, type="OSPA2", CV=0.3)
  res_faq$independent$gtt <- gmspat(h1, h2, method = faq_match, type="TT") 
  res_faq$independent$gospa1 <- gmspat(h1, h2, method = faq_match, type="OSPA1", CV=0.3) 
  res_faq$independent$gospa2 <- gmspat(h1, h2, method = faq_match, type="OSPA2", CV=0.3)
  # several of the expected results have perms that are non-unique; perm below is only returned in the
  # current sequence of random number generation with the above seed.
  res_expected <- list(perturbed =
      list(gtt = list(dist = 18.342097183295955, perm = c(2, 9, 12, 11, 7, 1, 10, 8, 5, 4, 3, 6, 14, 24, 22, 21, 19, 23, 16, 17, 13, 20, 18, 15)),
           gospa1 = list(dist = 0.32284780578403838, perm = c(1, 2, 6, 4, 5, 11, 7, 9, 8, 10, 3, 12)),
           gospa2 = list(dist = 0.32284780578403838, perm = c(1, 2, 6, 4, 5, 11, 7, 9, 8, 10, 3, 12))),
                       independent =
      list(gtt = list(dist = 15.435547318154443,  perm = c(4, 1, 16, 12, 3, 8, 6, 5, 7, 2, 19, 10, 13, 9, 17, 18, 14, 15, 11)),
           gospa1 = list(dist = 0.55706899032174517, perm = c(4, 1, 11, 6, 10, 8, 9, 5, 7, 2, 3)),
           gospa2 = list(dist = 0.61579682917869771, perm = c(6, 1, 11, 10, 8, 4, 9, 5, 3, 7, 2))))
  expect_equal(res_faq, res_expected)
})
# looks close enough



test_that("cross-section test auction, without TT", {
  res_auction <- list(perturbed=list(), independent=list())
  # res_auction$perturbed$gtt takes too long to get a first match  
  res_auction$perturbed$gospa1 <- gmspat(g1, g2, method = auction_match, type="OSPA1", CV=0.3)  
  res_auction$perturbed$gospa2 <- gmspat(g1, g2, method = auction_match, type="OSPA2", CV=0.3)
  # res_auction$independent$gtt takes too long to get a first match 
  res_auction$independent$gospa1 <- gmspat(h1, h2, method = auction_match, type="OSPA1", CV=0.3) 
  res_auction$independent$gospa2 <- gmspat(h1, h2, method = auction_match, type="OSPA2", CV=0.3)
  # several of the expected results have perms that are non-unique; 
  # for auction the resulting perm should be stable, but (due to small computational differences??)
  # $independent$gospa1 has slightly different result on all linuxes (tested ubuntu, debian, fedora).
  # dist = 0.559214105651555604, perm = c(4, 1, 3, 11, 10, 8, 6, 5, 7, 2, 9)
  # versus
  # dist = 0.54846412904483188, perm = c(6, 1, 10, 9, 4, 8, 11, 5, 3, 7, 2) 
  # on mac and windows
  # for all I was able to find out this should be due to a small (e.g. rounding) difference maybe in
  # one of the computations
  #   tolower(Sys.info()[["sysname"]])
  print(format(res_auction$independent$gospa1$dist, digits=18))
  print(res_auction$independent$gospa1$perm)
  res_expected <- list(perturbed =
      list(gospa1 = list(dist = 0.32804800621709695, perm = c(1, 2, 11, 4, 5, 6, 7, 9, 8, 10, 3, 12)),
           gospa2 = list(dist = 0.32804800621709695, perm = c(1, 2, 11, 4, 5, 6, 7, 9, 8, 10, 3, 12))), 
                       independent = 
      list(gospa1 = list(dist = 0.54846412904483188, perm = c(6, 1, 10, 9, 4, 8, 11, 5, 3, 7, 2)), 
           gospa2 = list(dist = 0.60477385370831649, perm = c(4, 3, 10, 9, 6, 8, 2, 5, 1, 7, 11)))) 
  expect_equal(res_auction$perturbed, res_expected$perturbed)
  expect_equal(res_auction$independent$gospa2, res_expected$independent$gospa2)
  expect_true(abs(res_auction$independent$gospa1$dist - res_expected$independent$gospa1$dist) < 0.02)
})
# looks close enough
