# small examples, verified by closeness (essentially non-neg. difference) to exact solution

skip_if_not_installed("ROI.plugin.cplex")
# because of cplex_match (test_that blocks below are hard to split up)
# in the long run do cplex tests separately or remove them here

testgraphs <- readRDS(file.path("_testdata/testgraphs.rds"))

test_that("gmspat for OSPA1, small, different size, all methods", {
  res_naive <- matrix(0,5,5)
  res_cplex <- matrix(0,5,5)
  res_faq <- matrix(0,5,5)
  res_auction <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res_naive[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = naive_match, type="OSPA1")$dist
      res_cplex[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = cplex_match, type="OSPA1")$dist
      suppressWarnings(res_faq[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = faq_match, type="OSPA1")$dist)
      res_auction[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = auction_match, type="OSPA1")$dist
    }
  }
  res_expected <- 
    structure(c(0, 0.94282933768246413, 1.1000716253872604, 1.1209574386007599, 
                1.1668574810569226, 0.94282933768246413, 0, 0.86698324634631729, 
                0.89839250596456122, 1.0975691613348677, 1.1000716253872604, 
                0.86698324634631729, 0, 0.77716500694471358, 0.90219736591598676, 
                1.1209574386007599, 0.89839250596456122, 0.77716500694471358, 
                0, 0.74003878645913068, 1.1668574810569226, 1.0975691613348677, 
                0.90219736591598676, 0.74003878645913068, 0), dim = c(5L, 5L))
  expect_equal(res_naive, res_expected)  
  expect_equal(res_cplex, res_expected)  
  res_expected_faq <-
    structure(c(0, 0.94282933768246413, 1.1000716253872604, 1.1209574386007599, 
                1.1673277336459376, 0.94282933768246413, 0, 0.88346515635821432, 
                0.89839250596456122, 1.0975691613348677, 1.1000716253872604, 
                0.88346515635821432, 0, 0.81630895563335959, 0.90219736591598687, 
                1.1209574386007599, 0.89839250596456122, 0.81630895563335959, 
                0, 0.74454262880514532, 1.1673277336459373, 1.0975691613348677, 
                0.90219736591598687, 0.74454262880514532, 0), dim = c(5L, 5L))
  expect_equal(res_faq, res_expected_faq)
  res_expected_auction <-
    structure(c(0, 0.94282933768246413, 1.1000716253872604, 1.1209574386007599, 
                1.1668574810569228, 0.94282933768246413, 0, 0.86698324634631729, 
                0.89839250596456122, 1.1011563224823757, 1.1000716253872604, 
                0.86698324634631729, 0, 0.81630895563335959, 0.90219736591598687, 
                1.1209574386007599, 0.91700794434276645, 0.81630895563335959, 
                0, 0.74003878645913079, 1.1668574810569228, 1.0975691613348677, 
                0.90219736591598687, 0.74003878645913079, 0), dim = c(5L, 5L))
  expect_equal(res_auction, res_expected_auction)  
  
  # only for naive 
  res <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = naive_match, type="OSPA1", CV=0.1)$dist
    }
  }
  res_expected <- 
    structure(c(0, 0.42406487648828706, 0.51666666666666661, 0.55000000000000004, 
                0.56596973002044515, 0.42406487648828706, 0, 0.34999999999999998, 
                0.42368673788766181, 0.49108667679551821, 0.51666666666666661, 
                0.34999999999999998, 0, 0.30000000000000004, 0.39692846742903121, 
                0.55000000000000004, 0.42368673788766181, 0.30000000000000004, 
                0, 0.32046096716942279, 0.56596973002044515, 0.49108667679551821, 
                0.39692846742903121, 0.32046096716942279, 0), dim = c(5L, 5L))
  expect_equal(res, res_expected)  
})


test_that("gmspat for OSPA2, small, different size, all methods", {
  res_naive <- matrix(0,5,5)
  res_cplex <- matrix(0,5,5)
  res_faq <- matrix(0,5,5)
  res_auction <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res_naive[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = naive_match, type="OSPA2")$dist
      res_cplex[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = cplex_match, type="OSPA2")$dist
      suppressWarnings(res_faq[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = faq_match, type="OSPA2")$dist)
      suppressWarnings(res_auction[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = auction_match, type="OSPA2")$dist)
    }
  }
  res_expected <- 
    structure(c(0, 0.9658405615683272, 1.2667382920539272, 1.47095743860076, 
                1.5335241477235895, 0.9658405615683272, 0, 0.86698324634631729, 
                1.0170079443427664, 1.3043251324630702, 1.2667382920539272, 0.86698324634631729, 
                0, 0.92187641706505963, 1.0021973659159868, 1.47095743860076, 
                1.0170079443427664, 0.92187641706505963, 0, 0.77337211979246412, 
                1.5335241477235895, 1.3043251324630702, 1.0021973659159868, 0.77337211979246412, 
                0), dim = c(5L, 5L))
  expect_equal(res_naive, res_expected)
  expect_equal(res_cplex, res_expected)
  res_expected_faq <-
    structure(c(0, 0.9658405615683272, 1.2667382920539272, 1.47095743860076, 
                1.5335241477235895, 0.9658405615683272, 0, 0.88346515635821432, 
                1.0170079443427664, 1.3062798179835664, 1.2667382920539272, 0.88346515635821432, 
                0, 0.9292422496407774, 1.0021973659159868, 1.47095743860076, 
                1.0170079443427664, 0.9292422496407774, 0, 0.77787596213847865, 
                1.5335241477235895, 1.3062798179835664, 1.0021973659159868, 0.77787596213847865, 
                0), dim = c(5L, 5L))
  expect_equal(res_faq, res_expected_faq)
  res_expected_auction <-
    structure(c(0, 0.9658405615683272, 1.2667382920539272, 1.47095743860076, 
                1.5335241477235895, 0.9658405615683272, 0, 0.86698324634631729, 
                1.0170079443427664, 1.3062798179835664, 1.2667382920539272, 0.86698324634631729, 
                0, 0.92187641706505952, 1.0317804094867886, 1.47095743860076, 
                1.0170079443427664, 0.9698380495935095, 0, 0.77787596213847865, 
                1.5339944003126043, 1.3174695727215544, 1.0021973659159868, 0.78396128483722072, 
                0), dim = c(5L, 5L))
  expect_equal(res_auction, res_expected_auction)
  
  # only for naive 
  res <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = naive_match, type="OSPA2", CV=0.1)$dist
    }
  }
  res_expected <-
    structure(c(0, 0.5907315431549538, 0.68333333333333346, 0.90000000000000013, 
                0.93263639668711174, 0.5907315431549538, 0, 0.35000000000000003, 
                0.62368673788766182, 0.72442001012885149, 0.68333333333333346, 
                0.35000000000000003, 0, 0.44474840931092691, 0.56359513409569784, 
                0.90000000000000013, 0.62368673788766182, 0.44474840931092691, 
                0, 0.35379430050275618, 0.93263639668711174, 0.72442001012885149, 
                0.56359513409569784, 0.35379430050275618, 0), dim = c(5L, 5L))
  expect_equal(res, res_expected)
})


test_that("gmspat for TT, small, different size, only naive", {
  skip_on_cran()  # takes too much time
  res <- matrix(0,3,3)
  for (i in 1:3) {  # takes still too much time otherwise
    for (j in 1:3) {
      res[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = naive_match, type="TT")$dist
    }
  }
  res_expected <- 
    structure(c(0, 2.3975216847049818, 3.7336198348823757, 2.3975216847049818, 
                0, 2.4679329853852692, 3.7336198348823757, 2.4679329853852692, 
                0), dim = c(3L, 3L))
  expect_equal(res, res_expected)
  
  res <- matrix(0,3,3)
  for (i in 1:3) {
    for (j in 1:3) {
      res[i,j] <- gmspat(g1=testgraphs[[i]], g2=testgraphs[[j]], method = naive_match, type="TT", vpen=0.1)$dist
    }
  }
  res_expected <-
    structure(c(0, 1.3721946294648613, 1.5992475902808845, 1.3721946294648613, 
                0, 1.429675403424447, 1.5992475902808845, 1.429675403424447, 
                0), dim = c(3L, 3L))
  expect_equal(res, res_expected)
})

