# older distcomp code
# very small examples, verified by manual computation

# basic setup for several examples: 
# xi: vertices as equilateral triangle, side length 2, center at (0, 1/sqrt(3)), edges <
g1 <- spatgraph(c(-1,1,0), c(0,0,sqrt(3)), from=c(1,1), to=c(2,3)) # equilateral triangle
# eta varies between three and zero points; for eta with three points:
# vertices equilateral triangle, side length 2*(1-lambda), center at (0, 1/sqrt(3)), edges >
# define function to generate the various scalings
scaled_g <- function(lambda) {
  x=(1-lambda)*c(-1,1,0) + lambda*rep(0,3)
  y=(1-lambda)*c(0,0,sqrt(3)) + lambda*rep(1/sqrt(3),3)
  spatgraph(x, y, from=c(1,2), to=c(2,3))
}


###### OSPA 1 ######

test_that("compdist for OSPA1, small and regular, same size", {
  # critical lambda-value is lam0 = 45/(8*sqrt(3)+48) = 0.7274913
  # below lam_0 identity matching is optimal, above lam_0 switching 1 and 2 is optimal for OSPA1
  lam <- 0.7275  # high
  g2 <- scaled_g(lam)
  d11 <- (2/sqrt(3))*lam     # dist from 1 to 1 (between graphs)
  d12 <- sqrt((4/3)*lam^2 - 4*lam + 4)  # dist from 1 to 2 (between graphs)
  d_exp <- d11/3 + 2*d12/3  # cost of switching 1 and 2
  # d11 + 1/3   # cost of identity permutation
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue")
  CV <- 1.34005 # so that d12 is just not cut 
  CE <- 1 # CE is always 1 in this block of tests
  res <- compdist(g1=g1, g2=g2, type="OSPA1", CV=CV, CE=CE)  # so that d12 is just not cut
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(2,1,3))
  
  lam <- 0.72748  # low
  g2 <- scaled_g(lam)
  d11 <- (2/sqrt(3))*lam     # dist from 1 to 1 (between graphs)
  d12 <- sqrt((4/3)*lam^2 - 4*lam + 4)  # dist from 1 to 2 (between graphs)
  d_exp <- d11 + 1/3   # cost of identity permutation
  # d11/3 + 2*d12/3  # cost of switching 1 and 2
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue")
  CV <- 1.34005 # so that d12 is just not cut 
  CE <- 1 # CE is always 1 in this block of tests
  res <- compdist(g1=g1, g2=g2, type="OSPA1", CV=CV, CE=CE)  # so that d12 is just not cut
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(1,2,3))
  
  CV <- 1.340023 # so that d12 is cut, but it just does not matter for the perm (and hence for the dist)
  # the critical value where it will start to matter is CV = (2/sqrt(3))*lam + 1/2 = 1.340021548 
  res <- compdist(g1=g1, g2=g2, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(1,2,3))  
  
  CV <- 1.340020 # d12 is cut and it does matter
  d_exp_cut <- (1/3)*d11 + (2/3)*CV  # 1 and 2 switched in optimal solution
  res <- compdist(g1=g1, g2=g2, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp_cut)
  expect_equal(res$perm, c(2,1,3))  
  
  CV <- 0.84 # just below d11, now all distances between vertices are 0.84
  res <- compdist(g1=g1, g2=g2, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, CV)
  expect_equal(res$perm[1], 2)
    # correct results are c(2,3,1) and c(2,1,3), though with the current code always the first one
  expect_equal(sort(res$perm), c(1,2,3))  # just to be super safe
})
  

test_that("compdist for OSPA1, small, different size", {
  # 2 vertices
  g2 <- spatgraph(c(-0.1,-0.1), c(sqrt(3),0), from=1, to=2)
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="darkseagreen")
  
  CV <- 1
  CE <- 1
  C <- CV + CE/2
  d_exp <- (1/3)*(C + min(0.9, CV) + min(0.1, CV)) + (1/12)*2*(0 + CE)  # Formel am Anfang von 8.2 (Stein-Paper)
  # (d_alt1 <- (1/3)*(C + min(0.9, CV) + min(0.1, CV)) + (1/6)*2*(0 + (1/2)*CE))  # hatten wir mal zwischendurch als Def. 2.4 (Metrik-Paper)
  # (d_alt2 <- (1/3)*(CV + min(0.9, CV) + min(0.1, CV)) + (1/12)*2*(0 + 2*CE))  # was wir m.E. in compdist rechnen
  # gibt das gleiche
  res <- compdist(g1=g1, g2=g2, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(2,3,1))
  # symmetry test
  res <- compdist(g1=g2, g2=g1, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(3,1,2))
  # hier evtl. noch weitere Tests mit anderem CV und evtl CE

  # 1 vertex
  g2 <- spatgraph(-0.1, sqrt(3), from=numeric(), to=numeric())
  CE <- 1
  CV <- 1
  C <- CV + CE/2
  d_exp <- (1/3)*(2*C + min(0.1, CV)) + (1/12)*(0 + 2*CE)  # Formel am Anfang von 8.2 (Stein-Paper)
  # (d_exp <- (1/3)*(2*C + min(0.1, CV)) + (1/6)*1*(0 + (1/2)*2*CE))  # hatten wir mal zwischendurch als Def. 2.4 (Metrik-Paper)
  # (d_alt2 <- (1/3)*(2*CV + min(0.1, CV)) + (1/12)*2*(0 + 3*CE))  # was wir m.E. in compdist rechnen
  # gibt alles das gleiche
  res <- compdist(g1=g1, g2=g2, type="OSPA1", CV=CV, CE=CE)  
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[3], 1)
  expect_equal(sort(res$perm), c(1,2,3))
  # symmetry test
  res <- compdist(g1=g2, g2=g1, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1], 3)
  expect_equal(sort(res$perm), c(1,2,3))
  
  # 0 vertices
  g0 <- spatgraph(numeric(), numeric(), from=numeric(), to=numeric())
  res <- compdist(g1=g1, g2=g0, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, 1.5)
  expect_equal(sort(res$perm), c(1,2,3))
  # symmetry test
  res <- compdist(g1=g0, g2=g1, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, 1.5)
  expect_equal(sort(res$perm), c(1,2,3))
  # 1 vs 0 points
  res <- compdist(g1=g2, g2=g0, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, 1.5)
  expect_equal(res$perm, 1)
  # symmetry test
  res <- compdist(g1=g0, g2=g2, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, 1.5)
  expect_equal(res$perm, 1)
  # 0 vs 0 points
  res <- compdist(g1=g0, g2=g0, type="OSPA1", CV=CV, CE=CE) 
  expect_equal(res$dist, 0)
  expect_equal(res$perm, integer(0))
})



###### OSPA 2 ######

test_that("compdist for OSPA2, small and regular, same size", {
  # Of course for same size, everything should be exactly the same as for OSPA1!
  # (still worth testing because there is separate code)
  # critical lambda-value is lam0 = 45/(8*sqrt(3)+48) = 0.7274913
  # below lam_0 identity matching is optimal, above lam_0 switching 1 and 2 is optimal for OSPA1
  lam <- 0.7275  # high
  g2 <- scaled_g(lam)
  d11 <- (2/sqrt(3))*lam     # dist from 1 to 1 (between graphs)
  d12 <- sqrt((4/3)*lam^2 - 4*lam + 4)  # dist from 1 to 2 (between graphs)
  d_exp <- d11/3 + 2*d12/3  # cost of switching 1 and 2
  # d11 + 1/3   # cost of identity permutation
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue")
  CV <- 1.34005 # so that d12 is just not cut 
  CE <- 1 # CE is always 1 in this block of tests
  res <- compdist(g1=g1, g2=g2, type="OSPA2", CV=CV, CE=CE)  # so that d12 is just not cut
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(2,1,3))
  
  lam <- 0.72748  # low
  g2 <- scaled_g(lam)
  d11 <- (2/sqrt(3))*lam     # dist from 1 to 1 (between graphs)
  d12 <- sqrt((4/3)*lam^2 - 4*lam + 4)  # dist from 1 to 2 (between graphs)
  d_exp <- d11 + 1/3   # cost of identity permutation
  # d11/3 + 2*d12/3  # cost of switching 1 and 2
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue")
  CV <- 1.34005 # so that d12 is just not cut 
  CE <- 1 # CE is always 1 in this block of tests
  res <- compdist(g1=g1, g2=g2, type="OSPA2", CV=CV, CE=CE)  # so that d12 is just not cut
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(1,2,3))
  
  CV <- 1.340023 # so that d12 is cut, but it just does not matter for the perm (and hence for the dist)
  # the critical value where it will start to matter is CV = (2/sqrt(3))*lam + 1/2 = 1.340021548 
  res <- compdist(g1=g1, g2=g2, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(1,2,3))  
  
  CV <- 1.340020 # d12 is cut and it does matter
  d_exp_cut <- (1/3)*d11 + (2/3)*CV  # 1 and 2 switched in optimal solution
  res <- compdist(g1=g1, g2=g2, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp_cut)
  expect_equal(res$perm, c(2,1,3))  
  
  CV <- 0.84 # just below d11, now all distances between vertices are 0.84
  res <- compdist(g1=g1, g2=g2, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, CV)
  expect_equal(res$perm[1], 2)
  # correct results are c(2,3,1) and c(2,1,3), though with the current code always the first one
  expect_equal(sort(res$perm), c(1,2,3))  # just to be super safe
})


test_that("compdist for OSPA2, small, different size", {
  # 2 vertices
  g2 <- spatgraph(c(-0.1,-0.1), c(sqrt(3),0), from=1, to=2)
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="darkseagreen")
  
  CV <- 1
  CE <- 1
  C2 <- CV + CE
  d_exp <- (1/3)*(C2 + min(0.9, CV) + min(0.1, CV)) + (1/12)*(0 + 1 + 1)  # Formel am Anfang von 8.4 (Stein-Paper)
  res <- compdist(g1=g1, g2=g2, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(2,3,1))
  # symmetry test
  res <- compdist(g1=g2, g2=g1, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm, c(3,1,2))
  # hier evtl. noch weitere Tests mit anderem CV und evtl CE
  
  # 1 vertex
  g2 <- spatgraph(-0.1, sqrt(3), from=numeric(), to=numeric())
  CE <- 1
  CV <- 1
  C2 <- CV + CE
  d_exp <- (1/3)*(2*C2 + min(0.1, CV)) + (1/12)*(0 + 2 + 2)  # Formel am Anfang von 8.4 (Stein-Paper)
  res <- compdist(g1=g1, g2=g2, type="OSPA2", CV=CV, CE=CE)  
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[3], 1)
  expect_equal(sort(res$perm), c(1,2,3))
  # symmetry test
  res <- compdist(g1=g2, g2=g1, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1], 3)
  expect_equal(sort(res$perm), c(1,2,3))
  
  # 0 vertices
  g0 <- spatgraph(numeric(), numeric(), from=numeric(), to=numeric())
  d_exp <- C2 + (1/12)*(0 + 2 + 2)
  res <- compdist(g1=g1, g2=g0, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(sort(res$perm), c(1,2,3))
  # symmetry test
  res <- compdist(g1=g0, g2=g1, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, d_exp)
  expect_equal(sort(res$perm), c(1,2,3))
  # 1 vs 0 points
  res <- compdist(g1=g2, g2=g0, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, 2)
  expect_equal(res$perm, 1)
  # symmetry test
  res <- compdist(g1=g0, g2=g2, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, 2)
  expect_equal(res$perm, 1)
  # 0 vs 0 points
  res <- compdist(g1=g0, g2=g0, type="OSPA2", CV=CV, CE=CE) 
  expect_equal(res$dist, 0)
  expect_equal(res$perm, integer(0))
})



###### TT ######

test_that("compdist for TT, small and regular, same size", {
  # For moderate to large vpen this is very similar as for OSPA, but since also the
  # inner normalizations change (Factor 1/(n-1)) computation of critical values are different.
  # critical lambda-value is lam0 = 0.75*sqrt(3)/(sqrt(3)+1) = 0.4754809
  # below lam_0 identity matching is optimal, above lam_0 switching 1 and 2 is optimal for OSPA1
  lam <- 0.47549  # high
  g2 <- scaled_g(lam)
  d11 <- (2/sqrt(3))*lam     # dist from 1 to 1 (between graphs)
  d12 <- sqrt((4/3)*lam^2 - 4*lam + 4)  # dist from 1 to 2 (between graphs)
  d_exp <- d11 + 2*d12  # cost of switching 1 and 2
  # 3*d11 + 2   # cost of identity permutation
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue")
  # The contribution of each point plus the half-edges attached to it under the switch1-2 matching
  # is always smaller than the cutoff (2*vpen+0.5*number of edges incident to the point and its pi-partner),
  # as long as vpen is > (d12-1)/2 = 0.2745151
  # MORE DETAILED, phase diagram in the lambda/vpen-space: there are two critical lines. L1: lambda = lam0 = 0.75*sqrt(3)/(sqrt(3)+1) = 0.4754809,
  # and L2: vpen = (1/sqrt(3)) * lambda (intersects with L1 at (0.4754809, 0.2745151)). Above L2 we have identity
  # as best permutation to the left of L1 and 213-Permutation (switch 1 und 2) to the right (I am pretty sure that
  # above L2 nothing else can happen). Below L2 we have locally 126-Permutation (i.e. 1-->2, 2-->2, 3-->6, rest does not matter)
  # to the left *and* right of L1 (and it seems pretty certain that more stuff is happening as we make vpen smaller).
  # A lesson to remember:
  # since dummy vertices do no come with edges (alternative choice probably "they always come with edges",
  # but seems unlikely that we can have arbitrary fixed weight for dummy edges), it is clearly not worthwhile
  # to tunnel 1 and 2 (as a rule there is only cheap tunneling for (almost) isolated points if edge cost is substantial)
  vpen <- 1 # vanilla
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1:3], c(2,1,3))
  vpen <- 0.2746 # just larger than critical value
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1:3], c(2,1,3))
  vpen <- 0.2744 # just smaller than critical value
  d_exp <- 2*d11 + 2*vpen + 2
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen) 
  expect_equal(res$dist, d_exp)        
  expect_equal(res$perm[1:2], c(1,2))
  expect_true(res$perm[3] %in% c(4,5,6))  # for current compdist it is 6
  expect_equal(sort(res$perm), 1:6)
  
  lam <- 0.47547  # low
  g2 <- scaled_g(lam)
  d11 <- (2/sqrt(3))*lam     # dist from 1 to 1 (between graphs)
  d12 <- sqrt((4/3)*lam^2 - 4*lam + 4)  # dist from 1 to 2 (between graphs)
  d_exp <- 3*d11 + 2   # cost of identity permutation
  # d11 + 2*d12  # cost of switching 1 and 2
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="lightblue")
  # see above
  vpen <- 1 # vanilla
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1:3], c(1,2,3))
  vpen <- 0.2746 # just larger than "critical" value (for lam low, there should be no critical value)
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen)
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1:3], c(1,2,3))
  vpen <- 0.2744 # just smaller than critical value (for lam low, there should be no critical value)
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen) 
  expect_false(isTRUE(all.equal(res$dist, d_exp)))          
  expect_false(isTRUE(all.equal(res$perm[1:3], c(1,2,3))))
})


test_that("compdist for TT, small, different size", {
  # 2 vertices
  g2 <- spatgraph(c(-0.1,-0.1), c(sqrt(3),0), from=1, to=2)
  # plot(g1, rescale=FALSE, asp=1, xlim=c(-1,1), ylim=c(0,sqrt(3)))  
  # plot(g2, add=TRUE, rescale=FALSE, vertex.color="darkseagreen")
  
  vpen <- 1
  d_exp <- (vpen + 0.9 + 0.1) + (0 + 0 + 1)  
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1:3], c(2,3,1))
  expect_equal(sort(res$perm), c(1,2,3,4,5))
  # symmetry test
  res <- compdist(g1=g2, g2=g1, type="TT", vpen=vpen) 
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1:3], c(3,1,2))
  expect_equal(sort(res$perm), c(1,2,3,4,5))
  # hier evtl. noch weitere Tests mit anderem vpen
  
  # 1 vertex
  g2 <- spatgraph(-0.1, sqrt(3), from=numeric(), to=numeric())
  vpen <- 1
  d_exp <- (2*vpen + 0.1) + (1+1+0)
  res <- compdist(g1=g1, g2=g2, type="TT", vpen=vpen)  
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[3], 1)
  expect_equal(sort(res$perm), c(1,2,3,4))
  # symmetry test
  res <- compdist(g1=g2, g2=g1, type="TT", vpen=vpen)
  expect_equal(res$dist, d_exp)
  expect_equal(res$perm[1], 3)
  expect_equal(sort(res$perm), c(1,2,3,4))
  
  # 0 vertices
  g0 <- spatgraph(numeric(), numeric(), from=numeric(), to=numeric())
  d_exp <- 3+2
  res <- compdist(g1=g1, g2=g0, type="TT", vpen=vpen)
  expect_equal(res$dist, d_exp)
  expect_equal(sort(res$perm), c(1,2,3))
  # symmetry test
  res <- compdist(g1=g0, g2=g1, type="TT", vpen=vpen)
  expect_equal(res$dist, d_exp)
  expect_equal(sort(res$perm), c(1,2,3))
  # 1 vs 0 points
  res <- compdist(g1=g2, g2=g0, type="TT", vpen=vpen)
  expect_equal(res$dist, 1)
  expect_equal(res$perm, 1)
  # symmetry test
  res <- compdist(g1=g0, g2=g2, type="TT", vpen=vpen)
  expect_equal(res$dist, 1)
  expect_equal(res$perm, 1)
  # 0 vs 0 points
  res <- compdist(g1=g0, g2=g0, type="TT", vpen=vpen)
  expect_equal(res$dist, 0)
  expect_equal(res$perm, integer(0))
})


# set.seed(221226)
# testgraphs <- list()
# par(mfrow=c(2,3), mai=rep(0,4))
# for (n in 1:5) {
#   testgraphs[[n]] <- rspatPER(n=n+1, p=0.5)
#   plot(testgraphs[[n]], rescale=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), vertex.size=8)
#   box()
# }
# saveRDS(testgraphs, "testgraphs3.rds", compress="xz")  # smallest

testgraphs <- readRDS(file.path("_testdata/testgraphs.rds"))

test_that("compdist for OSPA1, small, different size", {
  res <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res[i,j] <- compdist(g1=testgraphs[[i]], g2=testgraphs[[j]], type="OSPA1")$dist
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
  expect_equal(res, res_expected)
  
  res <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res[i,j] <- compdist(g1=testgraphs[[i]], g2=testgraphs[[j]], type="OSPA1", CV=0.1)$dist
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


test_that("compdist for OSPA2, small, different size", {
  res <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res[i,j] <- compdist(g1=testgraphs[[i]], g2=testgraphs[[j]], type="OSPA2")$dist
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
  expect_equal(res, res_expected)
  
  res <- matrix(0,5,5)
  for (i in 1:5) {
    for (j in 1:5) {
      res[i,j] <- compdist(g1=testgraphs[[i]], g2=testgraphs[[j]], type="OSPA2", CV=0.1)$dist
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


test_that("compdist for TT, small, different size", {
  skip_on_cran()  # takes too much time
  res <- matrix(0,3,3)
  for (i in 1:3) {  # takes still too much time otherwise
    for (j in 1:3) {
      res[i,j] <- compdist(g1=testgraphs[[i]], g2=testgraphs[[j]], type="TT")$dist
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
      res[i,j] <- compdist(g1=testgraphs[[i]], g2=testgraphs[[j]], type="TT", vpen=0.1)$dist
    }
  }
  res_expected <-
    structure(c(0, 1.3721946294648613, 1.5992475902808845, 1.3721946294648613, 
                0, 1.429675403424447, 1.5992475902808845, 1.429675403424447, 
                0), dim = c(3L, 3L))
  expect_equal(res, res_expected)
})