#include <Rcpp.h>
using namespace Rcpp;


// legacy function
// [[Rcpp::export]]
NumericMatrix crossdist2(NumericMatrix x, NumericMatrix y) {
  
  int nrow_x = x.nrow();
  int nrow_y = y.nrow();
  
  NumericMatrix res(nrow_x, nrow_y);
  
  for (int i = 0; i < nrow_x; i++) {
    for (int j = 0; j < nrow_y; j++) {
      double d = std::sqrt(std::pow(x(i, 0) - y(j, 0), 2) + std::pow(x(i, 1) - y(j, 1), 2));
      res(i,j) = d;
    }
  }

  return res;
}


// microbenchmark for k=2 does not give systematically worse performance than crossdist2
// Of course, if we needed really high speed crossdist is probably best computed with a lin alg package
// Compute Matrix of Distances Between Two Sequences of Points
//
// @param x,y k-column matrices specifying the k-dimensional coordinates of the points.
//
// @return a matrix of size \code{dim(x)[1]} times \code{dim(y)[1]}
// @export
//
// @examples
// x <- matrix(runif(8), 4, 2)
// y <- matrix(runif(10), 5, 2)
// crossdist(x, y)
//
// [[Rcpp::export]]
NumericMatrix crossdist(NumericMatrix x, NumericMatrix y) {
  
  int nrow_x = x.nrow();
  int nrow_y = y.nrow();
  int ncol = x.ncol();
  
  NumericMatrix res(nrow_x, nrow_y);
  
  for (int i = 0; i < nrow_x; i++) {
    for (int j = 0; j < nrow_y; j++) {
      double d = 0;
      for (int k = 0; k < ncol; k++) {
        d += std::pow(x(i, k) - y(j, k), 2);
      }
      res(i,j) = std::sqrt(d);
    }
  }
  
  return res;
}



// [[Rcpp::export]]
NumericVector edgedist(NumericVector e, NumericVector f, double maxedist = 1.0, bool bound = true, bool edge_weighted = false) {
  int n = e.length();
  NumericVector res(n);
  
  if (!bound) {
    res = abs(e - f);
  } else {
      // for (int i = 0; i < n; i++) {
      // if (!(e[i] == 0 || e[i] == 0.5 || e[i] == 1) || !(f[i] == 0 || f[i] == 0.5 || f[i] == 1) || e[i] * f[i] == 0.25) {
      //  stop("Invalid input values");
      //}
      // }
      if(edge_weighted){
        std::vector<double> res(n);
        std::transform (e.begin(), e.end(), f.begin(), std::back_inserter(res), [](double a, double b) {return std::min(a,b);}); 
      }
      else{
        res = maxedist * abs(e - f);
      }
  }
  return res;
}


// [[Rcpp::export]]
double edgedist1(double& e, double& f, double maxedist = 1.0, bool bound = true) {
  double res = 0.0;
  
  if (!bound) {
    res = abs(e - f);
  } else {
    // for (int i = 0; i < n; i++) {
    // if (!(e[i] == 0 || e[i] == 0.5 || e[i] == 1) || !(f[i] == 0 || f[i] == 0.5 || f[i] == 1) || e[i] * f[i] == 0.25) {
    //  stop("Invalid input values");
    //}
    // }
    res = maxedist * abs(e - f);
  }
  return res;
}


// This function still takes up 80% of the computation time of auction_match_Cpp
// Overall I think auction_match_Cpp could still run around 2-4 times as fast if we computed edgediff[j]
// only for one j and add/subtract the changes for other j. 
// The necessary investment in programming time currently does not seem worthwhile.
// For compute_pers_values2 (with compensation) both programming time and speed-up will probably increase
// [[Rcpp::export]]
NumericVector compute_pers_values(int i, int nn, IntegerVector pers_to_obj, IntegerVector obj_to_pers,
                                  double maxcost, NumericMatrix dmat, NumericMatrix sigma, NumericMatrix tau,
                                  NumericVector prices, double maxedist, double multe, bool bound, double vfac, double efac) {
  
  int j, k, l;
  NumericVector edgediff(nn);
  NumericVector persval(nn, maxcost);   // result
  IntegerVector pers_to_obj_hypo = clone(pers_to_obj);
  IntegerVector obj_to_pers_hypo = clone(obj_to_pers);
  //Rcout << pers_to_obj_hypo << "\n";
  //Rcout << obj_to_pers_hypo << "\n";
  int oldj = pers_to_obj[i];
  if (!IntegerVector::is_na(oldj)) {
    obj_to_pers_hypo[oldj] = NA_INTEGER; 
      // setting this happens every time we set i to some j, so we never have to undo this
      // unless if oldj = j but then it happens automatically (in that case oldi in NA now)
  }

  for (j = 0; j < nn; j++) {  // investigate values for matching person i to obj j
    // switch i to obj j, taking it away from previous bidder (if any)
    int oldi = obj_to_pers_hypo[j];  // _hypo is just for a tiny difference in efficiency (case j = oldj, saves some switching)
    bool takeaway_j = !IntegerVector::is_na(oldi);
    if (takeaway_j) {  // if obj j is already matched
      pers_to_obj_hypo[oldi] = NA_INTEGER;
    }
    pers_to_obj_hypo[i] = j;
    obj_to_pers_hypo[j] = i;

    double edgediff_temp = 0.0;
    double unass_sig_temp = 0.0;
    double unass_tau_temp = 0.0;
    for (k = 0; k < nn; k++) {
      if (IntegerVector::is_na(pers_to_obj_hypo[k])) {  // k is unassigned
        unass_sig_temp += sigma(k,i);  // sigma is col major
      } else {
        edgediff_temp += multe * abs(sigma(k,i) - tau(pers_to_obj_hypo[k],j));
      } 
    }
    for (l = 0; l < nn; l++) {
      if (IntegerVector::is_na(obj_to_pers_hypo[l])) {
        unass_tau_temp += tau(l,j);   // tau is col major
      }
    }
      
    edgediff[j] = edgediff_temp + multe * abs(unass_sig_temp - unass_tau_temp);

    // undo the switch (pers_to_obj_hypo[i] does not matter, obj_to_pers_hypo[oldj] is NA anyway)
    if (takeaway_j) {  // if obj j is already matched
      pers_to_obj_hypo[oldi] = j;
    }
    obj_to_pers_hypo[j] = oldi;
  }  
  
  persval = persval - vfac * dmat(i,_) - efac * edgediff - prices;

  return(persval);
}


// [[Rcpp::export]]
NumericVector compute_pers_values2(int i, int nn, IntegerVector pers_to_obj, IntegerVector obj_to_pers,
                                  double maxcost, NumericMatrix dmat, NumericMatrix sigma, NumericMatrix tau,
                                  NumericVector prices, double maxedist, double multe, bool bound, double vfac, double efac) {
  
  int j, k, l;
  NumericVector edgediff(nn);
  NumericVector transferpay(nn); 
  NumericVector persval(nn, maxcost);   // result
  IntegerVector pers_to_obj_hypo = clone(pers_to_obj);
  IntegerVector obj_to_pers_hypo = clone(obj_to_pers);
  //Rcout << pers_to_obj_hypo << "\n";
  //Rcout << obj_to_pers_hypo << "\n";
  int oldj = pers_to_obj[i];
  if (!IntegerVector::is_na(oldj)) {
    obj_to_pers_hypo[oldj] = NA_INTEGER; 
    // setting this happens every time we set i to some j, so we never have to undo this
    // unless if oldj = j but then it happens automatically (in that case oldi in NA now)
  }
  
  for (j = 0; j < nn; j++) {  // investigate values for matching person i to obj j
    // switch i to obj j, taking it away from previous bidder (if any)
    int oldi = obj_to_pers_hypo[j];  // _hypo is just for a tiny difference in efficiency (case j = oldj, saves some switching)
    bool takeaway_j = !IntegerVector::is_na(oldi);
    if (takeaway_j) {  // if obj j is already matched
      pers_to_obj_hypo[oldi] = NA_INTEGER;
    }
    pers_to_obj_hypo[i] = j;
    obj_to_pers_hypo[j] = i;
    
    double edgediff_temp = 0.0;
    double unass_sig_temp = 0.0;
    double unass_tau_temp = 0.0;
    for (k = 0; k < nn; k++) {
      if (IntegerVector::is_na(pers_to_obj_hypo[k])) {  // k is unassigned
        unass_sig_temp += sigma(k,i);  // sigma is col major
      } else {
        edgediff_temp += multe * abs(sigma(k,i) - tau(pers_to_obj_hypo[k],j));
      } 
    }

    for (l = 0; l < nn; l++) {
      if (IntegerVector::is_na(obj_to_pers_hypo[l])) {
        unass_tau_temp += tau(l,j);   // tau is col major
      }
    }
    
    edgediff_temp += multe * abs(unass_sig_temp - unass_tau_temp);
    
    //Rcout << edgediff_temp << "  edgediff only " << i << j << k << "\n";
    
    // externalities (we continue to add to edgediff_temp)
    // ----------------------------------------------------
    // there are two slightly different versions that only differ in the order they do things:
    // version 2 is typically quite a bit faster and gives the better results, which are
    // the same as in the R-version (also the order of things is almost exactly the same).
    // version 1 sometimes gives better solutions, especially if the graphs are of very unequal sizes
    // and we use OSPA1.
    
    // if (version == 1) {
    //   for (k = 0; k < nn; k++) {
    //     if (!IntegerVector::is_na(pers_to_obj_hypo[k]) && k != i) {
    //       if (!IntegerVector::is_na(oldj)) {
    //         // at least for OSPA2 and TT (where there are no 0.5 edges) this could be simplified
    //         // probably also for OSPA1 but less easily
    //         edgediff_temp += multe * (abs(sigma(k,i) - tau(pers_to_obj_hypo[k],j)) -
    //           abs(sigma(k,i) - tau(pers_to_obj_hypo[k],oldj)));
    //       }
    //       // Rcout << edgediff_temp << " edgediff + oldj " << i << j << k << "\n";
    //       if (!IntegerVector::is_na(oldi)) {
    //         edgediff_temp += multe * (abs(sigma(k,i) - tau(pers_to_obj_hypo[k],j)) -
    //           abs(sigma(k,oldi) - tau(pers_to_obj_hypo[k],j)));
    //       }
    //     }
    //   } 
    // }
    
    edgediff[j] = edgediff_temp;  
    
    // if (version == 2) {
    double firstsum = 0.0;
    if (!IntegerVector::is_na(oldj)) {
      for (k = 0; k < nn; k++) {
        if (!IntegerVector::is_na(pers_to_obj_hypo[k]) && k != i) {
          // at least for OSPA2 and TT (where there are no 0.5 edges) this could be simplified
          // probably also for OSPA1 but less easily
          firstsum += multe * (abs(sigma(k,i) - tau(pers_to_obj_hypo[k],j)) -
                               abs(sigma(k,i) - tau(pers_to_obj_hypo[k],oldj)));
        }
      }
    }
    double secondsum = 0.0;
    if (!IntegerVector::is_na(oldi)) {
      for (k = 0; k < nn; k++) {
        if (!IntegerVector::is_na(pers_to_obj_hypo[k]) && k != i) {
          secondsum += multe * (abs(sigma(k,i) - tau(pers_to_obj_hypo[k],j)) -
                                abs(sigma(k,oldi) - tau(pers_to_obj_hypo[k],j)));
        }
      }
    } 
    
    transferpay[j] = firstsum + secondsum;
    //}
    // Rcout << edgediff_temp << " edgediff + oldj + oldi " << i << j << k << "\n";
    // if (i==2 && j==2) stop("safety stop");

    // undo the switch (pers_to_obj_hypo[i] does not matter, obj_to_pers_hypo[oldj] is NA anyway)
    if (takeaway_j) {  // if obj j is already matched
      pers_to_obj_hypo[oldi] = j;
    }
    obj_to_pers_hypo[j] = oldi;
  }  
  
  persval = persval - vfac * dmat(i,_) - efac * edgediff - prices - efac * transferpay;
  
  return(persval);
}



// [[Rcpp::export]]
List auction_match_Cpp(NumericMatrix sigma, NumericMatrix tau, NumericMatrix dmat, int nn, double maxcost, bool compensate,
                       int numtype, double CE = 1.0, int stop_at = 1, int maxiter = 100, double eps = 0.01, int verbose = 0) {

  // Rcpp::Clock clock;
  
  // clock.tick("auction_match");
  IntegerVector best_assignment(nn, NA_INTEGER); // return value
  double best_cost = R_PosInf;                   // return value
  
  IntegerVector pers_to_obj(nn, NA_INTEGER);
  IntegerVector obj_to_pers(nn, NA_INTEGER);
  NumericVector prices(nn);
  NumericMatrix pers_values(nn, nn);
  
  // extra quantities depending on numtype (initialized for numtype=0)
  double maxedist = R_PosInf;  // names are unfortunate: this is where we cut off edges
  double multe = 1.0;  // this is by what we multiply difference of number of edges to unmatched 
  bool bound = false;
  double vfac = 1.0, efac = 0.5;  // division by 2 included in efac, saves a little time
  if(numtype == 1 || numtype == 2){
    maxedist = CE;
    multe = CE;
    bound = true;
    vfac = 1.0/nn;
    efac = 0.5/(nn*(nn-1));
  } 
  
  int iter = 1;
  int fullmatch_num = 0;
  bool fullmatch = false;
  std::list<std::vector<int> > assignments;  // Rcpp::List container has no push_back function
  NumericVector costs(0);  // Vector of length 0
  bool change = true;
  
  
  while (change && iter <= maxiter && fullmatch_num < stop_at) {
    if (verbose > 0) Rcout << "\n----------------------------------\n" << iter << "\n";
    change = false;
    
    for (int bidder = 0; bidder < nn; bidder++) {
      NumericVector my_values(nn);
      // in what follows we pass maxedist, bound, vfac, efac instead of numtype and CE to avoid repeating code
      // 2300424: it seems multe is now all we need.
      if (compensate) {
        my_values = pers_values(bidder, _) = compute_pers_values2(bidder, nn, pers_to_obj, obj_to_pers, maxcost,
                                                                  dmat, sigma, tau, prices, maxedist, multe, bound, vfac, efac);
      } else {
        // clock.tick("compute_pers_values");
        my_values = pers_values(bidder, _) = compute_pers_values(bidder, nn, pers_to_obj, obj_to_pers, maxcost,
                                                                 dmat, sigma, tau, prices, maxedist, multe, bound, vfac, efac);
        // clock.tock("compute_pers_values");
      }
      
      std::vector<double> v(my_values.begin(), my_values.end());
      std::nth_element(v.begin(), v.begin() + 1, v.end(), std::greater<double>{});
      // v[0] is largest element, v[1] second largest element  # these are references!!
      // std::nth_element(v.begin(), v.end() - 2, v.end()); is much easier, but according to cppreference
      // it is only guaranteed that the values to the left of v.end() - 2 are <= (and hence not that
      // the one value to the right of v.end() - 2 is >=; of course it will be automatically unless
      // there is a tie for the largest values)
      // std::nth_element(v.begin(), v.begin() + 1, v.end(), std::greater<NumericVector>{}); would also
      // be easier, but does not compile (it seems std::greater<NumericVector> is not implemented)
      
      // feels as if we (partially) sort twice, but this is probably just a pairwise comparison
      int bidfor = which_max(my_values);
      int oldbidfor = pers_to_obj[bidder];
      if (NumericVector::is_na(oldbidfor) || oldbidfor != bidfor) {
        change = true;
        if (!IntegerVector::is_na(oldbidfor)) {
          obj_to_pers[oldbidfor] = NA_INTEGER; 
          // IntegerVector::get_na() returns the NA for this Vector class; should be equivalent to = NA_INTEGER
        }
        double oldbidder = obj_to_pers[bidfor];
        if (!IntegerVector::is_na(oldbidder)) {
          pers_to_obj[oldbidder] = NA_INTEGER;
          // IntegerVector::get_na() returns the NA for this Vector class; should be equivalent to = NA_INTEGER
        }
        pers_to_obj[bidder] = bidfor;
        obj_to_pers[bidfor] = bidder;
        double bidincrement = v[0] - v[1] + eps;
        prices[bidfor] += bidincrement;
        
        int unassigned = sum(is_na(pers_to_obj));
        
        if (unassigned <= 1) {
          IntegerVector curpi = clone(pers_to_obj);
          if (unassigned == 1) {
            IntegerVector nnseq = seq(0,nn-1); // needed; inline in the next line does not work
            curpi[is_na(curpi)] = setdiff(nnseq, curpi);
            fullmatch = false;
          } else {
            if (!fullmatch) {
              fullmatch = true;
              fullmatch_num++;
            }
          }
          
          double tempcost = 0.0;
          NumericVector tempdvec(nn);
          NumericMatrix tempemat(nn);
          for (int i = 0; i < nn; i++) {
            tempdvec[i] = dmat(i, curpi[i]);
            for (int j = 0; j < nn; j++) {
              tempemat(i,j) = multe * abs(sigma(i,j) - tau(curpi[i],curpi[j]));
            }
          }

          tempcost = vfac * sum(tempdvec) + efac * sum(tempemat);
          // Rcout << curpi << "\n";
          // Rcout << tempdvec << "\n";
          // Rf_PrintValue(tempemat);
          // Rcout << vfac << ", " << efac << ", " << nn << "\n";
          // Rcout << tempcost << " = " << sum(tempdvec) << " + "<< sum(tempemat) << "\n";
          // stop("safety stop");
          
          std::vector<int> v = as<std::vector<int> >(curpi); 
          assignments.push_back(v);
          costs.push_back(tempcost);
          if (tempcost < best_cost) {
            best_cost = tempcost;
            best_assignment = clone(curpi);
          }
        } else {
          fullmatch = false;
        }
        
        if (verbose > 0) {
          Rcout << bidder << " for " << bidfor << "\n";
          if (verbose > 1) Rcout << pers_values;
          Rcout << pers_to_obj << "\n";
          Rcout << prices << "\n";
        }
      }
    }
    
    iter++;
  }
  
  if (verbose > 0) {
    Rcout << "--------------------------------------------------------------------------\n";
    Rcout << "Done at change = " << change << ", iter = " << iter << ", full matches = " << fullmatch_num << "\n";
    Rcout << "--------------------------------------------------------------------------\n";
    for (std::vector<int> aa : assignments) { // range based loop; C++11 feature
      for (int a : aa) {                      // range based loop; C++11 feature
        Rcout << a << " ";
      }
      Rcout << "\n";
    }
    Rcout << costs << "\n";
    Rcout << "best_assignment: " << best_assignment << "\n";
    Rcout << "best_cost: " << best_cost << "\n";
  }
  int nearmatch_num = costs.length();
  
  // clock.tock("auction_match");
  // clock.stop("auctiontimes");
  return List::create(Named("best_assignment") = best_assignment, Named("best_cost") = best_cost,
                      Named("change") = change, Named("iter") = iter, Named("fullmatch_num") = fullmatch_num,
                      Named("nearmatch_num") = nearmatch_num); 
}    

