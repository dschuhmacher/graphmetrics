#include <Rcpp.h>
using namespace Rcpp;


double compLstat(const NumericMatrix& dmat, const int& n1, const int& n2, const int& n, const int& N1, const int& N2, const int& N,
                 const bool& balanced, const IntegerVector& index1, const IntegerVector& index2) {
  double sum1 = 0; 
  double sum2 = 0;
  for (int i = 0; i < (n1-1); i++){
    for (int j = i+1; j < n1; j++)
    {
      sum1 += dmat(index1[i],index1[j]);
    }
  }
  for (int i = 0; i < (n2-1); i++){
    for (int j = i+1; j < n2; j++)
    {
      sum2 += dmat(index2[i],index2[j]);
    }
  }
  
  double mean1 = sum1/N1; 
  double mean2 = sum2/N2; 
  // Rcout << "mean1 " << mean1 << " mean2 " << mean2 << "\n";
  
  double enumerator = 0.0;
  if (balanced) {
    double meantot = (sum1 + sum2)/(N1+N2); 
    enumerator = n1*(mean1-meantot)*(mean1-meantot) + n2*(mean2-meantot)*(mean2-meantot);
  } else {
    enumerator = (n1*n2/n)*(mean1-mean2)*(mean1-mean2);
  }
  // Rcout << "enum " << enumerator << "\n";
  
  double denominator = 0;
  for (int i = 0; i < (n1-1); i++){
    for (int j = i+1; j < n1; j++)
    {
      double a = dmat(index1[i],index1[j])-mean1;
      denominator += a*a;
    }
  }
  for (int i = 0; i < (n2-1); i++){
    for (int j = i+1; j < n2; j++)
    {
      double a = dmat(index2[i],index2[j])-mean2;
      denominator += a*a;
    }
  }
  
  // Rcout << "denom " << denominator << "\n";
  return(enumerator/denominator);   // leaving away the factor (N-k)/(k-1) = N-2 !!
}


// [[Rcpp::export]]
List msm_levene_Cpp(NumericMatrix dmat, int n1, int n2, bool balanced = true, int nperm = 999, double alpha = 0.05) {
  // We would have to divide dmat by two to do exactly the same as in the paper
  // However this 1/2 already cancels on the level of Lstat (and then does not affect the test anyway)
  
  int n = n1+n2;
  int N1 = n1*(n1-1)/2;
  int N2 = n2*(n2-1)/2;
  int N = N1+N2;
  
  IntegerVector index1 = seq(0,n1-1);
  IntegerVector index2 = seq(n1,n-1);
  
  double Lstat = compLstat(dmat, n1, n2, n, N1, N2, N, balanced, index1, index2);
  // Rcout << "denom " << residuals << " Stat " << stat << \n";
  
  if (nperm == 0) {
    List res = List::create( Named("Lstat") = Lstat );
    return(res);
  } else {
    IntegerVector allindices = seq(0,n-1);
    NumericVector permstats(nperm);
    for (int k = 0; k < nperm; k++)
    {
      IntegerVector sampleindex = sample(allindices, n, false);
      IntegerVector newind1 = sampleindex[index1];
      IntegerVector newind2 = sampleindex[index2];

      permstats[k] = compLstat(dmat, n1, n2, n, N1, N2, N, balanced, newind1, newind2);
    }

    LogicalVector larger = (permstats >= Lstat); 
    double pvalue = (sum(larger) + 1.0)/(nperm + 1.0);
    bool reject = (pvalue <= alpha);
    Lstat *= (N-2);  //just for the output; irrelevant for the permutations
                     //should be also /(k-1), but k-1=1
    List res = List::create( Named("Lstat") = Lstat, Named("pval") = pvalue, Named("reject") = reject );
    return(res);
  }
}



