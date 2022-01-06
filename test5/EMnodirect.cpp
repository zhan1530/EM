#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector probLL_cpp(int dleft, int dright, int l, double mu, IntegerVector left, IntegerVector right, double a,double b,int readlength = 75)
{
  int n = left.size();
  NumericVector y(n);
  
  for(int i = 0;i<n;i++)
  {
    if(  ((dleft-readlength+2)<left[i]) &  (left[i]<(dright-2))){
      y[i]=0;
    }else if(  ((dleft-readlength+2)<right[i]) & (right[i]<(dright-2))){
      y[i]=0;
    }else{
      y[i]=1.0/((l-mu)*(1-(a+b)));
    }
  }
  return(y);
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

