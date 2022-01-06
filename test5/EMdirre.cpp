#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector lengthx_cpp1(int dleft, int dright, IntegerVector L_V4_1, IntegerVector L_V6m_1, NumericVector L_V4_2, IntegerVector L_V6m_2,IntegerVector R_V4_1, IntegerVector R_V6m_1, NumericVector R_V4_2, IntegerVector R_V6m_2, IntegerVector L_V6t, IntegerVector R_V6t){
  int n = L_V4_1.size();
  NumericVector x(n);
  for(int i = 0;i <n;i++)
  {
    
    if(  ( dleft< (L_V4_1[i]+L_V6m_1 [i]-2)) & (L_V4_1[i]<(dright-3)) ){
      x[i]=0;
    }
    else if( (!NumericVector::is_na(L_V4_2[i])) &( dleft< (L_V4_2[i]+L_V6m_2 [i]-2)) & (L_V4_2[i]< (dright-3))){
      x[i]=0;
    }
    else if (( dleft< (R_V4_1[i]+R_V6m_1 [i]-2)) & (R_V4_1[i]<(dright-3))  ){
      x[i]=0;
    }
    else if ((!NumericVector::is_na(R_V4_2[i]))&( (dleft) <( R_V4_2[i]+R_V6m_2 [i]-2)) & (R_V4_2[i]<(dright-3))){
      x[i]=0;
    }
//#####not supply but split########
    else if ( (NumericVector::is_na(L_V4_2[i])) &(L_V6m_1[i]!=L_V6t[i])&((dright-3)<=L_V4_1 [i])&(L_V4_1[i]<(dright+3))) {
      x[i]=R_V4_1[i]-L_V4_1[i]+R_V6t[i]+(L_V6t[i]-L_V6m_1[i]);
    }
    else if( (NumericVector::is_na(R_V4_2[i])) &(R_V6m_1[i]!=R_V6t[i]) & ((dright-3)<=R_V4_1[i])&(R_V4_1[i]<(dright+3))){
      x[i]=R_V4_1[i]-L_V4_1[i]-(dright-dleft+1)+R_V6m_1[i];
    }
    
    else if(L_V4_1[i]<=(dleft)&R_V4_1[i]>=(dright-3)){
      x[i]=R_V4_1[i]-L_V4_1[i]+R_V6t[i]-(dright-dleft+1);
    }
    else{
      x[i]=R_V4_1[i]-L_V4_1[i]+R_V6t[i];
    }
  }
  return(x);
}



// [[Rcpp::export]]
NumericVector probLL_cpp(int dleft, int dright, int l, int mu, IntegerVector L_V4_1, IntegerVector L_V6m_1, NumericVector L_V4_2, IntegerVector L_V6m_2,IntegerVector R_V4_1, IntegerVector R_V6m_1, NumericVector R_V4_2, IntegerVector R_V6m_2, IntegerVector L_V6t, IntegerVector R_V6t)
{
  int n = L_V4_1.size();
  NumericVector y(n);
  
  for(int i = 0;i<n;i++)
  { 
    if ((!NumericVector::is_na(L_V4_2[i])) &( (dleft-2)<(L_V4_1[i]+L_V6m_1[i]))&((L_V4_1[i]+L_V6m_1[i])< (dleft+2))
          &( (dright-2) <L_V4_2[i])&(L_V4_2[i]< (dright+2) ))
    {y[i]=1;}
    else if ((!NumericVector::is_na(R_V4_2[i])) &( (dleft-2)<(R_V4_1[i]+R_V6m_1[i]))&((R_V4_1[i]+R_V6m_1[i])<(dleft+2))
               &((dright-2)<R_V4_2[i])&(R_V4_2[i]<(dright+2)))
    {y[i]=1;}
//#####no supply L but split, report Left of left read########
    else if( (NumericVector::is_na(L_V4_2[i])) & (L_V6m_1[i]!=L_V6t[i])& ((dleft-2)<(L_V4_1[i]+L_V6m_1[i]))&((L_V4_1[i]+L_V6m_1[i])<(dleft+2))
               &(R_V4_1[i]>dright))
    {y[i]=1;}
//#####no supply L but split, report right of left read########
    else if((NumericVector::is_na(L_V4_2[i])) & (L_V6m_1[i]!=L_V6t[i])&((dright-2)<L_V4_1[i])&(L_V4_1[i]<(dright+2)) )
    {y[i]=1;}
//#####no supply R but split, report Left of right read########
    
    else if( (NumericVector::is_na(R_V4_2[i])) & (R_V6m_1[i]!=R_V6t[i])& ((dleft-2)<(R_V4_1[i]+R_V6m_1[i]))&((R_V4_1[i]+R_V6m_1[i])<(dleft+2))
               &(L_V4_1[i]<dleft))
    {y[i]=1;}
//#####no supply R but split, report right of right read########
    else if((NumericVector::is_na(R_V4_2[i]))&(R_V6m_1[i]!=R_V6t[i])&((dright-2)<R_V4_1[i])&(R_V4_1[i]<(dright+2)) )
    {y[i]=1;}
    
//#####start prob =0#########
    
    else if(  ( dleft < (L_V4_1[i]+L_V6m_1[i]-2)) & (L_V4_1[i]<(dright-3)) ){
      y[i]=0;
    }
    else if( (!NumericVector::is_na(L_V4_2[i]))&( dleft < (L_V4_2[i]+L_V6m_2[i]-2)) & (L_V4_2[i]< (dright-3))){
      y[i]=0;
    }
    else if (( dleft< (R_V4_1[i]+(R_V6m_1[i]-2))) & (R_V4_1[i]<(dright-3))  ){
      y[i]=0;
    }
    else if ((!NumericVector::is_na(R_V4_2[i])) & ( dleft < (R_V4_2[i]+R_V6m_2[i]-2)) & (R_V4_2[i]<(dright-3))){
      y[i]=0;
    }
    
    else 
    {y[i]=1.0/(l-mu);}
  }
  return(y);
}
