#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
bool str_eq(std::string s1, std::string s2){
  if((s1.compare(s2)) == 0) return true;
  return false;
}

// [[Rcpp::export]]
CharacterVector get_I1(IntegerVector V4, CharacterVector  V6, IntegerVector V8, NumericVector V9, int readlength = 75) {
  int n = V4.size();
  CharacterVector I1(n);
  for(int i = 0; i<n;i++){
    if((V8[i]-V4[i]+readlength== int(V9[i]))&(  str_eq(as<std::string>(V6[i]) , "75M") )){
      I1[i]="75M";
      }
    else if((V8[i]-V4[i]+readlength!= int(V9[i]) )&( str_eq(as<std::string>(V6[i]) , "75M") )){
      I1[i]="NO75MR";
      }
    else {I1[i]="NO75ML";}
  }

  return I1;
}



// [[Rcpp::export]]
NumericVector lengthx_cpp(int dleft, int dright, IntegerVector left, IntegerVector right, NumericVector t, CharacterVector I, int readlength=75){
  int n = t.size();
  NumericVector x(n);
  for(int i = 0;i <n;i++)
  {
    if(  ( left[i]<=(dleft))   &   (right[i]>=(dright-3) )  ) {
      x[i]=t[i]-(dright-dleft+1);
      
    }else if ( ((dright-3)<left[i]) &(left[i]<=dright+readlength)& ( str_eq(as<std::string>(I[i]) , "NO75ML") ) ){
        x[i]=t[i]+18;
      
    }else if( ( (dleft-readlength+1)<right[i])&(right[i]<=dleft)&( str_eq(as<std::string>(I[i]) , "NO75MR"))){
        x[i]=right[i]-left[i]+readlength;
      
    }else if ( ((dleft)<left[i]) &  (left[i]<(dright-3))  &   (  str_eq(as<std::string>(I[i]) , "75M")  | str_eq(as<std::string>(I[i]) , "NO75ML") )) {
      x[i]=0;
    }else if(  ((dleft)<right[i])  &  (right[i]<(dright-3))  &  (str_eq(as<std::string>(I[i]) , "75M") |str_eq(as<std::string>(I[i]) , "NO75MR") )  ){
      x[i]=0;
    } else{
      x[i]=t[i];
    }
  }
  return(x);
}

// [[Rcpp::export]]
NumericVector probLL_cpp(int dleft, int dright, int l, int mu, CharacterVector I, IntegerVector left, IntegerVector right, int readlength = 75)
{
  int n = I.size();
  NumericVector y(n);
  
  for(int i = 0;i<n;i++)
  {
    if(  ((dleft-readlength+2)<=left[i]) &  (left[i]<=(dleft)) & (right[i]>dright) & ( str_eq(as<std::string>(I[i]) , "NO75ML")  )){
      y[i]=1;
    }else if(  ((dleft-readlength+2)<=right[i])  &  (right[i]<=(dleft))  &  (str_eq(as<std::string>(I[i]) , "NO75MR") )  ){
      y[i]=1;
    }else if((dleft-readlength+3)<left[i]&left[i]<(dright-3)& (str_eq(as<std::string>(I[i]),"75M") | str_eq(as<std::string>(I[i]),"NO75ML") )  ){
      y[i]=0;
    }else if((dleft-readlength+3)<right[i]&right[i]<(dright-3)& (str_eq(as<std::string>(I[i]),"75M") | str_eq(as<std::string>(I[i]),"NO75MR") ) ){
      y[i]=0;
    }else if((((dright-readlength+2)<left[i]&left[i]<(dright+2))&(right[i]>=(dright+2)))&( str_eq(as<std::string>(I[i]),"NO75ML") )){ 
      y[i]=1;
    }else if(((dright-readlength+2)<right[i]&right[i]<(dright+2))&(str_eq(as<std::string>(I[i]),"NO75MR"))){
      y[i]=1;
    }else{
      y[i]=1.0/(l-mu);
    }
  }
  return(y);
}

// [[Rcpp::export]]
NumericVector compute_p(int iter, CharacterVector V6,
                        IntegerVector left, IntegerVector right, NumericVector t, 
                        double mu = 300,  double sd = 30,  
                        int dleft1=699, int dright1=1519, int dleft2=800,int dright2=1387, 
                        int dleft3=435, int dright3=1453, int dleft4=1000, int dright4=1049,
                        int l0 = 2060, int l1=1239, int l2=1472, int l3=1041, int l4=2010
                        ){
  
  CharacterVector I = get_I1(left, V6, right, t);
  NumericVector a =lengthx_cpp(dleft1,dright1,left, right, t,I);
  NumericVector b =lengthx_cpp(dleft2,dright2,left, right, t,I);
  NumericVector c =lengthx_cpp(dleft3,dright3,left, right, t,I);
  NumericVector d =lengthx_cpp(dleft4,dright4,left, right, t,I);
  
  double p_0 = 0.2;
  double p_1 = 0.2;
  double p_2 = 0.2;
  double p_3 = 0.2;
  double p_4 = 0.2;
  int i = 0;
  
  while(i < iter) {
    int n = I.size();
    NumericVector probL0(n);
    double initp0 = 1.0/(l0-mu);
    probL0.fill(initp0);
    NumericVector probL1=probLL_cpp(dleft1,dright1,l1, mu, I,left, right);
    NumericVector probL2=probLL_cpp(dleft2,dright2,l2, mu, I,left, right);
    NumericVector probL3=probLL_cpp(dleft3,dright3,l3, mu, I,left, right);
    NumericVector probL4=probLL_cpp(dleft4,dright4,l4, mu, I,left, right);
    
    //NumericVector d = dnorm(as<NumericVector>(t), double(mu), double(sd));
    NumericVector l_0 =  p_0 *probL0*dnorm(t, double(mu), double(sd));
    NumericVector l_1 =  p_1 *probL1*dnorm(a, double(mu), double(sd));
    NumericVector l_2 =  p_2 *probL2*dnorm(b, double(mu), double(sd));
    NumericVector l_3 =  p_3 *probL3*dnorm(c, double(mu), double(sd));
    NumericVector l_4 =  p_4 *probL4*dnorm(d, double(mu), double(sd)) ;
    
    NumericVector w_0 =l_0 / (l_0 +l_1+l_2+l_3+l_4);
    NumericVector w_1 =l_1 / (l_0 +l_1+l_2+l_3+l_4);
    NumericVector w_2 =l_2 / (l_0 +l_1+l_2+l_3+l_4);
    NumericVector w_3 =l_3 / (l_0 +l_1+l_2+l_3+l_4);
    NumericVector w_4 =l_4 / (l_0 +l_1+l_2+l_3+l_4);
    
    p_0=mean(w_0);
    p_1=mean(w_1);
    p_2=mean(w_2);
    p_3=mean(w_3);
    p_4=mean(w_4);
    mu=(sum(w_0*t)+sum(w_1*a)+sum(w_2*b)+sum(w_3*c)+sum(w_4*d))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4));
    sd=sqrt( (sum(w_0* pow( (t-mu),2) )+ sum(w_1* pow( (a-mu),2) )+sum(w_2*  pow( (d-mu),2) )+ 
      sum(w_3* pow( (c-mu),2)  )+sum(w_4*  pow( (d-mu),2) ))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4)));
  Rcout << "\n\n i  "  << i;
  i++;
  }
  
  return(NumericVector::create(p_0, p_1, p_2, p_3, p_4));
}

