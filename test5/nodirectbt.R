library(boot)
library(Rcpp)
sourceCpp("EMnodirect.cpp")

l0=2060
l1=1239
l2=1472
l3=1041
l4=2010


dleft1=699
dright1=1519
dleft2=800
dright2=1387
dleft3=435
dright3=1453
dleft4=1000
dright4=1049

readlength=75
data1=read.table("aln.sam", fill = T, col.names = paste("V", 1:19, sep = ""),skip = 2)

prob=function(data){
  perfectmactch=function(data){
    data0=data[,c(1,4,6,8,9)]
    datanew=subset(data0,V9>0)
    datanewnosplit0=subset(datanew,V6=="75M")#get rid of left end be splited
    datanewnosplit=subset(datanewnosplit0,datanewnosplit0$V8-
                            datanewnosplit0$V4+readlength==datanewnosplit0$V9)
    #get rid of right end be splited, do not have because deletion is 300-500
    return(datanewnosplit)
  }
  datamix=perfectmactch(data)
  
  left=datamix$V4
  right=datamix$V8
  t=datamix$V9
  x=rep(0,nrow(datamix))
  
  #calculate length x
  lengthx=function(dleft,dright){
    x=rep(0,nrow(datamix))
    for(i in 1:length(left))
    {
      if(left[i]<=(dleft-readlength+2)&right[i]>=(dright-2))
      {x[i]=right[i]-left[i]-(dright-dleft+1)+readlength}
      
      else if ((dleft-readlength+2)<left[i]&left[i]<(dright-2))
      {x[i]=0}
      
      else if((dleft-readlength+2)<right[i]&right[i]<(dright-2))
      {x[i]=0}
      
      else
      {x[i]=right[i]-left[i]+readlength}
    }
    return(x)
  }
  
  a=lengthx(dleft1,dright1)
  b=lengthx(dleft2,dright2)
  c=lengthx(dleft3,dright3)
  d=lengthx(dleft4,dright4)
  
  
  h=data.frame(t,a,b,c,d,left,right)
  #rightimpossible function
  rightimpossible=function(dleft,l,mu,sd){
    x=rep(0,dleft-readlength)
    for(i in 1:dleft-readlength)
    {
      x[i]=(pnorm(dleft-i+readlength,mu,sd)-pnorm(dleft-i,mu,sd))*(1/(l-mu))
    }
    return(sum(x))
  }
  #leftimpossible function
  leftimpossible=function(dright,l,mu,sd){
    y=rep(0,readlength)
    for(i in 1:readlength)
    {
      y[i]=(pnorm(l0-dright+readlength-i,mu,sd))*(1/(l-mu))
    }
    return(sum(y))
  }
  
  probL=function(dleft,dright,l,a,b){
    for(i in 1:nrow(datamix))
    {
      if((dleft-readlength+2)<left[i]&left[i]<(dright-2))
      {y[i]=0}
      
      else if((dleft-readlength+2)<right[i]&right[i]<(dright-2))
      {y[i]=0}
      
      else
      {y[i]=1/((l-mu)*(1-(a+b)))}
    }
    return(y)
  }
  
  t=h$t
  a=h$a
  b=h$b
  c=h$c
  d=h$d
  left=h$left
  right=h$right
  mu=300
  sd=30
  p_0 = 1/5
  p_1 = 1/5
  p_2 = 1/5
  p_3 = 1/5
  p_4 = 1/5
  
  y=rep(0,nrow(datamix))
  
  
  for( i in 1:100) {
    a1=rightimpossible(dleft1,l1,mu,sd)
    a2=leftimpossible(dright1,l1,mu,sd)
    b1=rightimpossible(dleft2,l2,mu,sd)
    b2=leftimpossible(dright2,l2,mu,sd)
    c1=rightimpossible(dleft3,l3,mu,sd)
    c2=leftimpossible(dright3,l3,mu,sd)
    d1=rightimpossible(dleft4,l4,mu,sd)
    d2=leftimpossible(dright4,l4,mu,sd)
    
    probL0=rep(1/(l0-mu),nrow(datamix))
    probL1=probLL_cpp(dleft1,dright1,l1, mu, left, right,a1,a2);
    probL2=probLL_cpp(dleft2,dright2,l2, mu, left, right,b1,b2);
    probL3=probLL_cpp(dleft3,dright3,l3, mu, left, right,c1,c2);
    probL4=probLL_cpp(dleft4,dright4,l4, mu, left, right,d1,d2);
    
    
    l_0 = p_0 *probL0* dnorm( h$t, mu, sd)
    l_1 = p_1 *probL1* dnorm( h$a, mu, sd)
    l_2 = p_2 *probL2* dnorm( h$b, mu, sd)
    l_3 = p_3 *probL3* dnorm( h$c, mu, sd)
    l_4 = p_4 *probL4* dnorm( h$d, mu, sd)
    
    
    w_0 =l_0 / (l_0 +l_1+l_2+l_3+l_4)
    w_1 =l_1 / (l_0 +l_1+l_2+l_3+l_4)
    w_2 =l_2 / (l_0 +l_1+l_2+l_3+l_4)
    w_3 =l_3 / (l_0 +l_1+l_2+l_3+l_4)
    w_4 =l_4 / (l_0 +l_1+l_2+l_3+l_4)
    
    
    p_0=mean(w_0)
    p_1=mean(w_1)
    p_2=mean(w_2)
    p_3=mean(w_3)
    p_4=mean(w_4)
    
    
    mu=(sum(w_0*h$t)+sum(w_1*h$a)+sum(w_2*h$b)+sum(w_3*h$c)+sum(w_4*h$d))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4))
    sd=sqrt((sum(w_0*(h$t-mu)^2)+sum(w_1*(h$a-mu)^2)+sum(w_2*(h$b-mu)^2)+sum(w_3*(h$c-mu)^2)+sum(w_4*(h$d-mu)^2))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4)))
   # print(i)
    #print(mu)
   # print(sd)
    #print(c(p_0,p_1,p_2,p_3,p_4))
  }
  original0=p_0*nrow(h)
  original1=p_1*nrow(h)/(1-(a1+a2))
  original2=p_2*nrow(h)/(1-(b1+b2))
  original3=p_3*nrow(h)/(1-(c1+c2))
  original4=p_4*nrow(h)/(1-(d1+d2))
  
  
  p0=original0/(original0+original1+original2+original3+original4)
  p1=original1/(original0+original1+original2+original3+original4)
  p2=original2/(original0+original1+original2+original3+original4)
  p3=original3/(original0+original1+original2+original3+original4)
  p4=original4/(original0+original1+original2+original3+original4)
  
  p_00=(p0/l0)/(p0/l0+p1/l1+p2/l2+p3/l3+p4/l4)
  p_11=(p1/l1)/(p0/l0+p1/l1+p2/l2+p3/l3+p4/l4)
  p_22=(p2/l2)/(p0/l0+p1/l1+p2/l2+p3/l3+p4/l4)
  p_33=(p3/l3)/(p0/l0+p1/l1+p2/l2+p3/l3+p4/l4)
  p_44=(p4/l4)/(p0/l0+p1/l1+p2/l2+p3/l3+p4/l4)
  result=cbind(p_00,p_11,p_22,p_33,p_44,mu,sd)
  return(result)
}
rsq=function(data, indices){
  d=data[indices,]
  return(prob(d))
}

h.boot=boot(data1,rsq,R=100)
bias=colMeans(h.boot$t)-h.boot$t0[1,]
citmean=matrix(NA,7,2)
citnorm=matrix(NA,7,2)
citbasic=matrix(NA,7,2)
citperc=matrix(NA,7,2)

sd=array(NA,7)
for( i in 1:7){
  sd[i]=sd(h.boot$t[,i]) 
  citmean[i,]=t.test(h.boot$t[,i])$conf.int
  citnorm[i,]=boot.ci(h.boot,index=i,type="norm")$norm[1,2:3]
  citbasic[i,]=boot.ci(h.boot,index=i,type="basic")$basic[1,4:5]
  citperc[i,]=boot.ci(h.boot,index=i,type="perc")$perc[1,4:5]
  
}
mean=colMeans(h.boot$t)
a=data.frame(h.boot$t0[1,],mean,bias,sd,citmean,citnorm,citbasic,citperc)
colnames(a)=c("original","mean","bias","sd","meanlow","meanhigh","nlow","nhigh","blow","bhigh","plow","phigh")
colnames(h.boot$t)=c("p0","p1","p2","p3","p4","mu","sd")
plot(h.boot)
write.csv(a, file = "nodirectbt.csv")
write.csv(h.boot$t, file = "nodirectbt_t.csv")
