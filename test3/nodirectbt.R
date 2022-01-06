library(boot)
library(Rcpp)
sourceCpp("EMnodirect.cpp")

l0=2060
l1=1860
l2=1560



dleft1=400
dright1=599
dleft2=800
dright2=1299

readlength=75

data1=read.table("aln.sam", fill = T, col.names = paste("V", 1:19, sep = ""),skip = 2)


#fliter data


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
  h=data.frame(t,a,b,left,right)
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
  
  t=h$t
  a=h$a
  b=h$b
  left=h$left
  right=h$right
  mu=300
  sd=30
  p_0 = 1/3
  p_1 = 1/3
  p_2 = 1/3
  y=rep(0,nrow(h))
  probL=function(dleft,dright,l,a,b){
    for(i in 1:nrow(h))
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
  
  
  
  for( i in 1:60) {
    a1=rightimpossible(dleft1,l1,mu,sd)
    a2=leftimpossible(dright1,l1,mu,sd)
    b1=rightimpossible(dleft2,l2,mu,sd)
    b2=leftimpossible(dright2,l2,mu,sd)
    
    
    probL0=rep(1/(l0-mu),nrow(datamix))
    probL1=probLL_cpp(dleft1,dright1,l1, mu, left, right,a1,a2);
    probL2=probLL_cpp(dleft2,dright2,l2, mu, left, right,b1,b2);
    #probL1=probL(dleft1,dright1,l1,a1,a2)
    #probL2=probL(dleft2,dright2,l2,b1,b2)
    
    
    
    l_0 = p_0 *probL0* dnorm( h$t, mu, sd)
    l_1 = p_1 *probL1* dnorm( h$a, mu, sd)
    l_2 = p_2 *probL2* dnorm( h$b, mu, sd)
    
    
    w_0 =l_0 / (l_0 +l_1+l_2)
    w_1 =l_1 / (l_0 +l_1+l_2)
    w_2 =l_2 / (l_0 +l_1+l_2)
    
    
    
    p_0=mean(w_0)
    p_1=mean(w_1)
    p_2=mean(w_2)
    
    
    
    mu=(sum(w_0*h$t)+sum(w_1*h$a)+sum(w_2*h$b))/(sum(w_0)+sum(w_1)+sum(w_2))
    sd=sqrt((sum(w_0*(h$t-mu)^2)+sum(w_1*(h$a-mu)^2)+sum(w_2*(h$b-mu)^2))/(sum(w_0)+sum(w_1)+sum(w_2)))
    #print(i)
    #print(mu)
    #print(sd)
    #print(c(p_0,p_1,p_2))
  }
  original0=p_0*nrow(h)
  original1=p_1*nrow(h)/(1-(a1+a2))
  original2=p_2*nrow(h)/(1-(b1+b2))
  
  
  
  p0=original0/(original0+original1+original2)
  p1=original1/(original0+original1+original2)
  p2=original2/(original0+original1+original2)
  
  p_00=(p0/l0)/(p0/l0+p1/l1+p2/l2)
  p_11=(p1/l1)/(p0/l0+p1/l1+p2/l2)
  p_22=(p2/l2)/(p0/l0+p1/l1+p2/l2)
  
  result=cbind(p_00,p_11,p_22,mu,sd)
  return(result)
}

rsq=function(data, indices){
  d=data[indices,]
  return(prob(d))
}

h.boot=boot(data1,rsq,R=100)
bias=colMeans(h.boot$t)-h.boot$t0[1,]
citmean=matrix(NA,5,2)
citnorm=matrix(NA,5,2)
citbasic=matrix(NA,5,2)
citperc=matrix(NA,5,2)
sd=array(NA,5)
for( i in 1:5){
  sd[i]=sd(h.boot$t[,i]) 
  citmean[i,]=t.test(h.boot$t[,i])$conf.int
  citnorm[i,]=boot.ci(h.boot,index=i,type="norm")$norm[1,2:3]
  citbasic[i,]=boot.ci(h.boot,index=i,type="basic")$basic[1,4:5]
  citperc[i,]=boot.ci(h.boot,index=i,type="perc")$perc[1,4:5]
}
mean=colMeans(h.boot$t)
a=data.frame(h.boot$t0[1,],mean,bias,sd,citmean,citnorm,citbasic,citperc)
colnames(a)=c("original","mean","bias","sd","meanlow","meanhigh","nlow","nhigh","blow","bhigh","plow","phigh")
colnames(h.boot$t)=c("p0","p1","p2","mu","sd")
plot(h.boot)
write.csv(a, file = "nodirectbt.csv")
write.csv(h.boot$t, file = "nobt_t.csv")




