library(boot)
library(Rcpp)
sourceCpp("EMdirre.cpp")
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

threepop=read.table("fill4all",header = TRUE)
prob=function(threepop){
  mu=500
  sd=50
  p_0 = 1/5
  p_1 = 1/5
  p_2 = 1/5
  p_3 = 1/5
  p_4 = 1/5
  L_V4_1=threepop$L_V4_1
  L_V6m_1=threepop$L_V6m_1
  L_V6_1=threepop$L_V6_1
  L_V4_2 =threepop$L_V4_2 
  L_V6m_2=threepop$L_V6m_2
  L_V6_2=threepop$L_V6_2 
  L_V2_s=threepop$L_V2_s
  L_V6t=threepop$L_V6t
  
  R_V4_1=threepop$R_V4_1
  R_V6m_1=threepop$R_V6m_1
  R_V6_1=threepop$R_V6_1
  R_V4_2 =threepop$R_V4_2 
  R_V6m_2=threepop$R_V6m_2
  R_V6_2=threepop$R_V6_2 
  R_V2_s=threepop$R_V2_s
  R_V6t=threepop$R_V6t
  V9=threepop$V9
  X0=R_V4_1-L_V4_1+R_V6t #even not accurate, but have insure mean not form POP1
  #X1=lengthx(dright1,dleft1)
  
  #X2=lengthx(dright2,dleft2)
  
  #X3=lengthx(dright3,dleft3)
  
  #X4=lengthx(dright4,dleft4)
  
  
  #X5=lengthx(dright5,dleft5)
  #X6=lengthx(dright6,dleft6)
  X1=lengthx_cpp1(dleft1,dright1,L_V4_1,L_V6m_1,L_V4_2,
                  L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                  L_V6t,R_V6t);
  X2=lengthx_cpp1(dleft2,dright2,L_V4_1,L_V6m_1,L_V4_2,
                  L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                  L_V6t,R_V6t);
  X3=lengthx_cpp1(dleft3,dright3,L_V4_1,L_V6m_1,L_V4_2,
                  L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                  L_V6t,R_V6t);
  X4=lengthx_cpp1(dleft4,dright4,L_V4_1,L_V6m_1,L_V4_2,
                  L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                  L_V6t,R_V6t);
  
  for( i in 1:80) {
    
    probL0=probL0=rep(1/(l0-mu),nrow(threepop))
    #probL1=probLL(dleft1,dright1,l1,mu)
    #probL2=probLL(dleft2,dright2,l2,mu)
    #probL3=probLL(dleft3,dright3,l3,mu)
    #probL4=probLL(dleft4,dright4,l4,mu)
    #probL5=probLL(dleft5,dright5,l5,mu)
    #probL6=probLL(dleft6,dright6,l6,mu)
    probL1=probLL_cpp(dleft1,dright1,l1,mu,L_V4_1,L_V6m_1,L_V4_2,
                      L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                      L_V6t,R_V6t);
    
    probL2=probLL_cpp(dleft2,dright2,l2,mu,L_V4_1,L_V6m_1,L_V4_2,
                      L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                      L_V6t,R_V6t);
    probL3=probLL_cpp(dleft3,dright3,l3,mu,L_V4_1,L_V6m_1,L_V4_2,
                      L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                      L_V6t,R_V6t);
    
    probL4=probLL_cpp(dleft4,dright4,l4,mu,L_V4_1,L_V6m_1,L_V4_2,
                      L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                      L_V6t,R_V6t);
    
    
    
    #probL0=probLL0(probL1,probL2,probL3,probL4,probL5,l0,mu)
    l_0 = p_0 *probL0* dnorm( X0, mu, sd)
    l_1 = p_1 *probL1* dnorm( X1, mu, sd)
    l_2 = p_2 *probL2* dnorm( X2, mu, sd)
    l_3 = p_3 *probL3* dnorm( X3, mu, sd)
    l_4 = p_4 *probL4* dnorm( X4, mu, sd)
    
    
    #4
    w_0 =l_0 / (l_0 +l_1+l_2+l_3 +l_4)
    w_1 =l_1 / (l_0 +l_1+l_2+l_3 +l_4)
    w_2 =l_2 / (l_0 +l_1+l_2+l_3 +l_4)
    w_3 =l_3 / (l_0 +l_1+l_2+l_3 +l_4)
    w_4 =l_4 / (l_0 +l_1+l_2+l_3 +l_4)
    
    p_0=mean(w_0)
    p_1=mean(w_1)
    p_2=mean(w_2)
    p_3=mean(w_3)
    p_4=mean(w_4)
    
    #p_00=median(w_0)
    #p_11=median(w_1)
    #p_22=median(w_2)
    #p_33=median(w_3)
    #p_44=median(w_4)
    
    
    mu=(sum(w_0*X0)+sum(w_1*X1)+sum(w_2*X2)+sum(w_3*X3)+sum(w_4*X4))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4))
    sd=sqrt((sum(w_0*(X0-mu)^2)+sum(w_1*(X1-mu)^2)+sum(w_2*(X2-mu)^2)+sum(w_3*(X3-mu)^2)+sum(w_4*(X4-mu)^2))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4)))
    
    #print(i)
    #print(mu)
    #print(sd)
    #print(c(p_0,p_1,p_2,p_3,p_4))
    
  }
  p_00=(p_0/l0)/(p_0/l0+p_1/l1+p_2/l2+p_3/l3+p_4/l4)
  p_11=(p_1/l1)/(p_0/l0+p_1/l1+p_2/l2+p_3/l3+p_4/l4)
  p_22=(p_2/l2)/(p_0/l0+p_1/l1+p_2/l2+p_3/l3+p_4/l4)
  p_33=(p_3/l3)/(p_0/l0+p_1/l1+p_2/l2+p_3/l3+p_4/l4)
  p_44=(p_4/l4)/(p_0/l0+p_1/l1+p_2/l2+p_3/l3+p_4/l4)
  result=cbind(p_00,p_11,p_22,p_33,p_44,mu,sd)
  return(result)
}

rsq=function(data, indices){
  d=data[indices,]
  return(prob(d))
}

h.boot=boot(threepop,rsq,R=100)

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
a=data.frame(h.boot$t0[1,],mean,bias,sd,citmean,citnorm)
a=data.frame(h.boot$t0[1,],mean,bias,sd,citmean,citnorm,citbasic,citperc)
colnames(a)=c("original","mean","bias","sd","meanlow","meanhigh","nlow","nhigh","blow","bhigh","plow","phigh")
plot(h.boot)
write.csv(a, file = "directbt.csv")
write.csv(h.boot$t, file = "directbt_t.csv")
  
  