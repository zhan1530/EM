library(boot)
library(Rcpp)
sourceCpp("EMdirLOG.cpp")
l0=1763#na
l3=1733#6
#l6=1601#5-7
l1=1535#5-8
l2=875#3-16
l4=1403#4-9
l5=875#1-14
l6=942#6-18


dleft3=499#6
dright3=528
#dleft6=433#5-7
#dright6=594
dleft1=433#5-8
dright1=660
dleft2=301#3-16
dright2=1188
dleft4=367#4-9
dright4=726
dleft5=169#1-14
dright5=1056
dleft6=499#6-18
dright6=1319

threepop1=read.table("fill5all",header = TRUE)
threepop=subset(threepop1,(!is.na(threepop1$L_V4_1))&(!is.na(threepop1$R_V4_1)))

prob=function(threepop){
  
  mu=5
  sd=1
  p_0 = 1/7
  p_1 = 1/7
  p_2 = 1/7
  p_3 = 1/7
  p_4 = 1/7
  p_5 = 1/7
  p_6 = 1/7
  
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
  X5=lengthx_cpp1(dleft5,dright5,L_V4_1,L_V6m_1,L_V4_2,
                  L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                  L_V6t,R_V6t);
  X6=lengthx_cpp1(dleft6,dright6,L_V4_1,L_V6m_1,L_V4_2,
                  L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                  L_V6t,R_V6t);
  #threepopnew=subset(threepop,(is.na(threepop$L_V4_2))&(is.na(threepop$R_V4_2)))
  #threepopnewnew=subset(threepopnew,(L_V6m_1!=L_V6t)&(R_V6m_1!=R_V6t))
  #threepopnewnew
  
  X00=log(X0)
  X11=log(X1)
  X22=log(X2)
  X33=log(X3)
  X44=log(X4)
  X55=log(X5)
  X66=log(X6)
  for( i in 1:100) {
    
    probL0=probL0=rep(1/(l0-exp(mu)),nrow(threepop))
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
    probL5=probLL_cpp(dleft5,dright5,l5,mu,L_V4_1,L_V6m_1,L_V4_2,
                      L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                      L_V6t,R_V6t);
    
    probL6=probLL_cpp(dleft6,dright6,l6,mu,L_V4_1,L_V6m_1,L_V4_2,
                      L_V6m_2,R_V4_1,R_V6m_1,R_V4_2,R_V6m_2,
                      L_V6t,R_V6t);
    #probL0=probLL0(probL1,probL2,probL3,probL4,probL5,l0,mu)
    l_0 = p_0 *probL0* dnorm( X00, mu, sd)
    l_1 = p_1 *probL1* dnorm( X11, mu, sd)
    l_2 = p_2 *probL2* dnorm( X22, mu, sd)
    l_3 = p_3 *probL3* dnorm( X33, mu, sd)
    l_4 = p_4 *probL4* dnorm( X44, mu, sd)
    l_5 = p_5 *probL5* dnorm( X55, mu, sd)
    l_6 = p_6 *probL6* dnorm( X66, mu, sd)
    #4
    #which(is.na(w_0))
    w_0 =l_0 / (l_0 +l_1+l_2+l_3 +l_4+l_5+l_6)
    w_1 =l_1 / (l_0 +l_1+l_2+l_3 +l_4+l_5+l_6)
    w_2 =l_2 / (l_0 +l_1+l_2+l_3 +l_4+l_5+l_6)
    w_3 =l_3 / (l_0 +l_1+l_2+l_3 +l_4+l_5+l_6)
    w_4 =l_4 / (l_0 +l_1+l_2+l_3 +l_4+l_5+l_6)
    w_5 =l_5 / (l_0 +l_1+l_2+l_3 +l_4+l_5+l_6)
    w_6 =l_6 / (l_0 +l_1+l_2+l_3 +l_4+l_5+l_6)
    
    p_0=mean(w_0)
    p_1=mean(w_1)
    p_2=mean(w_2)
    p_3=mean(w_3)
    p_4=mean(w_4)
    p_5=mean(w_5)
    p_6=mean(w_6)
    #p_00=median(w_0)
    #p_11=median(w_1)
    #p_22=median(w_2)
    #p_33=median(w_3)
    #p_44=median(w_4)
    
    
    mu=(sum(w_0*X00)+sum(w_1*X11)+sum(w_2*X22)+sum(w_3*X33)+sum(w_4*X44)+sum(w_5*X55)+sum(w_6*X66))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4)+sum(w_5)+sum(w_6))
    sd=sqrt((sum(w_0*(X00-mu)^2)+sum(w_1*(X11-mu)^2)+sum(w_2*(X22-mu)^2)+sum(w_3*(X33-mu)^2)+sum(w_4*(X44-mu)^2)+sum(w_5*(X55-mu)^2)+sum(w_6*(X66-mu)^2))/(sum(w_0)+sum(w_1)+sum(w_2)+sum(w_3)+sum(w_4)+sum(w_5)+sum(w_6)))
    
    print(i)
    print(mu)
    print(sd)
    print(c(p_0,p_1,p_2,p_3,p_4,p_5,p_6))
    
  }
  p00=(p_0/l0)/((p_0/l0)+(p_1/l1)+(p_2/l2)+(p_3/l3)+(p_4/l4)+(p_5/l5)+(p_6/l6))
  p11=(p_1/l1)/((p_0/l0)+(p_1/l1)+(p_2/l2)+(p_3/l3)+(p_4/l4)+(p_5/l5)+(p_6/l6))
  p22=(p_2/l2)/((p_0/l0)+(p_1/l1)+(p_2/l2)+(p_3/l3)+(p_4/l4)+(p_5/l5)+(p_6/l6))
  p33=(p_3/l3)/((p_0/l0)+(p_1/l1)+(p_2/l2)+(p_3/l3)+(p_4/l4)+(p_5/l5)+(p_6/l6))
  p44=(p_4/l4)/((p_0/l0)+(p_1/l1)+(p_2/l2)+(p_3/l3)+(p_4/l4)+(p_5/l5)+(p_6/l6))
  p55=(p_5/l5)/((p_0/l0)+(p_1/l1)+(p_2/l2)+(p_3/l3)+(p_4/l4)+(p_5/l5)+(p_6/l6))
  p66=(p_6/l6)/((p_0/l0)+(p_1/l1)+(p_2/l2)+(p_3/l3)+(p_4/l4)+(p_5/l5)+(p_6/l6))
  result=cbind(p00,p11,p22,p33,p44,p55,p66,mu,sd)
  return(result)
}

rsq=function(data, indices){
  d=data[indices,]
  return(prob(d))
}

h.boot=boot(threepop,rsq,R=100)

bias=colMeans(h.boot$t)-h.boot$t0[1,]
citmean=matrix(NA,9,2)
citnorm=matrix(NA,9,2)
citbasic=matrix(NA,9,2)
citperc=matrix(NA,9,2)
sd=array(NA,9)
for( i in 1:9){
  sd[i]=sd(h.boot$t[,i]) 
  citmean[i,]=t.test(h.boot$t[,i])$conf.int
  citnorm[i,]=boot.ci(h.boot,index=i,type="norm")$norm[1,2:3]
  citbasic[i,]=boot.ci(h.boot,index=i,type="basic")$basic[1,4:5]
  citperc[i,]=boot.ci(h.boot,index=i,type="perc")$perc[1,4:5]
}
mean=colMeans(h.boot$t)
a=data.frame(h.boot$t0[1,],mean,bias,sd,citmean,citnorm,citbasic,citperc)
colnames(a)=c("original","mean","bias","sd","meanlow","meanhigh","nlow","nhigh","blow","bhigh","plow","phigh")
colnames(h.boot$t)=c("p0","p1","p2","p3","p4","p5","p6","mu","sd")
write.csv(a, file = "btnore.csv")
write.csv(h.boot$t, file = "bt_tnore.csv")
