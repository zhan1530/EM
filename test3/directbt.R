library(boot)
library(Rcpp)
sourceCpp("gene_with_direct.cpp")

lengthx=function(dleft,dright,t,I){
  x=rep(0,nrow(datamix2))
  for(i in 1:length(left))
  {
    if(left[i]<=(dleft)&right[i]>=(dright-3))
    {x[i]=t[i]-(dright-dleft+1)}
    
    #can be more zhuque!!!!!!
    else if(((dright-3)<left[i])&(left[i]<=dright+3)&(I[i]=='NO75ML')){
      x[i]=t[i]+(readlength-as.numeric(gsub('\\d+[A-Z](\\d+)M','\\1',V6[i])))
    }
    
    else if(((dleft-readlength+1)<right[i])&(right[i]<=dleft)&(I[i]=='NO75MR')){
      x[i]=right[i]-left[i]+readlength
    }
    
    else if ((dleft)<left[i]&left[i]<(dright-3)&(I[i]=='75M'|I[i]=='NO75ML'))
    {x[i]=0}
    
    else if((dleft)<right[i]&right[i]<(dright-3)&(I[i]=='75M'|I[i]=='NO75MR'))
    {x[i]=0}
    
    else
    {x[i]=t[i]}
  }
  return(x)
}

l0=2060
l1=1860
l2=1560



dleft1=400
dright1=599
dleft2=800
dright2=1299
readlength=75

data=read.table("aln.sam", fill = T, col.names = paste("V", 1:19, sep = ""),skip = 2)
datamix1=subset(data,V2==97|V2==99)

#head(datamix1)
V4=datamix1$V4
V6=datamix1$V6
V8=datamix1$V8
V9=datamix1$V9
nrow(datamix1)
I1=rep(0,nrow(datamix1))
for(i in 1:nrow(datamix1)){
  if((V8[i]-V4[i]+readlength==V9[i])&(V6[i]=="75M")){
    I1[i]="75M"}
  else if((V8[i]-V4[i]+readlength!=V9[i])&(V6[i]=="75M")){
    I1[i]="NO75MR"}
  else {I1[i]="NO75ML"}
}
#table(I1)
#which(I1=="NO75ML")
#which(I1=="NO75MR")
datamix2=cbind(datamix1,I1)
#head(datamix2)
left=datamix2$V4
right=datamix2$V8
t=datamix2$V9
I1=datamix2$I1
a=lengthx(dleft1,dright1,t,I1)
b=lengthx(dleft2,dright2,t,I1)
#c=lengthx(dleft3,dright3,t,I1)
#d=lengthx(dleft4,dright4,t,I1)
h=data.frame(t,a,b,I1,left,right,V6,V9)

prob=function(h){
  mu=300
  sd=30
  p_0 = 1/3
  p_1 = 1/3
  p_2 = 1/3
  
  left=h$left
  right=h$right
  I1=h$I1
  for( i in 1:60) {
    probL0=rep(1/(l0-mu),nrow(h))
    #probL1=probLL(dleft1,dright1,l1,mu,I1)
    #probL2=probLL(dleft2,dright2,l2,mu,I1)
    #probL3=probLL(dleft3,dright3,l3,mu,I1)
    #probL4=probLL(dleft4,dright4,l4,mu,I1)
    probL1=probLL_cpp(dleft1,dright1,l1, mu, I1,left, right);
    probL2=probLL_cpp(dleft2,dright2,l2, mu, I1,left, right);
    
    
    l_0 = p_0 *probL0* dnorm( h$t, mu, sd)
    l_1 = p_1 *probL1* dnorm( h$a, mu, sd)
    l_2 = p_2 *probL2* dnorm( h$b, mu, sd)
    
    
    
    #4
    w_0 =l_0 / (l_0 +l_1+l_2)
    w_1 =l_1 / (l_0 +l_1+l_2)
    w_2 =l_2 / (l_0 +l_1+l_2)
    
    
    p_0=mean(w_0)
    p_1=mean(w_1)
    p_2=mean(w_2)
    
    #p_00=median(w_0)
    #p_11=median(w_1)
    #p_22=median(w_2)
    #p_33=median(w_3)
    #p_44=median(w_4)
    
    
    mu=(sum(w_0*h$t)+sum(w_1*h$a)+sum(w_2*h$b))/(sum(w_0)+sum(w_1)+sum(w_2))
    sd=sqrt((sum(w_0*(h$t-mu)^2)+sum(w_1*(h$a-mu)^2)+sum(w_2*(h$b-mu)^2))/(sum(w_0)+sum(w_1)+sum(w_2)))
    
    #print(i)
    #print(mu)
    #print(sd)
    #print(c(p_0,p_1,p_2))
    #print(c(p_00,p_11,p_22,p_33,p_44))
  }
  p_00=(p_0/l0)/(p_0/l0+p_1/l1+p_2/l2)
  p_11=(p_1/l1)/(p_0/l0+p_1/l1+p_2/l2)
  p_22=(p_2/l2)/(p_0/l0+p_1/l1+p_2/l2)
  result=cbind(p_00,p_11,p_22,mu,sd)
  return(result)
}
rsq=function(data, indices){
  d=data[indices,]
  return(prob(d))
}

h.boot=boot(h,rsq,R=100)
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
#plot(h.boot)
write.csv(a, file = "directbt.csv")
write.csv(h.boot$t, file = "directbt_t.csv")

