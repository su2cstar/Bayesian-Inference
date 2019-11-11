setwd("\\Users\\yonsei\\Desktop")
install.packages("MCMCpack");install.packages("statmod")
library(MCMCpack)
library(statmod)
library(bmsaplm)


rm(list=ls())
set.seed(123456)
data(betacaro)
RESULT.indicator   <-list(NA)
RESULT.quantile    <-list(NA)
RESULT.sigma       <-list(NA)
function.list      <-list(NA)
cord.y             <-list(NA)
DATA               <-list(NA)
basis              <-list(NA)
lambda.res         <-list(NA)

knot=function(x,knot){
  if(is.nan(abs((x-knot))^2*log(abs((x-knot))^2))==T){return(0)
  }else{
    return(abs((x-knot))^2*log(abs((x-knot))^2))}
  }

knot_radial            =function(x,knot){abs(x-knot)^3}
standardizing_function =function(x){(x-mean(x))/sd(x)}

#=======================================================#
#====================hyperparameter=====================#
#=======================================================#

k         =2
eta       =1/2
r         =1
delta     =10
v         =2
#=======================================================#
#=====================setting===========================#
#=======================================================#
iter          =25000     # iteration
burn          =25000      # burn-in period
iter          =iter+burn
knot_num      =15        # of knots
#betacaro      =betacaro[-62,]
n             =nrow(betacaro)
non.lin       =8
lin           =2
totl.x        =non.lin+lin
tau.res            <-vector("list",totl.x)
Beta.list          <-vector("list",totl.x)
#x setting -------------------------------------------------
X     <-list()
M.DATA<-list()
#-----------------------------------------------------------
X[[1]]<-betacaro$age
X[[2]]<-betacaro$bmi
X[[3]]<-betacaro$calories
X[[4]]<-betacaro$fat
X[[5]]<-betacaro$fiber
X[[6]]<-betacaro$alcohol
X[[7]]<-betacaro$chol
X[[8]]<-betacaro$betadiet
#X[[9]]<-betacaro$retdiet
#-------------------------
X[[9]]<-betacaro$sex
X[[10]]<-betacaro$smokestat
#X[[12]]<-betacaro$vituse

for(k in 1:non.lin){
#kn=seq(min(X[[k]]),max(X[[k]]),length.out = 11)[2:10]
kn=(max(X[[k]])-min(X[[k]]))*ppoints(knot_num)+min(X[[k]])
M.DATA[[k]]<-matrix(X[[k]],ncol=knot_num +1,nrow=n)
for(i in 1:knot_num )for(j in 1:n) M.DATA[[k]][j,i+1]<-knot(M.DATA[[k]][j,i+1],kn[i])
basis[[k]]=apply(M.DATA[[k]],2,standardizing_function)
}
basis[[9]]<-apply(cbind(X[[9]],X[[10]]),2,standardizing_function)

p<-knot_num+1  #non.linear
q<-2    #linear


#=====================================================
#=====================================================
#simulation===========================================
#=====================================================

#  for(i.n in index.num[1]:index.num[length(index.num)]){

#=====================================================
#save(gibbs)
#=====================================================


#lfdr.list  <-lapply(1:non.lin, function(x) matrix(0,nrow=iter,ncol=p))
#for(k in 1:non.lin) lfdr.list[[k]]<-matrix(0,nrow=iter,ncol=p)
#for(k in (totl.x-lin+1):totl.x) lfdr.list[[k]]<-matrix(0,nrow=iter,ncol=1)
#Beta.list  <-lapply(1:non.lin, function(x) matrix(0,nrow=iter,ncol=p))

for(k in 1:non.lin) Beta.list[[k]]<-matrix(0,nrow=iter,ncol=p)
for(k in (totl.x-lin+1):totl.x) Beta.list[[k]]<-matrix(0,nrow=iter,ncol=1)

for(k in 1:non.lin) tau.res[[k]]<-matrix(0.2^2,nrow=1,ncol=p)
for(k in (totl.x-lin+1):totl.x) tau.res[[k]]<-matrix(0.2^2,nrow=1,ncol=1)

lambda.res <-lapply(1:(non.lin+lin), function(x) matrix(0.2^2,nrow=1,ncol=1))

beta_set.1 <-c(rep(0,p))
beta_set.2 <-c(rep(0,p))
beta_set.3 <-c(rep(0,p))
beta_set.4 <-c(rep(0,p))
beta_set.5 <-c(rep(0,p))
beta_set.6 <-c(rep(0,p))
beta_set.7 <-c(rep(0,p))
beta_set.8 <-c(rep(0,p))
#beta_set.9 <-c(rep(0,p))
beta_set.10 <-0
beta_set.9 <-0
#beta_set.12 <-0

pi0         <-lapply(1:totl.x,function(x) matrix(0,nrow=1,ncol=1))

#v.matrix       <-diag(0.2^2,q)
sg             <-matrix(0.2^2,nrow=iter,ncol=1)             
mu             <-matrix(0,nrow=iter,ncol=1)                 
#linear.beta.mat<-matrix(0,nrow=iter,ncol=q)

# standardization=============================================#
#=============================================================#

ystar=log(betacaro$betacaro+1)
y    =ystar-mean(ystar)

par(mfrow=c(3,3),mar=c(4,2,1,1))
plot(X[[1]],y)   # centered data plot


#index
#================================================================#
#functions using in Gibbs========================================#
#================================================================#
local_fdr_ft=function(iter,data_set,beta_set,pi0,tau){
  Nu=sapply(1:p,function(index) t(y-data_set[,-index]%*%beta_set[-index])%*%data_set[,index])
  De=sapply(1:p,function(index) t(data_set[,index])%*%data_set[,index]+(1/tau[index]))
  local_fdr=sapply(1:p,function(index) pi0/(pi0+(1-pi0)*
                                              sqrt(1/(1+tau[index]*(t(data_set[,index])%*%data_set[,index])))*
                                              exp(Nu[index]^2/(2*De[index]*sg[iter-1,]))))
  return(local_fdr)
}
local_fdr_ft_2=function(iter,data_set,beta_set,pi0,tau){
  Nu=sapply(1:1,function(k) t(as.matrix(y))%*%as.matrix(data_set))
  De=sapply(1:1,function(k) t(as.matrix(data_set))%*%as.matrix(data_set)+(1/tau[k]))
  local_fdr=sapply(1:1,function(k) pi0/(pi0+(1-pi0)*
                                          sqrt(1/(1+tau[k]*(t(data_set)%*%data_set)))*
                                          exp(Nu[k]^2/(2*De[k]*sg[iter-1,]))))
  return(local_fdr)
}

beta_selection_ft=function(iter,index,data_set,beta_set,tau,z){
  Nu<-t(y-data_set[,-index]%*%beta_set[-index])%*%data_set[,index]
  De<-t(data_set[,index])%*%data_set[,index]+(1/tau[index])
  
  rn=sample(c(0,1),size=1,prob=c(z[index],1-z[index]))
  if(rn==0){
    beta_set[index]<-0}
  else{
    beta_set[index]<-rnorm(n=1,mean=Nu/De,sd=sqrt(sg[iter-1,1]/De))
  }
  return(beta_set[index])
} 
beta_selection_ft_2=function(iter,data_set,beta_set,tau,z){
  Nu<-t(y)%*%data_set
  De<-t(data_set)%*%data_set+(1/tau)
  
  rn=sample(c(0,1),size=1,prob=c(z,1-z))
  if(rn==0){
    beta_set<-0}
  else{
    beta_set<-rnorm(n=1,mean=Nu/De,sd=sqrt(sg[iter-1,1]/De))
  }
  return(beta_set)
} 

#================================================================#
# gibbs sampling=================================================#
#i is iteration                                                  
#j is the number of parameters                                   
#================================================================#
for(i in 2:iter){
  #non.linear
  local_fdr_1=local_fdr_ft(i,basis[[1]],beta_set.1,pi0[[1]],tau.res[[1]])
  local_fdr_2=local_fdr_ft(i,basis[[2]],beta_set.2,pi0[[2]],tau.res[[2]])
  local_fdr_3=local_fdr_ft(i,basis[[3]],beta_set.3,pi0[[3]],tau.res[[3]])
  local_fdr_4=local_fdr_ft(i,basis[[4]],beta_set.4,pi0[[4]],tau.res[[4]])
  local_fdr_5=local_fdr_ft(i,basis[[5]],beta_set.5,pi0[[5]],tau.res[[5]])
  local_fdr_6=local_fdr_ft(i,basis[[6]],beta_set.6,pi0[[6]],tau.res[[6]])
  local_fdr_7=local_fdr_ft(i,basis[[7]],beta_set.7,pi0[[7]],tau.res[[7]])
  local_fdr_8=local_fdr_ft(i,basis[[8]],beta_set.8,pi0[[8]],tau.res[[8]])
  #local_fdr_9=local_fdr_ft(i,basis[[9]],beta_set.9,pi0[[9]],tau.res[[9]])
  #---------------------------------------------------------------------------------  
  #linear
  local_fdr_9=local_fdr_ft_2(i,basis[[9]][,1],beta_set.9,pi0[[9]],tau.res[[9]])
  local_fdr_10=local_fdr_ft_2(i,basis[[9]][,2],beta_set.10,pi0[[10]],tau.res[[10]])
  #local_fdr_12=local_fdr_ft_2(i,basis[[10]][,3],beta_set.12,pi0[[12]],tau.res[[12]])
  #---------------------------------------------------------------------------------  
  
  for(j in 1:p){beta_set.1[j]<-beta_selection_ft(i,j,basis[[1]],beta_set.1,tau.res[[1]],local_fdr_1)}
  for(j in 1:p){beta_set.2[j]<-beta_selection_ft(i,j,basis[[2]],beta_set.2,tau.res[[2]],local_fdr_2)}
  for(j in 1:p){beta_set.3[j]<-beta_selection_ft(i,j,basis[[3]],beta_set.3,tau.res[[3]],local_fdr_3)}
  for(j in 1:p){beta_set.4[j]<-beta_selection_ft(i,j,basis[[4]],beta_set.4,tau.res[[4]],local_fdr_4)}
  for(j in 1:p){beta_set.5[j]<-beta_selection_ft(i,j,basis[[5]],beta_set.5,tau.res[[5]],local_fdr_5)}
  for(j in 1:p){beta_set.6[j]<-beta_selection_ft(i,j,basis[[6]],beta_set.6,tau.res[[6]],local_fdr_6)}
  for(j in 1:p){beta_set.7[j]<-beta_selection_ft(i,j,basis[[7]],beta_set.7,tau.res[[7]],local_fdr_7)}
  for(j in 1:p){beta_set.8[j]<-beta_selection_ft(i,j,basis[[8]],beta_set.8,tau.res[[8]],local_fdr_8)}
  #for(j in 1:p){beta_set.9[j]<-beta_selection_ft(i,j,basis[[9]],beta_set.9,tau.res[[9]],local_fdr_9)}
  #---------------------------------------------------------------------------------  
  for(j in 1:1){beta_set.9[j]<-beta_selection_ft_2(i,basis[[9]][,1],beta_set.9,tau.res[[9]],local_fdr_9)}
  for(j in 1:1){beta_set.10[j]<-beta_selection_ft_2(i,basis[[9]][,2],beta_set.10,tau.res[[10]],local_fdr_10)}
  #for(j in 1:1){beta_set.12[j]<-beta_selection_ft_2(i,basis[[10]][,3],beta_set.12,tau.res[[12]],local_fdr_12)}
  
  #---------------------------------------------------------------------------------  
  Beta.list[[1]][i,]<-beta_set.1
  Beta.list[[2]][i,]<-beta_set.2
  Beta.list[[3]][i,]<-beta_set.3
  Beta.list[[4]][i,]<-beta_set.4
  Beta.list[[5]][i,]<-beta_set.5
  Beta.list[[6]][i,]<-beta_set.6
  Beta.list[[7]][i,]<-beta_set.7
  Beta.list[[8]][i,]<-beta_set.8
 # Beta.list[[9]][i,]<-beta_set.9
  #---------------------------------------------------------------------------------  
  Beta.list[[9]][i,]<-beta_set.9
  Beta.list[[10]][i,]<-beta_set.10
 # Beta.list[[12]][i,]<-beta_set.12
  
  linear.beta=rbind(Beta.list[[9]][i,],Beta.list[[10]][i,])
  
  pi0=sapply(1:totl.x,function(k) rbeta(n=1,shape1=k*eta+length(which(Beta.list[[k]][i,]==0)),shape2=k*(1-eta)+p-length(which(Beta.list[[k]][i,]==0))),simplify = F)
  
  Dtau_1<-diag(c(tau.res[[1]]))
  Dtau_2<-diag(c(tau.res[[2]]))
  Dtau_3<-diag(c(tau.res[[3]]))
  Dtau_4<-diag(c(tau.res[[4]]))
  Dtau_5<-diag(c(tau.res[[5]]))
  Dtau_6<-diag(c(tau.res[[6]]))
  Dtau_7<-diag(c(tau.res[[7]]))
  Dtau_8<-diag(c(tau.res[[8]]))
#  Dtau_9<-diag(c(tau.res[[9]]))
  v.matrix<-diag(c(tau.res[[9]],tau.res[[10]]))
  
  
  #------------------------------------------------
  a1=length(which(Beta.list[[1]][i,]==0))
  a2=length(which(Beta.list[[2]][i,]==0))
  a3=length(which(Beta.list[[3]][i,]==0))
  a4=length(which(Beta.list[[4]][i,]==0))
  a5=length(which(Beta.list[[5]][i,]==0))
  a6=length(which(Beta.list[[6]][i,]==0))
  a7=length(which(Beta.list[[7]][i,]==0))
  a8=length(which(Beta.list[[8]][i,]==0))
  a9=length(which(Beta.list[[9]][i,]==0))
  a10=length(which(Beta.list[[10]][i,]==0))
  #a11=length(which(Beta.list[[11]][i,]==0))
 # a12=length(which(Beta.list[[12]][i,]==0))
  red<-(y-basis[[1]]%*%Beta.list[[1]][i,]
         -basis[[2]]%*%Beta.list[[2]][i,]
         -basis[[3]]%*%Beta.list[[3]][i,]
         -basis[[4]]%*%Beta.list[[4]][i,]
         -basis[[5]]%*%Beta.list[[5]][i,]
         -basis[[6]]%*%Beta.list[[6]][i,]
         -basis[[7]]%*%Beta.list[[7]][i,]
         -basis[[8]]%*%Beta.list[[8]][i,]
         #-basis[[9]]%*%Beta.list[[9]][i,]
         -basis[[9]]%*%linear.beta)
  
  #------------------------------------------------
  sg[i]<-rinvgamma(n=1,shape=(n+non.lin*p+q-1-sum(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10))/2,
                   scale=0.5*(t(red)%*%(red)
                        +(t(Beta.list[[1]][i,])%*%solve(Dtau_1)%*%Beta.list[[1]][i,]+
                          t(Beta.list[[2]][i,])%*%solve(Dtau_2)%*%Beta.list[[2]][i,]+
                          t(Beta.list[[3]][i,])%*%solve(Dtau_3)%*%Beta.list[[3]][i,]+
                          t(Beta.list[[4]][i,])%*%solve(Dtau_4)%*%Beta.list[[4]][i,]+
                          t(Beta.list[[5]][i,])%*%solve(Dtau_5)%*%Beta.list[[5]][i,]+
                          t(Beta.list[[6]][i,])%*%solve(Dtau_6)%*%Beta.list[[6]][i,]+
                          t(Beta.list[[7]][i,])%*%solve(Dtau_7)%*%Beta.list[[7]][i,]+
                          t(Beta.list[[8]][i,])%*%solve(Dtau_8)%*%Beta.list[[8]][i,]+
                         # t(Beta.list[[9]][i,])%*%solve(Dtau_9)%*%Beta.list[[9]][i,]+
                          t(linear.beta)%*%solve(v.matrix)%*%linear.beta)
  )
  )
  
  mu[i]<-rnorm(n=1,mean=mean(ystar),sd=sqrt(sg[i]/n))
  
  tau.res=lapply(1:non.lin,function(k)(Beta.list[[k]][i,]==0)*rexp(n=p,rate=lambda.res[[k]][1]/2)+
                   (Beta.list[[k]][i,]!=0)*1/rinvgauss(n=p,mean=sqrt((lambda.res[[k]][1]*sg[i,])/Beta.list[[k]][i,]^2),shape=lambda.res[[k]][1]))
  
  tau.res[[9]]<-(Beta.list[[9]][i,]==0)*rexp(n=1,rate=lambda.res[[9]][1]/2)+
    (Beta.list[[9]][i,]!=0)*1/rinvgauss(n=1,mean=sqrt((lambda.res[[9]][1]*sg[i,])/Beta.list[[9]][i,]^2),shape=lambda.res[[9]][1])
  tau.res[[10]]<-(Beta.list[[10]][i,]==0)*rexp(n=1,rate=lambda.res[[10]][1]/2)+
    (Beta.list[[10]][i,]!=0)*1/rinvgauss(n=1,mean=sqrt((lambda.res[[10]][1]*sg[i,])/Beta.list[[10]][i,]^2),shape=lambda.res[[10]][1])
  #tau.res[[12]]<-(Beta.list[[12]][i,]==0)*rexp(n=1,rate=lambda.res[[12]][1]/2)+
  #  (Beta.list[[12]][i,]!=0)*1/rinvgauss(n=1,mean=sqrt((lambda.res[[12]][1]*sg[i,])/Beta.list[[12]][i,]^2),shape=lambda.res[[12]][1])
  
  lambda.res      =lapply(1:non.lin,function(k) rgamma(n=1,shape=p+r,rate=sum(tau.res[[k]])/2+delta))
  lambda.res[[10]]=rgamma(n=1,shape=1+r,rate=sum(tau.res[[10]])/2+delta)
  lambda.res[[9]]=rgamma(n=1,shape=1+r,rate=sum(tau.res[[9]])/2+delta)
  #lambda.res[[12]]=rgamma(n=1,shape=1+r,rate=sum(tau.res[[12]])/2+delta)
  if(i%%1==0){cat("Iteration ", i, "\r",sep="")}
}

#================================================================#


data_t<-t<-seq(min(X[[1]]),max(X[[1]]),length.out=100)
kn=(max(X[[1]])-min(X[[1]]))*ppoints(knot_num)+min(X[[1]])
data_t=matrix(rep(data_t,knot_num+1),nrow=length(data_t),ncol=knot_num+1)
for(i in 1:knot_num)for(j in 1:100) data_t[j,i+1]=knot(data_t[j,i+1],kn[i])
data_t[data_t=="NaN"]=0
data_t=apply(data_t,2,standardizing_function)

#variable selection using FDR=========================================#

# R                 =list(NA)
# K_star            =list(NA)
# 
# for(j in 1:non.lin){R[[j]]=lapply(1:iter,function(i)cumsum(lfdr.list[[j]][i,c(order(lfdr.list[[j]][i,]))])/seq_along(lfdr.list[[j]][i,c(order(lfdr.list[[j]][i,]))]))}
# for(j in 1:non.lin){K_star[[j]]=lapply(1:iter,function(i)order(lfdr.list[[j]][i,])[which(R[[j]][[i]]<alpha)])}
#-----------------------------------------------------------------------------
alpha=0.35
order.signif       =order(unlist(indiplot),decreasing = T)
j.star             =if(length(which(cumsum(1-unlist(indiplot)[order.signif])/seq(1,non.lin*p+2,1)<alpha))==0){
  cat("need to reset alpha")}else{max(which(cumsum(1-unlist(indiplot)[order.signif])/seq(1,p*non.lin+2,1)<alpha))}

phi.alpha          =unlist(indiplot)[order.signif[j.star]] 
K.final            =which(unlist(indiplot)>phi.alpha)

#------------------------------------------------------------------- 
index =lapply(1:non.lin,function(j) seq((j-1)*p+1,j*p,1))
TT    =vector("list",8) 
RES.1<-c()
for(u in 1:length(K.final)){  
  if(K.final[u]<p+1){RES.1[u]<-1}
  else if(p<K.final[u]&K.final[u]<2*p+1){RES.1[u]<-2}
  else{RES.1[u]<-3}
}
RES=data.frame(RES.1,K.final-(RES.1-1)*p,K.final)
colnames(RES)<-c("a","b","c")

for(u in 1:non.lin){
  TT[[u]]=RES$b[which(RES$a==u)]
}
TT[[3]]<-TT[[3]][-5]
#-----------------------------------------------------------------------------
cord.x  =c(t,rev(t))
indiplot=vector("list",non.lin)
# plot of estimate of local fdr 
indicator.list=lapply(1:non.lin,function(x) matrix(0, nrow=iter, ncol=p))
for(j in 1:non.lin){
  for(i in 1:iter){
    indicator.list[[j]][i,]<-(Beta.list[[j]][i,]!=0)*1
  }
}
for(j in non.lin) indicator.list[[j]]<-indicator.list[[j]][-c(1:burn),]
for(k in 1:non.lin) indiplot[[k]]=apply(indicator.list[[k]],2,mean)
indiplot[[9]]<-mean(((Beta.list[[9]]!=0)*1)[-c(1:burn),])
indiplot[[10]]<-mean(((Beta.list[[10]]!=0)*1)[-c(1:burn),])
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(unlist(indiplot),lwd=2.5,type="h",xlab="",ylab="")
mtext(side=1,"j",2)
mtext(side=2,expression(p[j]),2)

for(kk in 1:non.lin) abline(v=16*kk+0.5,lty=2,col=1,lwd=1)
abline(v=16*8+1.5,lty=2,col=1,lwd=1)

abline(h=phi.alpha,lty=1,col=1)
hist(Beta.list[[10]][-c(1:burn)]/sd(X[[10]]))
hist(Beta.list[[9]][-c(1:burn)]/sd(X[[9]]))
hist(Beta.list[[12]][-c(1:burn)]/sd(X[[12]]))
ts.plot(sg)
#------------------------------------------------------------------------------------
# result without FDR based selection    
output.list.1=lapply(1:non.lin, function(x) matrix(NA,nrow=length(t),ncol=iter))
for(k in 1:non.lin){
  for(i in 1:length(t)){
    output.list.1[[k]][i,]=sapply(1:iter,function(j) (data_t[i,])%*%Beta.list[[k]][j,])
  } 
}

Quantile.list.1<-lapply(1:non.lin,function(x) matrix(NA,nrow=length(t),ncol=4))
for(j in 1:non.lin){
  for(i in 1:length(t)){
    Quantile.list.1[[j]][i,]<-c(quantile(output.list.1[[j]][i,-c(1:burn)],prob=c(0.05/2,0.5,1-(0.05/2))) ,mean(output.list.1[[j]][i,-c(1:burn)]))
  }
}
cord.y.1  =lapply(1:non.lin,function(i) c(Quantile.list.1[[i]][,1],rev(Quantile.list.1[[i]][,3])))
#--------------------------------------------------------
# result with FDR based selection 
#--------------------------------------------------------

output.list.2=lapply(1:non.lin,function(x) matrix(NA,nrow=length(t),ncol=iter))
for(k in 1:non.lin){
  for(i in 1:length(t)){
    output.list.2[[k]][i,]=sapply(1:iter,function(j) (data_t[i,TT[[k]]])%*%Beta.list[[k]][j,TT[[k]]])
  } 
}

Quantile.list.2=lapply(1:non.lin,function(x) matrix(NA,nrow=length(t),ncol=4))
for(j in 1:non.lin){
  for(i in 1:length(t)){
    Quantile.list.2[[j]][i,]<-c(quantile(output.list.2[[j]][i,-c(1:burn)],prob=c(0.05/2,0.5,1-(0.05/2))) ,mean(output.list.2[[j]][i,-c(1:burn)]))
  }
}

cord.y.2  =lapply(1:non.lin,function(i) c(Quantile.list.2[[i]][,1],rev(Quantile.list.2[[i]][,3])))




# RESULT.indicator[[7]]=apply(indicator.list[[1]][-(1:burn),],2,mean)
# RESULT.quantile [[7]]=Quantile_result
# RESULT.sigma    [[7]]=sqrt(sg[-(1:burn)])
# cord.y          [[7]]=c(RESULT.quantile[[7]][,1],rev(RESULT.quantile[[7]][,3]))

#  }

# MSE[seed]=sum((true_mu+function.list[[i.n]](t)-RESULT.quantile[[i.n]][,4])^2)

# cat(i.n,"th Simulation : ","# of data = ",n,",",seed," iteration","\n")
#}


# save data
save.image("simulation_real.RData")

ts.plot(sg[-c(1:burn)])
median(sg[-c(1:burn)])
#=============================================================================
#making plot==================================================================
#=============================================================================
postscript("model1_Z_250.eps",horizontal=F,onefile=F,print.it=F,height=6,width=8)

par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfcol=c(3,3))
nn=1
plot(c(seq(start_x,end_x,(end_x-start_x)/knot_num))[-1],indiplot[[nn]][-1],type="h",lwd=2.5,xlab="",xlim=range(start_x,end_x),ylim=c(0,1),ylab="")
points(start_x,indiplot[[nn]][1],type="h",lty=3)
points(start_x,indiplot[[nn]][1],pch=0) 
abline(h=phi.alpha,lty=2,col=2)
mtext("x", side=1, line=2) 
mtext("Probability", side=2, line=2.5,cex=0.9)

#---------------------------------------------------------------------
plot(x=c(min(X[[1]]),max(X[[1]])),y=c(min(cord.y.1[[1]]-0.5),max(cord.y.1[[1]])+0.5),type="n",xlab="",ylab="", main="")
cord.x<-c(t,rev(t))
polygon(cord.x,cord.y.1[[1]],col=gray(0:9/9)[8],border=NA)
lines(t,Quantile.list.1[[1]][,4])
mtext("x", side=1, line=2)
mtext(expression(italic(paste(f(x)))), side=2, line=2.5)

plot(x=c(start_x,end_x),y=c(min(cord.y.2[[nn]]-0.5),max(cord.y.2[[nn]])+0.5),type="n",xlab="",ylab="", main="")
mtext("x", side=1, line=2)
mtext(expression(italic(paste(f(x)))), side=2, line=2.5)
polygon(cord.x,cord.y.2[[nn]],col=gray(0:9/9)[8],border=NA)
lines(t,Quantile.list.2[[nn]][,2])
#---------------------------------------------------------------------
lines(t,true_mu+function.list[[1]](t)-mean(function.list[[1]](t)),lty=2)
lines(t,true_mu+function.list[[2]](t)-mean(function.list[[2]](t)),lty=2)
lines(t,true_mu+function.list[[6]](t)-mean(function.list[[6]](t)),lty=2)
#---------------------------------------------------------------------
plot(unlist(indiplot),type="h",lwd=3,ylab="")

#MSE-----------------------------------------------------------------
M11=mean((Quantile.list.1[[1]][,4]-(true_mu+function.list[[1]](t)+mean(a2+a3+b1+b2)))^2)
M21=mean((Quantile.list.2[[1]][,4]-(true_mu+function.list[[1]](t)+mean(a2+a3+b1+b2)))^2)

M12=mean((Quantile.list.1[[2]][,4]-(true_mu+function.list[[2]](t)+mean(a1+a3+b1+b2)))^2)
M22=mean((Quantile.list.2[[2]][,4]-(true_mu+function.list[[2]](t)+mean(a1+a3+b1+b2)))^2)

M13=mean((Quantile.list.1[[3]][,4]-(true_mu+function.list[[6]](t)+mean(a1+a2+b1+b2)))^2)
M23=mean((Quantile.list.2[[3]][,4]-(true_mu+function.list[[6]](t)+mean(a1+a2+b1+b2)))^2)

M11+M12+M13
M21+M22+M23
#---------------------------------------------------------------------
par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfcol=c(2,3))

hist(Beta.list[[4]][-c(1:burn),1]/sd(DATA[[4]]),xlab="",ylab="",main="",freq=F)
mtext(expression(paste(beta[1])), side=1, line=2.5,cex=0.9);abline(v=true.line.beta[1],lty=2,lwd=2,col=2)
mtext("Posterior density",side=2,2.2,cex=0.9)
ts.plot(Beta.list[[4]][c(10001:20000),1]/sd(DATA[[4]]),xlab="",ylab="",main="")

hist(Beta.list[[5]][-c(1:burn),1]/sd(DATA[[5]]),xlab="",ylab="",main="",freq=F)
mtext(expression(paste(beta[2])), side=1, line=2.5,cex=0.9);abline(v=true.line.beta[2],lty=2,lwd=2,col=2)
mtext("Posterior density",side=2,2.2,cex=0.9)
ts.plot(Beta.list[[5]][c(10001:20000),1]/sd(DATA[[5]]),xlab="",ylab="",main="")

hist(sg[-c(1:burn)],xlab="",ylab="",main="",freq=F);abline(v=true_sg^2,lty=2,lwd=2,col=2)
mtext(expression(paste(sigma^2)), side=1, line=2.5,cex =0.9)
mtext("Posterior density",side=2,2.2,cex=0.9)
ts.plot(sg[c(10001:20000)],xlab="",ylab="",main="")


#-------------------------------------------------------------------------------------------
par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfcol=c(1,1)) 
ts.plot(sg[-c(1:10000)])
hist(sg[-c(1:burn)],xlab="",ylab="",main="",freq=F);abline(v=true_sg^2,lty=2,lwd=2,col=2)
#abline(v=quantile(sg[-c(1:burn),1],prob=c(0.025,0.975)))


ts.plot(Beta.list[[5]][-c(1:burn),1]/sd(DATA[[5]]))


hist(sg[-c(1:10000)])
ts.plot(sg[-c(1:10000)])
ts.plot(sg)
abline(h=1,col=2)
mean(ystar)
# save(cord.x,file="cordx_Z_7_1000.RData")
# save(cord.y,file="cordy_Z_7_1000.RData")
# save(RESULT.quantile,file="quantile_Z_7_1000.RData")
# save(MSE,file="MSE_Z_7_1000.RData" )


#MSE-----------------------------------------------------------------------------
hist(MSE)
boxplot(MSE)	
round(mean(MSE),3)
round(sd(MSE),3)		
median(MSE)/250


#check --------------------------------------------------------------------------
# beta_1=beta_mat[-c(1:burn),1]
# 
# par(mfrow=c(1,2))
# hist(beta_1,main="",xlab="")
# mtext(expression(paste("histogram of  " ,beta[1])),side=3,line=1,cex=1.3)
# hist(beta_1/sd_data[1],xlab="",main="")
# mtext(expression(paste("histogram of  " , over(beta[1],sigma))),side=3,line=1,cex=1.3)
# abline(v=2,lwd=2)
# median(beta_1/sd_data[1])
# abline(lty=2,v=quantile(beta_1/sd_data[1],prob=c(alpha/2,1-alpha/2)))


#===========================================================
#parameters_sigma=============================================
#===========================================================
postscript("simplot3_sg_ex.eps",horizontal=F,onefile=F,print.it=F,height=2,width=8)

par(mex=0.1,mar=c(3.3,3.3,.1,.1)+1,mfrow=c(1,length(index.num)),cex=0.9,bty="n")
for(i.n in index.num[1]:index.num[length(index.num)]){
  
  h<-hist(RESULT.sigma[[i.n]],plot=F)
  plot(h$breaks,c(h$density,0),type="s",ylab="",xlab="")
  mtext(expression(sigma),side=1,line=2.1)
  mtext("Posterior pdf",side=2,line=2.2)
  abline(v=true_sg,lwd=2)
  abline(v=quantile(RESULT.sigma[[i.n]],prob=c(0.025,0.975)),lty=2)
}
dev.off()


# checking plot===========================================================

cairo_ps("noname.eps",onefile=F,height=4,width=8, fallback_resolution = 600)
par(mfrow=c(1,2))
plot(c(0,1),c(min(ystar)-0.5,max(ystar)+0.5),type="n",xlab="",ylab="",main=paste("N=",n,"iteration=",iter))
polygon(cord.x,cord.y[[i.n]],col="gray",border=NA)
lines(t,Quantile_result[,2])
lines(t,true_mu+function.list[[index.num]](t),lty=2,col="red")

plot(c(0,1),c(min(ystar)-0.5,max(ystar)+0.5),type="n",xlab="",ylab="",main=paste("N=",n,"iteration=",iter))
points(x,ystar,cex=0.3,col="gray40")
polygon(cord.x,cord.y[[i.n]],col=adjustcolor("#D3D3D3", alpha=0.8),border=NA)
lines(t,Quantile_result[,2])

dev.off()
#======================================================================


#======================================================================
# the number of selected betas================================================
K_starresult<-list(NA)
K_starresult=lapply(K_star[[3]],function(x){length(x)})
postscript("noname.eps",horizontal=F,onefile=F,print.it=F,height=4,width=8)
par(mfrow=c(1,1))
plot(ts(unlist(K_starresult)),ylab="number",main="The number of selected betas")
dev.off()
#=======================================================
summary(unlist(K_starresult)[-(1:burn)])
K_star[[1]][[26199]]

#parameters=============================================
# sigma-------------------------------------------------
par(mfrow=c(1,1))
hist(sqrt(sg[-c(1:burn)]),xlab="iteration",main="sigma")
abline(v=true_sg,lwd=2,col="red")
#=======================================================
# mu===================================================
plot(mu[-(1:10000)],type="l")
abline(h=true_mu+mean(ystar),col="red")
median(mu[-(1:1000)])



