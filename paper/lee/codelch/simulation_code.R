install.packages("MCMCpack");install.packages("statmod")
setwd("C:\\Users\\LeeChangHwan\\Google 드라이브\\SIMULATION\\1209")
load("simulation_1215.RData")
library(MCMCpack);library(statmod)

rm(list=ls())
seednum=115676634
#346269for(seed in 1:length(seednum)) {1213415
set.seed(seednum)

RESULT.indicator   <-list(NA)
RESULT.quantile    <-list(NA)
RESULT.sigma       <-list(NA)
function.list      <-list(NA)
cord.y             <-list(NA)
DATA               <-list(NA)
basis              <-list(NA)
lambda.res         <-list(NA)

knot                   =function(x,knot){abs((x-knot))^2*log(abs((x-knot))^2)}
knot_radial            =function(x,knot){abs(x-knot)^3}
standardizing_function =function(x){(x-mean(x))/sd(x)}
#--------------------------------------------------------
function.list<-
  c(
     simulation.plot.1=function(x){2*exp(-30*(x-0.2)^2)+exp(-50*(x-0.7)^2)}
    ,simulation.plot.2=function(x){sin(2*pi*x)}
    ,simulation.plot.3=function(x){x}
    ,simulation.plot.4=function(x){-0.75*x}
    ,simulation.plot.5=function(x){rep(0,length(x))}
    ,simulation.plot.6=function(x){rep(0,length(x))}
    ,simulation.plot.7=function(x){rep(0,length(x))}
    ,simulation.plot.8=function(x){rep(0,length(x))}
    )
    

#function setting----------------------------------------
#=======================================================#
#====================hyperparameter=====================#
#=======================================================#
k_pi      =2
eta       =1/2
r         =1
delta     =100

#=======================================================#
#=====================setting===========================#
#=======================================================#
n             =800             # of sample
iter          =90000           # iteration
burn          =10000          # burn-in period
iter          =iter+burn
knot_num      =15             # of knots
true_mu       =0
true_sg       =sqrt(1)       # sigma
true.line.beta=c(0.6,-1,rep(0,4))

non.lin       =8
lin           =6
totl.x        =non.lin+lin

tau.res            <-vector("list",totl.x)
lfdr.list          <-vector("list",totl.x)
Beta.list          <-vector("list",totl.x)
temp_beta_set      <-vector("list",totl.x)
pi0                <-vector("list",totl.x)
#x setting -------------------------------------------------
start_x       =0
end_x         =1
x             =runif(n,start_x,end_x)
kn            =ppoints(knot_num)
DATA          =lapply(1:non.lin,function(j) runif(n,start_x,end_x))

M.DATA        =lapply(1:non.lin,function(j) matrix(rep(DATA[[j]],knot_num+1),nrow=n,ncol=knot_num+1))
for(j in 1:non.lin){for(i in 1:knot_num){M.DATA[[j]][,i+1]<-knot_radial(M.DATA[[j]][,i+1],kn[i])}} # apply basis
for(j in 1:length(M.DATA)){basis[[j]]=apply(M.DATA[[j]],2,standardizing_function)}          # apply standardizing

basis.candid<-matrix(NA,nrow=n,ncol=lin)

#categorical data
for(j in c(1,3,5)){basis.candid[,j]<-sample(c(1,2),n,prob=c(0.5,0.5),replace=T)}
for(j in c(2,4,6)){basis.candid[,j]<-sample(c(1,2,3),n,prob=c(1/3,1/3,1/3),replace=T)}
basis[[length(M.DATA)+1]]<-apply(basis.candid,2,standardizing_function)
p<-knot_num+1                      #non.linear coef number
q<-ncol(basis[[length(M.DATA)+1]]) #linear number

#=====================================================
#simulation===========================================
#=====================================================
#=====================================================
#save(gibbs)
#=====================================================
#non.linear
for(k in 1:non.lin){
   lfdr.list[[k]]    <-matrix(0,nrow=iter,ncol=p)
   Beta.list[[k]]    <-matrix(1,nrow=iter,ncol=p)
   tau.res[[k]]      <-matrix(0.2^2,nrow=1,ncol=p)
   temp_beta_set[[k]]<-c(rep(0,p))
   
   }
#linear
for(k in (totl.x-lin+1):totl.x) {
    lfdr.list[[k]]    <-matrix(0,nrow=iter,ncol=1)
    Beta.list[[k]]    <-matrix(0,nrow=iter,ncol=1)
    tau.res[[k]]      <-matrix(0.2^2,nrow=1,ncol=1)
    temp_beta_set[[k]]<-0
    }

lambda.res     <-0.01
sg             <-matrix(0.2^2,nrow=iter,ncol=1)             
mu             <-matrix(0,nrow=iter,ncol=1)                 
pi0            <-lapply(1:totl.x, function(x) matrix(0,nrow=1,ncol=1))

# standardization=============================================#
#=============================================================#
yy1   = Reduce("+",lapply(1:non.lin,function(j) function.list[[j]](DATA[[j]])))
yy2   = basis.candid[,1]*true.line.beta[1]+true.line.beta[2]*basis.candid[,2]
ystar = yy1+yy2+rnorm(n,0,sd=true_sg)
y     = ystar-mean(ystar)


par(mfrow=c(2,4),mar=c(4,2,1,1))
for(j in 1:8) plot(DATA[[j]],y)   
# centered data plot
par(mfrow=c(2,3),mar=c(4,2,1,1))
for(j in 1:6) plot(basis.candid[,j],y) 


#================================================================#
#functions using in Gibbs========================================#
#================================================================#
local_fdr_ft       =function(iter,data_set,beta_set,pi0,tau){
  Nu=sapply(1:p,function(index) t(y-data_set[,-index]%*%beta_set[-index])%*%data_set[,index])
  De=sapply(1:p,function(index) t(data_set[,index])%*%data_set[,index]+(1/tau[index]))
  local_fdr=sapply(1:p,function(index) pi0/(pi0+(1-pi0)*
                                              sqrt(1/(1+tau[index]*(t(data_set[,index])%*%data_set[,index])))*
                                              exp(Nu[index]^2/(2*De[index]*sg[iter-1,]))))
  return(local_fdr)
}
local_fdr_ft_2     =function(iter,data_set,beta_set,pi0,tau){
  Nu=sapply(1:1,function(k) t(as.matrix(y))%*%as.matrix(data_set))
  De=sapply(1:1,function(k) t(as.matrix(data_set))%*%as.matrix(data_set)+(1/tau[k]))
  local_fdr=sapply(1:1,function(k) pi0/(pi0+(1-pi0)*
                                          sqrt(1/(1+tau[k]*(t(data_set)%*%data_set)))*
                                          exp(Nu[k]^2/(2*De[k]*sg[iter-1,]))))
  return(local_fdr)
}

beta_selection_ft  =function(iter,index,data_set,beta_set,tau,z){
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
  for(j in 1:non.lin){
  lfdr.list[[j]][i,]=local_fdr_ft(i,basis[[j]],temp_beta_set[[j]],pi0[[j]],tau.res[[j]])}
  #linear
  for(j in (totl.x-lin+1):totl.x){
  lfdr.list[[j]][i,]=local_fdr_ft_2(i,basis[[non.lin+1]][,j-non.lin],temp_beta_set[[j]],pi0[[j]],tau.res[[j]])
  }
  
  #---------------------------------------------------------------------------------  
  #non.linear
   for(k in 1:non.lin){
    for(j in 1:p){
      temp_beta_set[[k]][j]<-beta_selection_ft(i,j,basis[[k]],temp_beta_set[[k]],tau.res[[k]],lfdr.list[[k]][i,])
      }}
  
  #linear
  for(k in (totl.x-lin+1):totl.x){
    temp_beta_set[[k]]<-beta_selection_ft_2(i,basis[[non.lin+1]][,k-non.lin],temp_beta_set[[k]],tau.res[[k]],lfdr.list[[k]][i,])
  }
  #--------------------------------------------------------------------------------- 
  for(k in 1:totl.x) Beta.list[[k]][i,]<-temp_beta_set[[k]]
  
  #--------------------------------------------------------------------------------- 
  for(h in 1:non.lin){
  pi0[[h]]=rbeta(n=1,shape1=k_pi*eta+length(which(Beta.list[[h]][i,]==0)),shape2=k_pi*(1-eta)+p-length(which(Beta.list[[h]][i,]==0)))
  }
  for(h in (totl.x-lin+1):totl.x){
  pi0[[h]]=rbeta(n=1,shape1=k_pi*eta+length(which(Beta.list[[h]][i,]==0)),shape2=k_pi*(1-eta)+1-length(which(Beta.list[[h]][i,]==0)))
  }     
  #---------------------------------------------------------------------------------
  
  Dtau_nonlin = lapply(1:non.lin, function(k) diag(c(tau.res[[k]])))
  Dtau_lin<-c()
  for(k in 1:q) Dtau_lin[k]=tau.res[[k+non.lin]]
  Dtau_lin    = diag(Dtau_lin)
  #---------------------------------------------------------------------------------
  linear.beta<-c()
  for(k in 1:q) linear.beta[k]=Beta.list[[k+non.lin]][i,]
  
  ZZ       = Reduce("+",lapply(1:totl.x,function(k) length(which(Beta.list[[k]][i,]==0))))
  part1_sg = y-Reduce("+",lapply(1:non.lin,function(k) basis[[k]]%*%Beta.list[[k]][i,]))-basis[[non.lin+1]]%*%linear.beta
  part2_sg = Reduce("+",lapply(1:non.lin, function(k) t(Beta.list[[k]][i,])%*%solve(Dtau_nonlin[[k]])%*%Beta.list[[k]][i,]))
  
  sg[i]<-rinvgamma(n=1,shape=0.5*(n+non.lin*p+q-ZZ-1),
                   scale=0.5*(t(part1_sg)%*%part1_sg+part2_sg+t(linear.beta)%*%solve(Dtau_lin)%*%linear.beta)
  )

 #---------------------------------------------------------------------------------
  mu[i]<-rnorm(n=1,mean=mean(ystar),sd=sqrt(sg[i]/n))
  
  #---------------------------------------------------------------------------------
  #non.linear
  tau.res=lapply(1:non.lin,function(k) (Beta.list[[k]][i,]==0)*rexp(n=p,rate=lambda.res/2)+
                   (Beta.list[[k]][i,]!=0)*1/rinvgauss(n=p,mean=sqrt((lambda.res*sg[i,])/Beta.list[[k]][i,]^2),shape=lambda.res))
  #linear
  for(kk in (totl.x-lin+1):totl.x){
  tau.res[[kk]]<-(Beta.list[[kk]][i,]==0)*rexp(n=1,rate=lambda.res/2)+
    (Beta.list[[kk]][i,]!=0)*1/rinvgauss(n=1,mean=sqrt((lambda.res*sg[i,])/Beta.list[[kk]][i,]^2),shape=lambda.res)
  }
  #---------------------------------------------------------------------------------
  lambda.res=rgamma(n=1,shape=non.lin*p+lin+r,rate=sum(unlist(tau.res))/2+delta)
  if(i%%1==0){cat("Iteration ",i,"\r",sep="")}
  if(i%%10000==0){ts.plot(sg[-c(1:burn)])}
}

#=====================================================================#


data_t<-t<-seq(start_x,end_x,length.out=100)
data_t=matrix(rep(data_t,knot_num+1),nrow=length(data_t),ncol=knot_num+1)
for(i in 1:knot_num){data_t[,i+1]=knot_radial(data_t[,i+1],kn[i])}
data_t[data_t=="NaN"]=0
data_t=apply(data_t,2,standardizing_function)
#-----------------------------------------------------------------------------
alpha=0.15
order.signif       =order(unlist(indiplot),decreasing = T)
j.star             =if(length(which(cumsum(1-unlist(indiplot)[order.signif])/seq(1,non.lin*p+lin,1)<alpha))==0){
  cat("need to reset alpha")}else{max(which(cumsum(1-unlist(indiplot)[order.signif])/seq(1,non.lin*p+lin,1)<alpha))}


phi.alpha          =unlist(indiplot)[order.signif[j.star]] 
K.final            =which(unlist(indiplot)>phi.alpha)

#-----------------------------------------------------------------------------


index =lapply(1:non.lin,function(j) seq((j-1)*p+1,j*p,1))
TT    =vector("list",non.lin) 
RES.1<-c()
for(u in 1:length(K.final)){  
  if(K.final[u]<p+1){RES.1[u]<-1}
  else if(p<K.final[u]&K.final[u]<2*p+1){RES.1[u]<-2}
  else if(2*p<K.final[u]&K.final[u]<3*p+1){RES.1[u]<-3}
  else if(3*p<K.final[u]&K.final[u]<4*p+1){RES.1[u]<-4}
  else if(4*p<K.final[u]&K.final[u]<5*p+1){RES.1[u]<-5}
  else if(5*p<K.final[u]&K.final[u]<6*p+1){RES.1[u]<-6}
  else if(6*p<K.final[u]&K.final[u]<7*p+1){RES.1[u]<-7}
  else if(7*p<K.final[u]&K.final[u]<8*p+1){RES.1[u]<-8}
  else if(8*p<K.final[u]&K.final[u]<9*p+1){RES.1[u]<-9}
  else {RES.1[u]<-10}
  }
RES=data.frame(RES.1,K.final-(RES.1-1)*p,K.final)
colnames(RES)<-c("a","b","c")


for(u in 1:non.lin){
  TT[[u]]=RES$b[which(RES$a==u)]
}

#-----------------------------------------------------------------------------
burn=30000
indiplot=vector("list",totl.x)
# plot of estimate of local fdr 
indicator.list=lapply(1:non.lin,function(x) matrix(0, nrow=iter, ncol=p))
for(j in 1:non.lin){
  for(i in 1:iter){
    indicator.list[[j]][i,]<-(Beta.list[[j]][i,]!=0)*1
  }
}
for(j in 1:non.lin) indicator.list[[j]]<-indicator.list[[j]][-c(1:burn),]
for(k in 1:non.lin) indiplot[[k]]=apply(indicator.list[[k]],2,mean)
for(k in 9:14) indiplot[[k]]<-mean(((Beta.list[[k]]!=0)*1)[-c(1:burn),])

postscript("simlfdr1215_2.eps",horizontal=F,onefile=F,print.it=F,height=3,width=8)
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(unlist(indiplot)[61:126],lwd=2,type="h",xlab="",ylab="",ylim=c(0,1),xlim=c(1,64),xaxt="n")
mtext(side=2,expression(p[j]),2)
for(k in 1:4) abline(v=p*k+0.5,lty=2,col=1,lwd=1)
abline(h=phi.alpha,lty=1,col=1)
#points(K.final,unlist(indiplot)[K.final],pty=2,pch=1)
points(c(61,62),unlist(indiplot)[c(61+60,62+60)],pty=2,pch=1)
mtext(side=1,expression(f[5]),line=0.75,at=7.5);mtext(side=1,expression(f[6]),line=0.75,at=7.5+15)
mtext(side=1,expression(f[7]),line=0.75,at=7.5+30);mtext(side=1,expression(f[8]),line=0.75,at=7.5+45)
mtext(side=1,expression(beta),line=0.75,at=3.5+60)
dev.off()
#-----------------------------------------------------------
postscript("simsg1215.eps",horizontal=F,onefile=F,print.it=F,height=2,width=8)
par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfrow=c(1,3))
plot(0,xaxt="n",yaxt="n",bty="n",pch="",xlab="",ylab="")
hist(sg[-c(1:burn)],main="",xlab="",ylab="Posterior density",freq=F)
abline(v=true_sg^2,lty=2,lwd=2);mtext(side=1,expression(sigma^2),2.4)
dev.off()

#-----------------------------------------------------------
postscript("simpar1215.eps",horizontal=F,onefile=F,print.it=F,height=5,width=8)

layout(matrix(c(1,2,2,3,3,4,4,5,6,6,7,7,8,8,9,9),2,8,byrow=T))
par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfrow=c(3,2))
#plot(0,xaxt="n",yaxt="n",bty="n",pch="",xlab="",ylab="")

hist(Beta.list[[9]][-c(1:burn),1]/sd(basis.candid[,1]),main="",xlab="",ylab="",freq=F)
abline(v=true.line.beta[1],lty=2,lwd=2);mtext(side=1,expression(beta[1]),2.5);mtext(side=2,"Posterior density",cex=0.9,2.5)
hist(Beta.list[[10]][-c(1:burn),1]/sd(basis.candid[,2]),main="",xlab="",ylab="",freq=F)
abline(v=true.line.beta[2],lty=2,lwd=2);mtext(side=1,expression(beta[2]),2.5);mtext(side=2,"Posterior density",cex=0.9,2.5)
hist(Beta.list[[11]][-c(1:burn),1]/sd(basis.candid[,3]),main="",xlab="",ylab="",freq=F)
abline(v=true.line.beta[3],lty=2,lwd=2);mtext(side=1,expression(beta[3]),2.5);mtext(side=2,"Posterior density",cex=0.9,2.5)

#plot(0,xaxt="n",yaxt="n",bty="n",pch="",xlab="",ylab="")

hist(Beta.list[[12]][-c(1:burn),1]/sd(basis.candid[,4]),main="",xlab="",ylab="",freq=F)
abline(v=true.line.beta[4],lty=2,lwd=2);mtext(side=1,expression(beta[4]),2.5);mtext(side=2,"Posterior density",cex=0.9,2.5)
hist(Beta.list[[13]][-c(1:burn),1]/sd(basis.candid[,5]),main="",xlab="",ylab="",freq=F)
abline(v=true.line.beta[5],lty=2,lwd=2);mtext(side=1,expression(beta[5]),2.5);mtext(side=2,"Posterior density",cex=0.9,2.5)
hist(Beta.list[[14]][-c(1:burn),1]/sd(basis.candid[,6]),main="",xlab="",ylab="",freq=F)
abline(v=true.line.beta[6],lty=2,lwd=2);mtext(side=1,expression(beta[6]),2.5);mtext(side=2,"Posterior density",cex=0.9,2.5)

par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfrow=c(1,1))
hist(sg[-c(1:burn)],main="",xlab="",ylab="",freq=F);mtext(side=2,"Posterior density",cex=0.9,2.5)
abline(v=true_sg^2,lty=2,lwd=2);mtext(side=1,expression(sigma^2),2.4)


dev.off()

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

# save data
save.image("simulation_1215.RData")


#=============================================================================
#making plot==================================================================
#=============================================================================
postscript("sim1nonpar1215.eps",horizontal=F,onefile=F,print.it=F,height=8,width=8)
par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfcol=c(3,4))
cord.x  =c(t,rev(t))
for(nn in 1:4){
plot(c(seq(start_x,end_x,(end_x-start_x)/knot_num))[-1],indiplot[[nn]][-1],type="h",lwd=2.5,xlab="",xlim=range(start_x,end_x),ylim=c(0,1),ylab="")
points(start_x,indiplot[[nn]][1],type="h",lty=3)
points(start_x,indiplot[[nn]][1],pch=0) 
abline(h=phi.alpha,lty=2,col=1)
mtext("x", side=1, line=2) 
mtext("Probability", side=2, line=2.5,cex=0.9)

#---------------------------------------------------------------------
plot(x=c(start_x,end_x),y=c(min(cord.y.1[[nn]]-0.5),max(cord.y.1[[nn]])+0.5),type="n",xlab="",ylab="", main="")
polygon(cord.x,cord.y.1[[nn]],col=gray(0:9/9)[8],border=NA)
lines(t,Quantile.list.1[[nn]][,2])
lines(t,true_mu+function.list[[nn]](t)-mean(function.list[[nn]](t)),lty=2)
mtext("x", side=1, line=2)
mtext(expression(italic(paste(f(x)))), side=2, line=2.5)

plot(x=c(start_x,end_x),y=c(min(cord.y.2[[nn]]-0.5),max(cord.y.2[[nn]])+0.5),type="n",xlab="",ylab="", main="")
mtext("x", side=1, line=2)
mtext(expression(italic(paste(f(x)))), side=2, line=2.5)
polygon(cord.x,cord.y.2[[nn]],col=gray(0:9/9)[8],border=NA)
lines(t,Quantile.list.2[[nn]][,2])
lines(t,true_mu+function.list[[nn]](t)-mean(function.list[[nn]](t)),lty=2)
}
dev.off()


#---------------------------------------------------------------------
par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfcol=c(1,3))

hist(Beta.list[[4]][-c(1:burn),1]/sd(DATA[[4]]),xlab="",ylab="",main="",freq=F)
mtext(expression(paste(beta[1])), side=1, line=2.5,cex=0.9);abline(v=true.line.beta[1],lty=2,lwd=2,col=1)
mtext("Posterior density",side=2,2.2,cex=0.9)
ts.plot(Beta.list[[4]][c(10001:20000),1]/sd(DATA[[4]]),xlab="",ylab="",main="")

hist(Beta.list[[5]][-c(1:burn),1]/sd(DATA[[5]]),xlab="",ylab="",main="",freq=F)
mtext(expression(paste(beta[2])), side=1, line=2.5,cex=0.9);abline(v=true.line.beta[2],lty=2,lwd=2,col=1)
mtext("Posterior density",side=2,2.2,cex=0.9)
ts.plot(Beta.list[[5]][c(10001:20000),1]/sd(DATA[[5]]),xlab="",ylab="",main="")

hist(sg[-c(1:burn)],xlab="",ylab="",main="",freq=F);abline(v=true_sg^2,lty=2,lwd=2,col=1)
mtext(expression(paste(sigma^2)), side=1, line=2.5,cex =0.9)
mtext("Posterior density",side=2,2.2,cex=0.9)
ts.plot(sg[c(10001:20000)],xlab="",ylab="",main="")


#----------------------------------------------
library(glmnet)
knot_num=8
#knot_radial=function(x,knot){abs((x-knot))^2*log(abs((x-knot))^2)}
# knot_radial            =function(x,knot){abs(x-knot)^3}
# kn=(max(dataset$age)-min(dataset$age))*ppoints(knot_num)+min(dataset$age)
# M.DATA<-matrix(dataset$age,ncol=knot_num+1,nrow=314)
# for(i in 1:knot_num) M.DATA[,i+1]<-knot_radial(M.DATA[,i+1],kn[i])
# basis1=M.DATA
# 
# kn=(max(dataset$chol)-min(dataset$chol))*ppoints(knot_num)+min(dataset$chol)
# M.DATA<-matrix(dataset$chol,ncol=knot_num+1,nrow=314)
# for(i in 1:knot_num) M.DATA[,i+1]<-knot_radial(M.DATA[,i+1],kn[i])
# basis2=M.DATA

model_matrix=matrix(
  cbind(basis[[1]],basis[[2]],basis[[3]],basis[[4]],
    basis[[5]],basis[[6]],basis[[7]],basis[[8]],basis[[9]]),nrow=800)
warnings()
cv.lasso   = cv.glmnet(model_matrix,y,standardize =F,nfolds=20)
lasso.coef = predict(cv.lasso,type="coefficients", s=cv.lasso$lambda.min) # coefficients

length(which(lasso.coef!=0))

#
res.lasso.1=data.frame(DATA[[1]],basis[[1]]%*%lasso.coef[c(2:16),1])
res.lasso.1=res.lasso.1[order(res.lasso.1[,1]),]
#
res.lasso.2=data.frame(DATA[[2]],basis[[2]]%*%lasso.coef[c(17:31),1])
res.lasso.2=res.lasso.2[order(res.lasso.2[,1]),]
#
res.lasso.3=data.frame(DATA[[3]],basis[[3]]%*%lasso.coef[c(32:46),1])
res.lasso.3=res.lasso.3[order(res.lasso.3[,1]),]
#
res.lasso.4=data.frame(DATA[[4]],basis[[4]]%*%lasso.coef[c(47:61),1])
res.lasso.4=res.lasso.4[order(res.lasso.4[,1]),]
#
res.lasso.5=data.frame(DATA[[5]],basis[[5]]%*%lasso.coef[c(62:76),1])
res.lasso.5=res.lasso.5[order(res.lasso.5[,1]),]
#
res.lasso.6=data.frame(DATA[[6]],basis[[6]]%*%lasso.coef[c(77:91),1])
res.lasso.6=res.lasso.6[order(res.lasso.6[,1]),]
#
res.lasso.7=data.frame(DATA[[7]],basis[[7]]%*%lasso.coef[c(92:106),1])
res.lasso.7=res.lasso.7[order(res.lasso.7[,1]),]
#
res.lasso.8=data.frame(DATA[[8]],basis[[8]]%*%lasso.coef[c(107:121),1])
res.lasso.8=res.lasso.8[order(res.lasso.8[,1]),]
#
#Res.lasso
lasso.coef[122]/sd(basis.candid[,1])
lasso.coef[123]/sd(basis.candid[,2])
lasso.coef[124]/sd(basis.candid[,3])
lasso.coef[125]/sd(basis.candid[,4])
lasso.coef[126]/sd(basis.candid[,5])
lasso.coef[127]/sd(basis.candid[,6])
#-----------------------------------------------------------
library(ncvreg)
fit.scad<-cv.ncvreg(model_matrix,y,max.iter=100000,family="gaussian",penalty="SCAD")
summary(fit.scad)$min
which(colnames(fit.scad$fit$beta)=="0.02219")
fit.scad$fit$beta[,summary(fit.scad)$min]
length(which(fit.scad$fit$beta[,summary(fit.scad)$min]!=0))

res.scad.1=data.frame(DATA[[1]],basis[[1]]%*%fit.scad$fit$beta[c(2:16),53])
res.scad.1=res.scad.1[order(res.scad.1[,1]),]
#
res.scad.2=data.frame(DATA[[2]],basis[[2]]%*%fit.scad$fit$beta[c(17:31),53])
res.scad.2=res.scad.2[order(res.scad.2[,1]),]
#
res.scad.3=data.frame(DATA[[3]],basis[[3]]%*%fit.scad$fit$beta[c(32:46),53])
res.scad.3=res.scad.3[order(res.scad.3[,1]),]
#
res.scad.4=data.frame(DATA[[4]],basis[[4]]%*%fit.scad$fit$beta[c(47:61),53])
res.scad.4=res.scad.4[order(res.scad.4[,1]),]
#
res.scad.5=data.frame(DATA[[5]],basis[[5]]%*%fit.scad$fit$beta[c(62:76),53])
res.scad.5=res.scad.5[order(res.scad.5[,1]),]
#
res.scad.6=data.frame(DATA[[6]],basis[[6]]%*%fit.scad$fit$beta[c(77:91),53])
res.scad.6=res.scad.6[order(res.scad.6[,1]),]
#
res.scad.7=data.frame(DATA[[7]],basis[[7]]%*%fit.scad$fit$beta[c(92:106),53])
res.scad.7=res.scad.7[order(res.scad.7[,1]),]
#
res.scad.8=data.frame(DATA[[8]],basis[[8]]%*%fit.scad$fit$beta[c(107:121),53])
res.scad.8=res.scad.8[order(res.scad.8[,1]),]

fit.scad$fit$beta[122,53]/sd(basis.candid[,1])
fit.scad$fit$beta[123,53]/sd(basis.candid[,2])
fit.scad$fit$beta[124,53]/sd(basis.candid[,3])
fit.scad$fit$beta[125,53]/sd(basis.candid[,4])
fit.scad$fit$beta[126,53]/sd(basis.candid[,5])
fit.scad$fit$beta[127,53]/sd(basis.candid[,6])

#-----------------------------------------------------------
par(mex=0.1,mar=c(3,3,.1,.1)+1,cex=0.9,mfrow=c(2,2))
cord.x  =c(t,rev(t))
nn      =8
  plot(x=c(start_x,end_x),y=c(min(cord.y.2[[nn]]-0.5),max(cord.y.2[[nn]])+0.5),type="n",xlab="",ylab="", main="")
  polygon(cord.x,cord.y.2[[nn]],col=gray(0:9/9)[8],border=NA)
  lines(t,true_mu+function.list[[nn]](t)-mean(function.list[[nn]](t)),lty=1)
  lines(t,Quantile.list.2[[nn]][,2],lty=2)
  lines(res.lasso.8[,1],res.lasso.8[,2],lty=3) # lasso result
  lines(res.scad.8[,1],res.scad.8[,2],lty=4) #scad result
  mtext("x", side=1, line=2)
  mtext(expression(italic(paste(f[8](x)))), side=2, line=2.5)
  legend("bottom",cex=.8,legend=c("True","Proposed","Lasso","SCAD"),lty=c(1,2,3,4))
  
  