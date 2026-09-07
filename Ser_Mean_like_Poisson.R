#######################################################################

args <- commandArgs(trailingOnly = TRUE)

core <- as.numeric(args[1])

#######################################################################



#### We should find ka to control the censoring rate ##################

find.ka<-function(rate){
  
  set.seed(0)
  
  low=0.01; high=10; tol=1e-8
  
  n=50000
  
  calculate.rate<-function(ka){
    
    X=runif(n); Z1=rbinom(n,1,0.5); Z2=rnorm(n,0,1); 
    
    rateU=exp(gammaU[1]*X + gammaU[2]*Z1 + gammaU[3]*Z2); 
    
    U=pmin(5+rexp(n,rateU),rep(tau,n));C=pmin(5+ka*rexp(n,1),rep(tau,n))
    
    delta=ifelse(U<C,1,0)
    
    c.rate=1-mean(delta)
    
    return(c.rate)
    
  }
  
  c.rate=1
  
  while(abs(c.rate-rate)>tol&low<high){
    
    ka=(low+high)/2
    
    c.rate=calculate.rate(ka)
    
    if(c.rate<rate) high=ka   else low=ka
    
  }
  
  return(ka)
  
}

#######################################################################



################# Generate the Data Set ###############################

Gen.Data=function(n){
  
  X=runif(n); Z1=rbinom(n,1,0.5); Z2=rnorm(n,0,1);
  
  rateU=exp(gammaU[1]*X + gammaU[2]*Z1 + gammaU[3]*Z2);
  
  U=pmin(5+rexp(n,rateU),rep(tau,n));C=pmin(5+ka*rexp(n,1),rep(tau,n));Y=pmin(U,C);Delta=(Y==U)+0
  
  K=sample(1:Obs,n,TRUE);T=matrix(NA,n,Obs)
  
  for(i in 1:n){T[i,1:K[i]]=sort(runif(K[i]))*Y[i]}   
  
  
  
  
  
  rateN=exp(alpha[1]*Z1 + alpha[2]*Z2 + 2*X^2); #Linear: beta_0(x)=3x
  
  
  
  N=matrix(NA,n,Obs);
  
  # gam=rgamma(n,shape=2,scale=1/2)  
  
  gam=rep(1,n)
  
  if(group==1){for(i in 1:n){
    
    N[i,1]=rpois(1,gam[i]*(LAMBDA(U[i])-LAMBDA(U[i]-T[i,1]))*rateN[i])
    
    if(K[i]!=1){for(j in 2:K[i]){N[i,j]=rpois(1,gam[i]*rateN[i]*(LAMBDA(U[i]-T[i,j-1])-LAMBDA(U[i]-T[i,j])))+N[i,j-1]}}
    
  }}
  
  if(group==2){for(i in 1:n){
    
    N[i,1]=rnbinom(1,(LAMBDA(U[i])-LAMBDA(U[i]-T[i,1]))*rateN[i],0.5)
    
    if(K[i]!=1){for(j in 2:K[i]){N[i,j]=rnbinom(1,rateN[i]*(LAMBDA(U[i]-T[i,j-1])-LAMBDA(U[i]-T[i,j])),0.5)+N[i,j-1]}}
    
  }}
  
  dN=matrix(NA,length(K),Obs);dN[,1]=N[,1];for(j in 2:Obs){dN[,j]=N[,j]-N[,(j-1)]}
  
  
  
  
  
  return(list(Y=Y,Delta=Delta,N=N,T=T,K=K,X=X,Z1=Z1,Z2=Z2,dN=dN))
  
}

#######################################################################



########### The Lambda Function #######################################

LAMBDA1=function(s){L=s;return(L)}

# LAMBDA2=function(s){L=3*sqrt(s);return(L)}

#######################################################################



############# Calculate the knots of Spline ###########################

cal.knots=function(n){
  
  set.seed(80003)
  
  X=runif(n); Z1=rbinom(n,1,0.5); Z2=rnorm(n,0,1); 
  
  rateU=exp(gammaU[1]*X + gammaU[2]*Z1 + gammaU[3]*Z2);
  
  U=pmin(5+rexp(n,rateU),rep(tau,n));C=pmin(5+ka*rexp(n,1),rep(tau,n));Y=pmin(U,C);Delta=(Y==U)+0
  
  K=sample(1:Obs,n,TRUE);T=matrix(NA,n,Obs)
  
  for(i in 1:n){T[i,1:K[i]]=sort(runif(K[i]))*Y[i]}   
  
  ymt=matrix(NA,n,max(K));for(i in 1:n){for(j in 1:K[i]){ymt[i,j]=Y[i]-T[i,j]}}
  
  knots=c(quantile(ymt,seq(1/(nknots1+1),nknots1/(nknots1+1),1/(nknots1+1)),na.rm=TRUE,names=FALSE))
  
  return(knots)
  
}

#######################################################################



############# Calculate Spline at Each Time Point #####################

CalBB=function(Data,t){
  
  mm=length(t);nnn=length(Data$Y);T=Data$T;K=Data$K;T=cbind(rep(0,nnn),T);TT=na.omit(as.vector(t(T)))
  
  BB1=bSpline(c(0,pmax(rep(t,each=length(TT))-rep(TT,mm),0),tau),df=INN1,knots=knots1,degree=(INN1-1),intercept=FALSE)
  
  dBB1=array(NA,c(mm,nnn,max(K),INN1+nknots1-1));
  
  index=2
  
  for(m in 1:mm){for(i in 1:nnn){
    
    for(j in 1:K[i]){dBB1[m,i,j,]=BB1[index,]-BB1[(index+1),];index=index+1};index=index+1
    
  }}
  
  return(list(dBB1=dBB1))
  
}

#######################################################################



############# Calculate the LIKELIHOOD Loss Function ##################

Loss=function(theta){
  
  L=0;Y=Data$Y;Delta=Data$Delta;T=Data$T;K=Data$K;dN=Data$dN;mm=length(t);
  
  Z1=Data$Z1; Z2=Data$Z2; 
  
  
  
  len1 = INN1+nknots1-1
  
  len2 = INN2+nknots2-1
  
  id.xi1   = 1:len1
  
  id.xi2   = (len1+1) : (len1+len2)          
  
  id.alpha = (len1+len2+1) : length(theta)   
  
  
  
  
  
  expcoeff=exp(as.matrix(cbind(Z1, Z2))%*%theta[id.alpha] + BB2_X%*%theta[id.xi2])
  
  
  
  for(i in which(Delta==1)){
    
    for(j in 1:K[i]){
      
      mu = c(expcoeff[i]*(t(theta[id.xi1])%*%dBB1[rank[i],i,j,]))
      
      mu = max(mu, 1e-8)  
      
      L = L + (mu - dN[i,j] * log(mu))
      
    }
    
  }
  
  for(i in which(Delta==0)){
    
    if(sur[rank[i]]!=0&&(rank[i]!=mm)){
      
      temp=0;
      
      for(j in 1:K[i]){
        
        for(m in (rank[i]+1):mm){
          
          prob_diff = (sur[m-1]^Uexpcoeff[i] - sur[m]^Uexpcoeff[i])
          
          mu = c(expcoeff[i]*t(theta[id.xi1])%*%dBB1[m,i,j,])
          
          mu = max(mu, 1e-8)  
          
          temp = temp + prob_diff * (mu - dN[i,j] * log(mu))
          
        }
        
      }
      
      if((sur[rank[i]]^Uexpcoeff[i])>0.1^8){
        
        temp = temp / (sur[rank[i]]^Uexpcoeff[i])
        
        L = L + temp
        
      }
      
    }
    
  }
  
  return(L)
  
}

#######################################################################



######### Calculate 1st derivative (Speed Up the Optim) ###############

cal.DE <- function(theta) {
  
  gk = rep(0, length(theta)); 
  
  Y = Data$Y; Delta = Data$Delta; T = Data$T; K = Data$K; dN = Data$dN; mm = length(t);
  
  Z1 = Data$Z1; Z2 = Data$Z2; 
  
  
  
  len1 = INN1+nknots1-1; len2 = INN2+nknots2-1
  
  id.xi1   = 1:len1
  
  id.xi2   = (len1+1) : (len1+len2)
  
  id.alpha = (len1+len2+1) : length(theta)
  
  
  
  # 预先计算出 E_i
  
  expcoeff = exp(as.matrix(cbind(Z1, Z2)) %*% theta[id.alpha] + BB2_X %*% theta[id.xi2])
  
  
  
  # (1) 对于发生事件的样本 Delta == 1
  
  for(i in which(Delta==1)) {
    
    for(j in 1:K[i]) {
      
      mu = c(expcoeff[i] * (t(theta[id.xi1]) %*% dBB1[rank[i],i,j,]))
      
      mu = max(mu, 1e-8)  # 防御机制
      
      
      
      diff_mu = mu - dN[i,j]
      
      
      
      # 依据公式累加梯度
      
      gk[id.xi1]   = gk[id.xi1]   + (1 - dN[i,j]/mu) * expcoeff[i] * dBB1[rank[i],i,j,]
      
      gk[id.xi2]   = gk[id.xi2]   + diff_mu * BB2_X[i,]
      
      gk[id.alpha] = gk[id.alpha] + diff_mu * c(Z1[i], Z2[i])
      
    }
    
  }
  
  
  
  # (2) 对于删失的样本 Delta == 0
  
  for(i in which(Delta==0)) {
    
    if(sur[rank[i]] != 0 && (rank[i] != mm)) {
      
      temp_gk = rep(0, length(theta))
      
      for(j in 1:K[i]) {
        
        for(m in (rank[i]+1):mm) {
          
          prob_diff = (sur[m-1]^Uexpcoeff[i] - sur[m]^Uexpcoeff[i])
          
          
          
          mu = c(expcoeff[i] * (t(theta[id.xi1]) %*% dBB1[m,i,j,]))
          
          mu = max(mu, 1e-8)
          
          
          
          diff_mu = mu - dN[i,j]
          
          
          
          # 累加带概率权重的梯度
          
          temp_gk[id.xi1]   = temp_gk[id.xi1]   + prob_diff * (1 - dN[i,j]/mu) * expcoeff[i] * dBB1[m,i,j,]
          
          temp_gk[id.xi2]   = temp_gk[id.xi2]   + prob_diff * diff_mu * BB2_X[i,]
          
          temp_gk[id.alpha] = temp_gk[id.alpha] + prob_diff * diff_mu * c(Z1[i], Z2[i])
          
        }
        
      }
      
      # 除以分母的生存概率并累加到总梯度
      
      if((sur[rank[i]]^Uexpcoeff[i]) > 0.1^8) {
        
        gk = gk + temp_gk / (sur[rank[i]]^Uexpcoeff[i])
        
      }
      
    }
    
  }
  
  
  
  return(gk)
  
}

#######################################################################



################## Generate Bootstrap Data ############################

Genbots.Data=function(Data){
  
  nnn=length(Data$Y)
  
  id=sample(1:nnn,nnn,rep=T)
  
  Y=Data$Y[id];Delta=Data$Delta[id];N=Data$N[id,];T=Data$T[id,];K=Data$K[id]
  
  X=Data$X[id];Z1=Data$Z1[id];Z2=Data$Z2[id];dN=Data$dN[id,];
  
  return(list(Y=Y,Delta=Delta,N=N,T=T,K=K,X=X,Z1=Z1,Z2=Z2,dN=dN,id=id))
  
}

#######################################################################



############# Calculate Bootstrap LIKELIHOOD Loss Function ############

bots.Loss=function(theta){
  
  L=0;Y=bots.Data$Y;Delta=bots.Data$Delta;T=bots.Data$T;K=bots.Data$K;dN=bots.Data$dN;mm=length(bots.t);
  
  Z1=bots.Data$Z1; Z2=bots.Data$Z2;
  
  
  
  len1 = INN1+nknots1-1; len2 = INN2+nknots2-1
  
  id.xi1=1:len1; id.xi2=(len1+1):(len1+len2); id.alpha=(len1+len2+1):length(theta)
  
  
  
  #
  
  expcoeff=exp(as.matrix(cbind(Z1, Z2))%*%theta[id.alpha] + bots.BB2_X%*%theta[id.xi2])
  
  
  
  for(i in which(Delta==1)){
    
    for(j in 1:K[i]){
      
      mu = c(expcoeff[i]*(t(theta[id.xi1])%*%bots.dBB1[bots.rank[i],i,j,]))
      
      mu = max(mu, 1e-8)
      
      L = L + (mu - dN[i,j] * log(mu))
      
    }
    
  }
  
  for(i in which(Delta==0)){
    
    if(bots.sur[bots.rank[i]]!=0&&(bots.rank[i]!=mm)){
      
      temp=0;
      
      for(j in 1:K[i]){
        
        for(m in (bots.rank[i]+1):mm){
          
          prob_diff = (bots.sur[m-1]^bots.Uexpcoeff[i] - bots.sur[m]^bots.Uexpcoeff[i])
          
          mu = c(expcoeff[i]*t(theta[id.xi1])%*%bots.dBB1[m,i,j,])
          
          mu = max(mu, 1e-8)
          
          temp = temp + prob_diff * (mu - dN[i,j] * log(mu))
          
        }
        
      }
      
      if((bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])>0.1^8){
        
        temp = temp / (bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])
        
      } else {
        
        temp = 0
        
      }
      
      L = L + temp;
      
    }
    
  }
  
  return(L)
  
}

#######################################################################



######### Calculate Bootstrap 1st derivative ##########################

bots.DE <- function(theta) {
  
  gk = rep(0, length(theta));
  
  Y = bots.Data$Y; Delta = bots.Data$Delta; T = bots.Data$T; K = bots.Data$K; dN = bots.Data$dN; mm = length(bots.t);
  
  Z1 = bots.Data$Z1; Z2 = bots.Data$Z2;
  
  
  
  len1 = INN1+nknots1-1; len2 = INN2+nknots2-1
  
  id.xi1   = 1:len1
  
  id.xi2   = (len1+1) : (len1+len2)
  
  id.alpha = (len1+len2+1) : length(theta)
  
  
  
  expcoeff = exp(as.matrix(cbind(Z1, Z2)) %*% theta[id.alpha] + bots.BB2_X %*% theta[id.xi2])
  
  
  
  for(i in which(Delta==1)) {
    
    for(j in 1:K[i]) {
      
      mu = c(expcoeff[i] * (t(theta[id.xi1]) %*% bots.dBB1[bots.rank[i],i,j,]))
      
      mu = max(mu, 1e-8)
      
      diff_mu = mu - dN[i,j]
      
      
      
      gk[id.xi1]   = gk[id.xi1]   + (1 - dN[i,j]/mu) * expcoeff[i] * bots.dBB1[bots.rank[i],i,j,]
      
      gk[id.xi2]   = gk[id.xi2]   + diff_mu * bots.BB2_X[i,]
      
      gk[id.alpha] = gk[id.alpha] + diff_mu * c(Z1[i], Z2[i])
      
    }
    
  }
  
  
  
  for(i in which(Delta==0)) {
    
    if(bots.sur[bots.rank[i]] != 0 && (bots.rank[i] != mm)) {
      
      temp_gk = rep(0, length(theta))
      
      for(j in 1:K[i]) {
        
        for(m in (bots.rank[i]+1):mm) {
          
          prob_diff = (bots.sur[m-1]^bots.Uexpcoeff[i] - bots.sur[m]^bots.Uexpcoeff[i])
          
          
          
          mu = c(expcoeff[i] * (t(theta[id.xi1]) %*% bots.dBB1[m,i,j,]))
          
          mu = max(mu, 1e-8)
          
          diff_mu = mu - dN[i,j]
          
          
          
          temp_gk[id.xi1]   = temp_gk[id.xi1]   + prob_diff * (1 - dN[i,j]/mu) * expcoeff[i] * bots.dBB1[m,i,j,]
          
          temp_gk[id.xi2]   = temp_gk[id.xi2]   + prob_diff * diff_mu * bots.BB2_X[i,]
          
          temp_gk[id.alpha] = temp_gk[id.alpha] + prob_diff * diff_mu * c(Z1[i], Z2[i])
          
        }
        
      }
      
      if((bots.sur[bots.rank[i]]^bots.Uexpcoeff[i]) > 0.1^8) {
        
        gk = gk + temp_gk / (bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])
        
      }
      
    }
    
  }
  
  
  
  return(gk)
  
}

#######################################################################



library(survival);

library(splines2);

group=1;case=1;LAMBDA=LAMBDA1;n=200;tau=10;Obs=6;INN1=4;INN2=4;nknots1=round(2*n^(1/5));nknots2=round(n^(1/5));num.sim=2;num.bots=100;



gammaU=c(0.5, -1, 0.5); 

alpha=c(-2, 0.5); 

ka=find.ka(0.2);knots1=cal.knots(10000);knots2=seq(1/(nknots2+1),nknots2/(nknots2+1),1/(nknots2+1))



len1 = INN1+nknots1-1
len2 = INN2+nknots2-1
num_para = length(alpha)
total_len = len1 + len2 + num_para


theta.est=matrix(NA,num.sim,total_len); gammaU.est=matrix(NA,num.sim,length(gammaU));

bots.theta=array(NA,c(num.sim,num.bots,total_len));





## Constraint matrix for theta = c(xi1, xi2, alpha)
## 1) xi1[1] >= 0 and xi1 is nondecreasing
D1 = cbind(
  diag(-1, len1)[, 2:len1, drop = FALSE],
  rep(0, len1)
) + diag(1, len1)

AA.lambda = cbind(
  D1,
  matrix(0, nrow = len1, ncol = len2 + num_para)
)
BB.lambda = rep(0, len1)

## 2) alpha >= -5
AA.alpha.lower = cbind(
  matrix(0, nrow = num_para, ncol = len1 + len2),
  diag(1, num_para)
)

## 3) alpha <= 5, equivalently -alpha >= -5
AA.alpha.upper = cbind(
  matrix(0, nrow = num_para, ncol = len1 + len2),
  diag(-1, num_para)
)

AA = rbind(
  AA.lambda,
  AA.alpha.lower,
  AA.alpha.upper
)

BB = c(
  BB.lambda,
  rep(-5, num_para),
  rep(-5, num_para)
)

stopifnot(
  ncol(AA) == total_len,
  nrow(AA) == length(BB)
)


for(id.sim in 1:num.sim){
  
  i.sim=id.sim+(core-1)*num.sim
  
  set.seed(i.sim);
  
  Data=Gen.Data(n);
  
  
  
  Init=c(1:len1, rep(0, len2 + num_para))
  
  stopifnot(
    length(Init) == total_len,
    all(drop(AA %*% Init - BB) >= 0)
  )
  
  
  cox=coxph(Surv(Y, Delta) ~ X + Z1 + Z2, Data); 
  
  gammaUhat=cox$coefficients;sur=survfit(cox)$surv;t=survfit(cox)$time
  
  Uexpcoeff=exp(cbind(Data$X, Data$Z1, Data$Z2)%*%gammaUhat);rank=rep(NA,n); 
  
  for(i in 1:n){if(length(which(t==Data$Y[i]))==1){rank[i]=which(t==Data$Y[i])}else{rank[i]=which(abs(t-Data$Y[i])==min(abs(t-Data$Y[i])))}}
  
  dBB1=CalBB(Data,t)$dBB1;
  
  
  
  
  BB2_shared = bSpline(c(0, Data$X, 1), df=INN2, knots=knots2, degree=(INN2-1), intercept=FALSE)
  
  BB2_X = BB2_shared[ 2 : (n+1), ]
  
  
  
  optim.est=constrOptim(Init, Loss, grad=cal.DE, ui=AA, ci=BB, method="BFGS")
  
  Init=optim.est$par;theta.est[id.sim,]=optim.est$par;gammaU.est[id.sim,]=gammaUhat;
  
  
  
  ########## Bootstrap ############
  
  for(i.bots in 1:num.bots){
    
    set.seed(i.sim*num.bots+i.bots+1000);
    
    bots.Data=Genbots.Data(Data)
    
    
    
    
    
    bots.cox=coxph(Surv(Y, Delta) ~ X + Z1 + Z2, bots.Data)
    
    bots.gammaUhat=bots.cox$coefficients;bots.sur=survfit(bots.cox)$surv;bots.t=survfit(bots.cox)$time
    
    bots.Uexpcoeff=exp(cbind(bots.Data$X, bots.Data$Z1, bots.Data$Z2)%*%bots.gammaUhat)
    
    
    
    bots.rank=rep(NA,n);
    
    for(i in 1:n){
      
      if(length(which(bots.t==bots.Data$Y[i]))==1){bots.rank[i]=which(bots.t==bots.Data$Y[i])}
      
      else{bots.rank[i]=which(abs(bots.t-bots.Data$Y[i])==min(abs(bots.t-bots.Data$Y[i])))[1]}
      
    }
    
    bots.dBB1=CalBB(bots.Data,bots.t)$dBB1;
    
    
    
    
    bots.BB2_shared = bSpline(c(0, bots.Data$X, 1), df=INN2, knots=knots2, degree=(INN2-1), intercept=FALSE)
    
    bots.BB2_X = bots.BB2_shared[ 2 : (n+1), ]
    
    
    
    
    
    bots.est=constrOptim(Init, bots.Loss, grad=bots.DE, ui=AA, ci=BB, method="BFGS")
    
    bots.theta[id.sim,i.bots,]=bots.est$par
    
  }
  
  ##################################
  
  write.table(theta.est[id.sim,],file=paste0("n",n,"thetaest",i.sim,".txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
  
  write.table(bots.theta[id.sim,,],file=paste0("n",n,"botsTheta",i.sim,".txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
  
}

