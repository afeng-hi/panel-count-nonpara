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
    ### 🟢 X服从均匀分布, Z1服从二项分布, Z2服从标准正态分布
    X=runif(n); Z1=rbinom(n,1,0.5); Z2=rnorm(n,0,1);
    ### 🟢 对应 gammaU 的三个系数
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
  ### 🟢 数据生成
  X=runif(n); Z1=rbinom(n,1,0.5); Z2=rnorm(n,0,1);
  rateU=exp(gammaU[1]*X + gammaU[2]*Z1 + gammaU[3]*Z2);
  U=pmin(5+rexp(n,rateU),rep(tau,n));C=pmin(5+ka*rexp(n,1),rep(tau,n));Y=pmin(U,C);Delta=(Y==U)+0
  K=sample(1:Obs,n,TRUE);T=matrix(NA,n,Obs)
  for(i in 1:n){T[i,1:K[i]]=sort(runif(K[i]))*Y[i]}   
  
  ### 🟢 beta_0(x) =2*X^2，alpha =c(-2, 0.5), Linear:beta_0(x)=3X
  rateN=exp(alpha[1]*Z1 + alpha[2]*Z2 + 2*X^2);
  
  N=matrix(NA,n,Obs);
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
#######################################################################

############# Calculate the knots of Spline ###########################
cal.knots=function(n){
  set.seed(80003)
  ### 🟢 同步节点计算
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

############# Calculate the Loss Function #############################
Loss=function(theta){
  L=0;Y=Data$Y;Delta=Data$Delta;T=Data$T;K=Data$K;dN=Data$dN;mm=length(t)
  ### 🟢 提取 X, Z1, Z2
  X=Data$X; Z1=Data$Z1; Z2=Data$Z2;
  
  ### 🟢 回归到单一非参数项的索引
  len1 = INN1+nknots1-1
  len2 = INN2+nknots2-1
  id.xi1   = 1:len1
  id.xi2   = (len1+1) : (len1+len2)          
  id.alpha = (len1+len2+1) : length(theta)
  
  ### 🟢 线性预测包含 alpha_1*Z1 + alpha_2*Z2
  Zmat = cbind(Z1, Z2)
  expcoeff=exp(Zmat %*% theta[id.alpha] + BB2[2:(n+1),] %*% theta[id.xi2])
  
  for(i in which(Delta==1)){for(j in 1:K[i]){L=L+(dN[i,j]-c(expcoeff[i]*(t(theta[id.xi1])%*%dBB1[rank[i],i,j,])))^2}}
  for(i in which(Delta==0)){
    if(sur[rank[i]]!=0&&(rank[i]!=mm)){
      temp=0;
      for(j in 1:K[i]){for(m in (rank[i]+1):mm){
        temp=temp+(sur[m-1]^Uexpcoeff[i]-sur[m]^Uexpcoeff[i])*(dN[i,j]-c(expcoeff[i]*t(theta[id.xi1])%*%dBB1[m,i,j,]))^2
      }}
      if((sur[rank[i]]^Uexpcoeff[i])>0.1^8){temp=temp/(sur[rank[i]]^Uexpcoeff[i]);L=L+temp}
    }
  }
  return(L)
}
#######################################################################

######### calculate 1st derivative (Speed Up the Optim) ###############
cal.DE<-function(theta){
  gk=rep(0,length(theta));Y=Data$Y;Delta=Data$Delta;T=Data$T;K=Data$K;dN=Data$dN;mm=length(t)
  X=Data$X; Z1=Data$Z1; Z2=Data$Z2;
  
  len1 = INN1+nknots1-1; len2 = INN2+nknots2-1
  id.xi1=1:len1; id.xi2=(len1+1):(len1+len2); id.alpha=(len1+len2+1):length(theta)
  
  Zmat = cbind(Z1, Z2)
  expcoeff=exp(Zmat %*% theta[id.alpha] + BB2[2:(n+1),] %*% theta[id.xi2])
  Lambdatheta=array(NA,c(mm,n,max(K)));for(m in 1:mm){for(i in 1:n){for(j in 1:K[i]){Lambdatheta[m,i,j]=t(theta[id.xi1])%*%dBB1[m,i,j,]}}}
  
  for(i in which(Delta==1)){for(j in 1:K[i]){
    gk[id.xi1]=gk[id.xi1]+2*(c(expcoeff[i]*Lambdatheta[rank[i],i,j])-dN[i,j])*expcoeff[i]*dBB1[rank[i],i,j,]
    gk[id.xi2]=gk[id.xi2]+2*(c(expcoeff[i]*Lambdatheta[rank[i],i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[rank[i],i,j])*BB2[(i+1),]
    ### 🟢 梯度按向量 c(Z1[i], Z2[i]) 更新两个 alpha 参数
    gk[id.alpha]=gk[id.alpha]+2*(c(expcoeff[i]*Lambdatheta[rank[i],i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[rank[i],i,j])*c(Z1[i], Z2[i])
  }}
  
  for(i in which(Delta==0)){
    if(sur[rank[i]]!=0&&(rank[i]!=mm)){
      temp1=rep(0,len1); temp2=rep(0,len2); temp3=rep(0,length(alpha));
      for(j in 1:K[i]){for(m in (rank[i]+1):mm){
        temp1=temp1+(sur[m-1]^Uexpcoeff[i]-sur[m]^Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*expcoeff[i]*dBB1[m,i,j,]
        temp2=temp2+(sur[m-1]^Uexpcoeff[i]-sur[m]^Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[m,i,j])*BB2[(i+1),]
        temp3=temp3+(sur[m-1]^Uexpcoeff[i]-sur[m]^Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[m,i,j])*c(Z1[i], Z2[i])
      }}
      if((sur[rank[i]]^Uexpcoeff[i])>0.1^8){
        temp1=temp1/(sur[rank[i]]^Uexpcoeff[i]); gk[id.xi1]=gk[id.xi1]+temp1
        temp2=temp2/(sur[rank[i]]^Uexpcoeff[i]); gk[id.xi2]=gk[id.xi2]+temp2
        temp3=temp3/(sur[rank[i]]^Uexpcoeff[i]); gk[id.alpha]=gk[id.alpha]+temp3
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
  ### 🟢 提取 X, Z1, Z2
  X=Data$X[id]; Z1=Data$Z1[id]; Z2=Data$Z2[id]; dN=Data$dN[id,];
  return(list(Y=Y,Delta=Delta,N=N,T=T,K=K,X=X,Z1=Z1,Z2=Z2,dN=dN,id=id))
}
#######################################################################

############# Calculate Bootstrap Loss Function #######################
bots.Loss=function(theta){
  L=0;Y=bots.Data$Y;Delta=bots.Data$Delta;T=bots.Data$T;K=bots.Data$K;dN=bots.Data$dN;mm=length(bots.t)
  X=bots.Data$X; Z1=bots.Data$Z1; Z2=bots.Data$Z2;

  len1 = INN1+nknots1-1; len2 = INN2+nknots2-1
  id.xi1=1:len1; id.xi2=(len1+1):(len1+len2); id.alpha=(len1+len2+1):length(theta)

  Zmat = cbind(Z1, Z2)
  expcoeff=exp(Zmat %*% theta[id.alpha] + bots.BB2[2:(n+1),] %*% theta[id.xi2])

  for(i in which(Delta==1)){for(j in 1:K[i]){L=L+(dN[i,j]-c(expcoeff[i]*(t(theta[id.xi1])%*%bots.dBB1[bots.rank[i],i,j,])))^2}}
  for(i in which(Delta==0)){
    if(bots.sur[bots.rank[i]]!=0&&(bots.rank[i]!=mm)){
      temp=0;
      for(j in 1:K[i]){
        for(m in (bots.rank[i]+1):mm){
          temp=temp+(bots.sur[m-1]^bots.Uexpcoeff[i]-bots.sur[m]^bots.Uexpcoeff[i])*(dN[i,j]-c(expcoeff[i]*t(theta[id.xi1])%*%bots.dBB1[m,i,j,]))^2
        }
      }
      if((bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])>0.1^8){temp=temp/(bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])}else{temp=0}
      L=L+temp;
    }
  }
  return(L)
}
#######################################################################

######### calculate 1st derivative (Speed Up the Optim) ###############
bots.DE<-function(theta){
  gk=rep(0,length(theta));Y=bots.Data$Y;Delta=bots.Data$Delta;T=bots.Data$T;K=bots.Data$K;dN=bots.Data$dN;mm=length(bots.t)
  X=bots.Data$X; Z1=bots.Data$Z1; Z2=bots.Data$Z2;

  len1 = INN1+nknots1-1; len2 = INN2+nknots2-1
  id.xi1=1:len1; id.xi2=(len1+1):(len1+len2); id.alpha=(len1+len2+1):length(theta)

  Zmat = cbind(Z1, Z2)
  expcoeff=exp(Zmat %*% theta[id.alpha] + bots.BB2[2:(n+1),] %*% theta[id.xi2])
  Lambdatheta=array(NA,c(mm,n,max(K)));for(m in 1:mm){for(i in 1:n){for(j in 1:K[i]){Lambdatheta[m,i,j]=t(theta[id.xi1])%*%bots.dBB1[m,i,j,]}}}

  for(i in which(Delta==1)){for(j in 1:K[i]){
    gk[id.xi1]=gk[id.xi1]+2*(c(expcoeff[i]*Lambdatheta[bots.rank[i],i,j])-dN[i,j])*expcoeff[i]*bots.dBB1[bots.rank[i],i,j,]
    gk[id.xi2]=gk[id.xi2]+2*(c(expcoeff[i]*Lambdatheta[bots.rank[i],i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[bots.rank[i],i,j])*bots.BB2[(i+1),]
    gk[id.alpha]=gk[id.alpha]+2*(c(expcoeff[i]*Lambdatheta[bots.rank[i],i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[bots.rank[i],i,j])*c(Z1[i], Z2[i])
  }}

  for(i in which(Delta==0)){
    if(bots.sur[bots.rank[i]]!=0&&(bots.rank[i]!=mm)){
      temp1=rep(0,len1); temp2=rep(0,len2); temp3=rep(0,length(alpha));
      for(j in 1:K[i]){
        for(m in (bots.rank[i]+1):mm){
          temp1=temp1+(bots.sur[m-1]^bots.Uexpcoeff[i]-bots.sur[m]^bots.Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*expcoeff[i]*bots.dBB1[m,i,j,]
          temp2=temp2+(bots.sur[m-1]^bots.Uexpcoeff[i]-bots.sur[m]^bots.Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[m,i,j])*bots.BB2[(i+1),]
          temp3=temp3+(bots.sur[m-1]^bots.Uexpcoeff[i]-bots.sur[m]^bots.Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[m,i,j])*c(Z1[i], Z2[i])
        }
      }
      if((bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])>0.1^8){
        gk[id.xi1]=gk[id.xi1]+temp1/(bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])
        gk[id.xi2]=gk[id.xi2]+temp2/(bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])
        gk[id.alpha]=gk[id.alpha]+temp3/(bots.sur[bots.rank[i]]^bots.Uexpcoeff[i])
      }
    }
  }
  return(gk)
}
#######################################################################

#######################################################################
library(survival)
library(splines2)

group = 1
case = 1
LAMBDA = LAMBDA1

n = 200
tau = 10
Obs = 6

INN1 = 4
INN2 = 4

nknots1 = round(2 * n^(1/5))
nknots2 = round(n^(1/5))

num.sim = 2
num.bots = 100

## 真实参数
alpha = c(-2, 0.5)
gammaU = c(0.5, -1, 0.5)

ka = find.ka(0.2)

knots1 = cal.knots(10000)

knots2 = seq(
  1 / (nknots2 + 1),
  nknots2 / (nknots2 + 1),
  1 / (nknots2 + 1)
)

#######################################################################
## 参数维数
#######################################################################

len1 = INN1 + nknots1 - 1
len2 = INN2 + nknots2 - 1
num_para = length(alpha)

total_len = len1 + len2 + num_para

theta.est = matrix(
  NA,
  nrow = num.sim,
  ncol = total_len
)

gammaU.est = matrix(
  NA,
  nrow = num.sim,
  ncol = length(gammaU)
)

bots.theta = array(
  NA,
  dim = c(num.sim, num.bots, total_len)
)

#######################################################################
## 约束矩阵
##
## 参数排列：
## theta = c(xi1, xi2, alpha)
##
## xi1：Lambda_0(s) 的样条系数
## xi2：beta_0(x) 的样条系数
## alpha：Z1、Z2 的参数系数
#######################################################################

## Lambda_0(s) 的样条系数约束：
##
## xi1[1] >= 0
## xi1[2] - xi1[1] >= 0
## ...
## xi1[len1] - xi1[len1 - 1] >= 0

D1 = cbind(
  diag(-1, len1)[, 2:len1, drop = FALSE],
  rep(0, len1)
) + diag(1, len1)

AA.lambda = cbind(
  D1,
  matrix(
    0,
    nrow = len1,
    ncol = len2 + num_para
  )
)

BB.lambda = rep(0, len1)

## alpha >= -5
AA.alpha.lower = cbind(
  matrix(
    0,
    nrow = num_para,
    ncol = len1 + len2
  ),
  diag(1, num_para)
)

## alpha <= 5
## 等价于 -alpha >= -5
AA.alpha.upper = cbind(
  matrix(
    0,
    nrow = num_para,
    ncol = len1 + len2
  ),
  diag(-1, num_para)
)

## 合并所有约束
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

## 检查约束矩阵维数
stopifnot(
  ncol(AA) == total_len,
  nrow(AA) == length(BB)
)

#######################################################################
## Monte Carlo simulation
#######################################################################

for (id.sim in 1:num.sim) {
  
  i.sim = id.sim + (core - 1) * num.sim
  
  set.seed(i.sim)
  
  Data = Gen.Data(n)
  
  #####################################################################
  ## 初始值
  ##
  ## 参数排列：
  ## c(xi1, xi2, alpha)
  #####################################################################
  
  Init = c(
    1:len1,
    rep(0, len2),
    rep(0, num_para)
  )
  
  stopifnot(length(Init) == total_len)
  
  ## 检查初始值是否满足约束
  if (any(drop(AA %*% Init - BB) < 0)) {
    stop("Initial value does not satisfy the constraints.")
  }
  
  #####################################################################
  ## Cox model
  #####################################################################
  
  cox = coxph(
    Surv(Y, Delta) ~ X + Z1 + Z2,
    data = Data
  )
  
  gammaUhat = cox$coefficients
  
  cox.survfit = survfit(cox)
  
  sur = cox.survfit$surv
  t = cox.survfit$time
  
  Uexpcoeff = exp(
    cbind(
      Data$X,
      Data$Z1,
      Data$Z2
    ) %*% gammaUhat
  )
  
  rank = rep(NA, n)
  
  for (i in 1:n) {
    
    if (length(which(t == Data$Y[i])) == 1) {
      
      rank[i] = which(t == Data$Y[i])
      
    } else {
      
      rank[i] = which(
        abs(t - Data$Y[i]) ==
          min(abs(t - Data$Y[i]))
      )
    }
  }
  
  ## CalBB() 返回的对象名称是 dBB1
  dBB1 = CalBB(Data, t)$dBB1
  
  #####################################################################
  ## beta_0(x) 的 B-spline 基函数
  #####################################################################
  
  BB2 = bSpline(
    c(0, Data$X, 1),
    df = INN2,
    knots = knots2,
    degree = INN2 - 1,
    intercept = FALSE
  )
  
  #####################################################################
  ## 主样本估计
  #####################################################################
  
  optim.est = constrOptim(
    theta = Init,
    f = Loss,
    grad = cal.DE,
    ui = AA,
    ci = BB,
    method = "BFGS"
  )
  
  ## 检查优化结果是否满足约束
  constraint.value = drop(
    AA %*% optim.est$par - BB
  )
  
  if (min(constraint.value) < -1e-8) {
    warning(
      paste0(
        "Simulation ",
        i.sim,
        ": constrained estimate violates constraints."
      )
    )
  }
  
  ## Bootstrap 使用主样本估计作为初始值
  Init = optim.est$par
  
  theta.est[id.sim, ] = optim.est$par
  gammaU.est[id.sim, ] = gammaUhat
  
  #####################################################################
  ## Bootstrap
  #####################################################################
  
  for (i.bots in 1:num.bots) {
    
    set.seed(
      i.sim * num.bots + i.bots + 1000
    )
    
    bots.Data = Genbots.Data(Data)
    
    ###############################################################
    ## Bootstrap Cox model
    ###############################################################
    
    bots.cox = coxph(
      Surv(Y, Delta) ~ X + Z1 + Z2,
      data = bots.Data
    )
    
    bots.gammaUhat = bots.cox$coefficients
    
    bots.survfit = survfit(bots.cox)
    
    bots.sur = bots.survfit$surv
    bots.t = bots.survfit$time
    
    bots.Uexpcoeff = exp(
      cbind(
        bots.Data$X,
        bots.Data$Z1,
        bots.Data$Z2
      ) %*% bots.gammaUhat
    )
    
    bots.rank = rep(NA, n)
    
    for (i in 1:n) {
      
      if (
        length(
          which(bots.t == bots.Data$Y[i])
        ) == 1
      ) {
        
        bots.rank[i] = which(
          bots.t == bots.Data$Y[i]
        )
        
      } else {
        
        bots.rank[i] = which(
          abs(bots.t - bots.Data$Y[i]) ==
            min(abs(bots.t - bots.Data$Y[i]))
        )
      }
    }
    
    ## CalBB() 返回的对象名称是 dBB1
    bots.dBB1 = CalBB(
      bots.Data,
      bots.t
    )$dBB1
    
    ###############################################################
    ## Bootstrap beta_0(x) 的 B-spline 基函数
    ###############################################################
    
    bots.BB2 = bSpline(
      c(0, bots.Data$X, 1),
      df = INN2,
      knots = knots2,
      degree = INN2 - 1,
      intercept = FALSE
    )
    
    ###############################################################
    ## Bootstrap constrained optimization
    ###############################################################
    
    bots.est = constrOptim(
      theta = Init,
      f = bots.Loss,
      grad = bots.DE,
      ui = AA,
      ci = BB,
      method = "BFGS"
    )
    
    bots.theta[id.sim, i.bots, ] = bots.est$par
  }
  
  #####################################################################
  ## 保存结果
  #####################################################################
  
  write.table(
    theta.est[id.sim, ],
    file = paste0(
      "n",
      n,
      "thetaest",
      i.sim,
      ".txt"
    ),
    sep = "\t",
    append = TRUE,
    col.names = FALSE,
    row.names = FALSE
  )
  
  write.table(
    bots.theta[id.sim, , ],
    file = paste0(
      "n",
      n,
      "botsTheta",
      i.sim,
      ".txt"
    ),
    sep = "\t",
    append = TRUE,
    col.names = FALSE,
    row.names = FALSE
  )
}