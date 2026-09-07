############### Functions #######################################
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
##############################################################

################RUNNING SIMULATION IN PARALLEL ###############

library(parallel)

workerFunc <- function(core){
  system(paste("Rscript --vanilla Ser_Mean.R", core))
}

combinX=list()
no_cores <- 100

for(i in 1:no_cores)
{
  combinX[[i]]<- i
}

c1 <-makeCluster(no_cores)
parLapply(c1,combinX,workerFunc)
stopCluster(c1)

##############################################################

n=200;tau=10;Obs=6;INN1=4;INN2=4;nknots1=round(2*n^(1/5));nknots2=round(n^(1/5));num.sim=2;num.bots=100;


alpha=c(-2, 0.5)
gammaU=c(0.5, -1, 0.5)

ka=find.ka(0.2); 
knots1=cal.knots(10000); 
knots2=seq(1/(nknots2+1),nknots2/(nknots2+1),1/(nknots2+1))



len1 = INN1 + nknots1 - 1
len2 = INN2 + nknots2 - 1
total_len = len1 + len2 + length(alpha)


theta.est=matrix(NA, num.sim*no_cores, total_len);


bots.theta=array(NA,c(num.sim*no_cores, num.bots, total_len));
bots.sd=matrix(NA,num.sim*no_cores, total_len);

print("开始读取数据...")


  
for(i_core in 1:(num.sim*no_cores))
{
    #
  theta.est[i_core,] = as.matrix(read.table(file=paste0("n",n,"thetaest",i_core,".txt")))
    
  # #
  bots.data_temp = as.matrix(read.table(file=paste0("n",n,"botsTheta",i_core,".txt")))
  bots.theta[i_core,,] = matrix(bots.data_temp, nrow=num.bots, ncol=total_len)

  bots.sd[i_core,] = apply(bots.theta[i_core,,], 2, sd)
}
save.image(paste0("n",n,".RData"))

print("运行结束！")
