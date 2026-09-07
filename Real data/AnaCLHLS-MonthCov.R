################# Generate the Data Set ###############################
Gen.Data=function(CLHLS){
  Y=CLHLS[,10];Delta=CLHLS[,9];N=CLHLS[,11:16];T=CLHLS[,2:7];K=CLHLS[,20];Z=cbind(CLHLS[,c(17:18,8,26,23)])
  Z[which(Z[,1]==2),1]=0;Z[which(Z[,2]==2),2]=0; ## 1 M and U
  max.Z=matrix(apply(Z[,3:5],2,max),nr=n,nc=3,byrow=T)
  min.Z=matrix(apply(Z[,3:5],2,min),nr=n,nc=3,byrow=T)
  Z[,3:5]=(Z[,3:5]-min.Z)/(max.Z-min.Z)
  dN=matrix(NA,length(K),max(K));dN[,1]=N[,1];for(j in 2:max(K)){dN[,j]=N[,j]-N[,(j-1)]}
  return(list(Y=Y,Delta=Delta,N=N,T=T,K=K,Z=Z[,1:4],X=Z[,5],dN=dN))
}
#######################################################################

############# Calculate the knots of Spline ###########################
cal.knots=function(Data){
  Y=Data$Y;T=Data$T;K=Data$K
  ymt=matrix(NA,n,max(K));for(i in 1:n){for(j in 1:K[i]){ymt[i,j]=Y[i]-T[i,j]}}
  knots=c(quantile(ymt,seq(1/(nknots1+1),nknots1/(nknots1+1),1/(nknots1+1)),na.rm=TRUE,names=FALSE))
  knots[1]=knots[2]/2
  return(knots)
}
#######################################################################

#######################################################################
P_value <- function(cdf,x,paramet,side=0){
  n <-length(paramet)
  p <-switch(n+1,
             cdf(x),
             cdf(x,paramet),
             cdf(x,paramet[1],paramet[2]),
             cdf(x,paramet[1],paramet[2],paramet[3])
  )
  if(side <0) p
  else if(side >0) 1-p
  else
    if(p<1/2) 2*p
  else 2*(1-p)
}
#######################################################################

library(splines2);
load("CLHLS-MonthCovdata.RData")
n=nrow(CLHLS);Data=Gen.Data(CLHLS);tau=197;
INN1=4;INN2=4;nknots1=round(2*n^(1/5));nknots2=round(n^(1/5));knots1=cal.knots(Data);
ThetaEst=c(read.table(file="EstTheta.txt"))[[1]];Theta.05=Theta.95=rep(NA,length(ThetaEst));
num.bots=4;no_cores=25;ThetaBost=matrix(NA,num.bots*no_cores,length(ThetaEst));
for(i_core in 1:no_cores){
  ThetaBost[((i_core-1)*num.bots+1):(i_core*num.bots),]=as.matrix(read.table(file=paste0("BotsTheta_core",i_core,".txt")))
}
knots2=c(quantile(Data$X,seq(1/(nknots2+1),nknots2/(nknots2+1),1/(nknots2+1)),na.rm=TRUE,names=FALSE))
id.xi1=1:(INN1+nknots1-1);id.xi2=(INN1+nknots1):(INN1+nknots1+INN2+nknots2-2);
id.alpha=(length(ThetaEst)-4+1):length(ThetaEst)
tau=197;len=3000;s=0:len/len*tau;BB1=bSpline(s,df=INN1,knots=knots1,degree=(INN1-1),intercept=FALSE);
x=0:len/len;BB2=bSpline(x,df=INN2,knots=knots2,degree=(INN2-1),intercept=FALSE)
LambdaBost=BB1%*%t(ThetaBost[,id.xi1]);betaBost=BB2%*%t(ThetaBost[,id.xi2]);
id.Lambda=id.beta=rep(1,no_cores*num.bots)
id.Lambda[which(LambdaBost[1000,]<1)]=0
id.beta[which(abs(betaBost[2000,])<0.1)]=0
id.valid=which(id.Lambda*id.beta==1)
botsd=apply(ThetaBost[id.valid,],2,sd);
p_val=rep(NA,2);for(i in 1:4){p_val[i]=P_value(pnorm,ThetaEst[15+i],c(0,botsd[15+i]),side=0)};
##################### Calculate the Percentile ########################
#Theta.05=Theta.95=rep(NA,length(ThetaEst));Theta.05=ThetaEst-1.96*botsd;Theta.95=ThetaEst+1.96*botsd


# for(i in 1:length(ThetaEst)){Theta.05[i]=quantile(ThetaBost[id.valid,i],0.05);Theta.95[i]=quantile(ThetaBost[id.valid,i],0.95)}
# LambdaEst=BB1%*%ThetaEst[id.xi1];LAMBDA.05=BB1%*%Theta.05[id.xi1];LAMBDA.95=BB1%*%Theta.95[id.xi1];
# betaEst=BB2%*%ThetaEst[id.xi2];beta.05=BB2%*%Theta.05[id.xi2];beta.95=BB2%*%Theta.95[id.xi2];
# plot(s,LambdaEst,type="l",col="black",lty=1,ylab=expression(Lambda),xlab="s")
# plot(x,betaEst,type="l",col="black",lty=1,ylab=expression(beta),xlab="z")
# 
# plot(s,LambdaEst,type="l",col="black",lty=1,ylab=expression(Lambda),xlab="s",ylim=c(0,max(LAMBDA.95)))
# lines(s,LAMBDA.05,type="l",col="blue",lty=2);lines(s,LAMBDA.95,type="l",col="blue",lty=2)
# plot(x,betaEst,type="l",col="black",lty=1,ylab=expression(beta),xlab="z",ylim=c(min(beta.05),max(beta.95)))
# lines(x,beta.05,type="l",col="blue",lty=2);lines(x,beta.95,type="l",col="blue",lty=2)
# ThetaEst[id.alpha];botsd[id.alpha];p_val
for(i in 1:length(ThetaEst)){Theta.05[i]=quantile(ThetaBost[id.valid,i],0.05);Theta.95[i]=quantile(ThetaBost[id.valid,i],0.95)}
LambdaEst=BB1%*%ThetaEst[id.xi1];LAMBDA.05=BB1%*%Theta.05[id.xi1];LAMBDA.95=BB1%*%Theta.95[id.xi1];
betaEst=BB2%*%ThetaEst[id.xi2];beta.05=BB2%*%Theta.05[id.xi2];beta.95=BB2%*%Theta.95[id.xi2];
plot(s,LambdaEst,type="l",col="black",lty=1,ylab=expression(Lambda),xlab="s")
plot(x,betaEst,type="l",col="black",lty=1,ylab=expression(beta),xlab="z")

# plot(s,LambdaEst,type="l",col="black",lty=1,ylab=expression(Lambda),xlab="s",ylim=c(0,max(LAMBDA.95)))
# lines(s,LAMBDA.05,type="l",col="blue",lty=2);lines(s,LAMBDA.95,type="l",col="blue",lty=2)
# plot(x,betaEst,type="l",col="black",lty=1,ylab=expression(beta),xlab="z",ylim=c(min(beta.05),max(beta.95)))
# lines(x,beta.05,type="l",col="blue",lty=2);lines(x,beta.95,type="l",col="blue",lty=2)
ThetaEst[id.alpha];botsd[id.alpha];p_val

