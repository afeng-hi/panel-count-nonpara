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
  L=0;Y=Data$Y;Delta=Data$Delta;T=Data$T;K=Data$K;dN=Data$dN;mm=length(t);Z=Data$Z;X=Data$X;
  id.xi1=1:(INN1+nknots1-1);id.xi2=(INN1+nknots1):(INN1+nknots1+INN2+nknots2-2);id.alpha=(length(theta)-p+1):length(theta)
  expcoeff=exp(as.matrix(Z)%*%theta[id.alpha]+BB2[2:(n+1),]%*%theta[id.xi2])
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
  gk=rep(0,length(theta));Y=Data$Y;Delta=Data$Delta;T=Data$T;K=Data$K;dN=Data$dN;mm=length(t);Z=Data$Z;X=Data$X;
  id.xi1=1:(INN1+nknots1-1);id.xi2=(INN1+nknots1):(INN1+nknots1+INN2+nknots2-2);id.alpha=(length(theta)-p+1):length(theta)
  expcoeff=exp(as.matrix(Z)%*%theta[id.alpha]+BB2[2:(n+1),]%*%theta[id.xi2])
  Lambdatheta=array(NA,c(mm,n,max(K)));
  for(m in 1:mm){for(i in 1:n){for(j in 1:K[i]){Lambdatheta[m,i,j]=t(theta[id.xi1])%*%dBB1[m,i,j,]}}}
  for(i in which(Delta==1)){for(j in 1:K[i]){
    gk[id.xi1]=gk[id.xi1]+2*(c(expcoeff[i]*Lambdatheta[rank[i],i,j])-dN[i,j])*expcoeff[i]*dBB1[rank[i],i,j,]
    gk[id.xi2]=gk[id.xi2]+2*(c(expcoeff[i]*Lambdatheta[rank[i],i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[rank[i],i,j])*BB2[(i+1),]
    gk[id.alpha]=gk[id.alpha]+2*(c(expcoeff[i]*Lambdatheta[rank[i],i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[rank[i],i,j])*Z[i,]
  }}
  for(i in which(Delta==0)){
    if(sur[rank[i]]!=0&&(rank[i]!=mm)){
      temp1=rep(0,INN1+nknots1-1);temp2=rep(0,INN2+nknots2-1);temp3=rep(0,length(Data$Z[1,]));
      for(j in 1:K[i]){for(m in (rank[i]+1):mm){
        temp1=temp1+(sur[m-1]^Uexpcoeff[i]-sur[m]^Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*expcoeff[i]*dBB1[m,i,j,]
        temp2=temp2+(sur[m-1]^Uexpcoeff[i]-sur[m]^Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[m,i,j])*BB2[(i+1),]
        temp3=temp3+(sur[m-1]^Uexpcoeff[i]-sur[m]^Uexpcoeff[i])*2*(c(expcoeff[i]*Lambdatheta[m,i,j])-dN[i,j])*c(expcoeff[i]*Lambdatheta[m,i,j])*Z[i,]
      }}
      if((sur[rank[i]]^Uexpcoeff[i])>0.1^8){
        temp1=temp1/(sur[rank[i]]^Uexpcoeff[i]);gk[id.xi1]=gk[id.xi1]+temp1
        temp2=temp2/(sur[rank[i]]^Uexpcoeff[i]);gk[id.xi2]=gk[id.xi2]+temp2
        temp3=temp3/(sur[rank[i]]^Uexpcoeff[i]);gk[id.alpha]=gk[id.alpha]+temp3
      }
    }
  }
  return(gk)
}
#######################################################################

library(survival);
library(splines2);
load("CLHLS-MonthCovdata.RData")
n=nrow(CLHLS);Data=Gen.Data(CLHLS);p=ncol(Data$Z);q=1;tau=197;
INN1=4;INN2=4;nknots1=round(2*n^(1/5));nknots2=round(n^(1/5));knots1=cal.knots(Data);
knots2=c(quantile(Data$X,seq(1/(nknots2+1),nknots2/(nknots2+1),1/(nknots2+1)),na.rm=TRUE,names=FALSE))
Init=c(1:(INN1+nknots1-1),rep(0,p+q*(INN2+nknots2-1)))
cox=coxph(Surv(Y, Delta) ~ X + Z, Data);
gammaUhat=cox$coefficients;sur=survfit(cox)$surv;t=survfit(cox)$time
Uexpcoeff=exp(cbind(Data$X,Data$Z)%*%gammaUhat);rank=rep(NA,n);
for(i in 1:n){if(length(which(t==Data$Y[i]))==1){rank[i]=which(t==Data$Y[i])}else{rank[i]=which(abs(t-Data$Y[i])==min(abs(t-Data$Y[i])))}}
dBB1=CalBB(Data,t)$dBB;BB2=bSpline(c(0,Data$X,1),df=INN2,knots=knots2,degree=(INN2-1),intercept=FALSE)
AA=cbind(cbind(diag(-1,INN1+nknots1-1)[,2:(INN1+nknots1-1)],rep(0,(INN1+nknots1-1)))+diag(1,INN1+nknots1-1),matrix(0,INN1+nknots1-1,q*(INN2+nknots2-1)+p));
BB=rep(0,INN1+nknots1-1);
optim.est=constrOptim(Init,Loss,cal.DE,AA,BB,method = "BFGS")
EstTheta=optim.est$par
write.table(EstTheta,file=paste0("EstTheta.txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
save.image("CLHLS-Estimation.RData")

# Init=c(1:(INN1+nknots1-1),rep(0.1,p+q*(INN2+nknots2-1)))
# library(nlme)
# Loss(Init)
# cal.DE(Init)
# fdHess(Init,Loss)
