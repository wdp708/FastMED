logdis=function(x) if(s>0) p/s*log(sum((abs(x))^s)) else sum(log(abs(x)))
eucl=function(x) sqrt(sum(x^2))

library(pdist)
library(randtoolbox)
library(conf.design)
library(MaxPro)
library(SLHD)
library(cubature)
lower1=-100
upper1=150
lower2=-25
upper2=10
lf=function(para)
{
  theta1=lower1+(upper1-lower1)*para[1]
  theta2=lower2+(upper2-lower2)*para[2]+.03*theta1^2-3
  val=-.5*(theta1^2/100+theta2^2)
  return(val)
}
p=2
n=max(primes(100+5*p))

MED=function(lf){
  logdis=function(x) if(s>0) p/s*log(sum((abs(x))^s)) else sum(log(abs(x)))
  D=Lattice(n,p)
  #D=maximinSLHD(t=1,n,p)$StandDesign
  Dk=D
  nc1=round(n/2)
  D1=maximinSLHD(t=1,nc1,p)$StandDesign
  ini=D
  lfini=apply(ini,1,lf)
  K=ceiling(4*sqrt(p))
  gamma=1/K
  Sigmak=var(D)
  s=0
  A=med(ini,lfini,Sigmak,gamma,s)
  D=A[,1:p]
  lfD=A[,p+1]
  for(k in 2:K)
  {
    gamma=(k-1)/(K-1)
    Sigmak=(k-1)/k*Sigmak
    #s=round(2*(1-exp(-gamma*(max(lfD)-min(lfD)))),1)
    s=round(2*(1-exp(-gamma*(quantile(lfD,.9)-quantile(lfD,.1)))))
    D1=D1[,sample(1:p)]
    new=lf.new=Dk=lfDk=NULL
    R=numeric(n)
    A=as.matrix(dist(D))+10*p*diag(n)
    r=apply(A,1,min)
    J=apply(A, 1, function(x) order(x, decreasing=F)[2])
    for(i in 1:n)
      R[i]=A[i,J[i]]

    for(j in 1:n)
    {
      d=apply((ini-rep(1,dim(ini)[1])%*%t(D[j,])),1,eucl)
      cl=(order(d,1:dim(ini)[1]))[1:n]
      dD=apply((D-rep(1,n)%*%t(D[j,])),1,eucl)
      M=rep(1,nc1)%*%t(D[j,])+R[j]/sqrt(p)*2*(D1-.5)
      nc=round(nc1/2)
      clD=(order(dD,1:n))[1:(nc+1)]
      ru=c(.5,.5)
      if(j==1){ru[2]=0.45}#runif(1,.25,.75)}
      clD2=clD[-1][1:round(nc/2)]
      M1=(ru[1]*D[clD[-1],]+(1-ru[1])*rep(1,nc)%*%t(D[j,]))
      M2=(1+ru[2])*rep(1,round(nc/2))%*%t(D[j,])-ru[2]*D[clD2,]
      M0=rbind(M1,M2)
      inij=ini[cl,]
      l=apply(inij,2,min)
      u=apply(inij,2,max)
      delta=(gamma-1/(K-1))*(u-l)/6
      l=l-delta
      u=u+delta
      for(i in 1:p)
        M[,i]=pmin(pmax(M[,i],l[i]),u[i])
      M=rbind(M,M0)
      M=augeucldist(M,inij)

      afit=LKfit(inij,lfini[cl])
      obj=function(x){
        apred=LKpredict(x,inij,afit$coef,afit$denom,afit$dbar)
        min.d=apred$mind
        if(min.d<.1/n) penalty=log(min.d) else penalty=0
        if(j==1){
          val=gamma*(apred$pred)+penalty
        }  else {
          val=min(.5*gamma*lfDk+.5*gamma*(apred$pred)+apply((Dk-rep(1,dim(Dk)[1])%*%t(x)),1,logdis))+penalty
        }
        return(val)
      }
      vind=which.max(apply(M,1,obj))
      v=M[vind,]
      new=rbind(new,v)
      lfv=lf(v)
      lf.new=c(lf.new,lfv)
      if(j==1) {if(lfv<lfD[j]) {Dk=rbind(Dk,D[j,]); lfDk=c(lfDk,lfD[j])} else {Dk=rbind(Dk,v); lfDk=c(lfDk,lfv)}} else{
        ov=min(.5*gamma*lfDk+.5*gamma*lfv+apply((Dk-rep(1,dim(Dk)[1])%*%t(v)),1,logdis))
        oD=min(.5*gamma*lfDk+.5*gamma*lfD[j]+apply((Dk-rep(1,dim(Dk)[1])%*%t(D[j,])),1,logdis))
        if(ov<oD) {Dk=rbind(Dk,D[j,]); lfDk=c(lfDk,lfD[j])} else {Dk=rbind(Dk,v); lfDk=c(lfDk,lfv)}
      }
      ini=rbind(ini,matrix(v,ncol=p))
      lfini=c(lfini,lfv)
    }
    #Find MED
    A=med(ini,lfini,Sigmak,gamma,s)
    D=A[,1:p]
    lfD=A[,p+1]
    Sigmak=var(D)
  }
  return(list(D=D,lfD=lfD,ini=ini,lfini=lfini))
}

med=function(CAND,lfCAND,Sigmak,gamma,s){
  logdis=function(x) if(s>0) p/s*log(sum((abs(x))^s)) else sum(log(abs(x)))
  N=dim(CAND)[1]
  ind=which.max(lfCAND)
  v=CAND[ind,]
  M=matrix(v,ncol=p)
  val=NULL
  D.lf=lfCAND[ind]
  eig=eigen(Sigmak)
  sqrtSinv=eig$vec%*%diag(1/sqrt(eig$val))%*%t(eig$vectors)
  for(i in 2:n)
  {
    val=cbind(val,.5*gamma*lfCAND+.5*gamma*lfCAND[ind]+apply((CAND-rep(1,N)%*%t(v))%*%sqrtSinv,1,logdis))
    ind=which.max(apply(val,1,min))
    v=CAND[ind,]
    M=rbind(M,v)
    D.lf=c(D.lf,lfCAND[ind])
  }
  return(as.matrix(cbind(M,D.lf)))
}

augeucldist=function(CAND,exist){
  N1=dim(CAND)[1]
  N2=dim(exist)[1]
  DIST=matrix(0,nrow=N1,ncol=N2)
  # DIST=as.matrix(pdist(CAND,exist))
  DIST = test_fastpdist(CAND, exist)
  disCAND=as.matrix(dist(CAND))
  M=NULL
  for(i in 1:round(N1/2)){
    ind=which.max(apply(DIST,1,min))
    v=CAND[ind,]
    M=rbind(M,v)
    CAND=CAND[-ind,]
    DIST=DIST[-ind,]
    DISTv=disCAND[-ind,ind]
    disCAND=disCAND[-ind,-ind]
    DIST=cbind(DIST,DISTv)
  }
  return(M)
}
LKfit=function(S,y){
  N=dim(S)[1]
  basis=function(d)
    1/(1+(d/dbar)^2)
  dist.S=as.matrix(dist(S))
  quant=max(1-5*p/N,.5)
  dbar=p*2*quant*quantile(apply(dist.S+10*p*diag(N),1,min),quant)
  Rinv=solve(basis(dist.S)+10^(-6)*diag(N))
  coef=Rinv%*%y
  denom=Rinv%*%rep(1,N)
  return(list(coef=coef,denom=denom,dbar=dbar,Rinv=Rinv))
}
LKpredict=function(x,S,coef,denom,dbar){
  N=dim(S)[1]
  basis=function(d)
    1/(1+(d/dbar)^2)
  d=apply(S-rep(1,N)%*%t(x),1,eucl)
  rvec=basis(d)
  val=drop(t(rvec)%*%coef/t(rvec)%*%denom)
  return(list(mind=min(d),pred=val))
}
