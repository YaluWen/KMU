## Correlated Case ##
#load data and test
library(Matrix)
library(psych)
Beta.Weights<-function(MAF,weights.beta){
  n<-length(MAF)
  weights<-rep(0,n)
  IDX_0<-which(MAF == 0)
  if(length(IDX_0) == n){
    stop("No polymorphic SNPs")
  } else if( length(IDX_0) == 0){
    weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
  } else {
    weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
  }
  return(weights)
}

gen.ker<-function(covx,kernel.index)
{ 
  Lin0<-function(x) tcrossprod(x)/max(1,ncol(x))
  Quad1<-function(x) (tcrossprod(x)+1)^2/max(1,ncol(x))
  IBS<-function(x)  1 - as.matrix(dist(x, method='manhattan') * .5 /max(1, ncol(x)) )  ## dist does scaling in the presence of missing values
  Gau <-function(x,p){return(exp(-1*as.matrix(dist(x)^2)/p))}
  
  #covx:X; 
  #kernel.index=c("Gau","Lin","Quad","IBS"): index of kernel function;
  n <- nrow(covx)
  p <- ncol(covx)
  if (kernel.index=="Lin") ker=Lin0(covx);
  if (kernel.index=="Quad") ker=Quad1(covx);
  if (kernel.index=="IBS") ker=IBS(covx);
  if (kernel.index=="Gau") ker=Gau(covx,ncol(covx));

  ker0 <- ker
  diag(ker0) <- rep(0,n)
  J <- matrix(1,n,n)
  ker.cen <- ker-J%*%ker0/(n-1)-ker0%*%J/(n-1)+J%*%ker0%*%J/n/(n-1)
  v1.cen <- psych::tr(ker.cen)/n
  return(list(ker.cen=ker.cen,v1.cen=v1.cen))
}

gen.ker.hadamard<-function(K1,K2){
  n=nrow(K1);
  ker=matrixcalc::hadamard.prod(K1,K2);
  ker <- as.matrix(Matrix::forceSymmetric(ker))
  ker0 <- ker
  diag(ker0) <- rep(0,n)
  J <- matrix(1,n,n)
  ker.cen <- ker-J%*%ker0/(n-1)-ker0%*%J/(n-1)+J%*%ker0%*%J/n/(n-1)
  v1.cen <- psych::tr(ker.cen)/n
  return(list(ker.cen=ker.cen,v1.cen=v1.cen))
}

 
UPvalueMulti<-function(KAll,Yk)
{
  Kallnew=list();
  Kallnew[[1]]=KAll[[1]];num=1;
  if(length(KAll)>1)
  {
    for(j in 2:length(KAll)) {
      add=TRUE;
      for(i in 1:length(Kallnew)) if(all.equal(round(Kallnew[[i]]),KAll[[j]])=="TRUE"){add=FALSE;break;}
      if(add){Kallnew[[num+1]]=KAll[[j]]; num=num+1;}
    }
  }
  KAll=Kallnew;
  n=nrow(KAll[[1]])
  VarCov=matrix(NA,length(KAll),length(KAll));
  UStats=rep(NA,length(KAll));
  for(i in 1:length(KAll))
  {
    for(j in i:length(KAll))
    {
      K1=KAll[[i]];K1=K1/tr(K1)*10;
      K2=KAll[[j]];K2=K2/tr(K2)*10;
      if(i==j)
      {
        U=matrixcalc::hadamard.prod(K1,Yk);
        diag(U)=0;
        UU1=sum(U)/2;
        UStats[i]=UU1;
      }
      A1=matrixcalc::hadamard.prod(K1,K2);
      A2=matrixcalc::hadamard.prod(Yk,Yk);
      A3=matrixcalc::hadamard.prod(A1,A2)
      diag(A3)=0
      FF1=sum(A3)/n/(n-1);
      var1=(1/2*n^2*FF1);
      VarCov[i,j]=VarCov[j,i]=var1;
      }
  }
  factors=max(sqrt(1/diag(VarCov)))/sqrt(1000);
  VarCov=VarCov*factors*factors;
  UStats=UStats*factors;

  pobs=findpobs(UStats,VarCov)
  pvalue=main_UStatsPvalueMain(NoKernel=length(KAll), pobs=min(pobs$pobs), CovVarAll=VarCov, SimNumber=15000) 
  if(is.na(pvalue)) pvalue=0;
  1-pvalue;
}




findpobs<-function(UStats,VarCov)
{
  n=3^length(UStats)
#  print(UStats)
#  print(VarCov)
  UU=NULL;sds=NULL;
  
  for(i in 1:n)
  {
    seltmp=rep(0,length(UStats));
    inew=i
    for(j in 1:length(UStats))
    {
     # print(inew)
      if( inew <= (n/3^j) ) seltmp[j]=-1
      if( (inew > (n/3^j)) & (inew <= (n/3^j)*2) ) seltmp[j]=0
      if( inew > ((n/3^j)*2) ) seltmp[j]=1
      inew=inew-n/(3^j)*(seltmp[j]+1)
    }
    
    U=sum(seltmp*UStats);
    sd1=0;
    for(k in 1:length(UStats))
      for(kk in k:length(UStats))
      {
        if(k==kk) sd1=sd1+VarCov[k,k]*seltmp[k]*seltmp[k]
        if(k!=kk) sd1=sd1+2*seltmp[k]*seltmp[kk]*VarCov[k,kk]
      }
    #print(seltmp);
    UU=rbind(UU,U)
    #print(sd1)
    sd1=round(sd1,14);
    sds=rbind(sds,sqrt(sd1))
    #sds=rbind(sds,(sd1))
  }
  exc= which(sds==0);
  UU=UU[-exc];sds=sds[-exc];
  results=list();
  results$UU=UU;
  results$sds=sds;
  results$pobs=2*(1-pnorm(abs(UU/sds)))
  results;
}

eigenx<-function(x)
{
  ttt=eigen(cov(x))
  ttt$values[round(ttt$values,10)==0]=0;
  if (sum(ttt$values > 0) > 1)
  {
    As = matrix(0, sum(round(ttt$values,10) > 0), sum(round(ttt$values,10) > 0))
    diag(As) = 1 / sqrt(ttt$values[round(ttt$values,10) > 0])
    InvsqrtCovx = ttt$vectors[, round(ttt$values,10) > 0] %*% As %*% t(ttt$vectors[, round(ttt$values,10) >0])
    xnew = t(InvsqrtCovx %*% t(x))
    x = xnew
  }
  x
}
preprocess<-function(x,pred)
{
  x1=list();
  for(j in 1:ncol(x)) x[is.na(x[,j]),j]=mean(x[,j],na.rm=T)
  if(pred=="Cont")
  {
    if(ncol(x)>1) x=eigenx(x)
    x1$x=x;
  }
  if(pred!="Cont")
  {
    x[,apply(x,2,mean)>1]=2-x[,apply(x,2,mean)>1];
    MAF=apply(x,2,mean)/2
    weights <- Beta.Weights(MAF, weights.beta=c(1,25))
    if(length(MAF)>1) xscale=apply(x,2,scale);
    if(length(MAF)==1) xscale=scale(x);
    rarex05scale=t(t(xscale[,MAF<0.05]) * (weights[MAF<0.05]));
    if(length(MAF)>1) xscale=eigenx(xscale);
    if(sum(MAF<0.05)>1) rarex05scale=eigenx(rarex05scale)
    if(sum(MAF<0.05)<2) rarex05scale=rarex05scale;
    x1$common=xscale;
    x1$rare=rarex05scale
  }
  x1
}

getK<-function(x,y,z,pred)
{
  y=as.matrix(y);
  z=as.matrix(z);
  z=model.matrix(~.,data.frame(z));
  res=resid(lm(y~z));
  res=apply(res,2,rank);
  res=apply(res,2,scale);
  Y=gen.ker(covx=res,kernel.index="Lin");
  Y$ker.cen=Y$ker.cen/tr(Y$ker.cen)*nrow(y);
  x1=preprocess(x,pred);
  if(pred=="Cont")
  {
    G1=gen.ker(covx = x1$x,kernel.index = "Lin")
    G1$ker.cen=G1$ker.cen/(G1$v1.cen)
    G2=gen.ker(covx = x1$x,kernel.index = "Quad")
    G2$ker.cen=G2$ker.cen/(G2$v1.cen)
    G3=gen.ker(covx = x1$x,kernel.index = "Gau")
    G3$ker.cen=G3$ker.cen/(G3$v1.cen)
  }
  if(pred!="Cont")
  {
    G1=gen.ker(covx = x1$common,kernel.index = "Lin")
    G1$ker.cen=G1$ker.cen/(G1$v1.cen)
    G2=gen.ker(covx = x1$common,kernel.index = "Quad")
    G2$ker.cen=G2$ker.cen/(G2$v1.cen)
    G3=gen.ker(covx = x1$rare,kernel.index = "Lin")
    G3$ker.cen=G3$ker.cen/(G3$v1.cen)
  }
  Ks=list();
  Ks$Y=Y$ker.cen;
  Ks$K=list();
  Ks$K$G1=G1$ker.cen;
  Ks$K$G2=G2$ker.cen;
  Ks$K$G3=G3$ker.cen;
  Ks
}

#' A U-statistics method for association test between a set of predictors and multiple outcomes 
#'
#' @param y A n * p matrix of outcomes with n subjects and p phenotypes. 
#' @param x A n * p matrix of predictors with n subjects and p predictors. 
#' @param z A matrix of demographic variables (e.g., age and gender).
#' @param pred Specify the types of predictors ("Cont": Continuous predictors, where linear, quadratic and gaussian kernels are used. "Cat": categorical predictors, where linear, weighted linear with beta weights, and quadratic kernels are used.) 
#' @return The P-value of the test.
#' @export
#'
#' @examples
#' data(Example)
#' KMUTest(xcont,y,z,pred="Cont")
#' KMUTest(xcat,y,z,pred="Cat")
KMUTest<-function(x,y,z,pred="Cont")
{
  K=getK(x,y,z,pred)
  ps=UPvalueMulti(K$K,K$Y)
  ps
}
