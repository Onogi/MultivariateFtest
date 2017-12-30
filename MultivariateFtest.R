MultivariateFtest<-function(X,R,G,Y,A,Z=NULL,C=NULL,Ind=FALSE,MAFthr=NULL){  

  #Onogi (2018)
  
  #X is k x m matrix of marker genotypes (coded as 0, 1, 2)
  #Markers are tested one by one
  #R is d x d residual variance matrix
  #G is d x d genetic variance matrix
  #Y is (n x d) x 1 matrix or vector of phenotypes.
  #The first n elements are for variate 1, the second n elements are for variate 2, and so on.
  #A is k x k relationship matrix
  #Z is n x k design matrix connecting Y with X and A.
  #When Z is null, a diagonal matrix is used.
  #C is n x l matrix of covariates. These covariates are not tested.
  #C is shared among variates.
  #When C is null, the intercept is added automatically.
  #Ind is logical, test effects on all dimensions simultaneously (i.e., multivariate F-test) (FALSE) or independently (TRUE)
  #MAFthr is MAF threshold
  
  #NA in Y is allowed. But NA in X and C is not allowed.
  #Univariate F-test also can be peformed by setting d = 1.
  
  k<-nrow(X)
  stopifnot(nrow(A)==k&ncol(A)==k)
  if(is.vector(Y)) Y<-matrix(Y,nc=1)
  if(is.vector(R)&length(R)==1) R<-matrix(R,1,1)
  if(is.vector(G)&length(G)==1) G<-matrix(G,1,1)
  d<-nrow(R)
  stopifnot(ncol(R)==d&nrow(G)==d&ncol(G)==d)
  stopifnot(nrow(Y)%%d==0)
  n<-nrow(Y)/d
  m<-ncol(X)
  if(is.null(Z)){
    stopifnot(n==k)
    Z<-diag(n)
  }else{
    stopifnot(nrow(Z)==n&ncol(Z)==k)
  }
  if(!is.null(C)){
    stopifnot(nrow(C)==n)
    l<-ncol(C)
    dl<-d*l
  }

  if(is.null(MAFthr)) MAFthr<-0
  MAF<-numeric(d)

  cat("Processing data\n")
  use<-!is.na(Y)
  Y2<-Y[use,,drop=FALSE]
  nd<-nrow(Y2)
  #V<-G%x%(Z%*%A%*%t(Z))+R%x%diag(n)
  #V<-V[use,use]
  #sometimes V becomes too large
  
  use.d<-matrix(NA,nr=n,nc=d)
  for(trait in 1:d){
    use.d[,trait]<-!is.na(Y[((trait-1)*n+1):(trait*n),])
  }
  V<-matrix(0,nr=nd,nc=nd)
  for(trait in 1:d){
    if(trait==1) which.row<-1:sum(use.d[,1:trait]) else which.row<-(sum(use.d[,1:(trait-1)])+1):(sum(use.d[,1:trait]))
    for(j in trait:d){
      if(j==1) which.col<-1:sum(use.d[,1:j]) else which.col<-(sum(use.d[,1:(j-1)])+1):(sum(use.d[,1:j]))
      V[which.row,which.col]<-G[trait,j]*(Z[use.d[,trait],]%*%A%*%t(Z[use.d[,j],]))+R[trait,j]*diag(n)[use.d[,trait],use.d[,j]]
      if(trait!=j) V[which.col,which.row]<-t(V[which.row,which.col])
    }
  }
  
  cat("Calculate iV\n")
  count<-0
  repeat{
    iV<-try(solve(V),silent=T) 
    if(class(iV)!="try-error"){
      break
    }else{
      if(count==5){
        stop("V is singular")
      }
      diag(V)<-diag(V)+1e-3
      count<-count+1
    }
  }
  if(count>0) cat(count*1e-3,"was added to the diagonal of V\n")

  iVY<-iV%*%Y2
  YiVY<-t(Y2)%*%iVY
  
  if(is.null(C)){
    temp<-diag(d)
    temp<-temp[rep(1:d,each=n),]
    Xd<-cbind(temp,matrix(0,nr=n*d,nc=d))
    Xd<-Xd[use,]
    p<-d+d#number of fixed effects
    pr<-d#number of fixed effects not to be tested
  }else{
    Xd<-cbind(matrix(0,nr=n*d,nc=dl),matrix(0,nr=n*d,nc=d))
    for(trait in 1:d){
      for(covariate in 1:l){
        Xd[((trait-1)*n+1):(trait*n),(trait-1)*l+covariate]<-C[,covariate]
      }
    }
    Xd<-Xd[use,]
    
    Whichdim.col<-rep(1:d,each=l)
    Whichdim.row<-rep(1:d,each=n)
    Whichdim.row<-Whichdim.row[use]
    
    Remove<-colSums(abs(Xd[,1:dl]))==0#remove fixed effects with a single level
    if(any(Remove)){
      Xd<-Xd[,c(!Remove,rep(TRUE,d))]
      dl<-d*l-sum(Remove)
      Whichdim.col<-Whichdim.col[!Remove]
    }
    for(trait in 1:d){
      v<-which(Whichdim.col==trait)
      Intercept<-min(v)
      Xd[Whichdim.row==trait,Intercept]<-1
    }
    
    p<-dl+d#number of fixed effects
    pr<-dl#number of fixed effects not to be tested
  }
  cat(dim(Xd)," p",p," pr",pr,"\n")
  
  if(Ind) {
    K<-matrix(0,nr=1,nc=p)
    f<-1#number of fixed effects to be tested
    Pvalue<-matrix(NA,nr=m,nc=d)
  } else {
    K<-cbind(matrix(0,d,p-d),diag(d))
    f<-d#number of fixed effects to be tested
    Pvalue<-matrix(NA,nr=m,nc=1)
  }
  Beta<-matrix(NA,nr=m,nc=d)
  
  #calculate products between covariates and iV
  XiVX<-matrix(0,nr=p,nc=p)
  XiVX.int<-t(Xd[,1:pr,drop=F])%*%iV
  XiVX[1:pr,1:pr]<-XiVX.int%*%Xd[,1:pr,drop=F]
  
  cat("Start marker test\n")
  print(proc.time())
  for(marker in 1:m){
    if(marker%%100==0)cat(marker,"\n")
    
    for(trait in 1:d){
      if(trait==1) which.row<-1:sum(use.d[,1:trait]) else which.row<-(sum(use.d[,1:(trait-1)])+1):(sum(use.d[,1:trait]))
      Xd[which.row,pr+trait]<-Z[use.d[,trait],]%*%X[,marker,drop=F]
    } 
    TestedX<-Xd

    #check MAF
    for(trait in 1:d){
      target<-use.d
      target[,-trait]<-FALSE
      target<-as.vector(target)[as.vector(use.d)]
      MAF[trait]<-sum(TestedX[as.vector(target),pr+trait])/(2*sum(target))
    }
    MAF[MAF>0.5]<-1-MAF[MAF>0.5]
    
    if(!any(MAF<MAFthr)){
      XiVX[1:pr,(pr+1):p]<-XiVX.int%*%TestedX[,(pr+1):p]
      XiVX[(pr+1):p,1:pr]<-t(XiVX[1:pr,(pr+1):p])
      XiVX[(pr+1):p,(pr+1):p]<-t(TestedX[,(pr+1):p])%*%iV%*%TestedX[,(pr+1):p]
      iXiVX<-try(solve(XiVX))
      
      if(class(iXiVX)!="try-error"){
        XiVY<-t(TestedX)%*%iVY
        B<-iXiVX%*%XiVY
        SigmaE<-(YiVY-t(B)%*%XiVY)/(nd-p)
        Beta[marker,]<-as.vector(B)[-c(1:pr)]
        
        if(Ind){
          for(trait in 1:d){
            K2<-K
            K2[,pr+trait]<-1
            KB<-K2%*%B
            KXV<-K2%*%iXiVX%*%t(K2)
            Fvalue<-t(KB)%*%solve(KXV)%*%KB/f/SigmaE
            q<-(nd-p)/(nd-p+f*Fvalue)
            Pvalue[marker,trait]<--log10(pbeta(q,(nd-p)/2,f/2))       
          }
        }else{
          KB<-K%*%B
          KXV<-K%*%iXiVX%*%t(K)
          Fvalue<-t(KB)%*%solve(KXV)%*%KB/f/SigmaE
          q<-(nd-p)/(nd-p+f*Fvalue)
          Pvalue[marker,]<--log10(pbeta(q,(nd-p)/2,f/2))
        }      
      }
    }
  }
 
  cat("Finished\n")
  print(proc.time())
  rm(KXV,XiVY,XiVX,iXiVX,TestedX,Xd,iVY,YiVY,V,iV);gc();gc();gc();gc()
  cbind(Pvalue,Beta)
}