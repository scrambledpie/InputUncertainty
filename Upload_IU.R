##################################################################################################
##################################################################################################
library(MASS)
library(FastGP)
outer=function (X, Y,Z=1) {
  
  dX <- length(X)
  dY <- length(Y)
  # Y <- rep(Y, rep.int(dX, dY))
  # X <- rep(X, times = dY)
  # robj <- X-Y
  robj<- rep(Y,rep.int(dX,dY)) - rep(X,dY)
  dim(robj) <- c(dX, dY)
  robj
  # sapply(1:length(Y),function(i)Y[i]-X,USE.NAMES = FALSE)
}

##################################################################################################
##################################################################################################

# Input Uncertainty Parameter Distributions
# 
# DensGeneral=list(
#   rdens = function(n) return n random numbers,
#   ddens = function(x) return the probability density at x
#   rlhc  = function(n) return n random number according to latin hypercube
#   name  = the name of the distro for plots and saving results
#   Mean, Mode of the distro for use with varaib=tions of the EGO
# )

DensWedge=list(
  rdens = function(n){ 100*sqrt(runif(n))}, 
  ddens = function(x){  x*2e-4},
  rlhc  = function(n){  100*sqrt((1:n-runif(n))*(1/n))},
  name  = "Triangular",
  Mean = 200/3, Mode = 100 
)
  
DensUniform=list(
  rdens = function(n){ runif(n)*100},
  ddens = function(x){ rep(0.01,length(x))},
  rlhc  = function(n){  100*((1:n-runif(n))/n)},
  name  = "Uniform",
  Mean = 50, Mode=50
)
  
TrueHpars   = c(10,10,1)
TrueNoiseSD = 0.1

CreateTestFunction = function(Hpars=TrueHpars,seed=3){
  # returns a smooth random function on the domain [0,100]x[0,100]
  # the random function is the GP posterior mean of randomly generated data
  
  lX    =  Hpars[1];  ilX22=-0.5/lX^2
  lA    =  Hpars[2];  ilA22=-0.5/lA^2
  Sig0  =  Hpars[3]
  Sig2  =  Sig0*Sig0
  X0=seq(0,100,len=16);  nX=length(X0)
  A0=seq(0,100,len=15);  nA=length(A0)
  XA0=cbind(X0,rep(A0,each=nX))
  Sigma = Sig2 * exp( ilX22 * outer(XA0[,1],XA0[,1],"-")^2 + ilA22 * outer(XA0[,2],XA0[,2],"-")^2)
  diag(Sigma)=Sig2*1.001
  set.seed(seed)
  Y0 = mvrnorm(1,rep(0,nX*nA),Sigma)

  iKY0 = rcppeigen_invert_matrix(Sigma)%*%Y0
  
  
  GeneratedTestFun=function(xa,NoiseSD = TrueNoiseSD){
    if(length(xa)==2)xa=matrix(xa,1)
    ks=Sig2 * exp(ilX22 * (outer(xa[,1],XA0[,1],"-")^2 +  outer(xa[,2],XA0[,2],"-")^2))
    out=ks%*%iKY0
    out[1:length(out)] + rnorm(length(out),0,NoiseSD)
  }
  
  return(GeneratedTestFun)
}

RunBenchMark = function(TestFun=CreateTestFunction(),Distro=DensUniform, StartN=10, EndN=100, method = 2, verbose=T,seed=1){
  if(verbose)cat("starting Benchmark")
  
  ##################################################################################################
  ##################################################################################################
  # Gaussian Process Fitting Functions
    
  lX      = 10;  ilX22 = -0.5/lX^2
  lA      = 10;  ilA22 = -0.5/lA^2
  Sig0    = 1;   Sig2  = Sig0*Sig0
  noiseVar= 0.01
  
  invK=matrix(0,1,1)
  invKy=rep(0,1)
  
  UpdateGP=function(EstimateHpars = FALSE){
    if(EstimateHpars){
      
      LhoodOptimizer = function(xd=Data[,1:2],yd=Data[,3]){
        warning("Hyper-parameter estimation not implemented in this code, returning true values")
        return(c(TrueHpars,TrueNoiseSD^2))
      }
      MLHpars = LhoodOptimizer()
      
      lX       <<- MLHpars[1];  ilX22  <<-  -0.5/lX^2
      lA       <<- MLHpars[2];  ilA22  <<-  -0.5/lA^2
      Sig2     <<- MLHpars[3];  Sig0   <<-  sqrt(Sig2)
      noiseVar <<- MLHaprs[4]
    }
    
    invK       <<- Sig2*exp( ilX22*outer(Data[,1],Data[,1],"-")^2 + ilA22*outer(Data[,2],Data[,2],"-")^2 )
    diag(invK) <<- Sig2+noiseVar
    invK       <<- rcppeigen_invert_matrix(invK)
    invKy      <<- invK%*%Data[,3]
  }
  MU1=function(xa){
    Ks=Sig2*exp( ilX22 * (xa[1]- Data[,1])^2 +  ilA22*(xa[2]-Data[,2])^2) 
    sum(Ks*invKy)
  }
  MU=function(xai){
    if(is.null(dim(xai)))xai=matrix(xai,1)
    Ks=Sig2*exp( ilX22 * outer(xai[,1], Data[,1],"-")^2 +  ilA22*outer(xai[,2],Data[,2],"-")^2 ) 
    out=(Ks)%*%(invKy)
    out[1:length(out)]
  }
  COV1=function(xa1){
    Ks1=Sig2*exp( ilX22*(Data[,1]-xa1[1])^2 + ilA22*(Data[,2]-xa1[2])^2 )
    Kss=Sig2
    (Kss - t(Ks1) %*% invK %*% Ks1)[1]
  }
  COV=function(xa1,xa2){
    
    if(is.null(dim(xa1)))xa1=matrix(xa1,1)
    if(is.null(dim(xa2)))xa2=matrix(xa2,1)
    
    Ks1=Sig2*exp( ilX22*outer(Data[,1],xa1[,1],"-")^2 + ilA22*outer(Data[,2],xa1[,2],"-")^2)
    Ks2=Sig2*exp( ilX22*outer(Data[,1],xa2[,1],"-")^2 + ilA22*outer(Data[,2],xa2[,2],"-")^2)
    Kss=Sig2*exp( ilX22*outer(xa1[,1] ,xa2[,1],"-")^2 + ilA22*outer(xa1[,2] ,xa2[,2],"-")^2)
    
    Kss - t(Ks1) %*% invK %*% Ks2
    
  }
  SD=function(xa){
    if(!is.null(dim(xa))){
      A=apply(xa,1,function(xai)abs(COV(xai,xai)))
      A=A[1:length(A)]
    }else{
      A=abs(COV1(xa))[1]
    }
    A/sqrt(A+noiseVar)
  }
  
  ##################################################################################################
  ##################################################################################################
  # Sample allocations methods
  
  # method 1
  UNIdesign=function(N0){
    AA=100*(1:N0-runif(N0))/N0
    xa=cbind(Distro$rlhc(N0),sample(AA))
    y=TestFun(xa)
    Data<<-cbind(xa,y)
  }
  
  # method 2
  KG_MC_Input=function(Nx=nrow(Data)){
    KGCB = function(a,b){
      n=length(a)
      O=order(b)
      b=b[O]
      a=a[O]
      #sort (a,b) by ascending b
      
      A=1
      C=-Inf
      
      while(tail(A,1)<n){
        s=tail(A,1)
        si=(s+1):n
        Ci=-(a[s]-a[si])/(b[s]-b[si])
        bestsi=which.min(Ci)
        C=c(C,Ci[bestsi])
        A=c(A,si[bestsi])
      }
      C=c(C,Inf)
      
      pC=pnorm(C)
      dC=dnorm(C)
      
      sum(  a[A]*(pC[-1]-pC[-length(pC)]) - b[A]*(dC[-1]-dC[-length(pC)])  )-max(a)
      
    }
    
    Nx   = ceiling(Nx)
    
    Na   = nrow(Data)
    a0   = Data[,2]
    x0   = Distro$rlhc(Nx)
    
    xa0  = cbind(x0,rep(a0,each=Nx))
    MM   = apply(matrix(MU(xa0),nrow=Nx),2,sum)
    Ks0  = Sig2 * exp( ilX22*(outer(Data[,1],xa0[,1],"-")^2 + outer(Data[,2],xa0[,2],"-")^2 ) )
    iKKs0= invK%*%Ks0
    
    bestEVI = -10
    bestxa  = 0
    
    KG_IU=function(xa){
      if(abs(xa[1]-50)>50 | abs(xa[2]-50)>50 )return(-1000000)
      
      iysd = 1/sqrt(COV1(xa)+noiseVar)
      
      Kss  = Sig2 * exp( ilX22*( (xa0[,1]-xa[1])^2 + (xa0[,2]-xa[2])^2 ) )
      Ks1  = Sig2 * exp( ilX22*( (Data[,1] -xa[1])^2 + (Data[,2] -xa[2])^2 ) )
      sigt = matrix( Kss - Ks1 %*% iKKs0,nrow=Nx)
      
      newa = cbind(x0,xa[2])
      MM   = c(sum(MU(newa)),MM)
      sigt = cbind( COV(newa,xa)[1:Nx], sigt  )
      
      D    = Distro$ddens(xa[1])*100
      newx = cbind(xa[1],c(xa[2],a0))
      MM   = MU(newx)*D  + MM
      sigt = rbind( COV(newx,xa)[1:(Na+1)]*D , sigt  )
      
      sigt = sapply(1:ncol(sigt),function(i)sum(sigt[,i])) * iysd
      
      out  = KGCB(MM,sigt)
      if(out>bestEVI){
        bestEVI <<-out
        bestxa  <<-xa
      }
      return(out)
    }

    Nas = 15
    sX  = cbind( Distro$rlhc(Nas)[sample(1:Nas)],  100*(0:(Nas-1)+runif(Nas))/Nas)
    OO  = sapply(1:Nas,function(i)optim(sX[i,],KG_IU,method = "Nelder-Mead",control=list(maxit=100, fnscale=-1)))
    
    return(bestxa)

  }
  
  # method 3
  EGO_MC_Input=function(Nx=nrow(Data),Ns=15){
    
    bestEVI = -10
    bestxa  = c(-10,-10)
    Nx      = ceiling(Nx)
    Xi      = Distro$rlhc(Nx)
    
    obj     = function(a)sum(MU(cbind(Xi,a)))
    topA    = which.max(sapply(1:99,obj))
    out1    = optimise(obj,interval = c(topA-1,topA+1),maximum = T)
    topA    = out1$maximum
    topmu   = out1$objective
    
    EGO     = function(xa){
      if(abs(xa[1]-50)>50 | abs(xa[2]-50)>50 )return(-1000000)
      
      dXi     =  Distro$ddens(xa[1])*100
      topmi   =  topmu      + dXi*MU1(c(xa[1],topA))
      ami     =  obj(xa[2]) + dXi*MU1(xa)
      
      sigt1   = COV1(xa)
      sigt    = (sum(COV(xa,cbind(Xi,xa[2]))) + sigt1*dXi ) / sqrt(sigt1+noiseVar)
      
      z   = (topmi-ami)/sigt
      
      out = (ami-topmi)*pnorm(-z) + sigt*dnorm(z)
      
      if(out>bestEVI){
        bestEVI<<-out
        bestxa<<-xa
      }
      out
    }
    
    XAs     = cbind(sample(Distro$rlhc(Ns)), (1:Ns-runif(Ns))*100/Ns )
    A       = sapply(1:Ns,function(i)optim(par = XAs[i,],fn = EGO,method = "Nelder-Mead",control = list(maxit=100, fnscale=-1)))
    return(bestxa)
  
  }
  
  # methods 4, 5
  EGO_Single_Input=function(Xi=Distro$Mean){
    bestEVI=-10
    bestxa=c(-10,-10)
    
    
    obj=function(a)MU(c(Xi,a))
    topA=which.max(sapply(1:99,obj))
    out1=optimise(obj,interval = c(topA-1,topA+1),maximum = T)
    topA=out1$maximum
    topmu=out1$objective
    
    EGO=function(a){
      if(abs(a-50)>50)return(-1000000)
      
      ami     =  MU1(c(Xi,a))
      
      sigt1 = COV1(c(Xi,a))
      sigt  = sigt1  / sqrt(sigt1+noiseVar)
      
      z   = (topmu-ami)/sigt
      
      out = (ami-topmu)*pnorm(-z) + sigt*dnorm(z)

      
      if(out>bestEVI){
        bestEVI<<-out
        bestxa<<-c(Xi,a)
      }
      out
    }
    
    Ad=c(0,sort(Data[,2]),100)
    As=sort(runif(3+3*nrow(Data),Ad[-length(Ad)],Ad[-1]))
    Asi=which.max(sapply(As,EGO))+1
    As=c(0,As,100)
    
    A=optimise(EGO,interval = c(As[Asi-1],As[Asi+1]),maximum = T)
    
    
    return(bestxa)
  }
  
  ##################################################################################################
  ##################################################################################################
  # Quality measure/Opportunity Cost
  
  Test_Inputs       =  Distro$rlhc(1000)
  True_Quality      =  function(a)mean(TestFun(cbind(Test_Inputs,a),NoiseSD = 0))
  
  Recm_Inputs       =  Distro$rlhc(1000)
  Recommend_A_value =  function( Method = method ){
    
    if(Method==0){ # For finding the optimal A value of the true test function
      Pred_Quality = True_Quality
      
    }else if(Method==4){
      Pred_Quality = function(a) MU(c(Distro$Mean,a))
      
    }else if(Method==5){
      Pred_Quality = function(a) MU(c(Distro$Mode,a))
      
    }else{
      Pred_Quality = function(a) sum(MU(cbind(Recm_Inputs,a)))
      
    }
    
    bestA = which.max(sapply(1:99,Pred_Quality))
    out   = optimize(Pred_Quality,interval = c(bestA-1,bestA+1), maximum = T)
    bestA = out$maximum
    
    return(bestA)
    
  }
  
  TrueQ             =  True_Quality(a = Recommend_A_value(Method=0) )
  Get_Quality       =  function()True_Quality(a = Recommend_A_value(Method = method))
  
  
  ##################################################################################################
  ##################################################################################################
  # Initialize the benchmark and execute sampling from StartN samples to EndN samples
  
  set.seed(seed)
  UNIdesign(StartN)
  UpdateGP()
  Cost=data.frame(N=StartN, OC = TrueQ - Get_Quality())
  
  if(verbose)cat("Uni...",StartN,"   ", tail(Cost,1)[[2]],"\n")
  
  while(nrow(Data)< EndN){
    
    N0 = nrow(Data)
    if(method==1){
      if(verbose)cat("LHC...",N0+1,"   ")
      UNIdesign(nrow(Data)+1)
      
    }else if(method==2){
      if(verbose)cat("KG_IU...",N0+1,"   ")
      xanew = KG_MC_Input()
      Data  = rbind(Data,c(xanew,TestFun(xanew)))
      
    }else if(method==3){
      if(verbose)cat("EGO_IU...",N0+1,"   ")
      xanew = EGO_MC_Input()
      Data  = rbind(Data,c(xanew,TestFun(xanew)))
      
    }else if(method==4){
      if(verbose)cat("EGO Input Mean..",N0+1,"   ")
      xanew = EGO_Single_Input(Distro$Mean)
      Data  = rbind(Data,c(xanew,TestFun(xanew)))
      
    }else if(method==5){
      if(verbose)cat("EGO Input Mode..",N0+1,"   ")
      xanew = EGO_Single_Input(Distro$Mode)
      Data  = rbind(Data,c(xanew,TestFun(xanew)))
      
    }
  
    UpdateGP()
    Cost = rbind(Cost,c(nrow(Data),TrueQ-Get_Quality()))
    
    
    if(verbose){
      cat("..", tail(Cost,1)[[2]],"\n")
      
      par(mfrow=c(3,1),mar=c(4.1,4.1,0.1,2.1))
      plot(Cost[,1:2],log="y")
      lines(Cost)
      
      XX = seq(0,100,len=30)
      Y = matrix(MU(cbind(XX,rep(XX,each=length(XX)))),length(XX))
  
      contour(XX,XX,t(Y),xlab= "A", ylab="X")
      points(Data[,2:1])
      points(Data[nrow(Data),2],Data[nrow(Data),1],col="red",pch=19)
      
      SS = 2*sqrt(sapply(XX,function(ai)mean(COV(cbind(XX,ai),cbind(XX,ai)))))
      MM = apply(Y,2,mean)
      plot(range(XX),range(c(MM+SS,MM-SS)),xlab = "A",ylab = "MC",type = "l",col="white")
      lines(XX,MM)
      lines(XX,MM+SS)
      lines(XX,MM-SS)
      
      Sys.sleep(0.2)
    }
  }
  list(method = c("LHC","KG_IU","EGO_IU","EGO_Mean","EGO_Mode")[method],seed=seed,Cost = Cost,Distro = Distro$name)
}


