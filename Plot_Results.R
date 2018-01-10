Results = readRDS("~/Dropbox/PhD/ParREVI/TotalData733")

getRes = function( D,S,input = Results){
  
  mets = c("LHC","KG_IU","EGO_UI","EGO_Mean","EGO_Mode")
  
  DD = sapply(input,function(R)R$Distro)
  SS = sapply(input,function(R)R$method)
  
  O = input[DD==D&SS==mets[S]]
  
  OO = sapply(O,function(o)o$Cost$OC)
  
  M = apply(OO,1,mean)
  E = apply(OO,1,sd)*(1/sqrt(length(O)))
  cbind(M,E)
}
plotpoly=function(N,OC,se,col,low=-Inf){
  Xvals=c(N,rev(N))
  Yvals=c(OC+se,pmax(rev(OC-se),0.0001))
  col=adjustcolor(col,alpha.f = 0.4)
  # print(min(Yvals))
  Yvals = pmax(Yvals,low)
  
  polygon(Xvals,Yvals,col=col,border=NA)
  
} 

cols=c("pink","blue","green","red","purple")

par(mfrow=c(2,1))

OCU = sapply(1:5,function(s){
  getRes("Uniform",s)[,1]
})
ERU = sapply(1:5,function(s){
  getRes("Uniform",s)[,2]
})
plot(c(10,100),range(OCU),col="white",log="y",ylab = "OC",xlab = "N", main = "Uniform Uncertainty")
sapply(1:5,function(s){
  lines(10:100,OCU[,s],col=cols[s])
  plotpoly(10:100,OCU[,s],ERU[,s],cols[s])
})

OCT = sapply(1:5,function(s){
  getRes("Triangular",s)[,1]
})
ERT = sapply(1:5,function(s){
  getRes("Triangular",s)[,2]
})
plot(c(10,100),range(OCT),col="white",log="y",ylab = "OC",xlab = "N", main = "Triangular Uncertainty")
sapply(1:5,function(s){
  lines(10:100,OCT[,s],col=cols[s])
  plotpoly(10:100,OCT[,s],ERT[,s],cols[s])
})
