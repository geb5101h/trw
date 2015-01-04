########################
#Functions for calculating
#sufficient statistics
#
########################
legendre_Pl_array<-function(m,y){
  pf<-polynomial.functions(legendre.polynomials(m,normalized=F))
  lv<-sapply(X=1:(m+1),FUN=function(x){
    pf[[x]](y)
  })
  return(t(lv)*sqrt(2))
}

edge.means<-function(dat,edgelist,m,cores){
  n<-dim(dat)[1]
  nedge<- dim(edgelist)[1]
  if(nedge==0) return(array(dim=c(m,m,0)))
  me<-mclapply(
    1:nedge,
    function(e){
      ge<-edgelist[e,]
      s<-ge[1];t<-ge[2]
      tcrossprod(legendre_Pl_array(m,2*dat[,s]-1)[-1,,drop=F],
                 legendre_Pl_array(m,2*dat[,t]-1)[-1,,drop=F])/n
    },
    mc.cores=cores)
  return(array(unlist(me),dim=c(m,m,nedge)))
}

node.means<-function(dat,m,cores){
  if(is.vector(dat)){
    dat<-matrix(dat,length(dat),1)
  }
  d=dim(dat)[2]
  me<-mclapply(
    1:d,
    function(v){
      rowMeans(legendre_Pl_array(m,2*dat[,v]-1)[-1,,drop= F])
    },
    mc.cores=cores
  )
  return(matrix(unlist(me),nrow=d,byrow=T ))
}