#Sys.setenv("PKG_LIBS"="-lprofiler")
if(Sys.getenv("USER")=="ebjan"){
  setwd("~/trw")
}else{
  setwd("~/Dropbox/trw")
}
#Creates data used for tests

#source("~/Dropbox/trw/trw_run.r")
#source("~/Dropbox/trw/trw_fns.r")
source("suff_stats.r")
source("likelihoods.r")
Rcpp::sourceCpp("./trw.cpp")
source("trw_fit.r")
#source("~/Dropbox/trw/trw_global.r")

#setwd("./Experiments/Data")


libs=c(
  "Matrix",
  "orthopolynom",
  "igraph",
  "Rcpp",
  "inline",
  "RcppArmadillo",
  "parallel",
  "SnowballC",
  "ape",
  "mclust",
  "huge",
  "glasso",
  "abind")
if (sum(!(libs %in% .packages(all.available = TRUE))) > 0) {
  install.packages(libs[!(libs %in% .packages(all.available = TRUE))],
                   repos="http://cran.stat.ucla.edu")
}
for (i in 1:length(libs)) {
  library(libs[i],character.only = TRUE,quietly=TRUE)
}

library(RBGL)




res=128
seq.Out<-seq(from=0,to=1,length.out=res)
cores=detectCores()


legen.Vec<-legendre_Pl_array(m1,2*seq.Out-1)[-1,,drop=F]
#leg.mat<-legen.Vec%*%t(legen.Vec)
legen.array<-array(dim=c(res,res,m2^2))
k=1
for(i in 1:m2){
  for(j in 1:m2){
    legen.array[,,k]<-tcrossprod(legen.Vec[i,],legen.Vec[j,])
    k=k+1
  }
}

mixlik<-function(test,mclustout){
	d=dim(test)[2]
	out2=lapply(1:dim(test)[1],function(x){
		out=0
		for(i in 1:length(mclustout$pro)){
			out= out + mclustout$pro[i]*(det(mclustout$variance$sigma[,,i])^(-1/2))*
				exp(-.5*crossprod(test[x,]-mclustout$mean[,i],solve(mclustout$variance$sigma[,,i]))%*%(test[x,]-mclustout$mean[,i]))
		}
	return( log(out) - (d/2)*log(2*pi) )		
	})
	return(-mean(unlist(out2)))
}



