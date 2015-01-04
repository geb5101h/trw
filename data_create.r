rm(list=ls())
#Sys.setenv("PKG_LIBS"="-lprofiler")
setwd("~/Dropbox/trw")
#Creates data used for tests

#source("~/Dropbox/trw/trw_run.r")
#source("~/Dropbox/trw/trw_fns.r")
source("suff_stats.r")
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


d=20
n=10000
m=6
edgelist<- get.edgelist( graph.full(d,directed=F) )
g<-   erdos.renyi.game(d/2,.3); g<- g+g
#g<-graph.empty(d)
edgelist<- get.edgelist(g)
clust<- clusters(g)$membership
#edgelist<- get.edgelist( graph.tree(d))
#edgelist<- get.edgelist(graph.empty(d))
#edgelist<- cbind(edgelist[,2],edgelist[,1])
res=64
seq.Out<-seq(from=0,to=1,length.out=res)
weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])

###############

#Generate data
######################
h<-huge.generator(n,d=d,graph="random",prob=2/d)
h$data<-scale(h$data,center=T,scale=T)/10+0.5
#h$data<-qexp(pnorm(h$data),rate=12)
train<-h$data[1:floor(n*.5),] 
test<-h$data[ (floor(n*.5)+1) : n,]
#dat<-gauss.mix.sim(400,d=d)


train[1:1000,]<-apply(train[1:1000,],2,function(x){x+rnorm(d)*.08})
train[1001:2000,]<-apply(train[1001:2000,],2,function(x){x+rnorm(d)*.08})


#train[1100:2000,3:4]<-train[1100:2000,3:4]-.3
train[train>1] = .999 ;train[train<0] = .001
######################


cores<-parallel::detectCores()
edge_mean <- edge.means(train, edgelist, m, cores)
node_mean <- node.means(train , m,cores)
#gsl not working on ubuntu
legen.Vec<-legendre_Pl_array(m,2*seq.Out-1)[-1,]
#leg.mat<-legen.Vec%*%t(legen.Vec)
legen.array<-array(dim=c(res,res,m^2))
k=1
for(i in 1:m){
  for(j in 1:m){
    legen.array[,,k]<-tcrossprod(legen.Vec[i,],legen.Vec[j,])
    k=k+1
  }
}
if(dim(edgelist)[1]==0){
  edge_coef<- 
    array(dim=c(m,m,0))
}else{
  edge_coef<- sapply(1:dim(edgelist)[1],
                   function(x)matrix(tcrossprod(rnorm(m)*.000^(1:m) , rnorm(m)*.000^(1:m)),m,m),
                   simplify="array")
}
node_coef<- t(sapply(1:d,function(x)rnorm(m)*.00^(1:m)))







Rcpp::sourceCpp("./trw.cpp")

#parallelMatrixSqrt(matrix(10,10,10))

lambda=.1
alpha_start=1


out<-rwrapper(
  edge_coef,
  node_coef,
  edge_mean,
  node_mean,
  edgelist,
  weights,
  legen.Vec,
  legen.array, 
  lambda,
  alpha_start
)


par(mfcol=c(2,2))
for(i in 1:4){
  contour((out$edge_den[,,i]),nlevels=20)
  points(train[,(edgelist[i,])],cex=.01)
}


# r<-1:10
# #const_test(r)
# 
# microbenchmark::microbenchmark(
#   neg_ll(node_means,edge_means,node_coefs,edge_coefs,weights,
#          edgelist,legen.Vec,legen.array),
#   times=100
# )
# 
# microbenchmark(
# a<-beliefprop(edge_coefs, 
#               node_coefs, 
#               rep(1,dim(edgelist)[1]*2), 
#               rbind(edgelist, 
#               cbind(edgelist[,2],edgelist[,1])), 
#               legen.Vec,
#               legen.array,
#               cores=cores )
# )
# gq=get_quasi(a$edge_den, a$beliefs, edge_means,  node_means, legen.Vec)
# mi=mutual_info(a$edge_den,a$beliefs,edgelist)
# pf=partition_fn(mi$entropy,mi$mutual_info, node_coefs, edge_coefs,gq$vert_quasi,gq$edge_quasi, rep(1,dim(edgelist)[1]))
# 
# microbenchmark::microbenchmark(
# neg_ll(node_means,edge_means,node_coefs,edge_coefs,weights,
#        edgelist,legen.Vec,legen.array),
# times=1
# )
# 
# microbenchmark::microbenchmark(
#   get_quasi(a$edge_den, a$beliefs, edge_means,  node_means, legen.Vec),
#   get.Quasi(a$beliefs,aperm(a$edge_den,c(3,1,2)), node_means, aperm(edge_means, c(3,1,2))),
#   times=10
# )
# 
# microbenchmark::microbenchmark(
#   mutual_info(a$edge_den,a$beliefs,edgelist),
#   mut.info(a$beliefs,a$edge_den,edgelist),
#   times=5
# )
# microbenchmark::microbenchmark(
#   #beliefprop(edge_coefs, node_coefs, rep(1,dim(edgelist)[1]*2), rbind(edgelist, cbind(edgelist[,2],edgelist[,1])), legen.Vec, legen.array,cores=1 ),
#   beliefprop(edge_coefs, node_coefs, rep(1,dim(edgelist)[1]*2), rbind(edgelist, cbind(edgelist[,2],edgelist[,1])), legen.Vec, legen.array,cores=cores ),
#   times=10
# )
# 
# 
# pdf("dataplot.pdf")
# plot(h)
# dev.off()
# save(file="data.rdata",
#            list=c("d","n","h","train","test", "edge_means","node_means"))