rm(list=ls())
setwd("~/Dropbox/trw")
setwd("~/trw")

m1=2;m2=2;
m=m1
source("startup.r")
source("trw_fit.r")
####################
#generate data
#
####################
set.seed(11400)


d=100
n=400
h<-huge.generator(10000,d=d,graph="random",prob=2/d)
h$data<-scale(h$data,center=T,scale=T)/8+0.5
#h$data<- pbeta(h$data,2,1)
#h$data<-h$data-colMeans(h$data)
#h$data<-qexp(pnorm(h$data),rate=12)
h$data[h$data>=1] = .99
h$data[h$data<=0] = .001

test<-h$data[ 9000:9100,]
#dat<-gauss.mix.sim(400,d=d)


sparsity<- h$sparsity*d*(d-1)/2

g<-graph.adjacency(h$theta,mode="undirected")
edgelist<- get.edgelist(g)
clust<- clusters(g)$membership

nout<-round( seq(from=20,to=100,length.out=15) )


out<-list()
iter=0
for(n in nout){
  iter=iter+1
  train<-h$data[1:n,] 
  edge_mean <- edge.means(train, edgelist, m, cores)
  node_mean <- node.means(train , m,cores)
  edge_mean_test <- edge.means(test, edgelist, m, cores)
  node_mean_test <- node.means(test , m,cores)
  weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])
  
  
  edge_coef<- array(0,dim=c(m,m,dim(edgelist)[1]))
  node_coef<- matrix(0,d,m)
  
  
  
  
  tree_obj<- edge.Start(g,burn.in=10)
  out[[iter]]<-parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
                              edgelist,E(tree_obj)$weight,lambda=0,alpha_start=1,
                              g,cores=cores)
}

risk_path<-lapply(out,function(x){
  -trw_lik(node_mean_test,edge_mean_test,x$vert_coef,x$edge_coef,x$part_fn)
})

risk_path<-simplify2array(risk_path)
entropy<- (d/2)*(1+log(2*pi)) + (1/2)*log(det(cov(h$data)))
plot(nout,risk_path,type="l",ylim=range( c(risk_path,entropy) ))
abline(h=entropy,col="red",type="d")
dev.off()