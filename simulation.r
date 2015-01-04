rm(list=ls())
setwd("~/Dropbox/trw")
#setwd("~/trw")
source("startup.r")
source("trw_fit.r")
####################
#generate data
#
####################
set.seed(1140)


d=10
n=50
h<-huge.generator(n+1000,d=d,graph="random",prob=2/d)
h$data<-scale(h$data,center=T,scale=T)/8+0.5
#h$data<- pbeta(h$data,6,1)
#h$data<-h$data-colMeans(h$data)
#h$data<-qexp(pnorm(h$data),rate=12)
h$data[h$data>=1] = .99
h$data[h$data<=0] = .001
train<-h$data[1:n,] 
test<-h$data[ (n+1):(n+100),]
#dat<-gauss.mix.sim(400,d=d)


sparsity<- h$sparsity*d*(d-1)/2


####################
#fit path for glasso
#
########################################

cc<-cov(train)
rho_start<- max( abs(cc-diag(diag(cc))) )



run_seq <- exp( seq(from=log(rho_start*.03), to =log(rho_start*1.1),length.out=50 ) )
#run_seq<- seq(from=1e-5,to=rho_start,length.out=25)
test_error<- vector(length=length(run_seq))
ind=0


glasso_out<-glassopath(s=cc,
           rholist=run_seq)$wi

glasso_out_relax<-mcmapply(FUN=function(x){
    ge<-get.edgelist( 
      graph.adjacency( (glasso_out[,,x]==0),
                       mode="undirected") 
    )
    if(dim(ge)[1]==0) ge<-NULL
    glasso(s=cc,
             rho=1e-5*rho_start,
             zero=ge
      )$wi
    },
    x=1:length(run_seq),
    SIMPLIFY="array"
)


test_error<-mcmapply(FUN=function(x){
  print(x)
  -mean( apply(test,1,function(xx)gaussian_lik(xx-colMeans(test),glasso_out[,,x])))
  },
  x=1:length(run_seq),
  SIMPLIFY=T,
  mc.cores=cores
)


test_error_relax<-mcmapply(FUN=function(x){
  print(x)
  -mean( apply(test,1,function(xx)gaussian_lik(xx-colMeans(test),glasso_out_relax[,,x])))
},
x=1:length(run_seq),
SIMPLIFY=T,
mc.cores=cores
)

num_nonzero<-mapply(function(x)(sum(glasso_out[,,x]!=0) -d)/2,x=1:length(run_seq))
num_nonzero_relax<-mapply(function(x)(sum(glasso_out_relax[,,x]!=0) -d)/2,x=1:length(run_seq))


#pdf("sim.pdf")
plot(num_nonzero,test_error,type="l",ylim=range(test_error,test_error_relax))
points(num_nonzero,test_error_relax,col="blue",type="l")
abline(v=sparsity,col="red")
#dev.off()

which.min(test_error_relax)
which.min(test_error)
min(test_error);min(test_error_relax)

####################
#fit path for trw
#
########################################




#edgelist<- get.edgelist( graph.full(d,directed=F) )
#g<-   erdos.renyi.game(d,2/d); #  g<- g+g
#g<-graph.empty(d)
g<-graph.full(d,directed=F)
edgelist<- get.edgelist(g)
clust<- clusters(g)$membership




edge_mean <- edge.means(train, edgelist, m, cores)
node_mean <- node.means(train , m,cores)
edge_mean_test <- edge.means(test, edgelist, m, cores)
node_mean_test <- node.means(test , m,cores)
weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])


edge_coef<- array(.1,dim=c(m,m,dim(edgelist)[1]))
node_coef<- matrix(0,d,m)

weights<- E(edge.Start(g,burn.in=10))$weight

norm<-apply(edge_mean,3,function(x)norm(x,"F"))%>%sort(decreasing=T)

fit<-parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
               edgelist,weights,lambda=0,alpha_start=1,
               g,cores=cores)
fit$edge_quasi-edge_mean
fit$vert_quasi-node_mean


contour(fit$edge_den[,,1])
points(train[,edgelist[1,]])


lambda=.1
alpha_start=1

path<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
               edgelist,weights=NULL,alpha_start=10,
               g,cores=cores,limit=30,relax=F,edge.opt=T)
path_relax<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
               edgelist,weights=NULL,alpha_start=10,
               g,cores=cores,limit=30,relax=T,edge.opt=T)




risk_path<-lapply(path,function(x){
  -trw_lik(node_mean_test,edge_mean_test[,,x$eind,drop=F],x$vert_coef,x$edge_coef,x$part_fn)
})
risk_path_relax<-lapply(path_relax,function(x){
  -trw_lik(node_mean_test,edge_mean_test[,,x$eind,drop=F],x$vert_coef,x$edge_coef,x$part_fn)
})

risk_path<-simplify2array(risk_path)
risk_path_relax<-simplify2array(risk_path_relax)


plot(risk_path,type="l",ylim=range(risk_path,risk_path_relax,test_error_relax,test_error))
points(1:length(risk_path_relax),risk_path_relax,type="l",col="blue")
points(num_nonzero,test_error_relax,col="green",type="l")
points(num_nonzero,test_error,col="purple",type="l")
abline(v=sparsity,col="red")
legend("topright",c("trw","trw relax","glasso relax", "glasso"),col=c("red","blue","green","purple"))


par(mfcol=c(2,2))
for(i in sample(dim(path_relax[[20]]$edgelist)[1],4)){
contour(path_relax[[20]]$edge_den[,,i])
points(train[,path_relax[[20]]$edgelist[i,]])
}
########################################




summary(h$theta)
plot(train[,c(15,95)])
plot(train[,c(100,99)])
microbenchmark::microbenchmark(

path<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
         edgelist,weights,alpha_start=1,
         g,cores=cores,limit=25),
path<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
               edgelist,weights,alpha_start=1,
               g,cores=1,limit=25),
times=10

)



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


weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])
parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
               edgelist,weights,lambda=0,alpha_start=1,
               g,cores=1)$part_fn




lambda=.04;alpha_start=1
microbenchmark::microbenchmark(
  parallelTrwFit(edge_coef,node_coef,node_mean,edge_mean,
                 edgelist,weights,lambda=0,alpha_start=1,
                 g,cores=1),
  parallelTrwFit(edge_coef,node_coef,node_mean,edge_mean,
                 edgelist,weights,lambda=0,alpha_start=1,
                 g,cores=detectCores())
  ,times=10)








tree_obj<- edge.Start(g,burn.in=30)
parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
               edgelist,E(tree_obj)$weight,lambda=0,alpha_start=1,
               g,cores=1)$part_fn

