#Tree-reweighted belief propagation for
#continuous-valued density estimation

rm(list=ls())

if(Sys.info()["user"]=="eric") setwd("~/Dropbox/trw")
if(Sys.info()["user"]=="ebjan")setwd("~/trw")

libs=c("Matrix","orthopolynom","igraph","Rcpp","inline","RcppArmadillo","parallel","SnowballC","ape","mclust","glasso")
if (sum(!(libs %in% .packages(all.available = TRUE))) > 0) {
    install.packages(libs[!(libs %in% .packages(all.available = TRUE))],
    	repos="http://cran.stat.ucla.edu")
}
for (i in 1:length(libs)) {
    library(libs[i],character.only = TRUE,quietly=TRUE)
}
library("RBGL")



################
#GLOBAL VARS AND FNS
source("http://dl.dropbox.com/s/xk0gq4ae3faec52/trw_c.r")
source("http://dl.dropbox.com/s/b04iktw8bgosxki/trw_global.r")
source("http://dl.dropbox.com/s/mgreojd9wl2gmfr/trw_fns.r")


#####

#d=4 #num dimensions

# dat<-gauss.mix.sim(100)
# 
# #create graph object
# #graph.obj<-graph.full(d,directed=F)
# graph.obj<-graph.tree(d,mode="undirected")
# #graph.obj<-graph.lattice(c(5,3))
# graph.obj<-graph.create(graph.obj,dat=dat)
# graph.obj<-ISTA(graph.obj,alpha.start=1,lambda=.0,damp=1)
# 
# 
# 	
#graph.obj<-graph.full(d,directed=F)
set.seed(121455)	
n.set=floor(seq(from=10,to=1000,length.out=10))
d.set=c(3,4,5,6)

l=1
outfinal<-mclapply(mc.cores=3,X=d.set,FUN=function(d){
	graph.obj<-graph.full(d,directed=F)
	outp<-mclapply(mc.cores=8,X=1:25,FUN=function(x){
	graph.obj<-graph.full(d,directed=F)

		ss=sample(1:1000,1)
		k=1
		oo=vector()
		
		for(nn in n.set){
			test<-gauss.sim(1000,d=d,mu=.5,Sigma=diag(d)*.05)$dat
			print(paste("n: ",nn,"d: ",d))
			gt<-gauss.sim(nn,d=d,mu=.5,Sigma=diag(d)*.05)
			train<-gt$dat
			nll<-gt$nll
			graph.obj<-graph.create(graph.obj=graph.obj,dat=train,weight=2/d)
			graph.obj<-try(ISTA(graph.obj,alpha.start=1,lambda=0,damp=1))
			graph.obj<-graph.create(graph.obj=graph.obj,dat=test,weight=2/d)
			oo[k]<-neg.ll(graph.obj)$nll-nll
			k=k+1
		}
		oo	
	})
	outd<-rowMeans(
			matrix(unlist(outp),nrow=length(n.set))
	,na.rm=T)
	return(outd)	
})



save.image("simrun.rda")
# 
# ISTA(graph.obj,lambda=0,alpha.start=5)



# 
# 
# 
# 	
# #starting values for ML
# edge.Co1<-lapply(1:nedge.unique,function(x) matrix(0,m.set,m.set))
# names(edge.Co1)<-edges[which(less.edge(edges))]
# node.Co1<-lapply(1:d,function(x) rep(0,m.set))
# names(node.Co1)<-1:d
# co.start=c(unlist(node.Co1),unlist(edge.Co1))
# #ll<-lik.Pass(co.start,adj.Sim,mp.sim,edge.Weight=ew)
# #co.start<-c( rep( c(3/4,-1/4,rep(0,m.set-2)) , d ) , unlist(edge.Co1) )
# co.start<-c(unlist(node.Co1),unlist(edge.Co1))
# 
# 
# assign("MES.START",lapply(1:6,function(x) rep(1,res)),envir=lik.env)
# 
# ew1<-edge.Start(adj.Sim,burn.in=25)
# system.time(
# 	fista<-FISTA(co.start,x.Node.Mean,x.Edge.Mean,edge.Weight=ew1,lambda=lmax*2)
# )
# 
# 
# lik.Graph(fista[[1]],mean.Edge=x.Edge.Mean,mean.Node=x.Node.Mean,adj.Mat=adj.Sim,
# 				damp=1,edge.Weight=ew1,tree.optim=F)
# 
# q<-get("qua",lik.env)
# tt=tree.mix.lik(x.train,ew1,q$den)
# -mean(tt);fista[[2]]
# 
# 
# 
# system.time(
# 	ista<-ISTA2(co.start,x.Node.Mean,x.Edge.Mean,edge.Weight=ew1,lambda=2*0)
# )
# pgo<-list()
# i=0
# lam.seq<-exp(seq(from=log(lmax*.001),to=log(lmax*2),length.out=10))
# for(l in rev(lam.seq)){
# 	i=i+1
# 	print(paste("Step number",i))
#  pgo[[i]]<-ISTA2(co.start,x.Node.Mean,x.Edge.Mean,edge.Weight=ew1,lambda=l)
#  co.start<-pgo[[i]][[1]]
# }
# 
# 
# 
# 
# 
# 
# ist<-ISTA2(co.start,x.Node.Mean,x.Edge.Mean,edge.Weight=ew1,lambda=lmax+.05)
# 
# 
# 
# 
# ll=lik.Graph(co.start,mean.Edge=x.Edge.Mean,mean.Node=x.Node.Mean,adj.Mat=adj.Sim,
# 				damp=1,edge.Weight=ew1,tree.optim=F)
# q<-get("qua",lik.env)
# #pdf("testist.pdf")
# contour(q[[3]][[2]][[5]],nlevels=50)
# points(x.train[,c(3,2)],pch=16,cex=.5)
# dev.off()
# 
# 
# gr<-attr(ll,"grad")
# 
# co2<-co.start - 1*gr
# ll2=lik.Graph(co2,mean.Edge=x.Edge.Mean,mean.Node=x.Node.Mean,adj.Mat=adj.Sim,
# 				damp=1,edge.Weight=ew1,tree.optim=F)
# 				
# 				ll[1];ll2[1]
# 
# 	pdf("ist.pdf")
# 	plot(ist[[2]])
# 	dev.off()
# 	pdf("ist2.pdf")
# 	plot(attr(ll,"grad"))
# 	dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# pgo<-list()
# i=0
# for(l in rev(seq(from=0,to=lmax+.05,length.out=5))){
# 	i=i+1
# 	print(paste("Step number",i))
#  pgo[[i]]<-ISTA2(co.start,x.Node.Mean,x.Edge.Mean,edge.Weight=ew1,lambda=l)
#  co.start<-pgo[[i]][[1]]
# }
# 
# max.lik<-optim(fn=lik.Pass,gr=lik.Grad,par=co.start,damp=1,mp=mp.sim,
# 	edge.Weight=ew1,tree.optim=F,
# 	method="CG",control=list(trace=6))
# 
# save.image(file="save.RData")
# 
# #global variable to save edge weights	
# assign("EDGE.START",edge.Start(adj.Sim,burn.in=25),envir=lik.env)
# 
# max.lik.Edge<-optim(fn=lik.Pass,gr=lik.Grad,par=max.lik$par,damp=1,mp=mp.sim,
# 	edge.Weight=NULL,tree.optim=F,
# 	method="CG",control=list(trace=4,factr=1e12))
# save.image(file="save.RData")
# 
# 
#  node.Q<-lapply(1:d,function(x) max.lik$par[((x-1)*m.set+1):(x*m.set)])
#  	names(node.Q)<-1:d
#  edge.Q<-lapply(1:nedge.unique,function(x){
#  		matrix(max.lik$par[(d*m.set+(x-1)*m.set*m.set+1):(d*m.set+x*m.set*m.set)],m.set,m.set)
#  	})
#  	names(edge.Q)<-edges[which(less.edge(edges))]
# two<- tree.weight.optim(edge.Q,node.Q,ew.start=edge.Start(adj.Sim))
#   two$tree.list<-two$tree.list[order(two$weight.list,decreasing=T)]
#  two$weight.list<-sort(two$weight.list,decreasing=T)
#  par(mfcol=c(3,3))
#  for(i in 1:9){
#  	g<-graph.adjacency(two$tree.list[[i]],mode="undirected")
#  	#g$x<-rep(1:3,3);g$y=sort(rep(1:3,3))
#  	#attr(g,"x")<-rep(1:3,3);attr(g,"y")<-sort(rep(1:3,3))
#  	g$layout<-layout.circle
#  	plot(g,margin=0,main=paste(two$weight.list[i]))
#  }
# 
# a<- lik.Graph(max.lik.Edge$par,mean.Edge=test.Edge.Mean,test.Node.Mean,adj.Mat=adj.Sim,damp=1,
#  	edge.Weight=NULL,tree.optim=T)
# 	
# lik.Graph(max.lik$par,mean.Edge=test.Edge.Mean,test.Node.Mean,adj.Mat=adj.Sim,damp=1,
# 	edge.Weight=NULL,tree.optim=T)
	 
# lik.Graph(max.lik.Edge$par,mean.Edge=test.Edge.Mean,test.Node.Mean,adj.Mat=adj.Sim,damp=1,
# 	edge.Weight=NULL,tree.optim=T)
# 	
# lik.Graph(max.lik$par,mean.Edge=test.Edge.Mean,test.Node.Mean,adj.Mat=adj.Sim,damp=1,
# 	edge.Weight=NULL,tree.optim=T)
# 	
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###################
# #Preset coefs
# 
# 
# dec.co=.8
# dec<-dec.co^(1:m.set)%*%t(dec.co^(1:m.set))
# eco<-dec*matrix(rnorm(m.set^2,0,.1),m.set,m.set)
# edge.Co<-lapply(1:nedge.unique,function(x) dec*matrix(rnorm(m.set^2,0,.5),m.set,m.set))
# names(edge.Co)<-edges[which(less.edge(edges))]
# 
# #dxm matrix of node potential coefficients
# #node.Co<-lapply(1:d,function(x)rnorm(m.set,0,1)*dec.co^(1:m.set))
# node.Co<-lapply(1:d,function(x) c(1,rep(0,m.set-1)))
# 
# names(node.Co)<-1:d
# co1=c(unlist(node.Co),unlist(edge.Co))
# 
# 
# ew<-edge.Start(adj.Sim,burn.in=1)
# two=tree.weight.optim(edge.Co,node.Co,ew.start=ew)
# two$tree.list<-two$tree.list[order(two$weight.list,decreasing=T)]
# two$weight.list<-sort(two$weight.list,decreasing=T)
# par(mfcol=c(3,4))
# for(i in 1:9){
# 	g<-graph.adjacency(two$tree.list[[i]],mode="undirected")
# 	#g$x<-rep(1:3,3);g$y=sort(rep(1:3,3))
# 	#attr(g,"x")<-rep(1:3,3);attr(g,"y")<-sort(rep(1:3,3))
# 	g$layout<-layout.circle
# 	plot(g,margin=0,main=paste(two$weight.list[i]))
# }
# dev.off()
# ################
# assign("MES.START",lapply(edges,function(x) rep(1,res)),envir=lik.env)
# assign("EDGE.START",edge.Start(adj.Sim),envir=lik.env)
# 
# Rprof()
# 
# mp<-get.Quasi(m.set=m.set,edge.Co=edge.Co,node.Co=node.Co,damp=1,edge.Weight=edge.Start(adj.Sim)$ew)
# mut.info(mp[[3]])
# #	a<-series.bp(adj.Sim,edge.Co=edge.Co,
# #	node.Co=node.Co,bp.tol=1e-12,m=m.set,damp=.95)
# Rprof(NULL) 
# summaryRprof()
# 
# 
# xxx<-cbind(unlist(mut.info(mp[[3]])$mut.info),summary(two$ew)$x)
# cor(xxx)
# 
# 
# 	lik=0
# 	edge.Mat<-get("EDGE.START",lik.env)
# 	
# 	#put likelihood of maximum tree and mixture of trees
# 	
# 	
# node.Q<-lapply(1:d,function(x) max.lik$par[((x-1)*m.set+1):(x*m.set)])
# 	names(node.Q)<-1:d
# edge.Q<-lapply(1:nedge.unique,function(x){
# 		matrix(max.lik$par[(d*m.set+(x-1)*m.set*m.set+1):(d*m.set+x*m.set*m.set)],m.set,m.set)
# 	})
# 	names(edge.Q)<-edges[which(less.edge(edges))]
# tree.weight.optim(edge.Q,node.Q,ew.start=get("EDGE.START",lik.env))
# 
# 
# qm<-get.Quasi(m.set=m.set,edge.Co=edge.Q,
# 	node.Co=node.Q,damp=1,edge.Weight=edge.Mat)
# 	
# # a2<-series.bp(adj.Sim,edge.Co=edge.Q,
# # 	node.Co=node.Q,bp.tol=1e-12,m=m.set,damp=.95,edge.Weight=ew)
# 
# 
# 
# save.image(file="trw_out.RData")
# pdf("trw_plot2.pdf",width=14,height=8)
#  par(mfcol=c(2,length(edges)/4),oma=c(0,0,3,0))
#  i=1
#  for(e in edges[less.edge(edges)]){
#  f<-strsplit(e,".",fixed=T)[[1]]
#  t<-f[2];s<-f[1]
#  contour(a2[[2]][[i]],nlevels=15,xlab=paste(s),ylab=paste(t))
#  points(x[,as.numeric(s)],x[,as.numeric(t)],pch=16,cex=.5)
#  i=i+1
#  if(i>length(edges)/2)break
#  }
#  mtext(paste("TRW approx neg. log lik.: ",round(max.lik$value,3),"\n Exact mixture of Gaussian neg. log lik.: ",round(mix.lik(x),3)),outer=T,cex=1)
# 
# dev.off()
