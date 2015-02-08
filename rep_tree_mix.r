
library(huge)
library(clime)
library(glasso)
library(Rcpp)
library(clime)
library(parallel)
library(mclust)
setwd("~")
source("~/Forest Density/treeLoader.R");
loadTree()
set.seed(984)
########################

n<-100
ntest<-300
dlist<- c(30,50,100,120) #,120,150)
#dlist=30
reps=5

out_len=15
m1=m2=7
#err=syst=sparsity=array(0,dim=c(len,3,length(dlist),reps))
source("~/trw/startup.r")


err<-array(0,dim=c(6,length(dlist),reps))
roc<-array(0,dim=c(5,length(dlist),reps,2))
diter=0
for(d in dlist){
	limit=floor(d*2)

	diter=diter+1
	h1<-huge.generator(reps*(n+ntest),d=d,graph="scale-free")
	h2<-huge.generator(reps*(n+ntest),d=d,graph="scale-free")
	h3<-huge.generator(reps*(n+ntest),d=d,graph="scale-free")

	h1$data<-scale(h1$data,center=T,scale=T)/8+0.5
	h1$data<- sign(h1$data-.5)*(abs(h1$data-.5))^.6
    h1$data<-scale(h1$data,center=T,scale=T)/4+0.5

	h1$data[h1$data>=1] = .99
	h1$data[h1$data<=0] = .001

	h2$data<-scale(h2$data,center=T,scale=T)/8+0.5
	h2$data<- sign(h2$data-.5)*(abs(h2$data-.5))^.6
    h2$data<-scale(h2$data,center=T,scale=T)/4+0.5

	h2$data[h2$data>=1] = .99
	h2$data[h2$data<=0] = .001

	h3$data<-scale(h3$data,center=T,scale=T)/8+0.5
	h3$data<- sign(h3$data-.5)*(abs(h3$data-.5))^.6
    h3$data<-scale(h3$data,center=T,scale=T)/4+0.5

	h3$data[h3$data>=1] = .99
	h3$data[h3$data<=0] = .001

	h<-rbind(h1$data,h2$data,h3$data)
	h<-h[sample(1:dim(h)[1]),]
	for(r in 1:reps){
		dnum=0
		print(paste("dim : ",d))
		train<-h[((r-1)*n+1):(r*n),] 
		test<-h[(reps*n+((r-1)*ntest+1)) : (reps*n+(r*ntest)),]
		rho_start<- max( abs(cov(train)-diag(diag(cov(train)))) )
		sequ <- exp( seq(from=log(rho_start*.2), to =log(rho_start*1.1),length.out=20 ) )
		#sequ<-log(rev(seq(from=exp(.01),to=exp(0.3),length.out=len))) *(1/20)
		gl<-glassopath(cov(train),sequ)
		glasso_out_relax<-mcmapply(FUN=function(x){
		ge<-get.edgelist( 
			graph.adjacency( (gl$wi[,,x]==0),
				mode="undirected") 
			)
			if(dim(ge)[1]==0) ge<-NULL
			glasso(s=cov(train),
			rho=1e-4,
			zero=ge,
			penalize.diagonal=F
		)$wi
		},
		x=1:length(sequ),
		SIMPLIFY="array"
		)

		gl_seq<-apply(gl$wi,3,
			function(x){.5*sum(diag(x%*%cov(test)))-.5*log(det(x))+(d/2)*log(2*pi)}
		)
		gl_seq_relax<-apply(glasso_out_relax,3,
			function(x){.5*sum(diag(x%*%cov(test)))-.5*log(det(x))+(d/2)*log(2*pi)}
		)
		err[1,diter,r]<-min(gl_seq)
		err[2,diter,r]<-min(gl_seq_relax)
		# ind= row(h)<=(col(h)-1)
# 		spar=h$sparsity*d*(d-1)/2
# 
# 		roc[1,diter,r,1]<-1-sum( (h$theta[ind]!=0) & (gl$wi[,,which.min(gl_seq)][ind]==0) )/spar
# 		roc[1,diter,r,2]<- 1- sum( (gl$wi[,,which.min(gl_seq)][ind]!=0) & (h$theta[ind]==0) )/(d*(d-1)/2-spar)
# 		roc[2,diter,r,1]<-1-sum( (h$theta[ind]!=0) & (glasso_out_relax[,,which.min(gl_seq_relax)][ind]==0) )/spar
# 		roc[2,diter,r,2]<- 1- sum( (glasso_out_relax[,,which.min(gl_seq_relax)][ind]!=0) & (h$theta[ind]==0) )/(d*(d-1)/2-spar)



		g<-graph.full(d,directed=F)
		edgelist<- get.edgelist(g)
		clust<- clusters(g)$membership
		edge_mean <- edge.means(train, edgelist, m2, cores)
		node_mean <- node.means(train , m1,cores)
		edge_mean_test <- edge.means(test, edgelist, m2, cores)
		node_mean_test <- node.means(test , m1,cores)
		weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])
		edge_coef<- array(0,dim=c(m2,m2,dim(edgelist)[1]))
		node_coef<- matrix(0,d,m2)

		lambda=0
		alpha_start=1
		cores=detectCores()
		path<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
		edgelist,weights=NULL,alpha_start=alpha_start,
			g,cores=cores,limit=limit,relax=F,edge.opt=T,out_len=out_len)
		path_relax<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
		edgelist,weights=NULL,alpha_start=alpha_start,
			g,cores=cores,limit=limit,relax=T,edge.opt=T,out_len=out_len)

		risk_path<-lapply(path,function(x){
		-trw_lik(node_mean_test,edge_mean_test[,,x$eind,drop=F],x$vert_coef,x$edge_coef,x$part_fn)
		})
		risk_path_relax<-lapply(path_relax,function(x){
		-trw_lik(node_mean_test,edge_mean_test[,,x$eind,drop=F],x$vert_coef,x$edge_coef,x$part_fn)
		})
		err[3,diter,r]<-min(unlist(risk_path))
		err[4,diter,r]<-min(unlist(risk_path_relax))
		gr=graph.edgelist(path[[which.min(risk_path)]]$edgelist,directed=F)
		gr_relax=graph.edgelist(path_relax[[which.min(risk_path_relax)]]$edgelist,directed=F)
		# roc[3,diter,r,1]<-1-.5*sum( (h$theta!=0) & (gr[]==0) )/spar
# 		roc[3,diter,r,2]<- 1- .5*sum( (gr[]!=0) & (h$theta==0) )/(d*(d-1)/2-spar)
# 		roc[4,diter,r,1]<-1-.5*sum( (h$theta!=0) & (gr_relax[]==0) )/spar
# 		roc[4,diter,r,2]<- 1- .5*sum( (gr_relax[]!=0) & (h$theta==0) )/(d*(d-1)/2-spar)

		graphs = forestDensityEst(train, 2*d)
		likehd_out = computeLogLikehd(train, test, list(graphs));
		err[5,diter,r]<- -likehd_out$likehd_ls
		#roc[5,diter,r,1]<-1-sum( (h$theta[ind]!=0) & (graphs[][ind]==0) )/spar
		#roc[5,diter,r,2]<- 1- sum( (graphs[][ind]!=0) & (h$theta[ind]==0) )/(d*(d-1)/2-spar)
		mcl<-densityMclust(train)
		err[6,diter,r]<-mixlik(test,mcl$parameters)
		save(err,file="~/trw/Experiments/rep_tree_mix.rdata")
	}
}
save.image("~/trw/Experiments/rep_tree_mix.rdata")