bilinear.interp<-function(x,f){
	#interpolate a 2d point for fn f
	indx<-trunc(x[1]*(res-1)+1)
	indy<-trunc(x[2]*(res-1)+1)
	seq.Out<-seq(0,1,length.out=res)
	fr1<- (seq.Out[indx+1]-x[1])/(seq.Out[indx+1]-seq.Out[indx]) * f[indx,indy] + (x[1]-seq.Out[indx])/(seq.Out[indx+1]-seq.Out[indx]) * f[indx+1,indy]
	fr2<- (seq.Out[indx+1]-x[1])/(seq.Out[indx+1]-seq.Out[indx]) * f[indx,indy+1] + (x[1]-seq.Out[indx])/(seq.Out[indx+1]-seq.Out[indx]) * f[indx+1,indy+1] 
	fp<- (seq.Out[indy+1]-x[2])/(seq.Out[indy+1]-seq.Out[indy]) * fr1 + (x[2]-seq.Out[indy])/(seq.Out[indy+1]-seq.Out[indy]) * fr2
	return(fp)
}

bilinear.interp(c(.4,.5),matrix(rnorm(res^2),res,res))

univ.interp<-function(x,f){
	seq.Out<-seq(0,1,length.out=res)
	indx<-trunc(x*(res-1)+1)
	(seq.Out[indx+1]-x)/(seq.Out[indx+1]-seq.Out[indx]) * f[indx] + (x-seq.Out[indx])/(seq.Out[indx+1]-seq.Out[indx]) * f[indx+1]
}

univ.interp(.4,1:res)


tree.mix.lik<-function(x,edge_den,beliefs,edgelist,tree.list,weights,part_fn){
	e<- dim(edgelist)[1]
	d<-dim(beliefs)[1]
	uni.prod<-mapply(function(u){
		log(univ.interp(x[u],beliefs[u,]))
	},1:d)
	if(e==0){
		return(sum(uni.prod))
	}
		if(any(is.na(uni.prod))) print("help!")

	biv.prod<-mapply(
				function(l){
					s<- edgelist[l,1]
					t<- edgelist[l,2]
					log(bilinear.interp(c(x[s],x[t]) , edge_den[,,l]))  - (uni.prod[s]) - (uni.prod[t])
				},
				1:e
			)
	tree.prod<- lapply(tree.list,
					function(t){
						tree<- graph.adjacency(t,mode="undirected")
						el1<-apply(get.edgelist(tree), 1, paste, collapse="-")
						el2<-apply(get.edgelist(tree), 1, function(x)paste(rev(x),collapse="-"))
						el<-apply(edgelist, 1, paste, collapse="-")
						ekeep<- (el %in% el1) | (el %in% el2)
						sum( biv.prod[ekeep] )
					})
		tree.prod<- unlist(tree.prod)
		tmax<- max(tree.prod)		
		#tmax<-0
		#tree.weight<- exp(tree.prod-tmax)*weights
		#out<- sum( uni.prod ) + tmax + log( sum( tree.weight ) )
		out<- sum( uni.prod ) + weights * tree.prod
	return(out)
	}

# risk_path_tree<-lapply(path,function(x){
#   outvec<-lapply(1:dim(test)[1],
#   	function(i){
# 	  tree.mix.lik(test[i,],
# 		x$edge_den,
# 		x$beliefs,
# 		x$edgelist,
# 		x$tree_obj$tree.list,
# 		x$tree_obj$weight.list)
#   	})
#   	outvec<- unlist(outvec)
#   	return(- mean(outvec))
# })
# #risk_path_tree
# 
# # 
# # lik<-tree.mix.lik(test[1,],
# # 	fit$edge_den,
# # 	fit$beliefs,
# # 	edgelist,
# # 	tree_obj$tree.list,
# # 	tree_obj$weight.list
# # )