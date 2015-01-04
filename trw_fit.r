library(magrittr)

##############################
#Estimates the TRW model for fixed lambda
#Parallelizes by connected components
#
##############################

parallelTrwFit<-function(edge_coef,node_coef,edge_mean,node_mean,
                         edgelist,weights,lambda,alpha_start,
                         g,cores=1){
  clust<- clusters(g)$membership
  d<-dim(node_coef)[1]
  res<-dim(legen.Vec)[2]
  m<-dim(edge_mean)[1]
  out_list=list()
  out_list$beliefs<-matrix(nrow=d,ncol=res)
  out_list$edge_den<- array(dim=c(res,res,dim(edgelist)[1]))
  out_list$vert_coef<-matrix(nrow=d,ncol=m)
  out_list$vert_quasi<-matrix(nrow=d,ncol=m)
out_list$vert_gradient<-matrix(nrow=d,ncol=m)

  out_list$mutual_info<- vector(length=dim(edgelist)[1])
  
  out_list$edge_coef<- array(dim=c(m,m,dim(edgelist)[1]))
  out_list$edge_quasi<- array(dim=c(m,m,dim(edgelist)[1]))
  out_list$edge_gradient<- array(dim=c(m,m,dim(edgelist)[1]))
  out_list$mes <- array(dim=c(2*dim(edgelist)[1],res))
  clust_out<-mclapply(unique(clust),
                      function(c){
                        sub<- induced.subgraph(g,vids=(clust==c))
                        sub_edgelist<- get.edgelist(sub)
                        ekeep<- apply(edgelist,1,function(x)any( (clust==c)[x]) )
                        out<-rwrapper(
                          edge_coef[,,ekeep,drop=F],
                          node_coef[c==clust,,drop=F],
                          edge_mean[,,ekeep,drop=F],
                          node_mean[c==clust,,drop=F],
                          sub_edgelist,
                          weights[ekeep],
                          legen.Vec,
                          legen.array, 
                          lambda,
                          alpha_start
                        )
                        return(out)
                      },
                      mc.cores=cores,
                      mc.preschedule=F
  )
  out_list$part_fn<-0
  out_list$nll<-0
  for(c in unique(clust)){
    out_list$part_fn<-out_list$part_fn+clust_out[[c]]$part_fn
    out_list$nll<-out_list$nll+clust_out[[c]]$nll

    sub<- induced.subgraph(g,vids=(clust==c))
    sub_edgelist<- get.edgelist(sub)
    ekeep<- apply(edgelist,1,function(x)any( (clust==c)[x]) )
    out_list$beliefs[c==clust,]<-clust_out[[c]]$beliefs
    out_list$edge_den[,,ekeep]<- clust_out[[c]]$edge_den
    out_list$vert_coef[c==clust,]<-clust_out[[c]]$vert_coef
    out_list$edge_coef[,,ekeep]<- clust_out[[c]]$edge_coef
    out_list$mutual_info[ekeep]<- clust_out[[c]]$mutual_info
    out_list$vert_quasi[c==clust,]<- clust_out[[c]]$vert_quasi
    out_list$edge_quasi[,,ekeep]<- clust_out[[c]]$edge_quasi

    out_list$vert_gradient[c==clust,]<- clust_out[[c]]$vert_gradient
    out_list$edge_gradient[,,ekeep]<- clust_out[[c]]$edge_gradient
    out_list$mes[c(ekeep,ekeep),]<-clust_out[[c]]$mes
  }
  return(out_list)
}



##############################
#Fits TRW on path of lambdas
#using warm starts
##############################

trw_path<-function(edge_coef_start,node_coef_start,node_mean,edge_mean,
                   edgelist,weights=NULL,alpha_start=1,
                   g,cores=1,limit=100,relax=T,edge.opt=T,out_len=50,
                   lambda_range=NULL){
  #lambda_path<- apply(edge_mean,3,norm,"F")
  lambda_seq<- lapply(1:dim(edgelist)[1],
   	                    function(x){
   	                      norm(edge_mean[,,x]-tcrossprod(node_mean[edgelist[x,1],],node_mean[edgelist[x,2],]),"F")
   	                    })
  lambda_seq<-simplify2array(lambda_seq)
  if(!is.null(lambda_range)){
  	lambda_path <- exp( seq(from = log(max(lambda_seq)), to = log(max(lambda_seq)*.02), length.out = out_len) )
  }else{
  	lambda_path<- lambda_seq
  }
  edge_coef<-edge_coef_start
  node_coef<-node_coef_start
  out<-list()
  n=0
  seqlp<-round( seq(from=1,to=min(limit,length(lambda_path)),length.out=out_len) )
  lp<-(lambda_path%>%sort(decreasing=T))[seqlp]
  for(i in lp){
    lam<-ifelse(relax==T,.0001,i)
    n=n+1
    eind<- which(lambda_seq > i)
    elist_new<-edgelist[eind,,drop=F]
    sub<-subgraph.edges(g,eids=eind,delete.vertices=F)
    if(edge.opt==T){
	  tree_obj<-edge.Start(sub,burn.in=10)
      wt<-E(tree_obj)$weight
      if(is.null(wt))wt<-logical(0)
    }else{
      wt<- weights[eind]
    }
    out[[n]]<- parallelTrwFit(edge_coef[,,eind,drop=F],node_coef,edge_mean[,,eind,drop=F],node_mean,
                   elist_new,wt,lambda=lam,alpha_start,
                   sub,cores=cores)
    out[[n]]$edgelist<-elist_new
    out[[n]]$tree_obj<-tree_obj
    out[[n]]$eind<- eind
    out[[n]]$weight<- wt
    #node_coef<-out[[n]]$vert_coef
    #edge_coef[,,eind]<-out[[n]]$edge_coef
    #out[[n]]$edge_coef<-edge_coef
    out[[n]]$lambda<-i
    print(paste("DONE:",n))
  }
  return(out)
}


##############################
#Kruskal's algorithm
#
##############################
kruskal.min.span<-function (graph.obj) 
{
  #weights for Kruskal's algorithm in rand.weight
  #inputs a sparse adjacency matrix
  nv <- vcount(graph.obj)
  sum.mat<-get.edges(graph.obj,
                     E(graph.obj)
  )
  em <- t(sum.mat)
  ne <- dim(sum.mat)[1]
  eW <- E(graph.obj)$rand.weight
  ans <- .Call(
    "BGL_KMST_D", 
    as.integer(nv), 
    as.integer(ne), 
    as.integer(em - 1), 
    as.double(eW), 
    PACKAGE = "RBGL"
  )
  adj.out<-sparseMatrix(
    i=c(ans[[1]][1,],ans[[1]][2,])+1,
    j=c(ans[[1]][2,],ans[[1]][1,])+1,x=1,dims=c(nv,nv))
  
  #adj.out<-forceSymmetric(adj.out,uplo="L")
  
  g<-graph.adjacency(adjmatrix=adj.out)
  
  E(g)$weight<-1
  g
}


##############################
#Randomly generates a
#collection of spanning trees
#to create edge weights
##############################
edge.Start<-function(graph.obj,burn.in=5){
  #adj.Mat<-forceSymmetric(adj.Mat,uplo="L")
  g<-get.edgelist(graph.obj)
  
  dd<-vcount(graph.obj)
  ff<-ecount(graph.obj)
  ret<-sparseMatrix(i=c(),j=c(),dims=c(dd,dd),symmetric=T)
  ct<-0
  tree.list<-list()
  weight.list<-vector(length=0)
  while(TRUE){
    ct<-ct+1
    #edgeMatrix<-sparseMatrix(i=g[,1],j=g[,2],x=runif(ff,-1000,-999),dims=c(dd,dd),symmetric=T)
    #random edge weights
    E(graph.obj)$rand.weight<-lapply(E(graph.obj),function(e) runif(1,-1000,-999))
    
    new<-kruskal.min.span(graph.obj)[]
    tree.match=F
    if(length(tree.list)==0){
      #E(new)$weight<-1
      tree.list[[1]]<-new
      weight.list[1]<-1
      tree.match=T
    }else if(length(tree.list)>0){
      for(j in 1:length(tree.list)){
        if(all(tree.list[[j]][]==new[])){
          weight.list = weight.list*(ct-1)/ct
          weight.list[j] = weight.list[j] + 1/ct					
          tree.match=T
          break
        }
      }
    }
    if(tree.match==F){
      tree.list[[length(tree.list)+1]]<-new
      weight.list = weight.list*(ct-1)/ct
      weight.list[length(tree.list)]<-1/ct
    }
    
    ret<- ret + new
    if(all(sign(ret)==sign(graph.obj[])) & burn.in<=ct) break
  }
  E(graph.obj)$weight<-	summary(forceSymmetric(ret/ct))$x
  graph.obj$tree.list=tree.list
  graph.obj$weight.list=weight.list
  
  return(graph.obj)
}


##############################
#Add a new edge to graph correcting edge weights
#
##############################
edge.add<-function(graph.obj,edg){
  ge<-get.edges(graph.obj,edg)
  s=ge[1];t=ge[2]
  for(j in 1:length(graph.obj$tree.list)){
    g<-graph.adjacency(graph.obj$tree.list[[j]],mode="undirected",weighted=T)
    if(shortest.paths(g,s,t)==Inf){
      graph.obj$tree.list[[j]][s,t]=graph.obj$tree.list[[j]][t,s]=1
      E(graph.obj)[from(s) & to(t)]$weight<-E(graph.obj)[from(s) & to(t)]$weight+graph.obj$weight.list[[j]]
    }
    
  }
  return(graph.obj)
}


##############################
#Takes a single edge weight step
#
#
##############################
edge.weight.optim<-function(graph.obj,mutual_info,delta=.1){
  #Takes a single edge weight step
  
  E(graph.obj)$rand.weight<- -mutual_info
  
  tree.match=F
  new.tree<-	kruskal.min.span(graph.obj)
  for(j in 1:length(graph.obj$tree.list)){
    if(all(graph.obj$tree.list[[j]]==new.tree[])){
      graph.obj$weight.list = graph.obj$weight.list*(1-delta)
      graph.obj$weight.list[j] = graph.obj$weight.list[j] + delta				
      tree.match=T
      break
    }
  }
  
  if(tree.match==F){
    graph.obj$tree.list[[length(graph.obj$tree.list)+1]]<-new.tree[]
    graph.obj$weight.list = graph.obj$weight.list*(1-delta)
    graph.obj$weight.list[length(graph.obj$tree.list)]<-delta
  }
  new.weight<-graph.obj[]*(1-delta)+new.tree[]*delta
  E(graph.obj)$weight<-	E(graph.adjacency(new.weight,mode="undirected",weighted=T))$weight
  return(graph.obj)
  
}

##############################

new.edges<-function(enew, edge.coef,edge.mean, graph.obj)
{
  #For given edges, adds edge to current set of trees if
  #it doesn't create a cycle, then adds a new tree
  enew<- matrix(enew,ncol=2)
  #add to vertlist
  graph.obj<- add.edges(graph.obj, enew)
  #new edge means
  em.new<- edge.means(dat=data,vertlist=enew,m.set)
  #new edge coefs: all zero
  ec.new<- array(0,c(m.set,m.set,dim(enew)[1]))
  
  edge.coef<- abind(edge.coef, ec.new, along=1)
  edge.mean<- abind(edge.mean, em.new, along=1)
  
  #next: modify trees
  #E(graph.obj)$weight<-1
  #graph.obj<- edge.Start(graph.obj)
  E(graph.obj)$weight<-1
  if(is.null(graph.obj$tree.list)){
    graph.obj<-edge.Start(graph.obj)
  }else{
    lapply(
      1:dim(enew)[1],
      function(edg){
        lapply(
          1:length(graph.obj$tree.list),
          function(x){
            if( shortest.paths(
              graph.adjacency(graph.obj$tree.list[[x]]), 
              v=enew[edg,1],to=enew[edg,2]
            )==Inf
            )
            {
              graph.obj$tree.list[enew[edg,1],e[edg,2]] <<- 1 
              graph.obj$tree.list[enew[edg,2],e[edg,1]] <<- 1
            }
          }
        )
      }
    )
  }
  new<- Matrix(0,d,d)
  lapply(
    1:length(graph.obj$tree.list),
    function(x){
      new<<- new + graph.obj$tree.list[[x]]*graph.obj$weight.list[[x]]
    }
  )
 
  #add a random tree
  #not worked out yet
  #E(graph.obj)$rand.weight<-lapply(E(graph.obj),function(e) runif(1,-1000,-999))
  
  #new<-kruskal.min.span(graph.obj)[]
  
  
  E(graph.obj)$weight<- summary(forceSymmetric(new))$x
  
  return(list(edge.coef=edge.coef,edge.mean=edge.mean,graph.obj=graph.obj))
}



