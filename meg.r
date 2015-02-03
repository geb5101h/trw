rm(list=ls())
if(Sys.getenv("USER")=="ebjan"){
  setwd("~/trw")
}else{
  setwd("~/Dropbox/trw")
}

#m=3
m1=3
m2=3
source("startup.r")
source("trw_fit.r")
####################

####################
dat<-matrix( scan(file="MEG_art",what=double()) , ncol=122, byrow=F)
#save(file="meg.rdata",list=c("dat"))
#source("meg.rdata")
dat<-scale(dat,center=TRUE,scale=TRUE)
#dat[which(abs(dat)>6)]<- sign(dat[which(abs(dat)>6)])*6
data2<-dat
data2<- apply(data2,2,function(x)(x-min(x)+.001)/(max(x)-min(x)+.002))
train<-data2[2000:3000,]
test<-data2[3001:3500,]

d<-dim(train)[2]
####################

#MEG
#


cc<-cov(train)
rho_start<- max( abs(cc-diag(diag(cc))) )



run_seq <- exp( seq(from=log(rho_start*.2), to =log(rho_start*1.1),length.out=50 ) )
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




edge_mean <- edge.means(train, edgelist, m2, cores)
node_mean <- node.means(train , m1,cores)
edge_mean_test <- edge.means(test, edgelist, m2, cores)
node_mean_test <- node.means(test , m1,cores)
weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])


edge_coef<- array(0,dim=c(m2,m2,dim(edgelist)[1]))
node_coef<- matrix(0,d,m1)




lambda=.1
alpha_start=1
limit=30
out_len=10
path<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
               edgelist,weights=NULL,alpha_start=1,
               g,cores=cores,limit=limit,relax=F,edge.opt=T,out_len=out_len)
path_relax<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
                     edgelist,weights=NULL,alpha_start=1,
                     g,cores=cores,limit=limit,relax=T,edge.opt=T,out_len=out_len)




risk_path<-lapply(path,function(x){
  -trw_lik(node_mean_test,edge_mean_test[,,x$eind,drop=F],x$vert_coef,x$edge_coef,x$part_fn)
})
risk_path_relax<-lapply(path_relax,function(x){
  -trw_lik(node_mean_test,edge_mean_test[,,x$eind,drop=F],x$vert_coef,x$edge_coef,x$part_fn)
})

risk_path<-simplify2array(risk_path)
risk_path_relax<-simplify2array(risk_path_relax)


seqlm<-  seqlp<-round( seq(from=1,to=limit,length.out=out_len) )

setwd("./Experiments")


save.image("meg.rdata")



pdf("path_meg.pdf")
nf<- layout(matrix(c(1,1,1,2,3,4),nrow=2,ncol=3,byrow=T), respect = T,heights=c(8,6),widths=rep(4,3))

plot(seqlm,risk_path,type="l",ylim=range(risk_path,risk_path_relax,test_error,test_error_relax))
points(seqlm,risk_path_relax,type="l",col="blue")
#points(seqlm,risk_path_tree,type="l",col="orange")

points(num_nonzero,test_error_relax,col="green",type="l")
points(num_nonzero,test_error,col="purple",type="l")
legend("topright",c("trw","trw relax","glasso relax", "glasso"),
       pch=16,col=c("black","blue","green","purple"))
#dev.off()

########################################
#g1<- graph.adjacency(h$theta,mode="undirected")
g1<- graph.empty(122)
g2<- graph.edgelist(path_relax[[which.min(risk_path_relax)]]$edgelist,directed=F)
g3<-graph.adjacency(1*(glasso_out[,,33]!=0)-diag(d),
                mode="undirected")
l<-layout.auto(g1)


g12<- graph.union(g1,g2)
g13<- graph.union(g1,g3)

el1 <- apply(get.edgelist(g1), 1, paste, collapse="-")
el2 <- apply(get.edgelist(g2), 1, paste, collapse="-")
el3 <- apply(get.edgelist(g3), 1, paste, collapse="-")
el12 <- apply(get.edgelist(g12), 1, paste, collapse="-")
el13 <- apply(get.edgelist(g13), 1, paste, collapse="-")
E(g1)$color= "black"
E(g12)$color <- ifelse( !(el12 %in% el1)  , "red",
		ifelse(!(el12 %in% el2), "gray", "black") )
E(g13)$color <- ifelse( !(el13 %in% el1)  , "red",
		ifelse(!(el13 %in% el3), "gray", "black") )
edge.lty12 <- ifelse( (el12 %in% el2), "solid","solid")
edge.lty13 <- ifelse( (el13 %in% el3), "solid","solid")


#pdf("graph_gaussian.pdf")
#par(mfcol=c(1,3))
par(mar=rep(.5,4))  
# plot.igraph(g1,vertex.size=5,
#      layout=l,
#      vertex.label='',
#      edge.color="black",
#      edge.lty=1,
#      vertex.color="white");
#      title("True",line=-2) 
plot.igraph(g2,vertex.size=5,
     layout=l,
     vertex.label='',
     edge.lty=edge.lty12,
     vertex.color="white");
    title("TRW",line=-2) 
 
plot.igraph(g3,vertex.size=5,
     layout=l,
     vertex.label='',
     edge.lty=edge.lty13,
     vertex.color="white"
	)
          title("Glasso",line=-2) 

dev.off()


# g12<- graph.union(g1,g2)
# g13<- graph.union(g1,g3)
# 
# el1 <- apply(get.edgelist(g1), 1, paste, collapse="-")
# el2 <- apply(get.edgelist(g2), 1, paste, collapse="-")
# el3 <- apply(get.edgelist(g3), 1, paste, collapse="-")
# el12 <- apply(get.edgelist(g12), 1, paste, collapse="-")
# el13 <- apply(get.edgelist(g13), 1, paste, collapse="-")
# E(g1)$color= "black"
# E(g12)$color <- ifelse( !(el12 %in% el1)  , "red",
# 		ifelse(!(el12 %in% el2), "gray", "black") )
# E(g13)$color <- ifelse( !(el13 %in% el1)  , "red",
# 		ifelse(!(el13 %in% el3), "gray", "black") )
# edge.lty12 <- ifelse( (el12 %in% el2), "solid","solid")
# edge.lty13 <- ifelse( (el13 %in% el3), "solid","solid")
# 
# 
# pdf("graph_meg.pdf")
# par(mfcol=c(1,3))
# plot(g1,vertex.size=2,
#      layout=l,
#      vertex.label='',
#      edge.color="black",
#      edge.lty=1);
# plot(g12,vertex.size=2,
#      layout=l,
#      vertex.label='',
#      edge.lty=edge.lty12);
# plot(g13,vertex.size=2,
#      layout=l,
#      vertex.label='',
#      edge.lty=edge.lty13);
# dev.off()
# 
# 
# pdf("mi_weightng_meg.pdf")
# plot(path[[which.min(risk_path)]]$mutual_info,path[[which.min(risk_path)]]$weight,pch=16)
# dev.off()
