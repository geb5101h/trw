############
#Calculation of some likelihoods
#
############

gaussian_lik<-function(x,Omega){
  #Omega is precision matrix (inverse of covariance mat)
  d<-length(x)
  -log(2*pi)*d/2 + log( det(Omega) )/2 - crossprod(x,Omega)%*%x /2
}

trw_lik<-function(vert_means,edge_means,vert_coef,edge_coef,part_fn){
  sum( vert_means*vert_coef) + sum( edge_means*edge_coef ) - part_fn
}