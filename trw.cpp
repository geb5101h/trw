// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

class Coef{
  public:
    arma::mat vert_coef;
    arma::cube edge_coef;
    Coef(){};
    Coef(arma::mat&, arma::cube&);
    Coef(int d, int e, int m1, int m2);
    void proximal(arma::cube, double);
};
Coef::Coef(arma::mat& vc, arma::cube& ec){
  vert_coef = vc;
  edge_coef = ec;
}
Coef::Coef(int d, int e, int m1, int m2){
  vert_coef.zeros(d,m1);
  edge_coef.zeros(m2,m2,e);
}

class Mean{
  public:
  arma::mat vert_mean;
  arma::cube edge_mean;
  Mean(){};
  Mean(arma::mat, arma::cube);
};

Mean::Mean(arma::mat vm, arma::cube em){
  vert_mean = vm;
  edge_mean = em;
}

struct Gradient{
  public:
  arma::cube edge_gradient;
  arma::mat vert_gradient;
  arma::mat vert_quasi;
  arma::cube edge_quasi;
};

class Density{
  public:
  arma::mat beliefs;
  arma::mat messages;
  arma::cube edge_den;
  arma::vec mutual_info;
  arma::vec entropy;
  arma::mat messages_weighted;
  Density(){};
  Density(int d, int e, int res){
    beliefs.ones(d,res); 
    messages.ones(2*e,res);
    messages_weighted.ones(2*e,res);
    edge_den.ones(res,res,e);
    entropy.zeros(d);
    mutual_info.zeros(e);
    
    };
  Density(arma::mat vert_coef, arma::mat legen_vec);
  void beliefprop(
    const Coef&, 
    const arma::vec&, 
    const arma::mat&, 
    const arma::mat&, 
    const arma::cube&,
    int& dc_flag);
};


struct Nll_obj{
  public:
  double nll;
  double part_ret;
  Density density;
  Gradient gradient;
  Nll_obj(){}
};

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Bp : public Worker
{
   arma::mat& messages;
   arma::mat& messages_weighted;
   arma::mat& mesprod;
   arma::mat& messages_old;
   arma::cube& ten_save;
   arma::mat& edgelist2;
   arma::vec& weights2;
   
   Bp(arma::mat& messages, 
      arma::mat& messages_weighted, 
      arma::mat& mesprod, 
      arma::mat& messages_old,
      arma::cube& ten_save, 
      arma::mat& edgelist2, 
      arma::vec& weights2) : 
        messages(messages), 
        messages_weighted(messages_weighted),
        mesprod(mesprod), 
        messages_old(messages_old), 
        ten_save(ten_save),
        edgelist2(edgelist2), 
        weights2(weights2) {}
   
   void operator()(std::size_t begin, std::size_t end) {
     int s,t, flip;
     int e = messages.n_rows/2;
     for(int i=begin; i<end; ++i){
          s = edgelist2(i,0)-1;
          t = edgelist2(i,1)-1;
        if(i<e){
          flip = i+e;
          messages.row(i) = ( mesprod.row(s)/messages_old.row(flip) )  * ten_save.slice(i) ;
        }else{
          flip = i-e;
          messages.row(i) = ( mesprod.row(s)/messages_old.row(flip) )  * ten_save.slice(i-e).t() ;

        }
      messages.row(i) /= max(messages.row(i));
      //damping
      //messages.row(i) = messages.row(i)*.7 + messages_old.row(i)*.3;
      messages_weighted.row(i) = pow(messages.row(i) , weights2(i) );    
      }
   }
};

struct EdgeDenCreate : public Worker
{
   arma::cube& edge_den;
   const arma::mat& edgelist;
   arma::mat& beliefs;
   arma::mat& messages;
   arma::vec& mutual_info;
   arma::vec& entropy;
   arma::cube& ten_save;
   EdgeDenCreate(
      arma::cube& edge_den,
      const arma::mat& edgelist,
      arma::mat& beliefs,
      arma::mat& messages,
      arma::vec& mutual_info,
      arma::vec& entropy,
      arma::cube& ten_save) : 
        edge_den(edge_den),
        edgelist(edgelist),
        beliefs(beliefs),
        messages(messages),
        mutual_info(mutual_info),
        entropy(entropy),
        ten_save(ten_save) {}
   
   void operator()(std::size_t begin, std::size_t end) {
      int res = ten_save.n_cols;
   		int e = edgelist.n_rows;
      for(int i=begin; i<end; ++i){
    	  int s = edgelist(i,0)-1;
    	  int t = edgelist(i,1)-1;
   		edge_den.slice(i) = ten_save.slice(i);
        edge_den.slice(i) %= ( (beliefs.row(s)/messages.row(i+e)).t() * (beliefs.row(t)/messages.row(i)) ) ;
        edge_den.slice(i) /= accu(edge_den.slice(i))/(res*res);
        mutual_info(i) = accu( edge_den.slice(i)%log(edge_den.slice(i)) )/(res*res) + entropy(s) + entropy(t);
  		}
   }
};

struct tenHelper : public Worker
{
   const arma::cube& edge_coef;
   const arma::cube& legen;
   arma::vec weights;
   arma::cube& ten_out;

   tenHelper(
      const arma::cube& edge_coef,
      const arma::cube& legen,
      arma::vec weights,
      arma::cube& ten_out) :
        edge_coef(edge_coef),
        legen(legen),
        weights(weights),
        ten_out(ten_out) {}
   
   void operator()(std::size_t begin, std::size_t end) {
      int m2 = edge_coef.n_rows;
      int k;
      for(int ee=begin; ee< end ; ++ee ){
        k=0;
        for(int i=0;i<m2;i++){
          for(int j=0;j<m2;j++){
            ten_out.slice(ee) = ten_out.slice(ee) +  edge_coef(i,j,ee) * legen.slice(k);
            k = k+1;
         }
        }
        ten_out.slice(ee) /= weights(ee);
        ten_out.slice(ee) = exp(ten_out.slice(ee));
        int res = ten_out.n_cols;
        //ten_out.slice(ee) /= accu(ten_out.slice(ee))/(res*res);
       }
   }
};


/****************
Main function for conducting belief propagation
*
*
****************/
void Density::beliefprop(
    const Coef& coef, 
    const arma::vec& weights, 
    const arma::mat& edgelist, 
    const arma::mat& legen_vec, 
    const arma::cube& legen_array,
    int& dc_flag){
  int d = beliefs.n_rows;
  arma::mat node_pot = exp(coef.vert_coef * legen_vec );
  if(edgelist.n_rows == 0){
   for(int i=0; i<d; ++i){
      beliefs.row(i) = node_pot.row(i);
      beliefs.row(i) /= mean( node_pot.row(i));
    }
		entropy = - mean( beliefs%log(beliefs), 1);
        return;
   }
  int res = legen_array.n_cols;
  int nedge = edgelist.n_rows;
  arma::cube ten_save(res,res,nedge);
  ten_save.zeros();
  tenHelper th(
      coef.edge_coef,
      legen_array,
      weights,
      ten_save);
  parallelFor(0, nedge, th);
  
  /*pass messages until convergence
  *
  *
  */
  arma::uvec neighbors;
  double error;
  int drop=0;
  arma::vec weights2(join_cols(weights,weights));
  arma::mat edgelist2( join_cols(edgelist, fliplr(edgelist)) );
  
  messages_weighted.ones(2*nedge,res);
  messages.ones(2*nedge,res);
  arma::mat mesprod(d,res);
  
  arma::mat messages_old(2*nedge,res);
  arma::mat messages_weighted_old(2*nedge,res);
  while(true){
   drop++;
    if(drop>1000) 
    {
      Rcpp::Rcout << "didn't converge" << std::endl;
      dc_flag=1;
      break;
    }
    
    messages_old = messages;
    messages_weighted_old = messages_weighted;
    
    for(int t=0; t<d; ++t){
      neighbors = find(edgelist2.col(1) == t+1 );
      mesprod.row(t) = node_pot.row(t);
      mesprod.row(t) %= prod( messages_weighted_old.rows(neighbors) , 0);
    }
   
   Bp bp(messages,
        messages_weighted, 
        mesprod,
        messages_old,
        ten_save, 
        edgelist2, 
        weights2);
    parallelFor(0, messages.n_rows, bp);

	//messages = .2*messages + .8*messages_old;
    arma::vec inter =  mean( abs(messages - messages_old) , 1);
    error = inter.max();
    if( error <= 1.0e-6 ){
      //if(true){
      break;
    }
  }
    for(int t=0; t<d; ++t){
      neighbors = find(edgelist2.col(1) == t+1 );
      beliefs.row(t) = node_pot.row(t) % prod( messages_weighted.rows( neighbors ) , 0) ;
      beliefs.row(t) /= mean(beliefs.row(t));
    }
  entropy = - mean( beliefs%log(beliefs), 1);
 EdgeDenCreate edc(edge_den,
   edgelist,
   beliefs,
   messages,
   mutual_info,
   entropy,
   ten_save);
  parallelFor(0, edge_den.n_slices, edc);

// for(int i=0; i<e; ++i){
//     int s = edgelist(i,0)-1;
//     int t = edgelist(i,1)-1;
//     edge_den.slice(i) = ten_save.slice(i);
//     edge_den.slice(i) %= ( (beliefs.row(s)/messages.row(i+e)).t() * (beliefs.row(t)/messages.row(i)) ) ;
//     edge_den.slice(i) /= accu(edge_den.slice(i))/(res*res);
//     mutual_info(i) = accu( edge_den.slice(i)%log(edge_den.slice(i)) )/(res*res) + entropy(s) + entropy(t);
//   }

}

/****************
For given density and mean parameters
computes quasi-marginals and gradient
*****************/
// // [[Rcpp::export]]
Gradient get_quasi(
    const Density& density,   
    const Mean& mean, 
    const arma::mat& legen_vec){
  Gradient gradient;
  int e = mean.edge_mean.n_slices;
  int m2 = mean.edge_mean.n_cols;
  int res = density.beliefs.n_cols;
  gradient.vert_quasi = ( density.beliefs * legen_vec.t() )/res ;
  gradient.vert_gradient = gradient.vert_quasi; 
  gradient.vert_gradient -= mean.vert_mean;
  if(e == 0){
    return gradient;
  }
  
  arma::cube edge_quasi(m2,m2,e);
  for(int i=0; i<e; ++i){
    edge_quasi.slice(i) = ( legen_vec * density.edge_den.slice(i) * legen_vec.t() ) / (res*res);
}
	gradient.edge_quasi = edge_quasi;
  gradient.edge_gradient = gradient.edge_quasi; 
  gradient.edge_gradient -= mean.edge_mean;

  
  return gradient;
}


/***************
Computes negative log likelihood

****************/
// // [[Rcpp::export]]

Nll_obj neg_ll(
    const Coef& coef,
    const Mean& mean,  
    const arma::vec& weights,
    const arma::mat& edgelist,
    const arma::mat& legen_vec,
    const arma::cube& legen_array,
    int& dc_flag
    
    ){
 
      //IntegerVector edims = edge_coef.attr("dim");
      //arma::cube ec(edge_coef.begin(),edims[0],edims[1],edims[2],false);
      int e = edgelist.n_rows;
      Density density(mean.vert_mean.n_rows, e, legen_vec.n_cols);
      density.beliefprop(coef, weights,edgelist,legen_vec,legen_array,dc_flag);
      Nll_obj nll_obj;
      nll_obj.gradient = get_quasi(density, mean,legen_vec);
      nll_obj.density = density;
      double entropysum = accu( density.entropy );
      if(e==0){
      	nll_obj.part_ret = accu(coef.vert_coef % nll_obj.gradient.vert_quasi) + entropysum;
      	nll_obj.nll = nll_obj.part_ret - accu(coef.vert_coef % mean.vert_mean);
      }else{
      	double misum = accu(density.mutual_info % weights);

      	nll_obj.part_ret = accu(coef.vert_coef % nll_obj.gradient.vert_quasi) 
      						+ accu(coef.edge_coef % nll_obj.gradient.edge_quasi) 
      						+ entropysum-misum;
      	nll_obj.nll = nll_obj.part_ret - accu(coef.vert_coef % mean.vert_mean) - accu(coef.edge_coef % mean.edge_mean) ;
      }
      // 
//       double entropysum = accu( density.entropy );
//       //double vert_crossterm = accu(gradient.vert_quasi % coef.vert_coef);
//       double vert_crossterm = std::inner_product(gradient.vert_quasi.begin(), gradient.vert_quasi.end(), coef.vert_coef.begin(), 0.0);
//       nll_obj.part_ret = vert_crossterm + entropysum;
//       double vert_crossterm2 = std::inner_product(mean.vert_mean.begin(), mean.vert_mean.end(), coef.vert_coef.begin(), 0.0);
//       nll_obj.nll = nll_obj.part_ret- vert_crossterm2;
//       if(e==0){
//         return nll_obj;
//       }
//       double misum = accu( density.mutual_info % weights );
//       //double misum = std::inner_product(density.mutual_info.begin(), density.mutual_info.end(), weights.begin(), 0.0);
//       double edge_crossterm = accu(gradient.edge_quasi % coef.edge_coef);
//       //double edge_crossterm = std::inner_product(gradient.edge_quasi.begin(), gradient.edge_quasi.end(), coef.edge_coef.begin(), 0.0);
//       //nll_obj.part_ret += edge_crossterm - misum;
//       
// 
//       //Rcpp::Rcout << part_ret << std::endl;
//       //double part_fn = partition_fn(ent,mi,vert_coef,edge_coef,vert_quasi,edge_quasi)
//      
//       //IntegerVector emdims = edge_mean.attr("dim");
//       //arma::cube em(edge_mean.begin(),emdims[0],emdims[1],emdims[2],false);
//       //double vert_crossterm2 = accu( mean.vert_mean % coef.vert_coef);
//       double edge_crossterm2 = accu(mean.edge_mean % coef.edge_coef);
//       //double edge_crossterm2 = std::inner_product(mean.edge_mean.begin(), mean.edge_mean.end(), coef.edge_coef.begin(), 0.0);
//        //Rcpp::Rcout << "ec1: " << edge_crossterm << " ec2: " << edge_crossterm2 << std::endl
//      // << " vc1: " << vert_crossterm << " vc2: " << vert_crossterm2 << std::endl;
//       nll_obj.part_ret = vert_crossterm + edge_crossterm + entropysum - misum;
//       nll_obj.nll = nll_obj.part_ret - vert_crossterm2 - edge_crossterm2 ;

      return nll_obj;
    }

/*
Backtracking line search for ISTA

*/
double backtrack(
  const Coef& coef_new, 
  const Coef& coef_old, 
  const Nll_obj& nll_new, 
  const Nll_obj& nll_old, 
  double alpha,
  double& alpha_start){
  int e = coef_new.edge_coef.n_slices;
  //double armijo;
   //  if(e==0){
// 	armijo = (nll_new.nll - nll_old.nll
//   				-  arma::accu( (coef_new.vert_coef - coef_old.vert_coef) % nll_old.gradient.vert_gradient)
//   				-  ( 1/(2*alpha) ) * arma::accu( (coef_new.vert_coef-coef_old.vert_coef)%(coef_new.vert_coef-coef_old.vert_coef) ));
//   }else{
//   armijo = (nll_new.nll - nll_old.nll
//   				-  arma::accu( (coef_new.vert_coef - coef_old.vert_coef) % nll_old.gradient.vert_gradient)
//   				-  arma::accu( (coef_new.edge_coef - coef_old.edge_coef) % nll_old.gradient.edge_gradient)
//   				-  ( 1/(2*alpha) ) * arma::accu( (coef_new.vert_coef-coef_old.vert_coef)%(coef_new.vert_coef-coef_old.vert_coef) )
//   				-  ( 1/(2*alpha) ) * arma::accu( (coef_new.edge_coef-coef_old.edge_coef)%(coef_new.edge_coef-coef_old.edge_coef) ));
//   double sec = arma::accu( (coef_new.edge_coef - coef_old.edge_coef) % (nll_new.gradient.edge_gradient - nll_old.gradient.edge_gradient) ) ;
//   sec +=  arma::accu( (coef_new.vert_coef - coef_old.vert_coef) % (nll_new.gradient.vert_gradient - nll_old.gradient.vert_gradient) );
//   sec /= arma::accu( (coef_new.vert_coef - coef_old.vert_coef) % (coef_new.vert_coef - coef_old.vert_coef)  );
//   //alpha_start = sec;
//   }
  
  Coef coef_diff;
  //Gradient grad_diff;
  Nll_obj nll_diff;
  
  
  coef_diff.vert_coef = coef_new.vert_coef; 
  coef_diff.vert_coef -= coef_old.vert_coef;
  
  nll_diff.gradient.vert_gradient = nll_new.gradient.vert_gradient; 
  nll_diff.gradient.vert_gradient -= nll_old.gradient.vert_gradient;
    
  nll_diff.nll = nll_new.nll; 
  nll_diff.nll -= nll_old.nll;
  //double a = accu(coef_diff.edge_coef % coef_diff.edge_coef);
  //double b = accu(coef_diff.vert_coef % coef_diff.vert_coef);
  //double c = accu(coef_diff.vert_coef % (nll_diff.gradient.vert_gradient));
  //double d =accu(coef_diff.edge_coef % (nll_diff.gradient.edge_gradient));
  double b = std::inner_product(coef_diff.vert_coef.begin(), coef_diff.vert_coef.end(), coef_diff.vert_coef.begin(), 0.0);
  double c = std::inner_product(coef_diff.vert_coef.begin(), coef_diff.vert_coef.end(), nll_diff.gradient.vert_gradient.begin(), 0.0);
  double armijo = nll_diff.nll;
  armijo -= b/(2*alpha);
  armijo -= std::inner_product(nll_old.gradient.vert_gradient.begin(), nll_old.gradient.vert_gradient.end(), coef_diff.vert_coef.begin(), 0.0);
  if(e==0){
    //alpha_start = b/c;
    return armijo;
  }
  coef_diff.edge_coef = coef_new.edge_coef; 
  coef_diff.edge_coef -= coef_old.edge_coef;
  nll_diff.gradient.edge_gradient = nll_new.gradient.edge_gradient; 
  nll_diff.gradient.edge_gradient -= nll_old.gradient.edge_gradient;
  double a = std::inner_product(coef_diff.edge_coef.begin(), coef_diff.edge_coef.end(), coef_diff.edge_coef.begin(), 0.0);
  double d = std::inner_product(coef_diff.edge_coef.begin(), coef_diff.edge_coef.end(), nll_diff.gradient.edge_gradient.begin(), 0.0);;
   
  armijo -= a/(2*alpha);
  //armijo -= accu(nll_old.gradient.edge_gradient % coef_diff.edge_coef);
  armijo -= std::inner_product(nll_old.gradient.edge_gradient.begin(), nll_old.gradient.edge_gradient.end(), coef_diff.edge_coef.begin(), 0.0);
  //armijo -= accu(nll_old.gradient.vert_gradient % coef_diff.vert_coef);
  //alpha_start= (a+b)/(c+d);

  return armijo;
};

double softthresh(double lambda, double norm){
  if(lambda >= norm){
      return 0.0;
  }else{
      return 1.0 - lambda/norm;
  }
};

/*
Soft thresholding for proximal step
*/
void Coef::proximal(arma::cube ec, double lambda){
  double nrm;
  int e = ec.n_slices;
  if(e == 0){
    return;
  }
  for(int i=0; i<e; ++i){
    nrm = sqrt(arma::accu(ec.slice(i) % ec.slice(i)));
    //nrm = sqrt(std::inner_product(edge_coef_old.slice(i).begin(),edge_coef_old.slice(i).end(), edge_coef_old.slice(i).begin(), 0.0));
    edge_coef.slice(i) = softthresh(lambda, nrm) * ec.slice(i);
  }

};


double l1_norm(arma::cube& edge_coef){
  int e = edge_coef.n_slices;
  if(e==0){
    return 0.0;
  }
  double nrm = 0;
  for(int i=0; i<e; ++i){
    //nrm += sqrt(accu( edge_coef.slice(i) % edge_coef.slice(i)));
    nrm += sqrt(std::inner_product(edge_coef.slice(i).begin(), edge_coef.slice(i).end(),edge_coef.slice(i).begin(),0.0));
  }
  return nrm;
}
Nll_obj ista(Coef& coef,
          const Mean& mean, 
          const arma::vec& weights,
          const arma::mat& edgelist,
          const arma::mat& legen_vec,
          const arma::cube& legen_array,
          const double lambda,
          double alpha_start
          ){
  Coef coef_old = coef;
  int dc_flag=0;
  Nll_obj nll;
  Nll_obj nll_old = neg_ll(
        coef_old,
        mean,
        weights,
        edgelist,
        legen_vec,
        legen_array,
        dc_flag
      );
  int e = coef_old.edge_coef.n_slices;    
  double nrm_old = l1_norm(coef_old.edge_coef);
  double nrm,bt;
  int iter=0;
  
  while(true){
    double alpha = alpha_start;
    iter++;
    while(true){
      //Rcpp::Rcout << "alpha: " << alpha << std::endl;
      coef.vert_coef = coef_old.vert_coef; 
      coef.vert_coef -= alpha*nll_old.gradient.vert_gradient;
      if(e != 0){
        coef.proximal(coef_old.edge_coef - alpha*nll_old.gradient.edge_gradient ,lambda*alpha);
      }
      dc_flag=0;
      nll=neg_ll(
        coef,
        mean,
        weights,
        edgelist,
        legen_vec,
        legen_array,
        dc_flag
      );
      //return(nll);
      if(dc_flag==1){
      	alpha=alpha*.5;
      	//alpha_start=alpha_start*.8;
      	continue;
      }
      nrm= l1_norm(coef.edge_coef);
      bt = backtrack(coef, coef_old, nll, nll_old, alpha, alpha_start);
      //Rcpp::Rcout << " bt " <<  bt << std::endl;
      if(bt<=0.0){
        break;
      }
      alpha = alpha*0.5;
      
    }
    //Rcpp::Rcout << "pen nll : "<< nll.nll + lambda*nrm << std::endl;
     Rcpp::Rcout << "iteration #: " << iter << " pen nll: " << nll.nll + lambda*nrm << std::endl;

    //Rcpp::Rcout << "crit: "<< std::abs(nll_old.nll - nll.nll + nrm_old - nrm) << std::endl;
    if(std::abs(nll_old.nll - nll.nll + lambda*(nrm_old - nrm) ) < (1+e)*1e-4 | iter>1000){
      break;
    }
    coef_old = coef;
    nrm_old = nrm;
    nll_old = nll;
  }
	
  return(nll);
};
/*
Nll_obj fista(Coef& coef,
          double& part_fn,
          const Mean& mean, 
          const arma::vec& weights,
          const arma::mat& edgelist,
          const arma::mat& legen_vec,
          const arma::cube& legen_array,
          const double lambda,
          double alpha_start
          ){
  Coef coef_old = coef;
  Coef coef_y = coef;
  Nll_obj nll;
  Nll_obj nll_old = neg_ll(
        coef_y,
        mean,
        weights,
        edgelist,
        legen_vec,
        legen_array
      );
  int e = coef_old.edge_coef.n_slices;    
  double nrm_old = l1_norm(coef_old.edge_coef);
  //Rcpp::Rcout << nll_old.nll << std::endl;
  double nrm,bt;
  int iter=0;
  double t_old = 1;
  while(true){
    double alpha = alpha_start;
    iter++;
    //Rcpp::Rcout << "iteration #: " << iter << std::endl;
    while(true){
      //Rcpp::Rcout << "alpha" << alpha << std::endl;
      coef.vert_coef = coef_y.vert_coef; 
      coef.vert_coef -= alpha*nll_old.gradient.vert_gradient;
      if(e != 0){
        coef.proximal(coef_y.edge_coef - alpha*nll_old.gradient.edge_gradient ,lambda*alpha);
      }
      nll=neg_ll(
        coef,
        mean,
        weights,
        edgelist,
        legen_vec,
        legen_array
      );
      nrm= l1_norm(coef.edge_coef);
      bt = backtrack(coef, coef_y, nll, nll_old, alpha, alpha_start);

      if(bt<=0.0){
      //if(true){
        break;
      }
      alpha = alpha*0.5;
      
    }

    double t_new = (1+std::sqrt(1+4*t_old*t_old))/2;
    coef_y.edge_coef = coef.edge_coef + (t_old+1)/(t_new)*(coef.edge_coef-coef_old.edge_coef);
    coef_y.vert_coef = coef.vert_coef + (t_old+1)/(t_new)*(coef.vert_coef-coef_old.vert_coef);


    t_old=t_new;

    //Rcpp::Rcout << "pen nll : "<< nll.nll + lambda*nrm << std::endl;

    //Rcpp::Rcout << "crit: "<< std::abs(nll_old.nll - nll.nll + nrm_old - nrm) << std::endl;

    if(std::abs(nll_old.nll - nll.nll + lambda*(nrm_old - nrm) ) < (e+1)*1e-6 | iter>500){
      break;
    }
    
    coef_old = coef;
    nrm_old = nrm;
    nll_old=neg_ll(
        coef_y,
        mean,
        weights,
        edgelist,
        legen_vec,
        legen_array
      );
    //nll_old = nll;
  }
  part_fn = nll.part_ret;
  return(nll);
};
*/
// [[Rcpp::export]]
List rwrapper(
  const NumericVector& edge_coef,
  const arma::mat& vert_coef,
  const NumericVector& edge_mean,
  const arma::mat& vert_mean,
  const arma::mat& edgelist,
  const arma::vec& weights,
  const arma::mat& legen_vec,
  const NumericVector& legen_array,
  double lambda,
  double alpha_start,
  int alg=0
  )
{
  IntegerVector ldims = legen_array.attr("dim");
  arma::cube la(legen_array.begin(),ldims[0],ldims[1],ldims[2],false);
  IntegerVector ecdims = edge_coef.attr("dim");
  arma::cube ec(edge_coef.begin(),ecdims[0],ecdims[1],ecdims[2],false);
  IntegerVector emdims = edge_mean.attr("dim");
  arma::cube em(edge_mean.begin(),emdims[0],emdims[1],emdims[2],false);
  //int res = 128;  
  int d = vert_mean.n_rows;
  Coef coef(d, em.n_slices, em.n_rows, vert_mean.n_cols);
  Mean mean(vert_mean, em);
  Nll_obj nll;
  double part_fn=0;
  nll= ista(coef, mean, weights, edgelist,legen_vec,la, lambda,alpha_start);
  Gradient gr = get_quasi(nll.density,mean,legen_vec);

  List list;
  list["beliefs"]= nll.density.beliefs;
  list["edge_den"]= nll.density.edge_den;
  list["vert_coef"]= coef.vert_coef;
  list["edge_coef"]= coef.edge_coef;
  list["part_fn"]=nll.part_ret;
  list["nll"]=nll.nll;
  list["mutual_info"]=nll.density.mutual_info;
  list["edge_gradient"]=gr.edge_gradient;
  list["vert_gradient"]=gr.vert_gradient;
  list["edge_quasi"]=gr.edge_quasi;
  list["vert_quasi"]=gr.vert_quasi;
  list["mes"] = nll.density.messages;
  return list;
}

// [[Rcpp::export]]
List rwrapper2(
  const NumericVector& edge_coef,
  const arma::mat& vert_coef,
  const NumericVector& edge_mean,
  const arma::mat& vert_mean,
  const arma::mat& edgelist,
  const arma::vec& weights,
  const arma::mat& legen_vec,
  const NumericVector& legen_array,
  double lambda,
  double alpha_start,
  int alg=0
  )
{
  IntegerVector ldims = legen_array.attr("dim");
  arma::cube la(legen_array.begin(),ldims[0],ldims[1],ldims[2],false);
  IntegerVector ecdims = edge_coef.attr("dim");
  arma::cube ec(edge_coef.begin(),ecdims[0],ecdims[1],ecdims[2],false);
   IntegerVector emdims = edge_mean.attr("dim");
  arma::cube em(edge_mean.begin(),emdims[0],emdims[1],emdims[2],false);

  int dc_flag=0;
  //int res = 128;  
  int d = vert_coef.n_rows;
  Mean mean(vert_mean, em);
  Coef coef(d, ec.n_slices, ec.n_rows, vert_coef.n_cols);
  
   coef.edge_coef = ec;
   coef.vert_coef = vert_coef;
   Density density(vert_coef.n_rows, edgelist.n_rows, legen_vec.n_cols);
   density.beliefprop(coef, weights,edgelist,legen_vec,la,dc_flag);
   Nll_obj nll=neg_ll(
    coef,
    mean,  
    weights,
    edgelist,
    legen_vec,
    la,
    dc_flag);
  List list;
  list["beliefs"]= density.beliefs;
  list["edge_den"]= density.edge_den;
  list["vert_coef"]= coef.vert_coef;
  list["edge_coef"]= coef.edge_coef;
  list["mutual_info"]=density.mutual_info;
  list["mes"] = density.messages;
  list["part_fn"]=nll.part_ret;
  return list;
}


