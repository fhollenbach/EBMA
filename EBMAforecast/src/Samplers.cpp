#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector oneMultinomCalt(NumericVector probs) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  return(ans);
}

// [[Rcpp::export]]
NumericVector getRGamma(double shape) { 
  RNGScope scope;
  NumericVector x = rgamma(1, shape, 1 );
  return x;
}

// [[Rcpp::export]]
LogicalVector isNA(NumericVector x) {
  int n = x.size();
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::List emNorm(Rcpp::NumericVector outcome,
                  Rcpp::NumericMatrix prediction, Rcpp::NumericMatrix RSQ, Rcpp::NumericVector W, double tol, int maxIter
                    ,double wisdom, double sigma2) {
  double LL = 0.00;
  double LL_old = 0.00;
  int leng = prediction.nrow();
  int width = prediction.ncol();
  int iter = 0 ;
  double improv = 1.1;
  
  while((improv > tol) & (iter < maxIter)){
    Rcpp::NumericMatrix  znumerator(leng,width);
    Rcpp::NumericVector w_old(W);
    Rcpp::NumericVector zdenom(leng);
    Rcpp::NumericVector unnormalizedW(width);
    Rcpp::NumericMatrix Zs(leng,width); 
    Rcpp::NumericMatrix RZ(leng,width); 
    Rcpp::NumericMatrix  g(leng,width);
    Rcpp::NumericVector missZ(leng);
    Rcpp::NumericVector adjustConst(leng);
    for(int r = 0; r < leng; r++) {
      for (int c = 0; c < width; c++) {
        g(r,c) = R::dnorm(outcome[r],prediction(r,c),sqrt(sigma2),0);
      }
      for (int c = 0; c < width; c++) {
        znumerator(r,c) = g(r,c)*w_old(c);
        if(R_IsNA(znumerator(r,c))==TRUE){;
          znumerator(r,c) = 0;
        }
      }
      zdenom(r) = sum(znumerator(r,_));    
    } 
    for(int r = 0; r<leng; r++){;
      missZ(r) = 0;
      for(int c=0; c<width; c++){
        Zs(r,c) = znumerator(r,c)/zdenom(r);
        if(Zs(r,c)<1e-4){;
          Zs(r,c) = 0;
        }
      }
      LogicalVector na_test = is_na(Zs(r,_));
      for(int c=0; c<width; c++){
        if(na_test(c)==FALSE)
          missZ(r) +=1; 
      }
      adjustConst(r) = (wisdom*1.0)/missZ(r);
    }
    for(int r = 0; r<leng; r++){;
      for(int c=0; c<width; c++){
        Zs(r,c) = adjustConst(r) + ((1.0-wisdom)*Zs(r,c));
      }
      LogicalVector na_testZ = is_na(Zs(r,_));
      for (int c = 0; c < width; c++) {;
        if(na_testZ(c)==TRUE){;
          Zs(r,c)=0;
        }
      }
    }
    for(int c = 0; c<width; c++){;
      unnormalizedW[c] = (sum(Zs(_,c)));
      for(int r = 0; r < leng; r++){
        RZ(r,c) = Zs(r,c)*RSQ(r,c);
      }
    }
    for(int r = 0; r < leng; r++){
      LogicalVector na_testRZ = is_na(RZ(r,_));
      for (int c = 0; c < width; c++) {;
        if(na_testRZ(c)==TRUE){;
          RZ(r,c)=0;
        } 
      }
    }
    sigma2 = sum(RZ)/sum(Zs);
    
    W = unnormalizedW / double(sum(unnormalizedW));
    for(int c = 0; c<width; c++){;
      if(W(c)<1e-4){;
        W(c) = 0;
      }
    }
    LL = sum(log(zdenom));
    double top = double(fabs(LL_old-LL));
    double bot = (1.0+fabs(LL));
    improv  = top/bot;
    LL_old = LL;
    iter +=1;
  } 
  
  return Rcpp::List::create(Rcpp::Named("LL") = LL,Rcpp::Named("W") = W,Rcpp::Named("Sigma2") = sigma2,Rcpp::Named("Iterations") = iter, Rcpp::Named("improv") = improv);
}


// [[Rcpp::export]]
Rcpp::List GibbsNormal(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, Rcpp::NumericVector alpha, double sigma, int iterations, int burnin, int thin) {

  int length = prediction.nrow();
  int nmods = prediction.ncol();
  int outcount = 0;
  int output = round((iterations-burnin)/thin);
  Rcpp::List theta_post(iterations);
  Rcpp::NumericMatrix W_post(iterations,nmods);
  Rcpp::NumericMatrix W_out(output,nmods);
  Rcpp::NumericVector Sigma_post(iterations);
  Rcpp::NumericVector Sigma_out(output);
  
  
  for (int iterator = 0; iterator < iterations; iterator++){;
    Rcpp::NumericVector W_use(nmods);
    Rcpp::NumericMatrix evalsEach(length,nmods);
    Rcpp::NumericVector evalsAll(length);
    Rcpp::NumericMatrix theta(length,nmods);
    Rcpp::NumericMatrix T(length,nmods);
    Rcpp::NumericVector eta(nmods);
    Rcpp::NumericVector ssq(1);
    double temp;
    double sigma_use;
    Rcpp::NumericVector w_gamma(nmods);
    
    if(iterator == 0){;
      W_use = W;
      sigma_use = sigma;
    }
    if(iterator > 0){;
      W_use = W_post(iterator-1,_);
      sigma_use = Sigma_post(iterator-1);
    }    
    for(int m = 0; m < nmods; m++){;
      for(int i = 0; i < length; i++){;
        evalsEach(i,m) = W_use(m)*(R::dnorm(outcome(i),prediction(i,m),sigma_use,false));
      };
    };
    for(int i = 0; i < length; i++){
      evalsAll(i) = sum(evalsEach(i,_));
      theta(i,_) = evalsEach(i,_)/evalsAll(i);
    }
    
    for(int i=0; i<length; i++){;
      T(i,_) = oneMultinomCalt(theta(i,_));
    };
    
    for(int m=0; m<nmods; m++){
      eta(m) = alpha(m) + sum(T(_,m));
    }  
    
    Rcpp::NumericVector t_col = prediction.ncol();
    NumericVector pred_col = prediction.ncol();
    NumericVector test = prediction.ncol();
    NumericVector calc(1);
    for(int i=0; i<length; i++){;
      t_col = T(i,_);
      pred_col = prediction(i,_);
      ssq += pow((outcome(i) - as<NumericVector>(pred_col[t_col==1])),2.0);
    };
    
    
    double sample_sum = 0;
    
    for(int m = 0; m<nmods; m++){
      w_gamma(m) = as<double>(rgamma(1, eta(m), 1));
      sample_sum += w_gamma(m);
    }
    
    for(int m = 0; m<nmods; m++){
      W_post(iterator,m) =  w_gamma(m)/sample_sum;
      if(((iterator+1) % thin == 0) and (iterator+1 > burnin)){;
        W_out(outcount,m) =   W_post(iterator,m);
      }
    }
    temp  = (1/as<double>(Rcpp::rgamma(1,(length+1)/2, (1/(as<double>(ssq)/2)))));
    Sigma_post(iterator) = sqrt(temp);
    if(((iterator+1) % thin == 0) and (iterator+1 > burnin)){;
      Sigma_out(outcount) = Sigma_post(iterator);
    }
    
    if(((iterator+1) % thin == 0) and (iterator+1 > burnin)){;
      outcount += 1;
    }  
    if ((iterator+1) % 5000 == 0 ){;
      Rcpp::Rcout << "Iteration: " << iterator+1 << std::endl;
    }
  };   
  return Rcpp::List::create(Rcpp::Named("W") = W_out, Rcpp::Named("Sigma") = Sigma_out);
}


// [[Rcpp::export]]
Rcpp::List GibbsNormalMissing(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, Rcpp::NumericVector alpha, double sigma, int iterations, int burnin, int thin) {
  
  int length = prediction.nrow();
  int nmods = prediction.ncol();
  int outcount = 0;
  int output = round((iterations-burnin)/thin);
  Rcpp::NumericMatrix W_post(iterations,nmods);
  Rcpp::NumericMatrix W_out(output,nmods);
  Rcpp::NumericVector Sigma_post(iterations);
  Rcpp::NumericVector Sigma_out(output);
  Rcpp::NumericMatrix MissingInd(length,nmods);
  
  for(int i = 0; i < length; i++){;
    MissingInd(i,_) = isNA(prediction(i,_));
  }

  for (int iterator = 0; iterator < iterations; iterator++){;
    Rcpp::NumericVector W_use(nmods);
    Rcpp::NumericMatrix evalsEach(length,nmods);
    Rcpp::NumericVector evalsAll(length);
    Rcpp::NumericMatrix theta(length,nmods);
    Rcpp::NumericMatrix T(length,nmods);
    Rcpp::NumericVector eta(nmods);
    Rcpp::NumericVector ssq(1);
    double temp;
    double sigma_use;
    Rcpp::NumericVector w_gamma(nmods);
    
    if(iterator == 0){;
      W_use = W;
      sigma_use = sigma;
    }
    if(iterator > 0){;
      W_use = W_post(iterator-1,_);
      sigma_use = Sigma_post(iterator-1);
    }    
    for(int m = 0; m < nmods; m++){;
      for(int i = 0; i < length; i++){;
        if(MissingInd(i,m) == FALSE){;
          evalsEach(i,m) = W_use(m)*(R::dnorm(outcome(i),prediction(i,m),sigma_use,false));
        };
        if(MissingInd(i,m) == TRUE){;
          evalsEach(i,m) = 0;
        };
      };
    };
    for(int i = 0; i < length; i++){
      evalsAll(i) = sum(evalsEach(i,_));
      theta(i,_) = evalsEach(i,_)/evalsAll(i);
    }
    for(int i=0; i<length; i++){;
      T(i,_) = oneMultinomCalt(theta(i,_));
    };
    
    for(int m=0; m<nmods; m++){
      eta(m) = alpha(m) + sum(T(_,m));
    }  
    
    Rcpp::NumericVector t_col = prediction.ncol();
    NumericVector pred_col = prediction.ncol();
    NumericVector test = prediction.ncol();
    NumericVector calc(1);
    for(int i=0; i<length; i++){;
      t_col = T(i,_);
      pred_col = prediction(i,_);
      ssq += pow((outcome(i) - as<NumericVector>(pred_col[t_col==1])),2.0);
    };
    
    
    double sample_sum = 0;
    
    for(int m = 0; m<nmods; m++){
      w_gamma(m) = as<double>(rgamma(1, eta(m), 1));
      sample_sum += w_gamma(m);
    }
    
    for(int m = 0; m<nmods; m++){
      W_post(iterator,m) =  w_gamma(m)/sample_sum;
      if(((iterator+1) % thin == 0) and (iterator+1 > burnin)){;
        W_out(outcount,m) =   W_post(iterator,m);
      }
    }
    temp  =  (1/as<double>(Rcpp::rgamma(1,(length+1)/2, (1/(as<double>(ssq)/2)))));
    Sigma_post(iterator) = sqrt(temp);
    if(((iterator+1) % thin == 0) and (iterator+1 > burnin)){;
      Sigma_out(outcount) = Sigma_post(iterator);
    }
    
    if(((iterator+1) % thin == 0) and (iterator+1 > burnin)){;
      outcount += 1;
    }  
    if ((iterator+1) % 5000 == 0 ){;
      Rcpp::Rcout << "Iteration: " << iterator+1 << std::endl;
    }
  };   
  return Rcpp::List::create(Rcpp::Named("W") = W_out, Rcpp::Named("Sigma") = Sigma_out);
}


// [[Rcpp::export]]
Rcpp::List emLogit(Rcpp::NumericVector outcome,
                   Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, double tol, int maxIter,
                   double wisdom) {
  double LL = 0.00;
  double LL_old = 0.00;
  int leng = prediction.nrow();
  int width = prediction.ncol();
  int iter = 0 ;
  double improv = 1.1;
  
  while((improv > tol) & (iter < maxIter)){
    Rcpp::NumericMatrix  znumerator(leng,width);
    Rcpp::NumericMatrix Zs(leng,width); 
    Rcpp::NumericVector zdenom(leng);
    Rcpp::NumericVector unnormalizedW(width);
    Rcpp::NumericVector w_old(W);
    Rcpp::NumericVector missZ(leng);
    Rcpp::NumericVector adjustConst(leng);
    for (int r = 0; r < leng; r++) {; 
      if(outcome(r)==1){;
        for (int c = 0; c < width; c++) {
          znumerator(r,c) = prediction(r,c)*w_old(c);
        }
      } 
      if(outcome(r)==0){;
        for (int c = 0; c < width; c++) {
          znumerator(r,c) = (1-prediction(r,c))*(w_old(c));
        }
      }
      LogicalVector na_test = is_na(znumerator(r,_));
      for (int c = 0; c < width; c++) {;
        if(na_test(c)==TRUE or znumerator(r,c)< 1e-4 ){;
          znumerator(r,c)=0;
        }
      }
      zdenom(r) = sum(znumerator(r,_));
    }
    
    for(int c = 0; c<width; c++){;
      for(int r=0; r<leng; r++){
        Zs(r,c) = znumerator(r,c)/zdenom(r);
      }
      LogicalVector na_testz = is_na(Zs(_,c));
      for (int r = 0; r < leng; r++) {;
        if(na_testz(r)==TRUE or Zs(r,c)< 1e-4){;
          Zs(r,c)=0;
        }
      }
    }
    for(int r = 0; r<leng; r++){;
      missZ(r) = 0;
      LogicalVector na_test = is_na(Zs(r,_));
      for(int c=0; c<width; c++){
        if(na_test(c)==FALSE)
          missZ(r) +=1; 
      }
      adjustConst(r) = (wisdom*1.0)/missZ(r);
    }
    for(int r = 0; r<leng; r++){;
      for(int c=0; c<width; c++){
        Zs(r,c) = adjustConst(r) + ((1.0-wisdom)*Zs(r,c));
      }
    }  
    for(int r = 0; r<width; r++){;
      unnormalizedW[r] = (sum(Zs(_,r)));
    }
    W = unnormalizedW / double(sum(unnormalizedW));
    for(int r = 0; r<width; r++){;
      if(W(r)<1e-4){;
        W(r) = 0;
      }
    }
    LL = sum(log(zdenom));
    double top = double(fabs(LL_old-LL));
    double bot = (1.0+fabs(LL));
    improv  = top/bot;
    LL_old = LL;
    iter +=1;
  }
  return Rcpp::List::create(Rcpp::Named("LL") = LL,Rcpp::Named("W") = W,Rcpp::Named("Iterations") = iter, Rcpp::Named("improv") = improv);
}

// [[Rcpp::export]]
Rcpp::List GibbsLogit(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, int iterations, int burnin, int thin) {
  
  int leng = prediction.nrow();
  int width = prediction.ncol();
  int outcount = 0;
  int output = round((iterations-burnin)/thin);
  Rcpp::List theta_post(iterations);
  Rcpp::NumericMatrix W_post(iterations,width);
  Rcpp::NumericMatrix W_out(output,width);
  
  for (int iterator = 0; iterator < iterations; iterator++){;
    Rcpp::NumericVector W_use(width);
    Rcpp::NumericMatrix theta(leng,width);
    Rcpp::NumericMatrix numerator(leng,width);
    Rcpp::NumericVector denom(leng);
    Rcpp::NumericMatrix T(leng,width);
    Rcpp::NumericVector w_gamma(width);
    Rcpp::NumericVector nu(width);
    
    if(iterator == 0){;
      W_use = W;
    }
    if(iterator > 0){;
      W_use = W_post(iterator-1,_);
    }  
    for (int i =0; i < leng; i++){;
      if(outcome(i) == 1){;
        for(int m = 0; m < width; m++) {;
          numerator(i,m) = prediction(i,m)*(W_use(m));
        }
      }  
      if(outcome(i) == 0){;
        for(int m = 0; m < width; m++) {;
          numerator(i,m) = (1 - prediction(i,m))*(W_use(m));
        };
      };
      denom(i) = sum(numerator(i,_));
      theta(i,_) = numerator(i,_)/denom(i);
    };
    
    theta_post(iterator) = theta; 
    
    for(int i=0; i<leng; i++){;
      T(i,_) = oneMultinomCalt(theta(i,_));
    };
    
    for(int m=0; m<width; m++){
      nu(m) = 1 + sum(T(_,m));
    }  
    
    double sample_sum = 0;
    for(int m = 0; m<width; m++){
      w_gamma(m) = as<double>(rgamma(1, nu(m), 1));
      sample_sum += w_gamma(m);
    }
    for(int m = 0; m<width; m++){
      W_post(iterator,m) =  w_gamma(m)/sample_sum;
      if(((iterator+1) % 3 == 0) and (iterator+1 > burnin)){;
        W_out(outcount,m) =   W_post(iterator,m);
      }
    }
    if(((iterator+1) % thin == 0) and (iterator+1 > burnin)){;
      outcount += 1;
    }  
    if ((iterator+1) % 5000 == 0 ){;
      Rcpp::Rcout << "Iteration: " << iterator+1 << std::endl;
    }
  };   
  
  
  
  
  return Rcpp::List::create(Rcpp::Named("W_out") = W_out);
}

