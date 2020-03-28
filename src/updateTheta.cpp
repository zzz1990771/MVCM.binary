#include "ThetaA.h"
/// [[Rcpp::depends(RcppArmadillo)]] //only work for sourceCPP
// objective function is \sum||Y_i-(X_i^T \otimes I_q \Theta A^T B(T_i))||^2
//                  //   +\sum \lambda_j ||\Theta_j||

// Use Blockwise ISTA algorithm first

// In the following:
// Y=(Y_1,\cdots,Y_n): q \times n
// X=(X_1,\cdots,X_n): (p+1) \times n
// B=(B_1,\cdots,B_n): K \times n (plug in splines)
// A: K \times r
// \Theta: q(p+1) \times r
// \Theta_j: q \times r
// \lambda_j: tunning parameter for \Theta_j
// \tau: ISTA step-size

double logit(double x){
    return (1 / (1 + exp(-x))); 
} 

arma::mat Residual(const arma::mat *Y,const arma::mat *X,const arma::mat* B
                     ,const arma::mat* Theta,const arma::mat *A){
  int i;
  int j;
  // cube format of Theta:
  arma::cube Theta_cube((X->n_rows)*(Y->n_rows),A->n_cols,1);
  Theta_cube.slice(0)=*Theta;
  Theta_cube.reshape(Y->n_rows,X->n_rows,A->n_cols);

  // cube format for the X in order multiplication Theta*X
  arma::cube X_mass(Y->n_rows,Y->n_cols,A->n_cols);
  for(i=0;i<A->n_cols;++i) X_mass.slice(i) = Theta_cube.slice(i)*(*X);
  // X_mass.each_slice( [arma::mat r=X_temp] (arma::mat &K) { return K*r; }); // unsuccessful lambda function.

  // Times B->t()**A, for the reason that cube dimension matching.
  arma::cube beta(Y->n_cols,A->n_cols,1);
  beta.slice(0)= B->t()*(*A); // arma::trans(A->t()*(*B));
  beta.reshape(1,Y->n_cols,A->n_cols);
  for(i=0;i<Y->n_rows;++i) X_mass(i,0,0,size(beta)) = X_mass(i,0,0,size(beta)) % beta;

  arma::mat diff(Y->n_rows,Y->n_cols,arma::fill::zeros);
  for(i=0;i<A->n_cols;++i) diff+=X_mass.slice(i);
  
  //Y - g^-1(ita)
  for(i=0;i<Y->n_rows;++i){

    for(j=0;j<Y->n_cols;++j){
      diff(i,j) = (*Y)(i,j) - logit(diff(i,j));
    }
  }
  return diff;
}


arma::mat nLogLikelihood(const arma::mat *Y,const arma::mat *X,const arma::mat* B
                     ,const arma::mat* Theta,const arma::mat *A){
  int i;
  int j;
  // cube format of Theta:
  arma::cube Theta_cube((X->n_rows)*(Y->n_rows),A->n_cols,1);
  Theta_cube.slice(0)=*Theta;
  Theta_cube.reshape(Y->n_rows,X->n_rows,A->n_cols);

  // cube format for the X in order multiplication Theta*X
  arma::cube X_mass(Y->n_rows,Y->n_cols,A->n_cols);
  for(i=0;i<A->n_cols;++i) X_mass.slice(i) = Theta_cube.slice(i)*(*X);
  // X_mass.each_slice( [arma::mat r=X_temp] (arma::mat &K) { return K*r; }); // unsuccessful lambda function.

  // Times B->t()**A, for the reason that cube dimension matching.
  arma::cube beta(Y->n_cols,A->n_cols,1);
  beta.slice(0)= B->t()*(*A); // arma::trans(A->t()*(*B));
  beta.reshape(1,Y->n_cols,A->n_cols);
  for(i=0;i<Y->n_rows;++i) X_mass(i,0,0,size(beta)) = X_mass(i,0,0,size(beta)) % beta;

  arma::mat diff(Y->n_rows,Y->n_cols,arma::fill::zeros);
  for(i=0;i<A->n_cols;++i) diff+=X_mass.slice(i);
  
  //change to negative likelihood for binary logit link
  for(i=0;i<Y->n_rows;++i){

    for(j=0;j<Y->n_cols;++j){
      diff(i,j) = log(1+exp(diff(i,j)))-(*Y)(i,j)*(diff(i,j));
    }
  }
  //diff=(*Y-diff);
  return diff;
}

// This Proximal is for (adaptive) group lasso
arma::mat Proximal(const arma::mat &Z,const double &lambda){
  // Z and Z_next have the same size
  arma::mat Z_next=Z;
  double Z_norm=arma::norm(Z,"fro");
  if(lambda >=Z_norm) Z_next.zeros(Z.n_rows,Z.n_cols); // or Z_next.zeros();
  else Z_next=Z-lambda*(Z/Z_norm);
  return(Z_next);
}

// This ProximalScad is for scad group penalty
arma::mat ProximalScad(const arma::mat &Z, const double &lambda, const double &a){
  // Z_next is a block in whole question setting
  arma::mat Z_next=Z;
  double Z_norm=arma::norm(Z,"fro");
  if((2.0*lambda)>=Z_norm) Z_next=Proximal(Z,lambda);
  else if(a*lambda>=Z_norm) Z_next=(a-1.0)/(a-2.0)*(Proximal(Z,a*lambda/(a-1.0)));
  else { }
  return(Z_next);
}


arma::mat grad_H(const arma::mat *Y,const arma::mat *X,const arma::mat* B,const arma::mat* Theta,const arma::mat *A){
  arma::mat sum_gradH((Y->n_rows)*(X->n_rows),A->n_cols,arma::fill::zeros);
  arma::mat gradH_temp=sum_gradH ;
  int i,j,q=Y->n_rows;
  arma::mat diff=Residual(Y,X,B,Theta,A);
  arma::mat secondPart(q,A->n_cols,arma::fill::zeros);

  for(i=0;i<Y->n_cols;++i){
    secondPart=diff.col(i)*(B->col(i).t())*(*A);
    for(j=0;j < X->n_rows;++j) {
      gradH_temp.rows(j*q,(j+1)*q-1)=(*X)(j,i)*secondPart;
    }
    sum_gradH += -gradH_temp;
   // sum_gradH +=-2.0*arma::kron(X->col(i),I_q)*(Y->col(i)-(arma::kron((X->col(i)).t(),I_q))*(*Theta)*(A->t())*B->col(i))*B->col(i).t()*(*A);
  }
  return(sum_gradH);
}

double objective_Theta(const arma::mat *Y,const arma::mat *X,const arma::mat* B,const arma::mat* Theta,const arma::mat *A, const vector<double> *lambda,const int* m){
  double obj_Theta=0;
  int q=Y->n_rows;
  int p_plusOverm=(X->n_rows)/(*m);
  arma::mat I_q(Y->n_rows,Y->n_rows,arma::fill::eye);
  //construction error
  // for(i=0;i<Y->n_cols;++i) obj_Theta+=pow(arma::norm((Y->col(i)-(arma::kron((X->col(i)).t(),I_q))*(*Theta)*(A->t())*B->col(i)),"fro"),2.0);

  arma::mat diff=nLogLikelihood(Y,X,B,Theta,A);
  //get sqaure
  //diff %=diff;
  
  obj_Theta+=arma::accu(diff);


  //penalties
//   for(i=0;i<X->n_rows;++i){
//     obj_Theta+=lambda->at(i)*arma::norm(Theta->rows(i*q,(i+1)*q-1),"fro");
//   }

  arma::mat* penalty_temp=new arma::mat;
  *penalty_temp=Theta->t();

  penalty_temp->reshape(Theta->n_cols*q*(*m),p_plusOverm);  // reshaper to desired matrix dimensions
  *penalty_temp=arma::square(*penalty_temp);  // element-wise square

  arma::rowvec lambda_vec=arma::rowvec(*lambda);   // transfer vector<double> to arma::rowvec

  arma::rowvec norm_vec=arma::sum(*penalty_temp,0);
  norm_vec=arma::sqrt(norm_vec);                    // norm of desired columns.
  norm_vec.each_row() %= lambda_vec;     // each row vector element-wise multiplication; length of lambda is pOverm
  obj_Theta += accu(norm_vec);
  delete penalty_temp;

  return obj_Theta;
}


//[[Rcpp::export]]
 SEXP update_Theta(const SEXP Theta_r, const SEXP Y_r,const SEXP X_r,SEXP B_r,const SEXP A_r,
                             const SEXP lambda_r,const SEXP tol_r,const SEXP MaxIt_r,const SEXP tau_r,const SEXP m_r,
                             const SEXP method_r, const SEXP a_r){

  //Data and parameters transferred
  NumericMatrix B_s(B_r),A_s(A_r),Y_s(Y_r),Theta_s(Theta_r);
  NumericMatrix X_s(X_r);
  // NumericMatrix test_grad_s(test_grad_r);
  // arma::mat test_grad(test_grad_s.begin(),test_grad_s.nrow(),test_grad_s.ncol(),false);
  arma::mat Theta(Theta_s.begin(),Theta_s.nrow(),Theta_s.ncol(),false);
  arma::mat Y(Y_s.begin(),Y_s.nrow(),Y_s.ncol(),false);
  arma::mat X(X_s.begin(),X_s.nrow(),X_s.ncol(),false);
  arma::mat A(A_s.begin(),A_s.nrow(),A_s.ncol(),false);
  arma::mat B(B_s.begin(),B_s.nrow(),B_s.ncol(),false);
  arma::mat Theta_outer;

  vector<double> lambda=as<vector<double> >(lambda_r); //space is required for > >, length is pOverm
  // lambda always has length p_plus1, since we need to use this in the objetive_Theta function.
  double tau=as<double>(tau_r);
  double tol=as<double>(tol_r);
  int MaxIt=as<int>(MaxIt_r);
  int m=as<int>(m_r);
  // for sacd:
  int method=as<int>(method_r);
  double a=as<double>(a_r);

  int p_plus1;
  p_plus1=X.n_rows;
  int p_plusOverm=p_plus1/m; // in case for dummy group variable
  int q=Y.n_rows;
  int i=0,iter=0;
  double objValue;
  double tol_now,tol_outer,tol_outer2;
  objValue=objective_Theta(&Y,&X,&B,&Theta,&A,&lambda,&m);
  Rcout<<"The beginning objective function value is "<<objValue<<". "<<endl;
  arma::mat *Theta_now=new arma::mat;
  *Theta_now=Theta;
  arma::mat *gradH_now=new arma::mat;
  double tau_now, objValue_pre;

  //for saving and backtracking the scad proximal
  // the reason that lasso does not need this is because it always shrinks and satisfy the criterion
  // while A and scad may pass the "gold point" and stay the original point forever.
  vector<double> save_obj;
  vector<arma::mat> save_Theta;

  tol_outer=2.0*tol;
  do{   // iteratively solving Theta;
    Theta_outer=Theta;
    tau_now=tau;
    objValue_pre=objValue=objective_Theta(&Y,&X,&B,&Theta,&A,&lambda,&m);
    // begin to save the history
    save_obj.push_back(objValue_pre);
    save_Theta.push_back(Theta);

    *gradH_now=grad_H(&Y,&X,&B,&Theta,&A);
 //   Rcout<<"This is the iter "<<iter<<", And the tol is "<<tol_now<<endl;
   do{ // start to do backtracking
     if(method==1){
       for(i=0;i<p_plusOverm;++i){
         Theta_now->rows(i*m*q,(i+1)*m*q-1)=Proximal(Theta.rows(i*m*q,(i+1)*m*q-1)-tau_now*(gradH_now->rows(i*m*q,(i+1)*m*q-1)),tau_now*lambda[i]);
       }
     }
     if(method==2){
       if(tau_now<tol*0.00001) { // no need to keep doing, *Theta_now will stay the same as Theta, track the history.
//          for(i=0;i<p_plusOverm;++i){
//            Theta_now->rows(i*m*q,(i+1)*m*q-1)=save_Theta[distance(save_obj.begin(),min_element(save_obj.begin(),save_obj.end()))];
//          }
            *Theta_now=save_Theta[distance(save_obj.begin(),min_element(save_obj.begin(),save_obj.end()))];
         break;
       }
       // ordinary situation:
       for(i=0;i<p_plusOverm;++i){
         Theta_now->rows(i*m*q,(i+1)*m*q-1)=ProximalScad(Theta.rows(i*m*q,(i+1)*m*q-1)-tau_now*(gradH_now->rows(i*m*q,(i+1)*m*q-1)),tau_now*lambda[i],a);
       }
     }
      tol_now=arma::norm(*Theta_now-Theta,"fro");
      objValue_pre=objective_Theta(&Y,&X,&B,Theta_now,&A,&lambda,&m);
      tau_now=0.5*tau_now;
      save_obj.push_back(objValue_pre);
      save_Theta.push_back(*Theta_now);
 //    Rcout<<"This is the iter "<<iter<<", And the objValue_pre is "<<objValue_pre<<"reoutput the objvalue "<<objective_Theta(&Y,&X,&B,&Theta,&A,&lambda)<<" diff="<<0.1/(2.0*tau_now)*pow(tol_now,2.0)<<"tau_now="<<tau_now<<endl;

      if(tol_now<tol){ // two points are very close, the step-size is too small or close to optimal value;
        break;
      }
   }while(objValue_pre>=objValue-0.1/(2.0*tau_now)*pow(tol_now,2.0)); // out of backtracking
    Theta=*Theta_now;
    ++iter;
    tol_outer=arma::norm(Theta-Theta_outer,"fro");
    tol_outer2=abs(objValue-objValue_pre);
    // clean the history
    vector<double>().swap(save_obj);
    vector<arma::mat>().swap(save_Theta);
    // // test the gradient of H, updated in R
    // test_grad=*gradH_now;

  }while (iter<MaxIt && (tol_outer>tol && tol_outer2 > tol));

  for(i=0;i<p_plusOverm;++i){
    if(arma::norm(Theta.rows(i*m*q,(i+1)*m*q-1),"fro")<tol*m*q*Theta.n_cols) (Theta.rows(i*m*q,(i+1)*m*q-1)).zeros();
  }
  objValue=objective_Theta(&Y,&X,&B,&Theta,&A,&lambda,&m);
  Rcout<<"The current objective function value is "<<objValue<<". "<<endl;

  delete gradH_now;
  delete Theta_now;
  return wrap(Theta);

}
