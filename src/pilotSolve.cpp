#include "ThetaA.h"

double logitC(double x){
    return (1 / (1 + exp(-x))); 
} 

arma::mat Residual_forC(const arma::mat *Y, const arma::mat *X, const arma::mat *B, const arma::mat* C){
  int i,j,q=Y->n_rows;
  arma::mat diff=*Y;
  arma::mat diff_temp(q,1,arma::fill::zeros);
  arma::mat secondPart(C->n_rows,1);

  for(i=0;i<Y->n_cols;++i){
    secondPart=(*C)*B->col(i);
    for(j=0;j<X->n_rows;++j){
      diff_temp += (*X)(j,i)*secondPart.rows(j*q,(j+1)*q-1);
    }
    diff.col(i)=diff_temp;
    diff_temp.zeros();
  }
  for(i=0;i<Y->n_rows;++i){

    for(j=0;j<Y->n_cols;++j){
      //change to negative likelihood for binary logit link
      diff(i,j) = (*Y)(i,j) - logitC(diff(i,j));
    }
  }
  //diff=(*Y-diff);
  return diff;
}

arma::mat nLogLikelihood_forC(const arma::mat *Y, const arma::mat *X, const arma::mat *B, const arma::mat* C){
  int i,j,q=Y->n_rows;
  arma::mat diff=*Y;
  arma::mat diff_temp(q,1,arma::fill::zeros);
  arma::mat secondPart(C->n_rows,1);

  for(i=0;i<Y->n_cols;++i){
    secondPart=(*C)*B->col(i);
    for(j=0;j<X->n_rows;++j){
      diff_temp += (*X)(j,i)*secondPart.rows(j*q,(j+1)*q-1);
    }
    diff.col(i)=diff_temp;
    diff_temp.zeros();
  }
  for(i=0;i<Y->n_rows;++i){
    for(j=0;j<Y->n_cols;++j){
      //change to negative likelihood for binary logit link
      diff(i,j) = log(1+exp(diff(i,j)))-(*Y)(i,j)*(diff(i,j));
    }
  }
  //diff=(*Y-diff);
  return diff;
}


arma::mat pilotGrad(const arma::mat *Y, const arma::mat *X, const arma::mat* B,const arma::mat* C){
  int i,j,q=Y->n_rows;
  arma::mat secondPart(q,B->n_rows,arma::fill::zeros);
  arma::mat gradient((Y->n_rows)*(X->n_rows),B->n_rows,arma::fill::zeros);
  arma::mat gradient_temp=gradient;
  arma::mat diff=Residual_forC(Y,X,B,C);

  for(i=0;i<Y->n_cols;++i){
    secondPart=diff.col(i)*(B->col(i).t());
    for(j=0;j<X->n_rows;++j){
      gradient_temp.rows(j*q,(j+1)*q-1)=(*X)(j,i)*secondPart;
    }
    gradient += -gradient_temp;
   // gradient +=-2.0*arma::kron(X->col(i),I_q)*(Y->col(i)-(arma::kron((X->col(i)).t(),I_q))*(*C)*B->col(i))*B->col(i).t();
  }
  return(gradient);
}

double objPilot(const arma::mat *Y,const arma::mat *X,const arma::mat* B,const arma::mat* C){
  double obj_pilot=0.0;
  arma::mat diff=nLogLikelihood_forC(Y,X,B,C);
  //diff %=diff;
  obj_pilot=arma::accu(diff);
  return obj_pilot;
}

//// add projection and retraction on the fixed rank manifold on Jan, 2016.
void retractProjectPilot(const arma::mat *gradF, const arma::mat* C,arma::mat *grad,
                         arma::mat *retraction,const double* stepSize,const int *rank){
  arma::mat AB,AB_plus,A_plus_B,U,V,projU,projV,projU_plus,projV_plus;
  arma::vec s;
  svd(U,s,V,*C);
  U=U.cols(0,*rank-1);
  V=V.cols(0,*rank-1);
  s=s(arma::span(0,*rank-1));
  arma::mat I_p(C->n_rows,C->n_rows,arma::fill::eye);
  arma::mat I_q(C->n_cols,C->n_cols,arma::fill::eye);

  projU=U*U.t();
  projV=V*V.t();
  projU_plus=I_p-projU;
  projV_plus=I_q-projV;

  //// Projection on tangent space
  *grad=*gradF-projU_plus*(*gradF)*projV_plus;

  AB=projU*(*gradF)*projV*(*stepSize);
  AB_plus=projU*(*gradF)*projV_plus*(*stepSize);
  A_plus_B=projU_plus*(*gradF)*projV*(*stepSize);

  arma::vec s_inv(s.n_elem);
  int i;
  for(i=0;i<s.n_elem;++i){
    s_inv(i)=1.0/s(i);
  }
  arma::mat C_plus=V*arma::diagmat(s_inv)*U.t();

  arma::mat V1=*C-1.0/2.0*AB+A_plus_B-1.0/8.0*AB*C_plus*AB-1.0/2.0*A_plus_B*C_plus*AB;
  arma::mat V2=*C-1.0/2.0*AB+AB_plus-1.0/8.0*AB*C_plus*AB-1.0/2.0*AB*C_plus*AB_plus;

  *retraction=V1*C_plus*V2;
}

//[[Rcpp::export]]
SEXP pilotSolve(const SEXP Y_r, const SEXP X_r, const SEXP B_r,const SEXP C_r,
                           const SEXP tol_r, const SEXP Max_r,
                           const SEXP rank_r){
  //data input
  double tol=as<double>(tol_r);
  int Max=as<int>(Max_r);
  int rank=as<int>(rank_r);
  NumericMatrix Y_s(Y_r),X_s(X_r), B_s(B_r),C_s(C_r);
  arma::mat Y(Y_s.begin(),Y_s.nrow(),Y_s.ncol(),false);
  arma::mat X(X_s.begin(),X_s.nrow(),X_s.ncol(),false);
  arma::mat B(B_s.begin(),B_s.nrow(),B_s.ncol(),false);
  arma::mat C(C_s.begin(),C_s.nrow(),C_s.ncol(),false);


  // controls
  int iter=0,innerIter;
  double tol_now;
  arma::mat *C_temp=new arma::mat;
  arma::mat *steepD=new arma::mat;
  arma::mat *steepD_F=new arma::mat;
  double obj,obj_temp,stepSize,eAscent,obj_outer,obj_Des;
  obj=objPilot(&Y,&X,&B,&C);
  arma::mat C_outer;
  vector<double> save_obj;
  vector<arma::mat> save_C;
  *C_temp=C;
  bool flag=true;

  // Armijo parameter
  double beta=0.5;
  double sigma=0.05;
  int subMax=Max;

  // Steepest Ascent and Armijo back tracking
  do{
    // control of inner Armijo back tracking
    iter++;
    stepSize=1.0;
    obj_Des=-1.0;
    obj_outer=obj;
    C_outer=C;
    // negative gradient.
    *steepD_F=-1.0*pilotGrad(&Y,&X,&B,&C);
    retractProjectPilot(steepD_F,&C,steepD,C_temp,&stepSize,&rank);
    eAscent=sigma*arma::dot(*steepD,*steepD);

    innerIter=0;
    tol_now=-1.0;
    do{
      innerIter++;
      stepSize=beta*stepSize;
      eAscent=beta*eAscent;
      retractProjectPilot(steepD_F,&C,steepD,C_temp,&stepSize,&rank);
      obj_temp=objPilot(&Y,&X,&B,C_temp);
      save_obj.push_back(obj_temp);
      save_C.push_back(*C_temp);
    }while (innerIter<subMax && (obj-obj_temp)<eAscent ); //out of Armijo
    obj_Des=obj_outer-obj_temp;
    if(obj_Des<0.0) {
      obj=*min_element(save_obj.begin(),save_obj.end());
      C=save_C[distance(save_obj.begin(),min_element(save_obj.begin(),save_obj.end()))];
    } else{
      C=*C_temp;
      obj=objPilot(&Y,&X,&B,&C);
    }
    // clean the container
    vector<double>().swap(save_obj);
    vector<arma::mat>().swap(save_C);

    tol_now=arma::norm((C-C_outer),"fro");
    if(tol_now < tol) flag=false;

  } while (iter<Max && flag); // Out of iteration
  delete C_temp;
  delete steepD;
  delete steepD_F;
  return(wrap(C));
}
