#include "ThetaA.h"


arma::mat gradA(const arma::mat *Y,const arma::mat *X,const arma::mat *Theta,
                  const arma::mat *B, const arma::mat *A){
  int i,j,q=Y->n_rows;
  arma::mat diff=Residual(Y,X,B,Theta,A);
  arma::mat grad(B->n_rows,Theta->n_cols,arma::fill::zeros);   //K*r
  arma::mat grad_temp(B->n_rows,(q*X->n_rows),arma::fill::zeros);  // K *(pq)
  arma::mat grad_tempSum=grad_temp;
  arma::mat firstPart(A->n_rows,q,arma::fill::zeros);      //K*q

  for(i=0;i<Y->n_cols;++i){
    firstPart=B->col(i)*diff.col(i).t();

    for(j=0;j<X->n_rows;++j){
      grad_temp.cols(j*q,(j+1)*q-1)=firstPart*(*X)(j,i);
    }
    grad_tempSum += grad_temp;
    // grad=grad+(-2.0)*(B->col(i))*(Y->col(i)-arma::kron((X->col(i)).t(),I_q)*(*Theta)*A->t()*B->col(i)).t()*(arma::kron((X->col(i)).t(),I_q))*(*Theta);
  }
  grad=-grad_tempSum*(*Theta);
  return(grad);
}

arma::mat projectionA(arma::mat *gradF,const arma::mat *A){
  arma::mat grad_normal=A->t()*(*gradF);
  grad_normal=0.5*(*A)*(grad_normal+grad_normal.t());
  arma::mat grad=*gradF-grad_normal;
  return(grad);
}


arma::mat retractA(const arma::mat *direc,const double* stepSize,const arma::mat *A){
  arma::mat At=*A+*stepSize*(*direc);
  arma::mat retract_Q,retract_R;
  arma::qr_econ(retract_Q,retract_R,At);
  if(retract_R(0,0)<0){
    retract_Q=-retract_Q;
  }
  return(At=retract_Q);
}

double objA(const arma::mat *Y,const arma::mat *X,const arma::mat *Theta,
              const arma::mat *B, const arma::mat *A){
  double obj=0.0;
  arma::mat diff=nLogLikelihood(Y,X,B,Theta,A);
  //diff %=diff;
  obj=arma::accu(diff);
  return obj;
}


//[[Rcpp::export]]
SEXP updateA(const SEXP Y_r, const SEXP X_r, const SEXP Theta_r,
                           const SEXP B_r,const SEXP A_r,
                           const SEXP tol_r,
                           const SEXP Max_r,
                           const SEXP subMax_r,
                           SEXP beta_r,
                           SEXP alpha_r,
                           SEXP sigma_r){
  //data input
  double tol=as<double>(tol_r);
  int Max=as<int>(Max_r);
  int subMax=as<int>(subMax_r);
  double sigma=as<double>(sigma_r);
  double alpha=as<double>(alpha_r);
  double beta=as<double>(beta_r);
  NumericMatrix Y_s(Y_r),X_s(X_r), Theta_s(Theta_r),B_s(B_r),A_s(A_r);
  arma::mat Y(Y_s.begin(),Y_s.nrow(),Y_s.ncol(),false);
  arma::mat X(X_s.begin(),X_s.nrow(),X_s.ncol(),false);
  arma::mat Theta(Theta_s.begin(),Theta_s.nrow(),Theta_s.ncol(),false);
  arma::mat B(B_s.begin(),B_s.nrow(),B_s.ncol(),false);
  arma::mat A(A_s.begin(),A_s.nrow(),A_s.ncol(),false);


  // controls
  int iter=0,innerIter;
  double tol_now,tol_now2;
  arma::mat *A_temp=new arma::mat;
  arma::mat *steepD=new arma::mat;
  double obj,obj_temp,stepSize,eAscent,obj_outer,obj_Des;
  obj=objA(&Y,&X,&Theta,&B,&A);
  arma::mat A_outer;
  vector<double> save_obj;
  vector<arma::mat> save_A;
  *A_temp=A;
  bool flag=true;
  Rcout<<"Before Updating A, the objetive value is "<<obj<<endl;
  // Steepest Ascent and Armijo back tracking
  do{
    // control of inner Armijo back tracking
    iter++;
    stepSize=alpha;
    obj_Des=-1.0;
    obj_outer=obj;
    A_outer=A;
    // negative gradient.
    *steepD=-gradA(&Y,&X,&Theta,&B,&A);
    *steepD=projectionA(steepD,&A);
    eAscent=sigma*arma::dot(*steepD,*steepD);

    innerIter=0;
    tol_now=-1.0;
    do{
        innerIter++;
        stepSize=beta*stepSize;
        eAscent=beta*eAscent;
        *A_temp=retractA(steepD,&stepSize,&A);
        obj_temp=objA(&Y,&X,&Theta,&B,A_temp);
        save_obj.push_back(obj_temp);
        save_A.push_back(*A_temp);
      }while (innerIter<subMax && (obj-obj_temp)<eAscent ); //out of Armijo
    obj_Des=obj_outer-obj_temp;
    if(obj_Des<0.0) {
      obj=*min_element(save_obj.begin(),save_obj.end());
      A=save_A[distance(save_obj.begin(),min_element(save_obj.begin(),save_obj.end()))];
    } else{
      A=*A_temp;
      obj=objA(&Y,&X,&Theta,&B,&A);
    }
    // clean the container
    vector<double>().swap(save_obj);
    vector<arma::mat>().swap(save_A);

    tol_now=arma::norm((A-A_outer),"fro");
    tol_now2=abs(obj_outer-obj);
    if(tol_now < tol || tol_now2 < tol) flag=false;


  } while (iter<Max && flag); // Out of iteration
  delete A_temp;
  delete steepD;
  Rcout<<"After Updating A, the objetive value is "<<obj<<endl;
//   Rcout<<"After Updating A, the iteration is "<<iter<<endl;
//   Rcout<<"After Updating A, the sub iteration is "<<innerIter<<endl;
  return(wrap(A));
}


