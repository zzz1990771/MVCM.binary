### test the work of updateTheta
updateTheta<-function(Theta,Y,X,B,A,lambda,tol,MaxIt,tau=1,m=1,method=1,a=3.7){
  update_Theta(Theta,Y,X,B,A,lambda,tol,MaxIt,tau,m,method,a)
}



pilot_call<-function(Y,X,B,p,q,rank){
  ## set up pilot for adaptive lasso
  ## potentially use scope in CV to make it efficient.
  n<-dim(Y)[2]
  k<-dim(B)[1]  ## i.e., nbasis in model setting
  stopifnot(dim(X)[2]==n && dim(X)[1]==p && dim(Y)[1]==q)
## randomly choose the space where C begin with
  C=matrix(0,p*q,k)
  diag(C)[sample(1:(min(p*q,k)),rank)]=rep(1,rank)
  C_in_r<-pilotSolve(Y,X,B,C,1e-6,30,rank)
  return(C)
}


#' solve Theta and A
#'
#' @param ThetaStart Matrix. The startpoint of Theta, NULL is default.
#' @param Y A matrix.
#' A \eqn{q \times n} matrix, where the number of rows is the number of variables.
#' @param X A matrix. A \eqn{(p+1) \times n} matrix, where the number of rows is the number of variabels.
#' @param B A matrix. A \eqn{K \times r} matrix, where the \eqn{i}-th column is a plug-in spline basis at \eqn{T_i}.
#' @param tolTheta A number. \eqn{\epsilon_\Theta} is used in updating \eqn{\Theta} when \eqn{A} is fixed.
#' @param MaxItTheta A number. The maximum iteration for updating \eqn{\Theta} when \eqn{A} is fixed.
#' @param lambda A number. One of adaptive group lasso tuning parameters, \eqn{\lambda >0}.
#' @param gamma A number. One of adaptive group lasso tuning parameters, \eqn{\gamma=0} is equivalent to regular lasso.
#' @param rank A number. The number of principal components, used to give a random \eqn{\Theta}.
#' @param tolAll A number. A measure of two \eqn{\Theta A^T}, it is used to determine the overall
#' @param MaxItAll A number. The maximum iteration for updating \eqn{\Theta A^T}.
#' @param tolA A number. Tolerance number for updating A.
#' @param MaxItA A number. The maximum iteration for updating \eqn{A}.
#' @param seed A number. Seed value for randomness.
#' @param ... Other parameters. Including control parameter in Stiefel Manifold.
#'
#' @return A list including \eqn{\Theta} and \eqn{A}.
#' @export
#'
#' @examples
SolveThetaA <-
  function(ThetaStart=NULL,AStart=NULL,Y,X,B,tolTheta=1e-6,MaxItTheta=50,lambda,gamma = 0
           ,seed=0819,rank,tolAll=1e-6,MaxItAll=20,tolA=1e-6,MaxItA=50,tau=1,m=1,method="LASSO",a=3.7,...) {
    set.seed(seed) ## For generating ungiven ThetaStart
    n <- dim(Y)[2]
    q <- dim(Y)[1]
    p <- dim(X)[1]  # including the constant item, p+1 in data, p in computation
    pOverm <- p/m  ## must be integer, used for group dummy variable
    k <- dim(B)[1]
    r <- rank
    stopifnot(dim(X)[2] == n,r <= min(p*q,k),dim(B)[2]==n)
    if (is.null(ThetaStart))
      Theta = matrix(rnorm(n=p*q*rank),nrow = p*q)
    else {
      stopifnot(dim(ThetaStart) == c(p*q,r))
      Theta = ThetaStart
    }

    ## setup A's environment.
    if (is.null(AStart)){
      A<-matrix(0,k,r)
      diag(A)<-rep(1,r)
    } else {
      stopifnot(dim(AStart) == c(k,r))
      A = AStart
      A=qr.Q(qr(A))
    }
    if("iterSubMax" %in% names(list(...))){
      iterSubMax=list(...)[["iterSubMax"]]
    } else {
      iterSubMax=50
    }
    if("beta" %in% names(list(...))){
      beta=list(...)[["beta"]]
    } else {
      beta=min(0.5,10/n)
    }
    if("alpha" %in% names(list(...))){
      alpha=list(...)[["alpha"]]
    } else {
      alpha=5
    }
    if("sigma" %in% names(list(...))){
      sigma=list(...)[["sigma"]]
    } else {
      sigma=0.6
    }

    ## set up lambda and related:
    gamma_all <- gamma[1]
    lambda_all<-numeric(pOverm)
    # lambda_all <- numeric(p)
    lambda=as.double(lambda)
    if(lambda[1]<0){ ## for god sake class(seq(55,56,by=1))=="numeric" and class(55)=="numeric",but class(55:56)=="numeric" is false
      cat("Your input lambda is not valid, \n"
          , "It must be no less than 0.\n"
          , "reset to 10 this time.")
      lambda[1]<-10 ## this is very dangerous use in cluster, cannot get print in cluster use.
    }
    if(! toupper(method) %in% c("LASSO","SCAD")){
      cat("Choose one of (adaptive) group lasso or scad, \n"
          ,"The default is (adaptive) group lasso. \n")
      method<-"LASSO"
    }
    if(toupper(method)=="LASSO"){
      method_in_c<-1
        if (gamma_all < 0 || class(gamma_all) != "numeric") {
          cat(
            "Your input gamma is not valid, \n gamma=0 for regular lasso,\n
            and gamma>0 for adaptive. \n Set to regular lasso this time."
          )
          gamma_all = 0
        }
        if (gamma_all == 0) {
          lambda_all <- rep(lambda[1],pOverm) ## lambda always has length pOverm.
        } else {
          if(!("c_pilot" %in% names(list(...)))){
            ## c=Theta%*%t(A)
            c_pilot<-pilot_call(Y,X,B,p,q,rank)
          } else if(is.null(list(...)[["c_pilot"]])){
            c_pilot<-pilot_call(Y,X,B,p,q,rank)
          } else {
            c_pilot<-list(...)[["c_pilot"]]
          }
          ## adaptive lasso:
          for (i in 1:pOverm) {
            lambda_all[i] <- lambda[1] * ( sum(c_pilot[((i-1)*m*q+1):(i*m*q),]^2) ^ (-gamma_all/2))
          }
          # lambda_all=rep(lambda_short,each=m)
        }
    } else {
      method_in_c<-2
      if(a<=2){
        cat("a must be greater than 2. \n"
            ,"The default is 3.7")
        a<-3.7
      }
      lambda_all<-rep(lambda[1],pOverm)
    }

    this_tol <- 2*abs(tolAll)
    this_iter <- 0
    C<-Theta%*%t(A)
    ## start to do alternating direction:
    while(this_tol > tolAll && this_iter < MaxItAll){
      C_pre<-C
      ## update Theta:
      Theta_in_r<-update_Theta(Theta,Y,X,B,A,lambda_all,tolTheta,MaxItTheta,tau,m,method_in_c,a)

      ## update A:
      A_in_r<-updateA(Y,X,Theta,B,A,tolA,MaxItA,iterSubMax,beta,alpha,sigma)
      C<-Theta%*%t(A)
      this_tol<-sum((C-C_pre)^2,na.rm=T)
      this_iter<-this_iter+1
    }
    return(list(Theta=Theta,A=A,coefficients=Theta%*%t(A)%*%B,lambda=lambda_all))
  }


solveAll<-function(Y,X,Tpoints,nbasis=20,grid=100,ThetaStart=NULL,AStart=NULL,tolTheta=1e-6,MaxItTheta=50,lambda
                    ,gamma = 0,seed=0819,rank,tolAll=1e-6,MaxItAll=20,tolA=1e-6,MaxItA=50,tau=1,m=1,method="LASSO",a=3.7
                    ,plot=T,nplots=5,...){
  ## To extract sublist, we should use list[subnames] not [[]].
  paraToBspline<-list(...)[names(list(...))[names(list(...)) %in% names(formals(create.bspline.basis))]]
  splineInfo<-do.call("generateB",args=c(list(Tpoints=Tpoints,nbasis=nbasis,grid=grid)
                                         ,paraToBspline))
  B<-splineInfo$B
  ThetaA<-SolveThetaA(ThetaStart=ThetaStart,AStart=AStart,Y=Y,X=X,B=B,tolTheta=tolTheta
                       ,MaxItTheta=MaxItTheta,lambda=lambda,gamma = gamma,seed=seed
                       ,rank=rank,tolAll=tolAll,MaxItAll=MaxItAll,tolA=tolA,MaxItA=MaxItA
                       ,tau=tau,m=m,method=method,a=a,...)
  cat(range(splineInfo$gridlocation),"\n")
  if(plot==T){
    for(i in 1:nplots){
      plot(splineInfo$gridlocation
           ,((ThetaA$Theta)%*%t(ThetaA$A)%*%(splineInfo$Bgrid))[i,]
           ,xlab="Time",ylab="Function Values"
           ,main=paste0("The No.",i," Coefficients Functions"),type="l")
    }
  }
  ThetaA
}
