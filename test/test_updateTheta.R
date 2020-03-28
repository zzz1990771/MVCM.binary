#########################
## Test SolveThetaA  ####
#########################
library(fda)
## parameters:
set.seed(1225)
n<-200
p<-51  # including the constant covariate
p0<-11 # non-zero covariates, choose to be the top covariates
q<-15
r<-4
rank<-r
k<-30
nbasis<-k
lambda_e<-100000
gamma<-0
tol<-0.0000001
MaxIt<-100
tau<-1
rhoX<-0.3 # the base coefficients of X, i.e., cov*(x_j1,x_j2)=rho^|j1-j2|
sigmaX<-matrix(rep(0,(p-1)*(p-1)),nrow=p-1) ## for covariance matrix of X
for(i in 1:(p-1)){
  for(j in 1:(p-1)){
    sigmaX[i,j]<-rhoX^abs(i-j)
  }
}
sigma<-0.1 ## for the error terms
rangeval<-c(0,1)
grid<-100

## Theta are randomly generated as well:
Theta<-matrix(rnorm(p*q*r),nrow=p*q)

Theta_true<-Theta
Theta_true[-c(1:(p0*q)),]<-0


Tpoints<-runif(n)*(rangeval[2]-rangeval[1])+rangeval[1]
splineInfo<-generateB(Tpoints=Tpoints,nbasis=k,rangeval=rangeval)
B<-generateB(Tpoints=Tpoints,nbasis=k,rangeval=rangeval)$B # For pilot generating

## For spline produced functions:
# A<-matrix(rnorm(k*r),nrow=k)
# A_true<-qr.Q(qr(A))
# c_true<-Theta_true%*%t(A_true)

## For nonspline produced functions, use sin/exp to make normality constraint easy to get.
AtB_true<-matrix(rep(0,r*n),nrow=r) ## For rangeval=c(0,1) only
for(i in 1:(r/2)){
  for(j in 1:n)
    AtB_true[i,j]<-sin(2*pi*i*Tpoints[j])*sqrt(2)
  # AtB_true[i,j]<-exp(-1*i*Tpoints[j])/sqrt(1/(2*i)*(1-exp(-2*i)))
}
for(i in (r/2+1):r){
  for(j in 1:n){
    # AtB_true[i,j]<-exp(-1*i*Tpoints[j])/sqrt(1/(2*i)*(1-exp(-2*i)))
    AtB_true[i,j]<-cos(2*pi*i*Tpoints[j])*sqrt(2)
  }
}

## For spline produced functions:
# Y<-sapply(1:n, function(x){Y[,x]=kronecker(t(X[,x])
#                                            ,diag(q))%*%Theta_true%*%t(A_true)%*%B[,x]+error[,x]})


Y<-matrix(rep(0,q*n),nrow=q)
X<-t(MASS::mvrnorm(n=n,mu=rep(0,p-1),Sigma=sigmaX))
X<-rbind(rep(1,n),X)
## for error term:
 error<-matrix(rnorm(n*q,sd=sigma),nrow=q)
 error_pca = matrix(rnorm(n*r, sd = sigma), nrow = r)
# error<-matrix(0,q,n)

## For non-spline produced functions:
Y<-sapply(1:n,function(x){
  
  XThetaAtB <- kronecker(t(X[,x]),
            diag(q))%*%Theta_true%*%(AtB_true[,x]+error_pca[,x])+error[,x]
  mu <- apply(XThetaAtB,1,function(x) {1 / (1 + exp(-x))})
  Y[,x]=rbinom(n = q,size = 1,prob = mu)
    
    
})

# gradH<-function(Theta){
#   temp=matrix(rep(0,p*q*r),nrow=p*q)
#   for(x in 1:n){
#     temp=temp-2.0*(kronecker(X[,x],diag(q)))%*%(Y[,x]-kronecker(t(X[,x])
#                                                                 ,diag(q))%*%Theta%*%t(A)%*%B[,x])%*%t(B[,x])%*%A
#   }
#   return(temp)
# }
#
# gradH_atTheta<-gradH(Theta=Theta)

system.time(pilot<-MVCM.binary::pilot_call(Y=Y,X=X,B=B,p=p,q=q,rank=rank))

# anotherPilot<-solveAll2(ThetaStart=NULL,Y=Y,X=X,tolTheta=tol,MaxItTheta=MaxIt
#                        ,lambda=0,gamma = 0
#                        ,rank=r,tolAll=tol,MaxItAll=100,tolA=tol,MaxItA=MaxIt,tau=tau,seed =0819
#                        ,Tpoints=Tpoints,nbasis=k,rangeval=rangeval
#                        ,grid=grid, plot=T,nplots=1)
# pilot<-anotherPilot$Theta%*%t(anotherPilot$A)

# anotherPilot<-pilot_call2(Y=Y,X=X,B=B,p=p,q=q,rank=rank)

# sum((pilot-Theta_true%*%t(A_true))^2)

# result<-SolveThetaA(ThetaStart=Theta,Y=Y,X=X,B=B,tolTheta=tol,MaxItTheta=MaxIt,lambda=lambda_e/100,gamma = 1.5
#            ,rank=r,tolAll=tol,MaxItAll=MaxIt,tolA=tol,MaxItA=MaxIt,tau=tau,seed =0819,c_pilot=pilot)
#
# test_result<-result$Theta%*%t(result$A)
# sum((test_result-Theta_true%*%t(A_true))^2)
# result<-solveAll(ThetaStart=NULL,Y=Y,X=X,tolTheta=tol,MaxItTheta=MaxIt
#                ,lambda=lambda_e/10,gamma = 2.0
#               ,rank=r,tolAll=tol,MaxItAll=MaxIt,tolA=tol,MaxItA=MaxIt,tau=tau,seed =0819
#               ,c_pilot=pilot,iterSubMax=30,Tpoints=Tpoints,nbasis=k,rangeval=rangeval
#               ,grid=grid, plot=T,nplots=1)
# test_result<-result$Theta%*%t(result$A)
# sum((test_result-Theta_true%*%t(A_true))^2)

result1<-solveAll(ThetaStart=NULL,Y=Y,X=X,tolTheta=tol,MaxItTheta=MaxIt
                  ,lambda=10,gamma = 2.0
                  ,rank=r,tolAll=tol,MaxItAll=MaxIt,tolA=tol,MaxItA=MaxIt,tau=tau,seed =0819
                  ,c_pilot=pilot,Tpoints=Tpoints,nbasis=k,rangeval=rangeval
                  ,grid=grid, plot=T,nplots=1,method="lasso")
result1$Theta
result2<-solveAll(ThetaStart=NULL,Y=Y,X=X,tolTheta=tol,MaxItTheta=MaxIt
                 ,lambda=50,gamma = 2.0
                 ,rank=r,tolAll=tol,MaxItAll=MaxIt,tolA=tol,MaxItA=MaxIt,tau=tau,seed =0819
                 ,c_pilot=pilot,Tpoints=Tpoints,nbasis=k,rangeval=rangeval
                 ,grid=grid, plot=T,nplots=1,method="scad")
result2$Theta

plot(Tpoints,result1$coefficients[1,])
## CV Error:

system.time(thisER<-srrrVcmCvError(Kfolder = 5, Y=Y, X=X, Tpoints=Tpoints, nbasis=k
               , grid=grid, ThetaStart=NULL, tolTheta=tol
               , MaxItTheta=MaxIt
               , lambda=lambda_e/10, gamma=2.0, rank=r
               , tolAll=tol, MaxItAll=MaxIt,tolA=tol
               , MaxItA=MaxIt, tau=tau, seed =0819, plot=T
               , c_pilot=pilot
               , nplots=1,rangeval=rangeval,method="lasso"))
system.time(thisER2<-srrrVcmCvError(Kfolder = 5, Y=Y, X=X, Tpoints=Tpoints, nbasis=k
                                   , grid=grid, ThetaStart=NULL, tolTheta=tol
                                   , MaxItTheta=MaxIt
                                   , lambda=lambda_e/10, gamma=2.0, rank=r
                                   , tolAll=tol, MaxItAll=MaxIt,tolA=tol
                                   , MaxItA=MaxIt, tau=tau, seed =0819, plot=T
                                   , c_pilot=pilot
                                   , nplots=1,rangeval=rangeval,method="lasso",warmStart=T))

system.time(thisER3<-srrrVcmCvError(Kfolder = 5, Y=Y, X=X, Tpoints=Tpoints, nbasis=k
                                    , grid=grid, ThetaStart=matrix(0,p*q,r), tolTheta=tol
                                    , MaxItTheta=MaxIt
                                    , lambda=450, gamma=2.0, rank=r
                                    , tolAll=tol, MaxItAll=MaxIt,tolA=tol
                                    , MaxItA=MaxIt, tau=tau, seed =0819, plot=T
                                    , c_pilot=pilot
                                    , nplots=1,rangeval=rangeval,method="scad",warmStart=T))

system.time(bestER3<-srrrVcmCv(Kfolder = 5, Y=Y, X=X, Tpoints=Tpoints, nbasis=k
                                    , grid=grid, ThetaStart=NULL, tolTheta=tol
                                    , MaxItTheta=MaxIt
                                    , lambda=c(450,500,550), gamma=2.0, rank=3:5
                                    , tolAll=0.001, MaxItAll=MaxIt,tolA=tol
                                    , MaxItA=MaxIt, tau=tau, seed =0819, plot=T
                                    , c_pilot=pilot
                                    , nplots=1,rangeval=rangeval,method="scad",warmStart=T))
system.time(bestER4<-srrrVcmCv(Kfolder = 5, Y=Y, X=X, Tpoints=Tpoints, nbasis=k
                               , grid=grid, ThetaStart=NULL, tolTheta=tol
                               , MaxItTheta=MaxIt
                               , lambda=c(1000,2000), gamma=2.0, rank=4
                               , tolAll=1e-5, MaxItAll=MaxIt,tolA=tol
                               , MaxItA=MaxIt, tau=tau, seed =0819, plot=T
                               , c_pilot=pilot
                               , nplots=1,rangeval=rangeval,method="lasso",warmStart=T))
system.time(bestER5<-srrrVcmCv(Kfolder = 5, Y=Y, X=X, Tpoints=Tpoints
                               #, nbasis=k
                               , nbasis=15
                               , grid=grid, ThetaStart=NULL, tolTheta=tol
                               , MaxItTheta=MaxIt
                               , lambda=c(450,500,550), gamma=2.0, rank=3:5
                               , tolAll=0.001, MaxItAll=MaxIt,tolA=tol
                               , MaxItA=MaxIt, tau=tau, seed =0819, plot=T
                               , c_pilot=pilot
                               , nplots=1,rangeval=rangeval,method="scad",warmStart=T
                               , norder=3)) ## quadratic spline
## If it has rangeval, it must be included
