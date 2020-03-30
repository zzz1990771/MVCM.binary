
M=100
result_matrix <- matrix(0,ncol=6,nrow = M)
X_list <- list()
for(m in 1:M){
#set.seed(floor(runif(1,1,1000000)))
library(fda)
## parameters:
n<-200
p<-51  # including the constant covariate
p0<-11 # non-zero covariates, choose to be the top covariates
q<-15
r<-4
rank<-r
k<-30
nbasis<-k
lambda<-50
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

# For spline produced functions:
A<-matrix(rnorm(k*r),nrow=k)
A_true<-qr.Q(qr(A))
c_true<-Theta_true%*%t(A_true)

## For nonspline produced functions, use sin/exp to make normality constraint easy to get.
# AtB_true<-matrix(rep(0,r*n),nrow=r) ## For rangeval=c(0,1) only
# for(i in 1:(r/2)){
#   for(j in 1:n)
#     AtB_true[i,j]<-sin(2*pi*i*Tpoints[j])*sqrt(2)
#   # AtB_true[i,j]<-exp(-1*i*Tpoints[j])/sqrt(1/(2*i)*(1-exp(-2*i)))
# }
# for(i in (r/2+1):r){
#   for(j in 1:n){
#     # AtB_true[i,j]<-exp(-1*i*Tpoints[j])/sqrt(1/(2*i)*(1-exp(-2*i)))
#     AtB_true[i,j]<-cos(2*pi*i*Tpoints[j])*sqrt(2)
#   }
# }



# Y<-matrix(rep(0,q*n),nrow=q)
X<-t(MASS::mvrnorm(n=n,mu=rep(0,p-1),Sigma=sigmaX))
X<-rbind(rep(1,n),X)

# For spline produced functions:
Y<-sapply(1:n, function(x){
  XThetaAtB <- kronecker(t(X[,x]),diag(q))%*%Theta_true%*%t(A_true)%*%B[,x]
  mu <- apply(XThetaAtB,1,function(x) {1 / (1 + exp(-x))})
  rbinom(n = q,size = 1,prob = mu)
})



## for error term:
 #error<-matrix(rnorm(n*q,sd=sigma),nrow=q)
# error_pca = matrix(rnorm(n*r, sd = sigma), nrow = r)
# error<-matrix(0,q,n)

## For non-spline produced functions:
# Y<-sapply(1:n,function(x){
#   
#   XThetaAtB <- kronecker(t(X[,x]),
#             diag(q))%*%Theta_true%*%(AtB_true[,x]+error_pca[,x])
#   mu <- apply(XThetaAtB,1,function(x) {1 / (1 + exp(-x))})
#   Y[,x]=rbinom(n = q,size = 1,prob = mu)
#     
#     
# })


system.time(pilot<-MVCM.binary::pilot_call(Y=Y,X=X,B=B,p=p,q=q,rank=rank))



result2<-solveAll(ThetaStart=NULL,Y=Y,X=X,tolTheta=tol,MaxItTheta=MaxIt
                 ,lambda=45,gamma = 2.0
                 ,rank=r,tolAll=tol,MaxItAll=MaxIt,tolA=tol,MaxItA=MaxIt,tau=tau
                 ,c_pilot=pilot,Tpoints=Tpoints,nbasis=k,rangeval=rangeval
                 ,grid=grid, plot=T,nplots=1,method="scad")


Y1hat <- sapply(1:ncol(Y), function(x){
  uhat <- kronecker(t(X[,x]),
                   diag(nrow(Y)))%*%(result2$Theta)%*%t(result2$A)%*%B[,x]
  (1 / (1 + exp(-uhat)))
})

auc <- AUC::auc(AUC::roc(as.vector(Y1hat),as.factor(as.vector(Y))))




compare_theta <- function(Theta_fitted,Theta_true,tolerance = 0.05){
  N <- dim(Theta_true)[1]
  selected_index <- !(apply(abs(Theta_fitted),1,sum)<=tolerance)
  true_index <- !(apply(abs(Theta_true),1,sum)<=tolerance)
  TPR <- sum((selected_index==1)&(true_index==1))/N
  FPR <- sum((selected_index==1)&(true_index==0))/N
  TNR <- sum((selected_index==0)&(true_index==0))/N
  FNR <- sum((selected_index==0)&(true_index==1))/N
  return(c(TPR=TPR,FPR=FPR,TNR=TNR,FNR=FNR))
}

result_matrix[m,1:4] <- compare_theta(result2$Theta,Theta_true)

result_matrix[m,5] <- sum((result2$Theta-Theta_true)^2)
result_matrix[m,6] <- auc

}

apply(result_matrix,2,mean)


