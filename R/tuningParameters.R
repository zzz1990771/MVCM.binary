srrrVcmCvError <-
  function(Kfolder = 5,Y,X,Tpoints,nbasis=20,grid=100,ThetaStart=NULL, AStart=NULL
           ,tolTheta=1e-6,MaxItTheta=50,lambda=100
           ,gamma=0,seed=0819,rank,tolAll=1e-6,MaxItAll=20,tolA=1e-6,MaxItA=50
           ,tau=1,m=1,method="lasso",a=3.7,plot=F,nplots=0,warmStart=F,...) {
    ## Kfolder is the number of folds for cross-validation, usually 5 or 10.
    ## For all other paramters, see "solveAll" function.
    set.seed(seed)
    id <- sample(1:ncol(X)) # ncol is the # of samples
    X1 <- X[,id]; Y1 <- Y[,id];Tpoints1<-Tpoints[id]
    paraToBspline<-list(...)[names(list(...))[names(list(...)) %in% names(formals(create.bspline.basis))]]
    B1<-do.call("generateB",args=c(list(Tpoints=Tpoints1,nbasis=nbasis,grid=grid)
                                  ,paraToBspline))$B  ## Use Tpoints1 to interpolate.

    ## Group ID:
    group <- rep(1:Kfolder, ncol(X) / Kfolder + 1)[1:ncol(X)]
    real_folder <- unique(group) ## in case n is too small
    Y1hat <- matrix(0,nrow(Y),ncol(Y))
    ## calculate the first folder
    for (i in 1:1) {
      X_te <- X1[, group == i]
      X_tr <- X1[, group != i]
      Tpoints_te <- Tpoints1[group == i]
      Tpoints_tr <- Tpoints1[group != i]
      Y_tr <- Y1[, group != i]
      this_folder <-
        solveAll(
          Y=Y_tr, X=X_tr, Tpoints=Tpoints_tr, nbasis=nbasis, grid=grid
          , ThetaStart=ThetaStart, AStart=AStart, tolTheta=tolTheta, MaxItTheta=MaxItTheta
          ,lambda=lambda, gamma=gamma, seed =seed, rank=rank, tolAll=tolAll, MaxItAll=MaxItAll
          ,tolA=tolA, MaxItA=MaxItA, tau=tau, m=m,method=method,a=a, plot=plot, nplots=nplots, ...)
      Y1hat[,group==i] <- sapply(1:sum(group==i), function(x){
        Y1[,x]=kronecker(t(X_te[,x]),
                        diag(nrow(Y)))%*%(this_folder$Theta)%*%t(this_folder$A)%*%B1[,which(group==i)[x]]
      })
      cv_folder=sum((Y1hat[,group==i]-Y1[,group==i])^2)/(dim(Y1)[1]*sum(group==i))
    }
    if(warmStart==T){
      ThetaStart=this_folder$Theta
      AStart=this_folder$A
    }
    ### Then the rest folders
    for (i in 2:length(real_folder)) {
      X_te <- X1[, group == i]
      X_tr <- X1[, group != i]
      Tpoints_te <- Tpoints1[group == i]
      Tpoints_tr <- Tpoints1[group != i]
      Y_tr <- Y1[, group != i]
      this_folder <-
        solveAll(
          Y=Y_tr, X=X_tr, Tpoints=Tpoints_tr, nbasis=nbasis, grid=grid
          , ThetaStart=ThetaStart, AStart=AStart, tolTheta=tolTheta, MaxItTheta=MaxItTheta
          ,lambda=lambda, gamma=gamma, seed =seed, rank=rank, tolAll=tolAll, MaxItAll=MaxItAll
          ,tolA=tolA, MaxItA=MaxItA, tau=tau, m=m,method=method,a=a, plot=plot, nplots=nplots, ...)
      Y1hat[,group==i] <- sapply(1:sum(group==i), function(x){
        Y1[,x]=kronecker(t(X_te[,x]),
                         diag(nrow(Y)))%*%(this_folder$Theta)%*%t(this_folder$A)%*%B1[,which(group==i)[x]]
      })
      cv_folder=c(cv_folder,sum((Y1hat[,group==i]-Y1[,group==i])^2)/(dim(Y1)[1]*sum(group==i)))
    }
    cvError <- sum((Y1-Y1hat)^2)/(dim(Y)[1]*dim(Y)[2])

    return(c(lambda=lambda,gamma=gamma,rank=rank,cvError=cvError,cv_folder=cv_folder))
    ## CV error
  }


#### For lasso, we only export the group lasso, the adaptive group lasso is too slow
srrrVcmCv<-
  function(Kfolder = 5,Y,X,Tpoints
           ,lambda=c(100,1e3,1e4,1e5),gamma = c(0,1.5,2,3),rank=c(1,2,3,4,5)
           ,seed =0819,nbasis=20,grid=100,ThetaStart=NULL
           ,tolTheta=1e-6,MaxItTheta=50
           ,tolAll=1e-6,MaxItAll=20,tolA=1e-6,MaxItA=50
           ,tau=1,m=1,method="lasso",a=3.7,plot=F,nplots=0,warmStart=F,...) {
    if(! toupper(method) %in% c("LASSO","SCAD")){
      cat("Choose one of group 'lasso' or 'scad', \n"
          ,"The default is group 'lasso'. \n"
          ,"If the penalty is adaptive group lasso, \n"
          ,"it is suggested to use cluster to run 'srrrVcmCvError' parallelly. \n")
      method<-"LASSO"
    }
    ### currently, only rank and lambda, leave room to add tuning parameters in future
      lambda<-sort(lambda,decreasing=T)
     # gamma<-0 ## only do group lasso, no pilot
      gamma<-sort(gamma,decreasing=T)
      rank<-sort(rank,decreasing=T)
      i1 <- length(lambda)
      j1 <- length(gamma)
      r1 <- length(rank)
      best_err <- 1e8
      best_arg <- c(lambda[1],gamma[1],rank[1])

      ## pilot info
      c_pilot<-list()
      paraToBspline<-list(...)[names(list(...))[names(list(...)) %in% names(formals(create.bspline.basis))]]
      splineInfo<-do.call("generateB",args=c(list(Tpoints=Tpoints,nbasis=nbasis,grid=grid)
                                             ,paraToBspline))
      B<-splineInfo$B
      p<-dim(X)[1]
      q<-dim(Y)[1]
      k<-dim(B)[1]
      bestPilot=matrix(0,p*q,k)

      for(r in 1:r1){
        if(toupper(method)=="LASSO"){
          c_pilot[[r]]<-pilot_call(Y=Y,X=X,B=B,p=p,q=q,rank=rank[r])
        } else {
          c_pilot[[r]]<-matrix(0,p*q,k)
        }
        for(j in 1:j1){
          for(i in 1:i1){
            current_arg<-c(lambda[i],gamma[j],rank[r])
            current_err<-srrrVcmCvError(Kfolder = Kfolder, Y=Y, X=X, Tpoints=Tpoints, nbasis=nbasis
                                        , grid=grid, ThetaStart=ThetaStart, tolTheta=tolTheta
                                        , MaxItTheta=MaxItTheta
                                        , lambda=lambda[i], gamma=gamma[j], rank=rank[r]
                                        , seed =seed, tolAll=tolAll, MaxItAll=MaxItAll,tolA=tolA
                                        , MaxItA=MaxItA, tau=tau, m=m,method=method,a=a, plot=plot
                                        , nplots=nplots,warmStart=warmStart,c_pilot=c_pilot[[r]], ...)[4]
            if(current_err<best_err) {
              best_err<-current_err
              best_arg<-current_arg
              bestPilot=c_pilot[[r]]
            }
          }
        }
      }
    if(toupper(method)=="LASSO"){
      return(list(lambda=best_arg[1],gamma=best_arg[2],rank=best_arg[3],error=best_err,method=method,pilot=bestPilot))
    }
    ## The reason to add bestPilot in SCAD is that
    ## it may not converge in one algorithm circle, so add one circle to encourage convergence.
    ## But most cases, it does not work as we want.
    bestPilot=solveAll(Y=Y,X=X,Tpoints=Tpoints,nbasis=nbasis,grid=grid,ThetaStart=ThetaStart,tolTheta=tolTheta
                       ,MaxItTheta=MaxItTheta,lambda=best_arg[1], gamma = best_arg[2],rank=best_arg[3]
                       ,tolAll=tolAll,MaxItAll=MaxItAll,tolA=tolA,MaxItA=MaxItA,tau=tau,m=m,method="SCAD",a=3.7
                       ,plot=plot,nplots=nplots,...)
    return(list(lambda=best_arg[1],gamma=best_arg[2],rank=best_arg[3],error=best_err,method=method,pilot=bestPilot))
  }

