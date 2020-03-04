## spline preparation
generateB<-function(Tpoints,nbasis=20,grid=100,...){
  if(!"rangeval" %in% names(list(...))){
    rangeval<-range(Tpoints)
    bbasis<-create.bspline.basis(rangeval=rangeval,nbasis=nbasis,...)
  } else {
  rangeval<-list(...)[["rangeval"]]
  bbasis<-create.bspline.basis(nbasis=nbasis,...)
  }
  Bgrid<-eval.basis(seq(rangeval[1],rangeval[2],length.out=grid),bbasis)
  gridlocation<-seq(rangeval[1],rangeval[2],length.out=grid)
  ## Spline interpolations without normalization
  B<-t(eval.basis(Tpoints,bbasis))
  ## The normalization matrix:
  normalB<-sqrt(grid/(rangeval[2]-rangeval[1]))*t(solve(qr.R(qr(Bgrid))))
  B<-normalB%*%B
  Bgrid<-t(Bgrid) # cannot put this line go up.
  Bgrid<-normalB%*%Bgrid
  output<-list(B=B,Bgrid=Bgrid,gridlocation=gridlocation)
  return(output)
}
