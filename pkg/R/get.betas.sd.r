get.betas.sd <- function(MI.out,spline.list,eval.times){
  #browser()
  Q.tilde=eval.spline.list(spline.list,eval.times)
  sd.tilde=Q.tilde%*%MI.out$alpha.sd
  p=length(spline.list)
  beta.mat=matrix(as.numeric(sd.tilde),ncol=p,byrow=F)
  colnames(beta.mat)=MI.out$beta.names
  beta.mat
}
