get.betas <- function(MI.out,spline.list,eval.times){
  #browser()
  Q.tilde=eval.spline.list(spline.list,eval.times)
  beta.tilde=Q.tilde%*%MI.out$alpha
  p=length(spline.list)
  beta.mat=matrix(as.numeric(beta.tilde),ncol=p,byrow=F)
  colnames(beta.mat)=MI.out$beta.names
  beta.mat
}
