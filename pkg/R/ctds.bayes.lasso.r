ctds.bayes.lasso <- function(sim.obj,spline.list,stack.static,stack.grad,conspecifics.list=NULL,crw=TRUE,path.aug.interval=0,spline.period=0,intercept.start=0,intercept.tune=1,alpha.start,alpha.tune.mat,lambda=NA,lambda.prior.r=1,lambda.prior.q=2,n.mcmc=10){


  ##
  ## Preliminaries
  ##

  

  loglik.pois <- function(z,tau,Phi,alpha,intercept){
    xb=intercept+as.numeric(Phi%*%Matrix(as.numeric(alpha),ncol=1))
    sum(z*xb-tau*exp(xb),na.rm=T)
  }

  

  P=length(alpha.start)
  
  alpha=alpha.start
  alpha[which(alpha==0)] <- .000001
  #alpha=Matrix(0,ncol=1,nrow=ncol(out$Phi))

  intercept=intercept.start
  
  lambda.2=(lambda.prior.r/lambda.prior.q)^2
  if(!is.na(lambda)){
    lambda.2=lambda^2
  }

  D=Diagonal(length(alpha))
  one=Matrix(1,nrow=1,ncol=length(alpha))

  alpha.save=matrix(NA,nrow=n.mcmc,ncol=length(alpha))
  lambda.save=rep(NA,n.mcmc)
  intercept.save=rep(NA,n.mcmc)

  accept=0
  
  for(iter in 1:n.mcmc){

    cat(iter," ")

    #browser()

    
    ##
    ## Sample continuous path and convert to CTDS path
    ##

    out=make.Phi.mat(sim.obj,spline.list,stack.static,stack.grad,conspecifics.list,crw=crw,path.aug.interval=path.aug.interval,spline.period=spline.period)
    Phi=out$Phi
    XX=out$XX
    z=out$z
    QQ=out$QQ
    tau=out$tau
    crawl=out$samp.new

    
    ##
    ## Sample s2
    ##

    s2=rinvgauss(P,sqrt(lambda.2/alpha^2),lambda.2)
    
    alpha.prior.var=diag(as.numeric(s2))

    ##
    ## Sample lambda
    ##

    if(is.na(lambda)){
      lambda.2=rgamma(1,lambda.prior.r+P,,1/(lambda.prior.q+.5*sum(s2)))
    }

    ##
    ## Sample alphas
    ##

    alpha.star=rmvnorm(1,alpha,sigma=alpha.tune.mat)
    intercept.star=rnorm(1,intercept,sqrt(intercept.tune))

    mh1=loglik.pois(z,tau,Phi,alpha.star,intercept.star)+dmvnorm(alpha.star,alpha,alpha.prior.var,log=TRUE)
    mh2=loglik.pois(z,tau,Phi,alpha,intercept)+dmvnorm(alpha,alpha.star,alpha.prior.var,log=TRUE)

    if(runif(1)<exp(mh1-mh2)){
      alpha=alpha.star
      intercept=intercept.star
      accept=accept+1
    }

   
    

    alpha.save[iter,] <- alpha
    intercept.save[iter] <- intercept
    lambda.save[iter] <- sqrt(lambda.2)
  }
  alpha=apply(alpha.save,2,mean)
  alpha.sd=apply(alpha.save,2,sd)
  list(alpha=alpha,alpha.sd=alpha.sd,alpha.save=alpha.save,intercept.save=intercept.save,lambda.save=lambda.save,accept=accept)
}
