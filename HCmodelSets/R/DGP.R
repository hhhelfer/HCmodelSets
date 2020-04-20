DGP = function(s,a,sigStrength,rho,n,noise=NULL,var,d,intercept,type.response="N",DGP.seed=NULL, scale=NULL, shape = NULL, rate=NULL){

  if(type.response=="S" & is.null(scale)==TRUE){
    stop('You choose cox family! Therefore you must provide parameters scale, shape!')
  }
  if(type.response=="S" & is.null(shape)==TRUE){
    stop('You choose cox family! Therefore you must provide parameters scale, shape!')
  }
  if(type.response=="N" & is.null(shape)==FALSE){
    stop('Scale and shape parameters will not be used since type.response is gaussian!')
  }
  if(type.response=="N" & is.null(scale)==FALSE){
    stop('Scale and shape parameters will not be used since type.response is gaussian!')
  }
  if(type.response=="N" & is.null(rate)==FALSE){
    stop('Scale, shape and rate parameters will not be used since type.response is gaussian!')
  }
  if(type.response=="S" & !is.null(noise)==TRUE){
    stop('You choose cox family! Parameter noise is not use!')
  }
  if(!is.element(type.response,c("S","N"))){
    stop('Only supports gaussian (N) or survival data (S)!')
  }

  cov1=rho*rep(1,a+s)+(1-rho)*diag(a+s);
  covMatrixInit = rbind(cbind(cov1, matrix(0,s+a,d-(s+a))),cbind(matrix(0,d-(s+a),s+a),diag(d-(s+a))))
  covMatrix = diag(sqrt(var) * rep(1,d)) %*% covMatrixInit %*% diag(sqrt(var)* rep(1,d))
  trueBetaInit = c(sigStrength * rep(1,s) , rep(0,d-s))

  if(!is.null(DGP.seed)){
    set.seed(DGP.seed)
  }

  #### creating DGP
  permuteVec=sample(d);
  trueBeta=trueBetaInit[permuteVec] # permute the rows of the initial beta vector
  TRUE.idx = which(trueBeta!=0)
  I = diag(d)
  permMatrix=I[permuteVec,]
  covPerm=permMatrix%*%covMatrix%*%(solve(permMatrix)) # permute rows and columns of the covariance matrix accordingly.
  XAll = mvtnorm::rmvnorm(n,rep(0,d),covPerm)

  if(type.response!="S"){

    epsilon=rnorm(n,0,noise)
    YAll=intercept*rep(1,n)+XAll%*%trueBeta+epsilon;

    return(list("Y"=YAll,"X"=XAll,"TRUE.idx"=TRUE.idx))

  } else if(type.response=="S" & is.null(rate)==FALSE){
    v <- runif(n=nrow(XAll))
    Tlat <- (- log(v) / (scale * exp(XAll %*% trueBeta)))^(1 / shape)
    C <- rexp(n=nrow(XAll), rate=rate)
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C)

    return(list("Y"=time,"X"=XAll,"status"=status,"TRUE.idx"=TRUE.idx))

  } else if(type.response=="S" & is.null(rate)==TRUE){
    v <- runif(n=nrow(XAll))
    Tlat <- (- log(v) / (scale * exp(XAll %*% trueBeta)))^(1 / shape)
    status <- rep(0,length(Tlat))

    return(list("Y"=Tlat,"X"=XAll,"status"=status,"TRUE.idx"=TRUE.idx))

  }
}

