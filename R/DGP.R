
DGP = function(s,a,sigStrength,rho,n,noise,var,d,intercept,DGP.seed=NULL){

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
  epsilon=rnorm(n,0,noise)
  YAll=intercept*rep(1,n)+XAll%*%trueBeta+epsilon;


  return(list("Y"=YAll,"X"=XAll,"TRUE.idx"=TRUE.idx))

}

