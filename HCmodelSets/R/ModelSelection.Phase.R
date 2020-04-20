ModelSelection.Phase = function(X,Y, list.reduction, family=gaussian, signif=0.01, sq.terms=NULL, in.terms=NULL, modelSize=NULL, Cox.Hazard = FALSE){

  if(Cox.Hazard==TRUE){
    Y = as.matrix(Y)
    if(ncol(Y)!=2){
      stop('You choose cox family! Therefore you must provide the survival times in column 1 and a binary vector on column 2!')
    }
  } else{ Y = as.numeric(Y) }

  if(Cox.Hazard==FALSE & ncol(as.matrix(Y))==2){
    stop('You need to specify Cox.Hazard==TRUE to run cox model!')
  }

  if(!is.null(sq.terms)){
    X.SQ = X[,sq.terms]^2
    SQ.names = paste(sq.terms,"^2")
  }

  if(is.null(sq.terms)){
    X.SQ = SQ.names = NULL
  }

  if(!is.null(in.terms)){
    X.ITER = iter.names = NULL
    for(ii in 1:nrow(in.terms)){
      X.ITER = cbind(X.ITER,X[,in.terms[ii,1]]*X[,in.terms[ii,2]])
      iter.names = c(iter.names,paste(in.terms[ii,1],"*",in.terms[ii,2]))
    }
  }

  if(is.null(in.terms)){
    X.ITER = iter.names = NULL
  }

  X.full = cbind(X[,list.reduction],X.SQ,X.ITER)
  setSelected = c(list.reduction,SQ.names,iter.names)

  s=ncol(X.full)
  if(is.null(modelSize)){ modelSize=min(5,length(setSelected)) }

  ### Error message
  if(modelSize>7){
    stop('Sorry, this version only support model sizes<8')
  }

  n=nrow(X.full)

  if(Cox.Hazard==FALSE){
    dist.family = glm(Y~X.full, family = family)$family$family

    if(dist.family != "gaussian"){
      tag.family = "R"
    } else{
      tag.family = "G"
    }} else{
      tag.family = "R"
    }

  goodModels = list()

  if(tag.family == "G"){

    if(Cox.Hazard==TRUE){
      residBig = residuals(coxph(Surv(Y[,1], Y[,2])~X.full))
    } else{
      residBig = residuals(glm(Y~X.full, family = family))
    }
    RSSBig = t(residBig) %*%residBig ### same as dev comprehensive
    dfBig = n-(length(setSelected)+1)

    for (j in 1:modelSize){
      combinationMatrix=t(combn(1:length(setSelected),j))
      if(length(setSelected)>1){
        combinationMatrixNames=cbind(t(combn(setSelected,j)))
      } else{
        combinationMatrixNames=as.matrix(setSelected)
      }
      logicFitVectorF=matrix(0,nrow(combinationMatrix))
      for(l in 1:nrow(combinationMatrix)){
        XSelect = X.full[,combinationMatrix[l,]]
        if(Cox.Hazard==TRUE){
          residSmall = residuals(coxph(Surv(Y[,1], Y[,2])~XSelect))
        } else{
          residSmall = residuals(glm(Y~XSelect, family = family))
        }
        RSSSmall = t(residSmall) %*% residSmall
        fCrit = qf(1-signif,df1 = s - ncol(cbind(XSelect)),df2 = n-s)
        if(is.nan(fCrit)){
          logicFitVectorF[l]=1
        } else if (((RSSSmall - RSSBig)/(s - ncol(cbind(XSelect))))/(RSSBig/(n-ncol(cbind(XSelect))))<=fCrit){
          logicFitVectorF[l]=1
        }
      }
      if(j==1){goodModels[[paste('Model Size', j)]] = as.matrix(combinationMatrixNames[which(logicFitVectorF>0),])
      } else{
        goodModels[[paste('Model Size', j)]] = combinationMatrixNames[which(logicFitVectorF>0),]
      }
    }

  }

  if(tag.family == "R"){

    dfBig = s

    if(Cox.Hazard==TRUE){
      LBig = coxph(Surv(Y[,1], Y[,2])~X.full)$loglik[2]
    } else{
      LBig = logLik(glm(Y~X.full, family = family))
    }
    for (j in 1:modelSize){
      combinationMatrix=t(combn(1:length(setSelected),j))
      if(length(setSelected)>1){
        combinationMatrixNames=cbind(t(combn(setSelected,j)))
      } else{
        combinationMatrixNames=as.matrix(setSelected)
      }
      logicFitVectorF=matrix(0,nrow(combinationMatrix))

      for(l in 1:nrow(combinationMatrix)){
        XSelect = cbind(X.full[,combinationMatrix[l,]])
        dfSmall = ncol(XSelect)
        if(Cox.Hazard==TRUE){
          LSmall = coxph(Surv(Y[,1], Y[,2])~XSelect)$loglik[2]
        } else{
          LSmall = -glm(Y~XSelect, family = family)$deviance/2
        }
        chiCrit = qchisq(1-signif,df = dfBig - dfSmall)
        if(is.nan(chiCrit)){
          logicFitVectorF[l]=1
        } else if (-2*(LSmall-LBig)<=chiCrit){
          logicFitVectorF[l]=1
        }
      }
      if(j==1){goodModels[[paste('Model Size', j)]] = as.matrix(combinationMatrixNames[which(logicFitVectorF>0),])
      } else{
        goodModels[[paste('Model Size', j)]] = combinationMatrixNames[which(logicFitVectorF>0),]
      }
    }

  }

  return(list("goodModels"=goodModels))

}
