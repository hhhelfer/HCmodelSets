ModelSelection.Phase = function(X,Y, list.reduction, family=gaussian, signif=0.01, sq.terms=NULL, in.terms=NULL, modelSize=NULL, Cox.Hazard = FALSE){

  X = as.matrix(X)

  if(Cox.Hazard==TRUE){
    Y = as.matrix(Y)
    if(ncol(Y)!=2){
      stop('You choose cox family! Therefore you must provide the survival times in column 1 and a binary vector on column 2!')
    }
    if(length(unique(Y[,2]))!=2){
      stop('You choose cox family! Therefore your Y must contain binary responses!')
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
      iter.names = c(iter.names,paste(in.terms[ii,1],"-",in.terms[ii,2]))
    }
  }

  if(is.null(in.terms)){
    X.ITER = iter.names = NULL
  }

  X.full = cbind(X[,list.reduction],X.SQ,X.ITER)
  setSelected = c(list.reduction,SQ.names,iter.names)

  s=ncol(X.full)
  if(is.null(modelSize)){ modelSize=min(5,length(setSelected)) } ## check this with Heather

  ### Error message
  if(modelSize>7){
    stop('Sorry, this version only support model sizes<8')
  }

  n=nrow(X.full)
  if(Cox.Hazard==TRUE){
    residBig = residuals(coxph(Surv(Y[,1], Y[,2])~X.full))
  } else{
    residBig = residuals(glm(Y~X.full, family = family))
  }
  RSSBig = t(residBig) %*%residBig
  dfBig = n-(length(setSelected)+1)

  goodModels = list()

  for (j in 1:modelSize){
    combinationMatrix=t(combn(1:length(setSelected),j))
    combinationMatrixNames=cbind(t(combn(setSelected,j)))
    logicFitVectorF=matrix(0,nrow(combinationMatrix))
    for(l in 1:nrow(combinationMatrix)){
      XSelect = X.full[,combinationMatrix[l,]]
      if(Cox.Hazard==TRUE){
        residSmall = residuals(coxph(Surv(Y[,1], Y[,2])~XSelect))
      } else{
        residSmall = residuals(glm(Y~XSelect, family = family))
      }

      RSSSmall = t(residSmall) %*% residSmall
      diffRSS = RSSSmall - RSSBig
      dfSmall = n-(ncol(cbind(XSelect))+1)
      fCrit = qf(1-signif,df1 = (dfSmall-dfBig),df2 = dfBig)
      if(is.nan(fCrit)){
        logicFitVectorF[l]=1
      } else if((diffRSS/(dfSmall-dfBig))/(RSSBig/dfBig)<=fCrit){
        logicFitVectorF[l]=1
      }
    }
    if(j==1){goodModels[[paste('Model Size', j)]] = as.matrix(combinationMatrixNames[which(logicFitVectorF>0),])
    } else{
      goodModels[[paste('Model Size', j)]] = combinationMatrixNames[which(logicFitVectorF>0),]
    }
  }

  return(list("goodModels"=goodModels))

}
