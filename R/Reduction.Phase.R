Reduction.Phase = function(X,Y,family=gaussian, dmHC=NULL,vector.signif=NULL,seed.HC = NULL, Cox.Hazard = FALSE){

  glmRoutine = function(subsetX,Y,family=gaussian,intercept=TRUE,significance=NULL, Cox.Hazard=FALSE){

    n = nrow(subsetX)

    if(intercept==TRUE){
      intercept=1
    } else{
      intercept=0
    }
    count = 0

    if(intercept==0){
      if(Cox.Hazard==TRUE){
        stats = coxph(Surv(Y[,1], Y[,2]) ~ subsetX)
        pVals = coef(summary(stats))[,5]
      } else{
        stats = glm(Y~subsetX-1, family = family)
        pVals = coef(summary(stats))[,4]
      }
    } else{
      if(Cox.Hazard==TRUE){
        stats = coxph(Surv(Y[,1], Y[,2]) ~ subsetX)
        pVals = coef(summary(stats))[,5]
      } else{
        stats = glm(Y~subsetX, family = family)
        pVals = coef(summary(stats))[,4]
        pVals = pVals[-1]
      }
    }

    if(length(which(pVals>0.99))==length(pVals) | length(which(pVals==0))==length(pVals)){
      sTmp=ncol(subsetX)
      sB1=floor(sTmp/2)
      block1=subsetX[,1:sB1]
      block2=subsetX[,(sB1+1):ncol(subsetX)]
      stats = glm(Y~block1, family = family)
      pVals1 = coef(summary(stats))[,4]
      pVals1 = pVals1[-1]
      if(length(which(pVals1>0.99))==length(pVals1)){
        pVals1 = 0.0001*rep(1,length(pVals1))
        count = count+1
      }
      stats = glm(Y~block2, family = family)
      pVals2=coef(summary(stats))[,4]
      pVals2 = pVals2[-1]
      if(length(which(pVals2>0.99))==length(pVals2)){
        pVals2 = 0.0001*rep(1,length(pVals2))
        count = count+1
      }
      if(sTmp>1){pVals=c(pVals1,pVals2)} else{pVals=c(pVals1)}
    }

    if(significance==2){
      idxSelected=which(pVals %in% sort(pVals)[1:2])
    } else if(significance==1){
      idxSelected=which(pVals %in% sort(pVals)[1])
    }  else{
      idxSelected=which(pVals<significance)
    }

    return(list("idxSelected"=idxSelected,"count"=count))
  }

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

  n = nrow(X)
  d = ncol(X)

  ### Hypercube dimension Not specified by the user
  if(is.null(dmHC)){
    v=2:5
    dmHC = v[ceiling(d^(1/v)) <= 15 & ceiling(d^(1/v))>=10]
    if(dmHC==0){
      stop('Hypercube dimension was not specified and conditioned to the number of variables in the design matrix could not the calculated.','\n',
           'Consider delete/include/transform some variables...')
    }
  }

  ### Error messages
  if(dmHC>5){
    stop('Sorry, this version only support cube with maximum 5 dimensions! More dimensions will be avaliable in the next version...')
  }

  if(!is.null(dmHC) & !is.null(vector.signif) & (length(vector.signif)+1)!=dmHC){
    stop(paste('The lenght of vector.signif is not the same as the HC dimension',dmHC,'.'))
  }

  ### Signif vectors
  if(is.null(vector.signif)){
    signif.Default =TRUE
  } else{
    signif.Default =FALSE
    vector.signif = c(NA,vector.signif)
  }

  highest.dmHC = dmHC

  ### Outputs
  Matrix.Selection = list()
  List.Selection = list()

  aux.dmHC5 = aux.dmHC4 = aux.dmHC3 = aux.dmHC2 = 'N'

  ########## case in which dmHC=5 ##########
  if(dmHC==5){

    if(signif.Default==TRUE & highest.dmHC==5){
      aux.signif = 2
    } else{
      aux.signif = vector.signif[5]
    }

    dimHC5 = ceiling(d^(1/5))
    nearestHypercube5=dimHC5^5
    remainderHC5=nearestHypercube5-d

    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      hypercube5=array(sample(c((1:d),rep(0,remainderHC5))),dim=c(dimHC5,dimHC5,dimHC5,dimHC5,dimHC5))
    } else{
      hypercube5=array(c((1:d),rep(0,remainderHC5)),dim=c(dimHC5,dimHC5,dimHC5,dimHC5,dimHC5))
    }

    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC5,dimHC5,dimHC5,dimHC5,dimHC5))

    for(ind5 in 1:dim(hypercube5)[5]){
      for(ind4 in 1:dim(hypercube5)[4]){
        for(ind3 in 1:dim(hypercube5)[3]){
          for(indR in 1:dim(hypercube5)[1]){
            if(length(which(hypercube5[indR,,ind3,ind4,ind5]!=0))>0){
              idx.aux = which(hypercube5[indR,,ind3,ind4,ind5]!=0)
              subsetX = cbind(X[,hypercube5[indR,idx.aux,ind3,ind4,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[indR,idx.aux[idxSelected$idxSelected],ind3,ind4,ind5]=hypercubeSelect[indR,idx.aux[idxSelected$idxSelected],ind3,ind4,ind5]+1
              }
            }
          } # indR
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[,indC,ind3,ind4,ind5]!=0))>0){
              idx.aux = which(hypercube5[,indC,ind3,ind4,ind5]!=0)
              subsetX = cbind(X[,hypercube5[idx.aux,indC,ind3,ind4,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[idx.aux[idxSelected$idxSelected],indC,ind3,ind4,ind5]=hypercubeSelect[idx.aux[idxSelected$idxSelected],indC,ind3,ind4,ind5]+1
              }
            }
          } # indC
        } # ind3
        for(indR in 1:dim(hypercube5)[1]){
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[indR,indC,,ind4,ind5]!=0))>0){
              idx.aux = which(hypercube5[indR,indC,,ind4,ind5]!=0)
              subsetX = cbind(X[,hypercube5[indR,indC,idx.aux,ind4,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[indR,indC,idx.aux[idxSelected$idxSelected],ind4,ind5]=hypercubeSelect[indR,indC,idx.aux[idxSelected$idxSelected],ind4,ind5]+1
              }
            }
          } # indR
        } # indC
      } # ind4
      for(ind3 in 1:dim(hypercube5)[3]){
        for(indR in 1:dim(hypercube5)[1]){
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[indR,indC,ind3,,ind5]!=0))>0){
              idx.aux = which(hypercube5[indR,indC,ind3,,ind5]!=0)
              subsetX = cbind(X[,hypercube5[indR,indC,ind3,idx.aux,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[indR,indC,ind3,idx.aux[idxSelected$idxSelected],ind5]=hypercubeSelect[indR,indC,ind3,idx.aux[idxSelected$idxSelected],ind5]+1
              }
            }
          } # indC
        } # indR
      } # ind3
    } # ind5
    # now traverse in the 5th dimension
    for(ind4 in 1:dim(hypercube5)[4]){
      for(ind3 in 1:dim(hypercube5)[3]){
        for(indR in 1:dim(hypercube5)[1]){
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[indR,indC,ind3,ind4,]!=0))>0){
              idx.aux = which(hypercube5[indR,indC,ind3,ind4,]!=0)
              subsetX = cbind(X[,hypercube5[indR,indC,ind3,ind4,idx.aux]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[indR,indC,ind3,ind4,idx.aux[idxSelected$idxSelected]]=hypercubeSelect[indR,indC,ind3,ind4,idx.aux[idxSelected$idxSelected]]+1
              }
            }
          } # indC
        } # indR
      } # ind3
    }

    setSelected5Times5 = hypercube5[which(hypercubeSelect>4, arr.ind = TRUE)]
    setSelected4Times5 = hypercube5[which(hypercubeSelect>3, arr.ind = TRUE)]
    setSelected3Times5 = hypercube5[which(hypercubeSelect>2, arr.ind = TRUE)]
    setSelected2Times5 = hypercube5[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times5 = hypercube5[which(hypercubeSelect>0, arr.ind = TRUE)]

    numSelected5times5=length(setSelected5Times5)
    numSelected4times5=length(setSelected4Times5)
    numSelected3times5=length(setSelected3Times5)
    numSelected2times5=length(setSelected2Times5)
    numSelected1times5=length(setSelected1Times5)

    Matrix.Selection[[paste('Hypercube with dim',dmHC)]] = c(numSelected1times5,numSelected2times5,numSelected3times5,numSelected4times5,numSelected5times5)
    names(Matrix.Selection[[paste('Hypercube with dim',dmHC)]]) = c('numSelected1','numSelected2','numSelected3','numSelected4','numSelected5')

    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected1')]] = setSelected1Times5
    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected2')]] = setSelected2Times5
    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected3')]] = setSelected3Times5
    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected4')]] = setSelected4Times5
    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected5')]] = setSelected5Times5

    aux.dmHC5 <- 'Y'

  }

  ########## case in which dmHC=4 ##########
  if(dmHC==4 | aux.dmHC5=='Y'){

    if(signif.Default==TRUE & highest.dmHC==4){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>4){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[4]
    }

    if(dmHC==4){
      dimHC = ceiling(d^(1/4))
      nearestHypercube=dimHC^4
      remainderHC=nearestHypercube-d
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC,dimHC))
      } else{
        hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC,dimHC))
      }} else if(dmHC>4){
        dimHC = ceiling(length(setSelected3Times5)^(1/4))
        nearestHypercube=dimHC^4
        remainderHC=nearestHypercube-length(setSelected3Times5)
        if(!is.null(seed.HC)){
          set.seed(seed.HC)
          hypercube=array(sample(c(setSelected3Times5,rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC,dimHC))
        } else{
          hypercube=array(c(setSelected3Times5,rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC,dimHC))
        }
      }

    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC,dimHC))

    for(ind4 in 1:dim(hypercube)[4]){
      for(ind3 in 1:dim(hypercube)[3]){
        for(indR in 1:dim(hypercube)[1]){
          if(length(which(hypercube[indR,,ind3,ind4]!=0))>0){
            idx.aux = which(hypercube[indR,,ind3,ind4]!=0)
            subsetX = cbind(X[,hypercube[indR,idx.aux,ind3,ind4]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,idx.aux[idxSelected$idxSelected],ind3,ind4]=hypercubeSelect[indR,idx.aux[idxSelected$idxSelected],ind3,ind4]+1
            }
          }
        } # indR
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[,indC,ind3,ind4]!=0))>0){
            idx.aux = which(hypercube[,indC,ind3,ind4]!=0)
            subsetX = cbind(X[,hypercube[idx.aux,indC,ind3,ind4]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[idx.aux[idxSelected$idxSelected],indC,ind3,ind4]=hypercubeSelect[idx.aux[idxSelected$idxSelected],indC,ind3,ind4]+1
            }
          }
        } # indC
      } # ind3
      for(indR in 1:dim(hypercube)[1]){
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[indR,indC,,ind4]!=0))>0){
            idx.aux = which(hypercube[indR,indC,,ind4]!=0)
            subsetX = cbind(X[,hypercube[indR,indC,idx.aux,ind4]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,indC,idx.aux[idxSelected$idxSelected],ind4]=hypercubeSelect[indR,indC,idx.aux[idxSelected$idxSelected],ind4]+1
            }
          }
        } # indC
      } # indR
    } # ind4
    # now traverse in the 4th dimension
    for(ind3 in 1:dim(hypercube)[3]){
      for(indR in 1:dim(hypercube)[1]){
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[indR,indC,ind3,]!=0))>0){
            idx.aux = which(hypercube[indR,indC,ind3,]!=0)
            subsetX = cbind(X[,hypercube[indR,indC,ind3,idx.aux]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,indC,ind3,idx.aux[idxSelected$idxSelected]]=hypercubeSelect[indR,indC,ind3,idx.aux[idxSelected$idxSelected]]+1
            }
          }
        } # indC
      } # indR
    } # ind3

    setSelected4Times = hypercube[which(hypercubeSelect>3, arr.ind = TRUE)]
    setSelected3Times = hypercube[which(hypercubeSelect>2, arr.ind = TRUE)]
    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]

    numSelected4times5=length(setSelected4Times)
    numSelected3times5=length(setSelected3Times)
    numSelected2times5=length(setSelected2Times)
    numSelected1times5=length(setSelected1Times)

    Matrix.Selection[[paste('Hypercube with dim',4)]] = c(numSelected1times5,numSelected2times5,numSelected3times5,numSelected4times5)
    names(Matrix.Selection[[paste('Hypercube with dim',4)]]) = c('numSelected1','numSelected2','numSelected3','numSelected4')

    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected1')]] = setSelected1Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected2')]] = setSelected2Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected3')]] = setSelected3Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected4')]] = setSelected4Times


    aux.dmHC4 <- 'Y' #readline(cat("Reduction of dimension 4 done!", "\n", length(setSelected4Times5),
    #        "Variables selected at least 4 times","\n",length(setSelected3Times5),
    #        "Variables selected at least 3 times","\n",length(setSelected2Times5),
    #       "Variables selected at least 2 times","\n", "Wanna proceed with reduction?[Y/N]"))
  }

  ########## case in which dmHC=3 ##########
  if(dmHC==3 | aux.dmHC4=='Y'){

    if(signif.Default==TRUE & highest.dmHC==3){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>3){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[3]
    }

    if(dmHC==3){
      dimHC = ceiling(d^(1/3))
      nearestHypercube=dimHC^3
      remainderHC=nearestHypercube-d
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC))
      } else{
        hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC))
      }} else if(dmHC>3){
        dimHC = ceiling(length(setSelected2Times)^(1/3))
        nearestHypercube=dimHC^3
        remainderHC=nearestHypercube-length(setSelected2Times)
        if(!is.null(seed.HC)){
          set.seed(seed.HC)
          hypercube=array(sample(c(setSelected2Times,rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC))
        } else{
          hypercube=array(c(setSelected2Times,rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC))
        }
      }

    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC))

    for(indL in 1:dim(hypercube)[3]){
      for(indR in 1:dim(hypercube)[1]){
        if(length(which(hypercube[indR,,indL]!=0))>0){
          idx.aux = which(hypercube[indR,,indL]!=0)
          subsetX = cbind(X[,hypercube[indR,idx.aux,indL]])
          idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,idx.aux[idxSelected$idxSelected],indL]=hypercubeSelect[indR,idx.aux[idxSelected$idxSelected],indL]+1
          }
        }
      } # indR
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[,indC,indL]!=0))>0){
          idx.aux = which(hypercube[,indC,indL]!=0)
          subsetX = cbind(X[,hypercube[idx.aux,indC,indL]])
          idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[idx.aux[idxSelected$idxSelected],indC,indL]=hypercubeSelect[idx.aux[idxSelected$idxSelected],indC,indL]+1
          }
        }
      } # indC
    } # indL
    for(indR in 1:dim(hypercube)[1]){
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[indR,indC,]!=0))>0){
          idx.aux = which(hypercube[indR,indC,]!=0)
          subsetX = cbind(X[,hypercube[indR,indC,idx.aux]])
          idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,indC,idx.aux[idxSelected$idxSelected]]=hypercubeSelect[indR,indC,idx.aux[idxSelected$idxSelected]]+1
          }
        }
      } # indC
    } # indR

    setSelected3Times = hypercube[which(hypercubeSelect>2, arr.ind = TRUE)]
    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]

    numSelected3times=length(setSelected3Times)
    numSelected2times=length(setSelected2Times)
    numSelected1times=length(setSelected1Times)

    Matrix.Selection[[paste('Hypercube with dim',3)]] = c(numSelected1times,numSelected2times,numSelected3times)
    names(Matrix.Selection[[paste('Hypercube with dim',3)]]) = c('numSelected1','numSelected2','numSelected3')

    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected1')]] = setSelected1Times
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected2')]] = setSelected2Times
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected3')]] = setSelected3Times

    aux.dmHC3 <- 'Y'
  }

  ########## case in which dmHC=2 ##########
  if(dmHC==2 | aux.dmHC3=='Y'){

    if(signif.Default==TRUE & highest.dmHC==2){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>2){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[2]
    }

    if(dmHC==2){
      dimHC = ceiling(d^(1/2))
      nearestHypercube=dimHC^2
      remainderHC=nearestHypercube-d
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC))
      } else{
        hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC))
      }} else if(dmHC>2){
        dimHC = ceiling(numSelected1times^(1/2))
        nearestHypercube=dimHC^2
        remainderHC=nearestHypercube-length(setSelected2Times)
        if(!is.null(seed.HC)){
          set.seed(seed.HC)
          hypercube=array(sample(c(setSelected2Times,rep(0,remainderHC))),dim=c(dimHC,dimHC))
        } else{
          hypercube=array(c(setSelected2Times,rep(0,remainderHC)),dim=c(dimHC,dimHC))
        }
      }

    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC))

    for(indR in 1:dim(hypercube)[1]){
      if(length(which(hypercube[indR,]!=0))>0){
        idx.aux = which(hypercube[indR,]!=0)
        subsetX = cbind(X[,hypercube[indR,idx.aux]])
        idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[indR,idx.aux[idxSelected$idxSelected]]=hypercubeSelect[indR,idx.aux[idxSelected$idxSelected]]+1
        }
      }
    } # indR
    for(indC in 1:dim(hypercube)[2]){
      if(length(which(hypercube[,indC]!=0))>0){
        idx.aux = which(hypercube[,indC]!=0)
        subsetX = cbind(X[,hypercube[idx.aux,indC]])
        idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif, Cox.Hazard = Cox.Hazard)

        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[idx.aux[idxSelected$idxSelected],indC]=hypercubeSelect[idx.aux[idxSelected$idxSelected],indC]+1
        }
      }
    } # indC

    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]

    numSelected2times=length(setSelected2Times)
    numSelected1times=length(setSelected1Times)

    Matrix.Selection[[paste('Hypercube with dim',2)]] = c(numSelected1times,numSelected2times)
    names(Matrix.Selection[[paste('Hypercube with dim',2)]]) = c('numSelected1','numSelected2')

    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected1')]] = setSelected1Times
    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected2')]] = setSelected2Times

    aux.dmHC2 <- 'Y' #readline(cat("Reduction of dimension 2 done!", "\n", length(setSelected2Times),
    #"Variables selected at least 2 times","\n",length(setSelected1Times),
    #"Variables selected at least 1 times","\n", "Wanna proceed with reduction?[Y/N]"))
  }

  return(list("Matrix.Selection" = Matrix.Selection, "List.Selection" = List.Selection))
}
