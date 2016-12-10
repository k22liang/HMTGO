HMT<-function(treedata,pval,method="Both",control=list())
{
#############check tree structure##########################

  if(length(treedata$childrenCount)!=treedata$nTree)
  {
    stop("childrenCount did not match the tree node number!")
  }
  if(ncol(treedata$childrenIndexMat)!=treedata$nTree)
  {
    stop("childrenIndexMat did not match the tree node number!")
  }
  if(ncol(treedata$childrenIndexMat)!=treedata$nTree)
  {
    stop("childrenIndexMat did not match the tree node number!")
  }
  if(length(treedata$depth)!=treedata$nTree)
  {
    stop("depth didnot match the tree node number")
  }
  if(length(treedata$parentIndex)!=treedata$nTree)
  {
    stop("parentIndex didnot match the tree node number")
  }
  if(length(treedata$propInParent)!=treedata$nTree)
  {
    stop("propInParent match the tree node number")
  }
  if(ncol(treedata$componMappingMat)!=treedata$nTree)
  {
    stop("componMappingMat didnot match the tree node number")
  }
  if(nrow(treedata$componMappingMat)!=treedata$nDAG)
  {
    stop("componMappingMat didnot match the Merge number")
  }

################check method#######################
  if( !((method=="DA")|(method=="RS") | (method=="Both")) )
  {
    stop("wrong method, should be DA, RS or Both")
  }
################check control parameter############
  if (is.null(control$penalty))
  {
    penalty=FALSE
  }else if((control$penalty==T) | (control$penalty==F))
  {
    penalty=control$penalty
  }else
  {
    stop("wrong control$penalty, should be Ture or False")
  }

  if (is.null(control$null_dist))
  {
    null_dist="Mixture"
  }else if((control$null_dist == "Mixture" |control$null_dist == "Beta") | (control$null_dist == "Unif"))
  {
    null_dist=control$null_dist
  }else
  {
    stop("wrong control$null_dist, should be Mixture,Unif or Beta")
  }

  if (is.null(control$num_method))
  {
    num_method= "Nelder-Mead"
  }else if((control$num_method== "Nelder-Mead")|(control$num_method== "BFGS")|(control$num_method== "CG")|
          (control$num_method=="L-BFGS-B")|(control$num_method== "SANN")|(control$num_method=="Brent"))
  {
    num_method=control$num_method
  }else
  {
    stop("wrong control$num_method, should be Nelder-Mead, BFGS, CG, L-BFGS-B, SANN or Brent")
  }

  if (is.null(control$iternum_final))
  {
    iternum_final=1000
  }else if(is.numeric(control$iternum_final) & (control$iternum_final>0))
  {
    if (control$iternum_final%%1==0)
    {
      iternum_final=control$iternum_final
    }else
    {
        stop("wrong control$iternum_final, should be a positive integer")
    }
  }else
  {
    stop("wrong control$iternum_final, should be a positive integer")
  }

  if (is.null(control$iternum_inner))
  {
    iternum_inner=1000
  }else if(is.numeric(control$iternum_inner) & (control$iternum_inner>0))
  {
    if (control$iternum_inner%%1==0)
    {
      iternum_inner=control$iternum_inner
    }else
    {
      stop("wrong control$iternum_inner, should be a positive integer")
    }
  }else
  {
    stop("wrong control$iternum_inner, should be a positive integer")
  }

  if (is.null(control$iternum_hot))
  {
    iternum_hot=10
  }else if(is.numeric(control$iternum_hot) & (control$iternum_hot>0))
  {
    if (control$iternum_hot%%1==0)
    {
      iternum_hot=control$iternum_hot
    }else
    {
      stop("wrong control$iternum_hot, should be a positive integer")
    }
  }else
  {
    stop("wrong control$iternum_hot, should be a positive integer")
  }


  if (is.null(control$iternum_start))
  {
    iternum_start=5
  }else if(is.numeric(control$iternum_start) & (control$iternum_start>1))
  {
    if (control$iternum_start%%1==0)
    {
      iternum_start=control$iternum_start
    }else
    {
      stop("wrong control$iternum_start, should be a positive integer and greater than 1")
    }
  }else
  {
    stop("wrong control$iternum_start, should be a positive integer and greater than 1")
  }

  if (is.null(control$RS_threshold))
  {
    RS_threshold=0.01
  }else if(is.numeric(control$RS_threshold) & (control$RS_threshold>0)&(control$RS_threshold<1))
  {
    RS_threshold=control$RS_threshold
  }else
  {
    stop("wrong control$RS_threshold, should be a positive number less than 1")
  }
  if (is.null(control$RS_num))
  {
    RS_num=30
  }else if(is.numeric(control$RS_num) & (control$RS_num>0))
  {
    if (control$RS_num%%1==0)
    {
      RS_num=control$RS_num
    }else
    {
      stop("wrong control$RS_num, should be an positive integer")
    }
  }else
  {
    stop("wrong control$RS_num, should be an integer")
  }

  if (is.null(control$conv_threshold))
  {
    conv_threshold=1e-8
  }else if(is.numeric(control$conv_threshold) & (control$conv_threshold>0))
  {
    conv_threshold=control$conv_threshold
  }else
  {
    stop("wrong control$conv_threshold??should be a positive number")
  }

  if (is.null(control$schedule))
  {
    schedule <- c(1e-9, 1e-6, 0.001, 0.01, .02, .04, .1, seq(.2, .9, by=.1), .95, 1)
  }else if(is.numeric(control$schedule))
  {
    if(control$schedule[length(control$schedule)]!=1)
    {
      stop("wrong control$schedule??should be end with 1")
    }
    if(control$schedule[0]<=0)
    {
      stop("wrong control$schedule??should start with a very small positive number")
    }

    if(sum(sort(control$schedule)!=control$schedule)!=0)
    {
      stop("wrong control$schedule??should be in an ascending order")
    }
    schedule<-control$schedule
  }else
  {
    stop("wrong control$schedule??should be a numeric vector")
  }

  if (is.null(control$init_origin))
  {
    init_origin <- c(0.3, 10, 0.5, 0.5, 1.1, 1.1, 0.5)
  }else if(is.numeric(control$init_origin))
  {
    init_origin<-control$init_origin
    if( (init_origin[5]<1) |(init_origin[5]>20))
    {
      stop("wrong control$init_origin, alpha_null should be locate between 1 and 20")
    }
    if((init_origin[6]<1 )|(init_origin[6]>20))
    {
      stop("wrong control$init_origin, beta_null should be locate between 1 and 20")
    }
    if((init_origin[7]<0) |(init_origin[5]>1))
    {
      stop("wrong control$init_origin, lambda should be locate between 0 and 1")
    }
  }else
  {
    stop("wrong control$init_origin??should be a numeric vector")
  }

######################################################

  if(method=="DA")
  {
    result<-HMT_DAEMC(pval, treedata, optim,getNegAltLlikeC,getNegNullLlikeC,schedule,init_origin, penalty = penalty, null_dist = null_dist,
                     num_method = num_method,iternum_final = iternum_final,
                     iternum_hot = iternum_hot, iternum_inner=iternum_inner, conv_threshold = conv_threshold)
    return(result)
  }else if(method=="RS")
  {
      return(HMT_RSEM(pval, treedata, RS_num, null_dist, penalty, num_method, iternum_start, iternum_final))
  }else if(method=="Both")
  {
    #DA
    result_DA<-HMT_DAEMC(pval, treedata, optim,getNegAltLlikeC,getNegNullLlikeC,schedule,init_origin, penalty = penalty, null_dist = null_dist,
                      num_method = num_method,iternum_final = iternum_final,
                      iternum_hot = iternum_hot, iternum_inner=iternum_inner, conv_threshold = conv_threshold)

    #RS
    result_RS <- HMT_RSEM(pval, treedata, RS_num, null_dist, penalty, num_method, iternum_start, iternum_final)
    
    if(result_DA$LL>result_RS$LL){
      result<-result_DA
      result$Method<-"DA"
    }else{
      result<-result_RS
      result$Method<-"RS"
    }
    result$LL_dif <- abs(result_DA$LL-result_RS$LL)

    return(result)
  }else
  {
    stop("undefind method!")
  }
}

HMT_RSEM <- function(pval, treedata, RS_num, null_dist, penalty, num_method, iternum_start, iternum_final){
    initial_value<-sapply(1:RS_num, function(x){
        result.temp<-list(LL=-Inf)
        while(result.temp$LL==-Inf){
            initial_value_temp<-c(runif(1, 0.01, 0.8),runif(1, 1.01, 15),runif(1, 0.01, 0.99),runif(1, 0.88, 0.99),
                                  runif(1, 1.01, 5),runif(1, 1.01, 9),runif(1, 0.01, 0.99))
            result.temp<-tryCatch(
                HMTEMC(optim,getNegAltLlikeC,getNegNullLlikeC,initial_value_temp,pval,
                       treedata$depth,treedata$leaf,treedata$nonleaf,treedata$childrenIndexMat,
                       treedata$propInParent,treedata$childrenCount,treedata$parentIndex,
                       treedata$nDAG,null_dist = null_dist,penalty = penalty,
                       num_method = num_method,max_it = iternum_start,innerit = 200,conv_threshold = 1e-8)
                ,error=function(y){list(LL=-Inf)})
        }
        return(c(result.temp$LL,result.temp$par))
    }
    )
    LL0<-initial_value[1,]
    initial_value<-t(initial_value[-1,])
    
    result_RS <- HMTEMC(optim,getNegAltLlikeC, getNegNullLlikeC, initial_value[which.max(LL0),], pval,
                        treedata$depth,treedata$leaf,treedata$nonleaf,treedata$childrenIndexMat,
                        treedata$propInParent,treedata$childrenCount,treedata$parentIndex,
                        treedata$nDAG,null_dist = null_dist,penalty = penalty,
                        num_method = num_method,max_it = iternum_final,innerit = 1000,conv_threshold = 1e-8)
    
    result_RS$PDE <- getPDEC(result_RS$cp, treedata$parentIndex, treedata$componMappingMat, treedata$sepLocation, treedata$formula)
    return(result_RS)
}

    
    

