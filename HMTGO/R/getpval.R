getpval=function(tData,treedata,comparelevel,method="globaltest")
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


pval=rep(0,treedata$nTree)

if(method=="globaltest")
{
  Y=as.factor(comparelevel)
  for(i in 1:treedata$nTree){
    testdata <- tData[treedata$geneSet[[i]],]
    if(is.null(nrow(testdata))){
      testdata<-t(as.matrix(testdata))
      rownames(testdata)<-treedata$geneSet[[i]]

    }
    pval[i]<-p.value(gt(Y,t(testdata)))
  }
  return(pval)
}
else
{
  stop("unrecognized method")
}
}
