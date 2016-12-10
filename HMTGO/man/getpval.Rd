\name{getpval}
\alias{getpval}
\title{Computing p-value for tree structure}
\usage{
getpval(tData,treedata,comparelevel,method="globaltest")
}

\description{
This function could compute p-value of a tree structure based on testing data. Only support globaltest for now. Further test method would be added if needed. It will return a vector of p-value.
}

\arguments{
  \item{tData}{Gene expression matrix, gene by subject.}
  \item{treedata}{ A list of tree structure data, created by \link{buildTree} or use \link{treedata} built in this package.}
  \item{method}{Method for computing p-value, only support "globaltest" in this version.}
  \item{comparelevel}{A vector indicating the group labels.}
}

\details{
\code{tData} is the testing data. It should be a matrix with genes names being row names and each column represents a subject.
}

\value{
A vector of p-values for each tree node.
}

\examples{
# to run the following example, install R package ALL from Bioconductor
# library(ALL)
# data(ALL)
# data(treedata)
# pval<-getpval(exprs(ALL), treedata, comparelevel, method="globaltest")

}