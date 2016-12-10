\name{buildTree}
\alias{buildTree}
\title{
Build the tree structure
}
\description{
The tree transformation step for the hidden Markov tree model.}
\usage{
buildTree(root = "GO:0008150", DAG, geneSet)
}
\arguments{
  \item{root}{
The root node of the tree.
}
  \item{DAG}{
The DAG structure returned by getDAG.
}
  \item{geneSet}{
The gene sets returned by getGeneSet.
}
}
\details{
This function returns a list of data to represent tree structure.
}
\value{
  Tree data created by this function contains following argument:
  \item{nDAG}{The number of original GO DAG nodes}
  \item{nTree}{The number of tree nodes}
  \item{parentIndex}{A vector that records parent-child relationship. The ith value is the parent index of the ith tree node.}
  \item{childrenIndexMat}{A matrix that records parent-child relationship. The ith column contains the index of children of the ith tree node.}
  \item{childrenCount}{A vector of the number of children of tree nodes.}
  \item{componMappingMat}{A nDAG by nTree matrix. Each index of row represent a node in orginal graph and each index of column represent a node of tree structure. The value in each cell should be either 0 or 1. 1 means the relative tree nodes is a subset of corresponding graph node.}
  \item{leaf}{A vector that record the index of all leaf nodes}
  \item{nonleaf}{A vector that record the index of all non-leaf nodes}
  \item{formula, sepLocation}{Lists that used to generate PDE}
  \item{geneSet}{List that records all genes in each tree node}

}

\seealso{
See Also \code{\link{getDAG}}, \code{\link{getGeneSet}}.
}
\examples{
library(GO.db)
dag <- getDAG(listChildren=GO.db::GOBPCHILDREN)

# to run the following example, install R package hgu95av2.db from Bioconductor
# library(hgu95av2.db)
# gs <- getGeneSet(hgu95av2.db::hgu95av2GO2ALLPROBES, min.genes=1, filter="IEA")

# treedata <- buildTree(root="GO:0008150", DAG=dag, geneSet=gs)

}
