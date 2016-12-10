\name{getGeneSet}
\alias{getGeneSet}
\title{
Obtain gene sets
}
\description{
Obtain a list of gene sets corresponding to GO terms.
}
\usage{
getGeneSet(listGene, min.genes = 1, filter = "IEA", available.genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listGene}{
A list of gene sets with each gene annotated by its evidence code.
}
  \item{min.genes}{
Minimum number of genes each set to be considered a set
}
  \item{filter}{
Optionly gene annotation filter critetia. By default, we filter out evidence code of 'IEA', which is not annotated by human curators.
}
  \item{available.genes}{
Optionly gene set to indicate available genes.
}
}
\details{
For an existing list of gene sets, performs selection of size (min.genes), annotation (filter) and gene (available.genes).
}
\value{
a list of gene sets indexed by GO term IDs.
}

\seealso{
See Also \code{\link{getDAG}}, \code{\link{buildTree}}.
}
\examples{
# to run the following example, install R package hgu95av2.db from Bioconductor
# library(hgu95av2.db)
# gs <- getGeneSet(hgu95av2.db::hgu95av2GO2ALLPROBES, min.genes=1, filter="IEA")

}

