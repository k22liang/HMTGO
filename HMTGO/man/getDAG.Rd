\name{getDAG}
\alias{getDAG}
\title{Obtain GO DAG}
\usage{
getDAG(listChildren=GO.db::GOBPCHILDREN)
}

\description{
Obtain GO DAG.
}

\arguments{
  \item{listChildren}{
  A list of children nodes.
  }
}

\details{
  This function return a list of children nodes to represent GO DAG structure. For proper subset relationship, we only keep the 'part_of' and 'is_a' relationships.
}

\value{
  A list of children nodes.
}
