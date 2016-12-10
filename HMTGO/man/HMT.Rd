\name{HMT}
\alias{HMT}
\alias{HMTEMC}
\alias{HMT_DAEMC}
\alias{HMT_RSEMC}
\title{EM Step of HMT approach}
\usage{
HMT(treedata, pval, method="Both", control=list())
}

\description{
This is the last step for the hidden Markov tree model, EM algorithm. Function returns a list of final Likelihood, Parameter, probability of differential expression (PDE).
}

\arguments{
  \item{treedata}{
  Tree structure. Created by \link{buildTree}.
  }
  \item{pval}{
  p-value for all tree structure nodes. Could be created by \link{getpval}
  }
  \item{method}{
  Method for running EM algorithm. Could be set as "RS" (Random Start) or "DA" Deterministic Annealing.
  }
  \item{control}{
A list including all other parameters for running EM algorithm.
  }
}

\details{
Note that the names of following arguments in control list must be matched exactly but the order does not a matter.

list(schedule, init_origin, penalty, null_dist, num_method, iternum_final, iternum_hot, iternum_start, iternum_inner, conv_threshold).

schedule. This argument is the temprature schedule for deterministic annealing approach. it should be a numeric array. Starting from a little number that close to but great than 0 and end with 1. Default setting is c(1e-9, 1e-6, 0.001, 0.01, 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)

init_origin is the initial value of all parameter of p-value distribution. (alpha, beta, q, q0, alpha_null, beta_null, lambda). Default setting is c(0.3, 10, 0.5, 0.5, 1.1, 1.1, 0.5). Note that alpha_null should be locate between 1 and 20, beta_null should be locate between 1 and 20, lambda should be locate between 0 and 1.

penalty is used in the application for numerical stability. Should be True or False.

null_dist indicate the distribution under under null hypothesis. Should be one of "Mixture", "Beta" and "Unif".

num_method is the name of optimization method, which should be one of "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" or "Brent". Default value is "Nelder-Mead". Note that this should be a character.

iternum_final is the iteration number of the final/last run for both Deterministic Annealing and Random Start approaches. Default value is 1000.

iternum_hot is the iteration number of non-final runs for Deterministic Annealing approach. Default value is 100.

iternum_start is the iteration number of starting points for Random Start approach. Default value is 5.

iternum_inner is the iteration number for optimization function. Default value is 1000.

conv_threshold is the threshold to determine converge. Default value is 1e-8.

RS_num is the number of random start points for the Random Start approach. Default value is 30.
}

\value{
  \item{LL}{likelihood value of final EM computation}
  \item{par}{Estimated Parameters of hidden Markov tree models}
  \item{cp}{probability of differential expression (PDE) of graph nodes}
  \item{iter}{iteration number of final run. If less than 1000(default setting) or final iteration parameter value, then the EM algorithm converge in the end}
}
